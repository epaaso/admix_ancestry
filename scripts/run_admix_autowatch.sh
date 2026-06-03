#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STUDY_GLOB="${STUDY_GLOB:-/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz}"
STUDY_DIR="${STUDY_DIR:-/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names}"
PRIMARY_OUT="${PRIMARY_OUT:-/datos/migccl/leukemia_exoma/outs-results/005-admix_whole}"
FALLBACK_OUT="${FALLBACK_OUT:-/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_no_hm3}"
REF_POPS="${REF_POPS:-EUR:/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz;AFR:/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz}"
HM3="${HM3:-/datos/home/epaaso/ancestry/hm3_combined_ids_and_loci.txt}"
EXCLUDE_REGIONS="${EXCLUDE_REGIONS:-/datos/home/epaaso/ancestry/leukemia_genes_hg38.bed}"
GRCH38="${GRCH38:-/datos/migccl/arriba_refs/GRCh38.primary_assembly.genome.fa}"
MAIL_TO="${MAIL_TO:-ernesto.paas@ciencias.unam.mx}"
AGY_TIMEOUT="15m"
AGY_WALL_TIMEOUT="20m"
STATUS_INTERVAL=60
MIN_PRUNED_MARKERS=100
MAX_TOTAL_REPAIRS=5
MAX_SAME_ERROR=3
RUN_DIR="/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch"
LOG_DIR="${RUN_DIR}/logs"
STATUS="${RUN_DIR}/status.tsv"
LAST_STATUS="${RUN_DIR}/last_status.md"
FINGERPRINTS="${RUN_DIR}/error_fingerprints.tsv"
AGY_LOCK_DIR="${RUN_DIR}/agy_agent.lock"
PRECHECK_ONLY=0
TEST_MAIL=0
STATUS_ONLY=0
NEXTFLOW_PID=""
AGY_PID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --preflight-only) PRECHECK_ONLY=1; shift ;;
    --test-mail) TEST_MAIL=1; shift ;;
    --status-only) STATUS_ONLY=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "$RUN_DIR" "$LOG_DIR"
touch "$STATUS" "$FINGERPRINTS"

cleanup_children() {
  local rc=$?
  if [[ -n "${AGY_PID:-}" ]] && kill -0 "$AGY_PID" >/dev/null 2>&1; then
    kill "$AGY_PID" >/dev/null 2>&1 || true
  fi
  if [[ -n "${NEXTFLOW_PID:-}" ]] && kill -0 "$NEXTFLOW_PID" >/dev/null 2>&1; then
    kill "$NEXTFLOW_PID" >/dev/null 2>&1 || true
  fi
  rmdir "$AGY_LOCK_DIR" >/dev/null 2>&1 || true
  exit "$rc"
}

trap cleanup_children INT TERM

log() {
  printf '[%s] %s\n' "$(date -Is)" "$*" | tee -a "${LOG_DIR}/autowatch.log"
}

latest_nextflow_workdir() {
  local log_file="$1"
  { grep -oE 'workDir: [^]]+|workDir=[^, ]+' "$log_file" 2>/dev/null || true; } |
    sed -E 's/^workDir[:=] //' |
    tail -n 1
}

latest_nextflow_task() {
  local log_file="$1"
  grep -E 'Submitted process >|TaskHandler\[|executor local > tasks to be completed|status:|ERROR|WARN|failed|terminated' "$log_file" 2>/dev/null |
    tail -n 12 || true
}

file_stat_line() {
  local path="$1"
  if [[ -e "$path" ]]; then
    stat -c '%s bytes; modified %y; path=%n' "$path"
  else
    printf 'missing; path=%s\n' "$path"
  fi
}

live_snapshot() {
  local mode="${1:-unknown}"
  local outdir="${2:-unknown}"
  local log_file="${3:-}"
  local workdir=""
  [[ -n "$log_file" ]] && workdir="$(latest_nextflow_workdir "$log_file")"

  {
    printf 'Mode: %s\n' "$mode"
    printf 'Output directory: %s\n' "$outdir"
    printf 'Autowatch run directory: %s\n' "$RUN_DIR"
    printf 'Nextflow log: %s\n' "${log_file:-unknown}"
    printf 'Tmux session/window: 0:admix-auto\n'
    printf 'Repair attempts seen: %s\n' "$(wc -l < "$FINGERPRINTS" 2>/dev/null || echo 0)"
    printf 'mail command available: '
    if have mail; then printf 'yes\n'; else printf 'no\n'; fi
    printf '\n'

    printf 'Live processes:\n'
    ps -eo pid,ppid,stat,etime,pcpu,pmem,cmd |
      awk '/run_admix_autowatch|nextflow run main.nf|java.*nextflow|bcftools|plink|admixture|agy --dangerously/ && !/awk/ && !/--status-only/ {print}' || true
    printf '\n'

    if [[ -n "$log_file" && -s "$log_file" ]]; then
      printf 'Latest Nextflow task lines:\n'
      latest_nextflow_task "$log_file"
      printf '\n'
    fi

    if [[ -n "$workdir" ]]; then
      printf 'Active/last work directory: %s\n' "$workdir"
      printf 'merged_raw.vcf.gz: '
      file_stat_line "${workdir}/merged_raw.vcf.gz"
      printf '.command.log: '
      file_stat_line "${workdir}/.command.log"
      printf '\n'
    fi

    printf 'Expected final outputs:\n'
    printf 'ancestry.tsv: '
    file_stat_line "${outdir}/ancestry/ancestry.tsv"
    printf 'ancestry_plot.png: '
    file_stat_line "${outdir}/plots/ancestry_plot.png"
    printf 'pca_all_samples.png: '
    file_stat_line "${outdir}/plots/pca_all_samples.png"
  }
}

write_status() {
  local state="$1"
  local message="$2"
  {
    printf '# admix_whole autowatch status\n\n'
    printf 'State: %s\n\n' "$state"
    printf 'Updated: %s\n\n' "$(date -Is)"
    printf '%s\n' "$message"
  } > "$LAST_STATUS"
  printf '%s\t%s\t%s\n' "$(date -Is)" "$state" "$message" >> "$STATUS"
}

write_live_status() {
  local state="$1"
  local message="$2"
  local mode="$3"
  local outdir="$4"
  local log_file="$5"
  {
    printf '# admix_whole autowatch status\n\n'
    printf 'State: %s\n\n' "$state"
    printf 'Updated: %s\n\n' "$(date -Is)"
    printf '%s\n\n' "$message"
    printf '## Live Snapshot\n\n'
    live_snapshot "$mode" "$outdir" "$log_file"
  } > "$LAST_STATUS"
  printf '%s\t%s\t%s\n' "$(date -Is)" "$state" "$message" >> "$STATUS"
}

have() {
  command -v "$1" >/dev/null 2>&1
}

ensure_mail() {
  if have mail; then
    return 0
  fi
  log "mail command not found; attempting noninteractive mailutils install"
  if have sudo && sudo -n true >/dev/null 2>&1; then
    sudo apt-get update >>"${LOG_DIR}/mailutils_install.log" 2>&1
    sudo apt-get install -y mailutils >>"${LOG_DIR}/mailutils_install.log" 2>&1
  else
    log "sudo without password is unavailable; cannot install mailutils noninteractively"
  fi
  have mail
}

notify_blocked() {
  local subject="$1"
  local body="$2"
  write_status "BLOCKED" "$body"
  if ensure_mail; then
    printf '%s\n' "$body" | mail -s "$subject" "$MAIL_TO" || true
  elif have wall; then
    printf '%s\n' "$subject: $body" | wall || true
  fi
}

cleanup_orphaned_repair_agents() {
  local reason="$1"
  local killed=0
  local pid cmd

  while read -r pid cmd; do
    [[ -n "${pid:-}" ]] || continue
    # Only kill repair agents launched by this watcher. Do not touch unrelated agy sessions.
    if [[ "$cmd" == *"agy --dangerously-skip-permissions"* &&
          "$cmd" == *"autonomous repair agent for the admix_whole"* ]]; then
      if [[ -n "${AGY_PID:-}" && "$pid" == "$AGY_PID" ]]; then
        continue
      fi
      log "Stopping orphaned admix_whole agy repair agent pid=${pid}; reason=${reason}"
      kill "$pid" >/dev/null 2>&1 || true
      killed=$((killed + 1))
    fi
  done < <(ps -eo pid=,args=)

  if [[ "$killed" -gt 0 ]]; then
    sleep 3
    while read -r pid cmd; do
      [[ -n "${pid:-}" ]] || continue
      if [[ "$cmd" == *"agy --dangerously-skip-permissions"* &&
            "$cmd" == *"autonomous repair agent for the admix_whole"* ]]; then
        log "Force-stopping stubborn admix_whole agy repair agent pid=${pid}; reason=${reason}"
        kill -9 "$pid" >/dev/null 2>&1 || true
      fi
    done < <(ps -eo pid=,args=)
  fi
}

count_glob() {
  local glob="$1"
  compgen -G "$glob" | wc -l
}

preflight() {
  log "Running preflight"
  cd "$ROOT"

  local study_count
  study_count="$(find "$STUDY_DIR" -maxdepth 1 -name '*.vcf.gz' | wc -l)"
  [[ "$study_count" -eq 50 ]] || { echo "Expected 50 study VCFs, found $study_count" >&2; return 1; }

  local ref_count
  ref_count="$(count_glob "/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz")"; [[ "$ref_count" -gt 0 ]] || return 1
  ref_count="$(count_glob "/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz")"; [[ "$ref_count" -gt 0 ]] || return 1
  ref_count="$(count_glob "/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz")"; [[ "$ref_count" -gt 0 ]] || return 1
  ref_count="$(count_glob "/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz")"; [[ "$ref_count" -gt 0 ]] || return 1

  test -s "$HM3"
  test -s "$EXCLUDE_REGIONS"
  test -s "$GRCH38"
  test -s "${GRCH38}.fai"
  test -s "${GRCH38%.fa}.dict" || test -s "${GRCH38%.fasta}.dict"

  have agy
  have nextflow
  have bcftools
  PLINK_BIN="${PLINK_BIN:-/home/epaaso/bin/plink}"
  have plink || test -x "$PLINK_BIN"
  have admixture

  ensure_mail || log "mail remains unavailable; wall/log fallback will be used for blocked notifications"
}

index_study_vcfs() {
  log "Ensuring adjacent indexes for study VCF symlinks"
  local indexed=0
  local skipped=0
  local vcf
  for vcf in "$STUDY_DIR"/*.vcf.gz; do
    if [[ -s "${vcf}.tbi" || -s "${vcf}.csi" ]]; then
      skipped=$((skipped + 1))
    else
      bcftools index -t "$vcf" >>"${LOG_DIR}/index_study_vcfs.log" 2>&1
      indexed=$((indexed + 1))
    fi
  done
  log "Study VCF indexes created: $indexed; already present: $skipped"
}

run_nextflow() {
  local mode="$1"
  local outdir="$2"
  local log_file="${LOG_DIR}/nextflow_${mode}.log"
  local rc=0
  local common_args=(
    --input "$STUDY_GLOB"
    --ref_pops "$REF_POPS"
    --k 3
    --exclude_regions "$EXCLUDE_REGIONS"
    --maf 0.05
    --geno 0.05
    --mind 0.8
    --relaxed_geno 0.99
    --ld_window 50
    --ld_step 5
    --ld_r2 0.2
    --max_cpus 18
    --max_memory "120 GB"
  )
  mkdir -p "$outdir"
  cleanup_orphaned_repair_agents "entering_nextflow_${mode}"
  log "Starting Nextflow mode=${mode} outdir=${outdir}"
  write_live_status "RUNNING" "mode=${mode}; outdir=${outdir}; log=${log_file}" "$mode" "$outdir" "$log_file"
  cd "$ROOT"

  if [[ "$mode" == "primary_hm3" ]]; then
    NXF_LOG_FILE="$log_file" nextflow run main.nf -resume \
      --outdir "$outdir" \
      --extract_snps "$HM3" \
      "${common_args[@]}" >>"$log_file" 2>&1 &
  else
    NXF_LOG_FILE="$log_file" nextflow run main.nf -resume \
      --outdir "$outdir" \
      "${common_args[@]}" >>"$log_file" 2>&1 &
  fi
  NEXTFLOW_PID=$!

  while kill -0 "$NEXTFLOW_PID" >/dev/null 2>&1; do
    write_live_status "RUNNING" "mode=${mode}; outdir=${outdir}; nextflow_pid=${NEXTFLOW_PID}; log=${log_file}" "$mode" "$outdir" "$log_file"
    sleep "$STATUS_INTERVAL"
  done

  set +e
  wait "$NEXTFLOW_PID"
  rc=$?
  set -e
  NEXTFLOW_PID=""

  if [[ "$rc" -eq 0 ]]; then
    write_live_status "NEXTFLOW_EXITED" "mode=${mode}; exit_code=0; validating required outputs" "$mode" "$outdir" "$log_file"
  else
    write_live_status "NEXTFLOW_FAILED" "mode=${mode}; exit_code=${rc}; log=${log_file}" "$mode" "$outdir" "$log_file"
  fi
  return "$rc"
}

validate_outputs() {
  local outdir="$1"
  test -s "${outdir}/ancestry/ancestry.tsv"
  test -s "${outdir}/plots/ancestry_plot.png"
  test -s "${outdir}/plots/pca_all_samples.png"
}

pruned_marker_count() {
  local outdir="$1"
  if [[ -s "${outdir}/plink/merged_pruned.bim" ]]; then
    wc -l < "${outdir}/plink/merged_pruned.bim"
  else
    echo 0
  fi
}

latest_error_context() {
  local mode="$1"
  local log_file="${LOG_DIR}/nextflow_${mode}.log"
  {
    echo "Mode: $mode"
    echo "Primary log: $log_file"
    echo
    echo "Tail of Nextflow log:"
    tail -n 160 "$log_file" 2>/dev/null || true
    echo
    echo "Recent .nextflow.log tail:"
    tail -n 120 "${ROOT}/.nextflow.log" 2>/dev/null || true
  }
}

fingerprint_error() {
  local mode="$1"
  local fp
  fp="$(
    latest_error_context "$mode" |
      awk '
        /ERROR ~|Caused by:|Process .* terminated|Command error|No such file|ERROR:/ {
          gsub(/[0-9a-f]{2}\/[0-9a-f]{6,}/, "[work]", $0)
          print
        }' |
      head -n 12 |
      sha256sum |
      awk '{print $1}'
  )"
  [[ -n "$fp" ]] || fp="unknown"
  echo "$fp"
}

fingerprint_count() {
  local fp="$1"
  awk -F'\t' -v fp="$fp" '$1 == fp {n++} END {print n+0}' "$FINGERPRINTS"
}

acquire_agy_lock() {
  local waited=0
  while ! mkdir "$AGY_LOCK_DIR" >/dev/null 2>&1; do
    if [[ "$waited" -ge 300 ]]; then
      log "agy lock held for ${waited}s; removing stale lock ${AGY_LOCK_DIR}"
      rm -rf "$AGY_LOCK_DIR"
      continue
    fi
    log "Waiting for existing agy repair agent lock: ${AGY_LOCK_DIR}"
    sleep 30
    waited=$((waited + 30))
  done
}

release_agy_lock() {
  rmdir "$AGY_LOCK_DIR" >/dev/null 2>&1 || true
}

launch_repair_agent() {
  local mode="$1"
  local outdir="$2"
  local fp="$3"
  local attempt="$4"
  local prompt_file="${RUN_DIR}/agy_prompt_${mode}_${attempt}.md"
  local agy_log="${LOG_DIR}/agy_${mode}_${attempt}.log"

  acquire_agy_lock
  latest_error_context "$mode" > "${RUN_DIR}/last_error_context.txt"
  cat > "$prompt_file" <<EOF
You are an autonomous repair agent for the admix_whole Nextflow ancestry pipeline.

Read these first:
- ${ROOT}/AGENTS.md
- ${ROOT}/SKILLS.md

The pipeline failed or became inconsistent.

Mode: ${mode}
Output directory: ${outdir}
Error fingerprint: ${fp}
Attempt: ${attempt} of ${MAX_TOTAL_REPAIRS}

Instructions:
- Inspect the logs and failed work directories.
- Make only targeted fixes needed for this failure.
- Preserve the study input glob: ${STUDY_GLOB}
- Preserve reference paths and default ancestry params documented in AGENTS.md/SKILLS.md.
- Do not delete FASTQs, reference files, old outputs, final VCF symlinks, or run outputs.
- If you edit code, explain what changed.
- Rerun with -resume only if needed to verify the fix.
- Return a concise final status.

Error context:

\`\`\`
$(cat "${RUN_DIR}/last_error_context.txt")
\`\`\`
EOF

  log "Launching agy repair agent attempt=${attempt} mode=${mode} fingerprint=${fp}"
  set +e
  timeout --kill-after=2m "$AGY_WALL_TIMEOUT" agy --dangerously-skip-permissions \
    --print \
    --print-timeout "$AGY_TIMEOUT" \
    --add-dir /datos/migccl/leukemia_exoma \
    --add-dir /datos/migccl/ancestry_refs \
    "$(cat "$prompt_file")" >"$agy_log" 2>&1 &
  AGY_PID=$!
  wait "$AGY_PID"
  local agy_rc=$?
  AGY_PID=""
  set -e
  release_agy_lock
  cleanup_orphaned_repair_agents "repair_attempt_${attempt}_finished"
  if [[ "$agy_rc" -eq 124 ]]; then
    log "agy repair agent timed out after ${AGY_WALL_TIMEOUT}; log=${agy_log}"
  else
    log "agy repair agent exited rc=${agy_rc}; log=${agy_log}"
  fi
}

run_with_repair() {
  local mode="$1"
  local outdir="$2"
  local attempt=0
  while true; do
    if run_nextflow "$mode" "$outdir"; then
      if validate_outputs "$outdir"; then
        log "Nextflow mode=${mode} produced required outputs"
        return 0
      fi
      log "Nextflow mode=${mode} exited 0 but required outputs are missing"
    else
      log "Nextflow mode=${mode} failed"
    fi

    attempt=$((attempt + 1))
    local fp
    fp="$(fingerprint_error "$mode")"
    printf '%s\t%s\t%s\t%s\n' "$fp" "$(date -Is)" "$mode" "$attempt" >> "$FINGERPRINTS"
    local same_count
    same_count="$(fingerprint_count "$fp")"

    if [[ "$same_count" -ge "$MAX_SAME_ERROR" || "$attempt" -ge "$MAX_TOTAL_REPAIRS" ]]; then
      notify_blocked "admix_whole autowatch blocked" \
        "Blocked in mode=${mode}. fingerprint=${fp}; same_count=${same_count}; attempts=${attempt}. See ${RUN_DIR}."
      return 1
    fi

    write_live_status "REPAIRING" "mode=${mode}; attempt=${attempt}; fingerprint=${fp}; same_count=${same_count}; agy_lock=${AGY_LOCK_DIR}" "$mode" "$outdir" "${LOG_DIR}/nextflow_${mode}.log"
    launch_repair_agent "$mode" "$outdir" "$fp" "$attempt"
  done
}

main() {
  if [[ "$STATUS_ONLY" -eq 1 ]]; then
    write_live_status "LIVE_STATUS" "manual status refresh; no pipeline action taken" \
      "primary_hm3" "$PRIMARY_OUT" "${LOG_DIR}/nextflow_primary_hm3.log"
    exit 0
  fi

  write_status "STARTING" "Autowatch starting"
  preflight
  index_study_vcfs

  if [[ "$TEST_MAIL" -eq 1 ]]; then
    if ensure_mail; then
      printf 'admix_whole autowatch test mail at %s\n' "$(date -Is)" | mail -s "admix_whole autowatch test" "$MAIL_TO"
      log "Sent test mail to $MAIL_TO"
    else
      echo "mail unavailable" >&2
      exit 1
    fi
  fi

  if [[ "$PRECHECK_ONLY" -eq 1 ]]; then
    write_status "PREFLIGHT_OK" "Preflight completed without launching Nextflow"
    exit 0
  fi

  if run_with_repair "primary_hm3" "$PRIMARY_OUT"; then
    local markers
    markers="$(pruned_marker_count "$PRIMARY_OUT")"
    if [[ "$markers" -ge "$MIN_PRUNED_MARKERS" ]]; then
      write_status "SUCCESS" "primary_hm3 completed; pruned_markers=${markers}; outdir=${PRIMARY_OUT}"
      exit 0
    fi
    log "Primary completed but marker retention is poor (${markers} < ${MIN_PRUNED_MARKERS}); running fallback"
  else
    exit 1
  fi

  if run_with_repair "fallback_no_hm3" "$FALLBACK_OUT"; then
    local markers
    markers="$(pruned_marker_count "$FALLBACK_OUT")"
    write_status "SUCCESS" "fallback_no_hm3 completed; pruned_markers=${markers}; outdir=${FALLBACK_OUT}"
    exit 0
  fi

  exit 1
}

main "$@"
