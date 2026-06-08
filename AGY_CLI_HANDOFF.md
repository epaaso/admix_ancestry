# AGY CLI Handoff: leukemia exome ancestry run

This file captures the current state of the autonomous `admix_whole` run so a future
`agy` CLI session can continue without rediscovering the failure history.

## Stop State

The current autowatch/Nextflow run was intentionally stopped by the user on:

```text
2026-06-03 21:07 UTC
```

At stop time the active task was:

```text
MERGE_ALL (1)
```

Active work directory:

```text
/datos/home/epaaso/ancestry/admix_whole/work/87/e667176c986e61df6377227cde0e47
```

Partial file at stop time:

```text
/datos/home/epaaso/ancestry/admix_whole/work/87/e667176c986e61df6377227cde0e47/merged_raw.vcf.gz
253231104 bytes
```

Do not treat this partial `merged_raw.vcf.gz` as resumable. `bcftools merge` writes a
single bgzip stream and a killed partial output is not a safe restart point. Nextflow
`-resume` can reuse upstream cached tasks, but `MERGE_ALL` itself must rerun.

No watcher, Nextflow, `bcftools`, `plink`, `admixture`, or `agy` repair process should
be running after this stop.

Check with:

```bash
ps -eo pid,ppid,stat,etime,%cpu,%mem,cmd --sort=start_time |
  rg 'run_admix_autowatch|nextflow run main.nf|java.*nextflow|bcftools merge|bcftools annotate|plink|admixture|agy --dangerously' |
  rg -v rg || true
```

## Fixed Analysis Inputs

Study VCFs:

```text
/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz
```

Expected count: `50`.

Primary output:

```text
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole
```

Fallback no-HM3 output:

```text
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_no_hm3
```

Autowatch run/log directory:

```text
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch
```

Important logs:

```text
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch/logs/autowatch.log
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch/logs/nextflow_primary_hm3.log
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch/logs/agy_primary_hm3_1.log
```

Primary run parameters:

```text
--k 3
--extract_snps /datos/home/epaaso/ancestry/hm3_combined_ids_and_loci.txt
--exclude_regions /datos/home/epaaso/ancestry/leukemia_genes_hg38.bed
--maf 0.05
--geno 0.05
--mind 0.8
--relaxed_geno 0.99
--ld_window 50
--ld_step 5
--ld_r2 0.2
--max_cpus 18
--max_memory "120 GB"
```

Reference populations:

```text
EUR:/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz;AFR:/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz
```

## Failure History

The run repeatedly reached `MERGE_ALL`, then failed or became impractically slow.

Observed failures:

1. Initial primary HM3 run failed on `2026-06-02 14:36 UTC`.
2. Later run failed on `2026-06-03 00:46 UTC`.
3. Current restarted run failed again on `2026-06-03 12:54 UTC`, then `agy` repaired it
   and restarted at `2026-06-03 12:57 UTC`.

The important direct error from the threaded `bcftools merge` attempts was:

```text
bcftools: thread_pool.c:676: wake_next_worker: Assertion `p->njobs >= q->n_input' failed.
Aborted (core dumped)
```

Interpretation:

- This is not a normal resource exhaustion error.
- It is an internal `bcftools`/`htslib` thread-pool assertion triggered by threaded
  merge behavior.
- Running other user processes may slow the job, but does not explain this assertion.
- Do not simply restore `bcftools merge --threads 14` unless also changing the merge
  strategy or validating with a small controlled test.

## Current Code State

In `main.nf`, `MERGE_ALL` currently:

- caps `bcftools_threads` at 14 for annotate/view steps,
- runs the central `bcftools merge` without `--threads`,
- indexes `merged_raw.vcf.gz`,
- then filters it to `merged.vcf.gz`.

This single-threaded merge was intentionally introduced by an `agy` repair agent to
avoid the `wake_next_worker` assertion, but it is very slow. At around `2026-06-03
20:57 UTC`, the partial merge was growing at roughly 30-35 MB/hour.

## Recommended Continuation

Best next engineering step:

Implement a chromosome-sharded `MERGE_ALL` strategy:

1. Use `study_sites.tsv` to define per-chromosome target site files.
2. Run one `bcftools merge` per chromosome or contig, each single-threaded or with a
   very small thread count.
3. Run chromosome merges in parallel through Nextflow task parallelism or a bounded
   shell background worker pool.
4. Index each per-chromosome merged VCF.
5. Concatenate with `bcftools concat -Oz -o merged.vcf.gz`.
6. Index `merged.vcf.gz`.

This uses available cores without relying on the failure-prone single giant
`bcftools merge --threads N` path.

Conservative fallback:

Let the current single-threaded merge run to completion, but expect it may take many
hours to days depending on final output size.

Do not delete:

- study VCF symlinks,
- VCF indexes,
- FASTQs,
- reference VCFs,
- previous output directories,
- Nextflow `work/` cache.

## Commands For Agy CLI

Start from:

```bash
cd /datos/home/epaaso/ancestry/admix_whole
```

Suggested `agy` CLI invocation:

```bash
agy --dangerously-skip-permissions \
  --add-dir /datos/migccl/leukemia_exoma \
  --add-dir /datos/migccl/ancestry_refs \
  "Read AGENTS.md, SKILLS.md, and AGY_CLI_HANDOFF.md. Continue the leukemia exome ancestry analysis. The previous autowatch run was stopped intentionally on 2026-06-03 after MERGE_ALL became too slow. Do not delete input VCFs, references, old outputs, or the Nextflow work cache. Inspect the logs under /datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch/logs. Implement a targeted fix for MERGE_ALL, preferably chromosome-sharded merge plus bcftools concat, avoiding a single giant bcftools merge --threads run because it previously crashed with thread_pool.c wake_next_worker. Then rerun with nextflow -resume using the primary HM3 parameters documented in SKILLS.md. Validate that ancestry/ancestry.tsv, plots/ancestry_plot.png, and plots/pca_all_samples.png exist."
```

If using non-interactive print mode:

```bash
agy --dangerously-skip-permissions --print --print-timeout 15m \
  --add-dir /datos/migccl/leukemia_exoma \
  --add-dir /datos/migccl/ancestry_refs \
  "Read AGENTS.md, SKILLS.md, and AGY_CLI_HANDOFF.md. Continue the leukemia exome ancestry analysis as described in AGY_CLI_HANDOFF.md."
```

## Relaunch Command After Fix

Primary HM3 run:

```bash
nextflow run main.nf -resume \
  --input "/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz" \
  --outdir /datos/migccl/leukemia_exoma/outs-results/005-admix_whole \
  --ref_pops "EUR:/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz;AFR:/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz" \
  --k 3 \
  --extract_snps /datos/home/epaaso/ancestry/hm3_combined_ids_and_loci.txt \
  --exclude_regions /datos/home/epaaso/ancestry/leukemia_genes_hg38.bed \
  --maf 0.05 --geno 0.05 --mind 0.8 --relaxed_geno 0.99 \
  --ld_window 50 --ld_step 5 --ld_r2 0.2 \
  --max_cpus 18 --max_memory "120 GB"
```

Use tmux for long runs:

```bash
tmux new-window -t 0 -n admix-agy
tmux send-keys -t 0:admix-agy 'cd /datos/home/epaaso/ancestry/admix_whole' C-m
```

## Validation

```bash
outdir=/datos/migccl/leukemia_exoma/outs-results/005-admix_whole
test -s "$outdir/ancestry/ancestry.tsv"
test -s "$outdir/plots/ancestry_plot.png"
test -s "$outdir/plots/pca_all_samples.png"
wc -l "$outdir/ancestry/ancestry.tsv"
```

If primary HM3 completes but marker retention is too low, run the no-HM3 fallback
documented in `SKILLS.md`.
