nextflow.enable.dsl=2

/*
 * Pipeline: Whole-VCF supervised ADMIXTURE ancestry inference
 *
 * Merges study VCFs with reference population VCFs, runs PLINK QC,
 * then ADMIXTURE in supervised mode (--supervised / .pop file) to fix
 * reference populations and infer ancestry for study samples.
 *
 * Reference input modes:
 *   A) --ref_pops "Pop1:/path/*.vcf.gz;Pop2:/path/*.vcf.gz;..."
 *   B) --ref_vcf_dir "/path/*.vcf.gz" --ref_pop_map pop_map.tsv
 *   C) --ref_vcf_dir "/path/*.vcf.gz" --ref_manifests "Pop1:manifest.tsv;Pop2:manifest.tsv"
 *
 * Tools required: bcftools, plink2, plink 1.9, admixture
 */

// ── Parameters ──────────────────────────────────────────────────────────────
params.input          = "/datos/migccl/leukemia_vcfs/*.vcf.gz"           // Study VCFs glob (bgzipped)
params.outdir         = "results"
params.k              = 3              // Number of ancestral populations

// Reference mode A: per-population VCF globs
params.ref_pops       = "EUR:/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz;AFR:/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz"  // "Pop1:glob;Pop2:glob;..."

// Reference mode B: shared VCF dir + 2-col pop map (sample_id\tpopulation)
params.ref_vcf_dir    = null
params.ref_pop_map    = null

// Reference mode C: shared VCF dir + per-pop manifests (fastq_manifest.tsv format)
params.ref_manifests  = null           // "Pop1:manifest.tsv;Pop2:manifest.tsv"

// QC / PLINK
params.maf            = 0.05
params.geno           = 0.05
params.mind           = 0.8
params.ld_window      = 50
params.ld_step        = 5
params.ld_r2          = 0.2
params.seed           = 42
params.plink2         = "plink2"
params.plink          = "/home/epaaso/bin/plink"
params.admixture      = "admixture"
params.relaxed_geno   = 0.99
params.max_cpus       = 18
params.max_memory     = '120 GB'

// ── Helper: build reference entries ─────────────────────────────────────────
// Returns list of [sample_id, population, vcf_path]
def buildRefEntries() {
    def entries = []

    if (params.ref_pops) {
        // Mode A: "Pop1:glob1;Pop2:glob2"
        params.ref_pops.toString().split(';').each { entry ->
            def parts = entry.trim().split(':', 2)
            if (parts.size() != 2) error "Invalid --ref_pops entry: '${entry}'. Use 'PopName:/path/*.vcf.gz'"
            def pop  = parts[0].trim()
            def glob = parts[1].trim()
            def vcfs = file(glob)
            if (!vcfs) error "No VCFs matched for population '${pop}' with glob '${glob}'"
            vcfs.each { vcf ->
                if (vcf.name.endsWith('.tbi') || vcf.name.endsWith('.csi')) return
                def sid = vcf.name
                    .replaceAll(/\.hard-filtered\.vcf\.gz$/, '')
                    .replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
                entries << [sid, pop, vcf]
            }
        }
    } else if (params.ref_vcf_dir && params.ref_pop_map) {
        // Mode B: ref_vcf_dir + 2-col pop_map TSV
        def pop_map = [:]
        file(params.ref_pop_map).eachLine { line ->
            if (line.trim() && !line.startsWith('#')) {
                def parts = line.trim().split('\t')
                if (parts.size() >= 2) pop_map[parts[0].trim()] = parts[1].trim()
            }
        }
        file(params.ref_vcf_dir).each { vcf ->
            if (vcf.name.endsWith('.tbi') || vcf.name.endsWith('.csi')) return
            def sid = vcf.name
                .replaceAll(/\.hard-filtered\.vcf\.gz$/, '')
                .replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
            if (pop_map.containsKey(sid)) {
                entries << [sid, pop_map[sid], vcf]
            }
        }
    } else if (params.ref_vcf_dir && params.ref_manifests) {
        // Mode C: ref_vcf_dir + per-pop manifests (fastq_manifest.tsv format)
        def pop_map = [:]
        params.ref_manifests.toString().split(';').each { entry ->
            def parts = entry.trim().split(':', 2)
            if (parts.size() != 2) error "Invalid --ref_manifests entry: '${entry}'. Use 'PopName:/path/manifest.tsv'"
            def pop      = parts[0].trim()
            def manifest = file(parts[1].trim())
            manifest.eachLine { line ->
                if (line.trim()) {
                    def cols = line.trim().split('\t')
                    // Column 1 may be semicolon-separated sample IDs
                    cols[0].split(';').each { sid ->
                        if (sid.trim()) pop_map[sid.trim()] = pop
                    }
                }
            }
        }
        file(params.ref_vcf_dir).each { vcf ->
            if (vcf.name.endsWith('.tbi') || vcf.name.endsWith('.csi')) return
            def sid = vcf.name
                .replaceAll(/\.hard-filtered\.vcf\.gz$/, '')
                .replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
            if (pop_map.containsKey(sid)) {
                entries << [sid, pop_map[sid], vcf]
            }
        }
    } else {
        error """Provide reference populations via one of:
  --ref_pops "Pop1:glob;Pop2:glob;..."
  --ref_vcf_dir <glob> --ref_pop_map <tsv>
  --ref_vcf_dir <glob> --ref_manifests "Pop1:manifest;Pop2:manifest"
"""
    }

    if (!entries) error "No reference VCFs found matching population configuration"

    // Deduplicate by sample_id (keep first)
    def seen = [] as Set
    entries = entries.findAll { e ->
        if (seen.contains(e[0])) return false
        seen << e[0]
        return true
    }

    def pops = entries.collect { it[1] }.unique().sort()
    log.info "Reference: ${entries.size()} samples across ${pops.size()} populations: ${pops.join(', ')}"
    return entries
}

// ── Workflow ─────────────────────────────────────────────────────────────────
workflow {
    if (!params.input) error "Provide --input with study VCF glob"

    // 1. Build reference pop map and collect VCFs
    def ref_entries  = buildRefEntries()
    def populations  = ref_entries.collect { it[1] }.unique().sort()
    def actual_k     = populations.size()

    if (actual_k != params.k as int) {
        log.warn "Detected ${actual_k} reference populations (${populations}), but --k=${params.k}. Using K=${actual_k}."
    }

    // Write pop_map.tsv (sample_id\tpopulation) for downstream processes
    def pop_map_file = file("${workDir}/pop_map.tsv")
    pop_map_file.text = ref_entries.collect { "${it[0]}\t${it[1]}" }.join('\n') + '\n'
    def pop_map_ch = Channel.value(pop_map_file)

    // 2. Channels
    Channel.fromPath(params.input, checkIfExists: true)
        .set { study_vcfs }

    Channel.from(ref_entries.collect { it[2] }.unique())
        .set { ref_vcfs }

    // 3. Index all VCFs
    study_vcfs
        .mix(ref_vcfs)
        .map { vcf ->
            def id = vcf.name
                .replaceAll(/\.hard-filtered\.vcf\.gz$/, '')
                .replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
            tuple(id, vcf)
        }
        .set { all_vcfs }

    INDEX_VCF(all_vcfs)

    // 3b. Collect study sample sites (for intersection with ref)
    Channel.fromPath(params.input, checkIfExists: true)
        .map { it.toString() }
        .collectFile(name: 'study_vcfs.list', newLine: true)
        .set { study_vcf_list_for_sites }

    COLLECT_STUDY_SITES(study_vcf_list_for_sites)

    // 4. Merge all VCFs (restricted to study sites = intersection)
    INDEX_VCF.out
        .flatMap { id, vcf, tbi -> [vcf, tbi] }
        .collect()
        .set { all_indexed_files }

    MERGE_ALL(all_indexed_files, COLLECT_STUDY_SITES.out)

    // 5. PLINK QC
    PREPARE_PLINK(MERGE_ALL.out, pop_map_ch, Channel.value(actual_k))

    // 6. Create .pop file for supervised ADMIXTURE
    CREATE_POP_FILE(PREPARE_PLINK.out.fam, pop_map_ch)

    // 7. Run supervised ADMIXTURE
    RUN_ADMIXTURE_SUPERVISED(
        PREPARE_PLINK.out.pruned,
        CREATE_POP_FILE.out,
        Channel.value(actual_k)
    )

    // 8. Summarize Q
    SUMMARIZE_Q(RUN_ADMIXTURE_SUPERVISED.out, pop_map_ch)

    // 9. Plot
    PLOT_ANCESTRY(SUMMARIZE_Q.out, pop_map_ch)
}

// ── Processes ────────────────────────────────────────────────────────────────

process INDEX_VCF {
    tag { sample_id }

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path(vcf), path("${vcf}.tbi")

    script:
    """
    set -euo pipefail
    if [ ! -f "${vcf}.tbi" ]; then
        bcftools index -t -f ${vcf}
    fi
    """
}

process COLLECT_STUDY_SITES {
    input:
    path study_vcf_list

    output:
    path "study_sites.tsv"

    script:
    """
    set -euo pipefail

    n_vcfs=\$(wc -l < "${study_vcf_list}")
    echo "Collecting study sites from \$n_vcfs study VCFs" >&2

    while IFS= read -r vcf; do
        [ -n "\$vcf" ] || continue
        bcftools query -f '%CHROM\\t%POS\\n' "\$vcf"
    done < "${study_vcf_list}" | LC_ALL=C sort -u > study_sites.tsv

    n=\$(wc -l < study_sites.tsv)
    echo "Collected \$n unique study sites" >&2
    """
}

process MERGE_ALL {
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path all_files
    path study_sites

    output:
    tuple path("merged.vcf.gz"), path("merged.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    # List only VCF files for merge
    # Ignore possible leftovers from prior failed attempts in the same work dir.
    (ls *.vcf.gz | grep -Ev '^merged(_raw)?\\.vcf\\.gz\$' || true) | sort -V > all_files.txt
    n_all=\$(wc -l < all_files.txt)
    echo "Found \$n_all VCF files for merge" >&2

    # Drop fully redundant files with duplicated sample IDs.
    # If a file has a partial overlap (some duplicated + some new samples), fail fast.
    : > file_list.txt
    declare -A seen_samples=()
    n_kept=0
    n_skipped=0

    while IFS= read -r vcf; do
        [ -n "\$vcf" ] || continue

        mapfile -t samples < <(bcftools query -l "\$vcf")
        if [ "\${#samples[@]}" -eq 0 ]; then
            echo "ERROR: No samples found in \$vcf" >&2
            exit 1
        fi

        dup_samples=()
        new_count=0
        for s in "\${samples[@]}"; do
            if [[ -v "seen_samples[\$s]" ]]; then
                dup_samples+=("\$s")
            else
                new_count=\$((new_count + 1))
            fi
        done

        if [ "\${#dup_samples[@]}" -gt 0 ] && [ "\$new_count" -gt 0 ]; then
            echo "ERROR: Input VCF '\$vcf' mixes duplicate and new samples." >&2
            echo "Duplicate samples: \${dup_samples[*]}" >&2
            echo "Please remove overlapping multi-sample VCFs from the merge input." >&2
            exit 1
        fi

        if [ "\$new_count" -eq 0 ]; then
            n_skipped=\$((n_skipped + 1))
            echo "Skipping redundant VCF \$vcf (all samples already present)" >&2
            continue
        fi

        printf '%s\n' "\$vcf" >> file_list.txt
        n_kept=\$((n_kept + 1))
        for s in "\${samples[@]}"; do
            seen_samples["\$s"]=1
        done
    done < all_files.txt

    n_files=\$(wc -l < file_list.txt)
    echo "Merging \$n_files VCF files (skipped \$n_skipped redundant files)..." >&2

    if [ "\$n_files" -eq 0 ]; then
        echo "ERROR: No non-redundant VCF files available after sample deduplication." >&2
        exit 1
    fi

    # Harmonize chromosome naming across cohorts (e.g., "1" vs "chr1").
    # Use study sites naming style as the target convention.
    : > merge_list.txt
    first_site_chr=\$(awk 'NR==1 {print \$1}' ${study_sites})
    if [ -z "\$first_site_chr" ]; then
        echo "ERROR: Study site list is empty." >&2
        exit 1
    fi

    target_style="nochr"
    if [[ "\$first_site_chr" == chr* ]]; then
        target_style="chr"
        cat > chr_map.tsv <<'EOF'
1	chr1
2	chr2
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
10	chr10
11	chr11
12	chr12
13	chr13
14	chr14
15	chr15
16	chr16
17	chr17
18	chr18
19	chr19
20	chr20
21	chr21
22	chr22
X	chrX
Y	chrY
MT	chrM
M	chrM
EOF
    else
        cat > chr_map.tsv <<'EOF'
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr10	10
chr11	11
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr19	19
chr20	20
chr21	21
chr22	22
chrX	X
chrY	Y
chrM	MT
chrMT	MT
EOF
    fi

    n_renamed=0
    while IFS= read -r vcf; do
        [ -n "\$vcf" ] || continue
        first_chr=\$(bcftools query -f '%CHROM\n' "\$vcf" | head -n 1 || true)
        if [ -z "\$first_chr" ]; then
            printf '%s\n' "\$vcf" >> merge_list.txt
            continue
        fi

        needs_rename=0
        if [ "\$target_style" = "chr" ] && [[ ! "\$first_chr" == chr* ]]; then
            needs_rename=1
        fi
        if [ "\$target_style" = "nochr" ] && [[ "\$first_chr" == chr* ]]; then
            needs_rename=1
        fi

        if [ "\$needs_rename" -eq 1 ]; then
            harmonized_vcf="\${vcf%.vcf.gz}.harmonized.vcf.gz"
            bcftools annotate --threads ${task.cpus} --rename-chrs chr_map.tsv -Oz -o "\$harmonized_vcf" "\$vcf"
            bcftools index -t "\$harmonized_vcf"
            printf '%s\n' "\$harmonized_vcf" >> merge_list.txt
            n_renamed=\$((n_renamed + 1))
        else
            printf '%s\n' "\$vcf" >> merge_list.txt
        fi
    done < file_list.txt
    echo "Contig harmonization: renamed \$n_renamed VCFs to '\$target_style' style" >&2

    if [ "\$n_files" -eq 1 ]; then
        bcftools view -T ${study_sites} \$(cat merge_list.txt) -Oz -o merged.vcf.gz
    else
        # Merge only study target sites to avoid off-target REF conflicts across cohorts
        bcftools merge --threads ${task.cpus} -R ${study_sites} -l merge_list.txt -Oz -o merged_raw.vcf.gz
        bcftools index -t merged_raw.vcf.gz

        bcftools view --threads ${task.cpus} -T ${study_sites} merged_raw.vcf.gz -Oz -o merged.vcf.gz
        rm -f merged_raw.vcf.gz merged_raw.vcf.gz.tbi
    fi

    # Cleanup temporary harmonized VCFs
    grep -F '.harmonized.vcf.gz' merge_list.txt | while IFS= read -r f; do
        rm -f "\$f" "\$f.tbi"
    done

    bcftools index -t merged.vcf.gz

    n_sites=\$(bcftools query -f '%CHROM\\n' merged.vcf.gz | wc -l)
    echo "Merged VCF: \$n_sites sites" >&2
    """
}

process PREPARE_PLINK {
    publishDir "${params.outdir}/plink", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path pop_map
    val expected_k

    output:
    tuple path("merged_pruned.bed"),
          path("merged_pruned.bim"),
          path("merged_pruned.fam"), emit: pruned
    path "merged_pruned.fam",        emit: fam

    script:
    """
    set -euo pipefail
    plink=\${PLINK2:-${params.plink2}}
    mem_mb=${task.memory.toMega()}

    # VCF -> PLINK2 pgen (autosomes only, sort variants)
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --vcf ${vcf} --max-alleles 2 --autosome --make-pgen --sort-vars --out sorted

    # pgen -> BED
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --pfile sorted --make-bed --out base

    # Keep biallelic SNPs only
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile base --snps-only just-acgt --max-alleles 2 --make-bed --out snps

    # QC (strict first): MAF, missingness, mind
    strict_ok=1
    if ! \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile snps --maf ${params.maf} --geno ${params.geno} --mind ${params.mind} --make-bed --out clean_strict
    then
        strict_ok=0
        echo "Strict QC failed; switching to relaxed QC." >&2
    fi

    # Ensure strict QC retains reference samples from all expected populations.
    # If not, rerun with relaxed missingness filter to preserve supervision anchors.
    if [ "\$strict_ok" -eq 1 ]; then
        if python3 - <<'PY'
import sys

pop_map = {}
with open("${pop_map}", "r") as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            parts = line.split("\\t")
            if len(parts) >= 2:
                pop_map[parts[0]] = parts[1]

expected = sorted(set(pop_map.values()))
present = set()
with open("clean_strict.fam", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2 and parts[1] in pop_map:
            present.add(pop_map[parts[1]])

print(f"Strict QC reference populations retained: {sorted(present)}", file=sys.stderr)
if len(present) < int("${expected_k}"):
    sys.exit(1)
PY
        then
            mv clean_strict.bed clean.bed
            mv clean_strict.bim clean.bim
            mv clean_strict.fam clean.fam
        else
            strict_ok=0
            echo "Strict QC dropped reference populations; switching to relaxed QC." >&2
        fi
    fi

    if [ "\$strict_ok" -eq 0 ]; then
        rm -f clean_strict.bed clean_strict.bim clean_strict.fam
        \$plink --threads ${task.cpus} --memory \$mem_mb \\
            --bfile snps --maf ${params.maf} --geno ${params.relaxed_geno} --make-bed --out clean
    fi

    # Deduplicate
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile clean --set-all-var-ids @:# --rm-dup exclude-all --make-bed --out dedup

    # LD pruning
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile dedup --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} --bad-ld --out prune

    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile dedup --extract prune.prune.in --make-bed --out merged_pruned
    """
}

process CREATE_POP_FILE {
    publishDir "${params.outdir}/admixture", mode: 'copy'

    input:
    path fam
    path pop_map

    output:
    path "merged_pruned.pop"

    script:
    """
    python3 - <<'PY'
import sys

# Read pop map (sample_id -> population)
pop_map = {}
with open("${pop_map}", 'r') as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split('\\t')
            if len(parts) >= 2:
                pop_map[parts[0]] = parts[1]

# Read FAM and write .pop file
n_ref = 0
n_study = 0
with open("${fam}", 'r') as fin, open("merged_pruned.pop", 'w') as fout:
    for line in fin:
        parts = line.strip().split()
        if len(parts) >= 2:
            iid = parts[1]
            if iid in pop_map:
                fout.write(pop_map[iid] + '\\n')
                n_ref += 1
            else:
                fout.write('-\\n')
                n_study += 1

print(f"Pop file created: {n_ref} reference, {n_study} study samples", file=sys.stderr)
if n_ref == 0:
    sys.exit("ERROR: No reference samples matched FAM IIDs. Check pop_map sample IDs.")
if n_study == 0:
    print("WARNING: All individuals are reference — no study samples to infer.", file=sys.stderr)
PY
    """
}

process RUN_ADMIXTURE_SUPERVISED {
    publishDir "${params.outdir}/admixture", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)
    path pop_file
    val k

    output:
    tuple path("${bed.baseName}.${k}.Q"),
          path("${bed.baseName}.${k}.P"),
          path(fam),
          path(pop_file)

    script:
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${task.cpus}

    # ADMIXTURE requires .pop to share the .bed base name
    if [ "${pop_file}" != "${bed.baseName}.pop" ]; then
        cp ${pop_file} ${bed.baseName}.pop
    fi

    ${params.admixture} -j${task.cpus} --seed=${params.seed} --cv --supervised \\
        ${bed} ${k} | tee admixture.log
    """
}

process SUMMARIZE_Q {
    publishDir "${params.outdir}/ancestry", mode: 'copy'

    input:
    tuple path(qfile), path(pfile), path(fam), path(pop_file)
    path pop_map

    output:
    path "ancestry.tsv"

    script:
    """
    python3 - <<'PY'
import csv, sys, itertools

qfile       = "${qfile}"
famfile     = "${fam}"
pop_map_f   = "${pop_map}"
outfile     = "ancestry.tsv"

# Read pop map
pop_map = {}
with open(pop_map_f, 'r') as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split('\\t')
            if len(parts) >= 2:
                pop_map[parts[0]] = parts[1]

populations = sorted(set(pop_map.values()))

# Read Q
q_rows = []
with open(qfile, 'r') as f:
    for line in f:
        line = line.strip()
        if not line: continue
        q_rows.append([float(x) for x in line.split()])

# Read FAM IIDs
iids = []
with open(famfile, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            iids.append(parts[1])

if len(q_rows) != len(iids):
    sys.exit(f"Q rows ({len(q_rows)}) != FAM IIDs ({len(iids)})")

num_k = len(q_rows[0]) if q_rows else 0
if num_k != len(populations):
    sys.exit(f"K={num_k} columns in Q but {len(populations)} populations detected")

# Map Q columns -> populations via reference sample means
col_means  = {p: [0.0] * num_k for p in populations}
col_counts = {p: 0 for p in populations}
for idx, iid in enumerate(iids):
    if iid in pop_map:
        pop = pop_map[iid]
        for k in range(num_k):
            col_means[pop][k] += q_rows[idx][k]
        col_counts[pop] += 1

for p in populations:
    if col_counts[p] > 0:
        col_means[p] = [x / col_counts[p] for x in col_means[p]]
    else:
        sys.exit(f"No reference samples from '{p}' found in FAM")

# Best 1:1 assignment (max score)
best_score = None
best_assign = None
for perm in itertools.permutations(range(num_k), len(populations)):
    score = sum(col_means[p][ci] for p, ci in zip(populations, perm))
    if best_score is None or score > best_score:
        best_score = score
        best_assign = {p: ci for p, ci in zip(populations, perm)}

print(f"Column assignment: {best_assign}", file=sys.stderr)

# Write
with open(outfile, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\\t')
    w.writerow(["IID", "is_reference"] + populations)
    for iid, q in zip(iids, q_rows):
        is_ref = "ref" if iid in pop_map else "study"
        vals   = [q[best_assign[p]] for p in populations]
        w.writerow([iid, is_ref] + [f"{v:.6f}" for v in vals])

print(f"Wrote ancestry for {len(iids)} individuals", file=sys.stderr)
PY
    """
}

process PLOT_ANCESTRY {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path ancestry_tsv
    path pop_map

    output:
    path "ancestry_plot.png"
    path "ancestry_study_only.tsv"

    script:
    """
    python3 - <<'PY'
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

df = pd.read_csv("${ancestry_tsv}", sep='\\t')
study_df = df[df['is_reference'] == 'study'].copy()
study_df.to_csv("ancestry_study_only.tsv", sep='\\t', index=False)

if study_df.empty:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, 'No study samples', ha='center', va='center')
    plt.savefig("ancestry_plot.png")
else:
    pop_cols = [c for c in study_df.columns if c not in ('IID', 'is_reference')]
    plot_df  = study_df.set_index('IID')[pop_cols].astype(float)
    plot_df  = plot_df.sort_values(by=pop_cols, ascending=False)

    width = max(10, min(200, len(plot_df) * 0.3))
    ax = plot_df.plot(kind='bar', stacked=True, figsize=(width, 6), width=1.0)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title("Supervised ADMIXTURE Ancestry (K=${params.k})")
    plt.xlabel("Sample")
    plt.ylabel("Proportion")
    plt.tight_layout()
    plt.savefig("ancestry_plot.png", dpi=300)
PY
    """
}
