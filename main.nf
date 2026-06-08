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
params.extract_snps   = null           // File containing list of SNPs to keep (e.g., AIMs)
params.exclude_snps   = null           // File containing list of SNPs to exclude
params.exclude_regions = null          // BED file containing regions to exclude (e.g., leukemia genes)
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

// Liftover (GRCh37 -> GRCh38)
params.grch38_ref     = "/datos/migccl/arriba_refs/GRCh38.primary_assembly.genome.fa"
params.chain_file     = null           // Path to hg19ToHg38.over.chain.gz (auto-downloaded if null)
params.picard_version = "3.2.0"        // Used only when --picard_jar is not provided
params.picard_jar     = null           // Optional explicit path to picard.jar

// Infer sample ID from VCF filename, removing common technical suffixes.
def inferSampleIdFromVcfName(String filename) {
    def sid = filename
        .replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
        .replaceAll(/\.hard-filtered$/, '')
        .replaceAll(/\.annotated\.nh$/, '')
        .replaceAll(/\.annotated$/, '')
        .replaceAll(/\.nh$/, '')
    return sid
}

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
                def sid = inferSampleIdFromVcfName(vcf.name)
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
            def sid = inferSampleIdFromVcfName(vcf.name)
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
            def sid = inferSampleIdFromVcfName(vcf.name)
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

    // 2b. Prepare chain file for liftover (download if not provided)
    def chain_ch
    if (params.chain_file) {
        chain_ch = Channel.value(file(params.chain_file))
    } else {
        chain_ch = DOWNLOAD_CHAIN_FILE().chain
    }

    // 2c. Prepare Picard jar for liftover (download if not provided)
    def picard_jar_ch
    if (params.picard_jar) {
        picard_jar_ch = Channel.value(file(params.picard_jar))
    } else {
        picard_jar_ch = DOWNLOAD_PICARD_JAR().picard
    }

    def grch38_ref_ch  = Channel.value(file(params.grch38_ref))
    def grch38_fai_ch  = Channel.value(file("${params.grch38_ref}.fai"))
    def grch38_dict_path = params.grch38_ref.replaceAll(/\.fa(sta)?$/, '.dict')
    def grch38_dict_ch = Channel.value(file(grch38_dict_path))

    // 3. Detect genome build and liftover GRCh37 VCFs
    //    Study VCFs are assumed GRCh38; only ref VCFs are checked.
    ref_vcfs
        .map { vcf ->
            def id = inferSampleIdFromVcfName(vcf.name)
            tuple(id, vcf)
        }
        .set { ref_vcfs_tuples }

    DETECT_AND_LIFTOVER(ref_vcfs_tuples, chain_ch, grch38_ref_ch, grch38_fai_ch, grch38_dict_ch, picard_jar_ch)

    // 4. Index all VCFs (study + harmonised refs)
    study_vcfs
        .map { vcf ->
            def id = inferSampleIdFromVcfName(vcf.name)
            tuple(id, vcf)
        }
        .mix(DETECT_AND_LIFTOVER.out.vcf)
        .set { all_vcfs }

    INDEX_VCF(all_vcfs)

    // 4b. Collect study sample sites (for intersection with ref)
    Channel.fromPath(params.input, checkIfExists: true)
        .map { it.toString() }
        .collectFile(name: 'study_vcfs.list', newLine: true)
        .set { study_vcf_list_for_sites }

    COLLECT_STUDY_SITES(study_vcf_list_for_sites)

    // 5. Merge all VCFs (restricted to study sites = intersection)
    HARMONIZE_CONTIGS(INDEX_VCF.out, COLLECT_STUDY_SITES.out.first())

    HARMONIZE_CONTIGS.out.harmonized
        .flatMap { id, vcf, tbi -> [vcf, tbi] }
        .collect()
        .set { all_harmonized_files }

    SPLIT_STUDY_SITES(COLLECT_STUDY_SITES.out)

    SPLIT_STUDY_SITES.out.sites_files
        .flatten()
        .map { file ->
            def chr = file.name.replaceAll(/^sites_/, "").replaceAll(/\.tsv$/, "")
            tuple(chr, file)
        }
        .set { split_sites_ch }

    MERGE_CHROMOSOME(split_sites_ch, all_harmonized_files)

    MERGE_CHROMOSOME.out.merged_chr
        .flatMap { chr, vcf, tbi -> [vcf, tbi] }
        .collect()
        .set { all_merged_chrs }

    CONCAT_VCFS(all_merged_chrs, SPLIT_STUDY_SITES.out.chr_list, COLLECT_STUDY_SITES.out)

    // 6. PLINK QC
    PREPARE_PLINK(CONCAT_VCFS.out, pop_map_ch, Channel.value(actual_k))

    // 7. Create .pop file for supervised ADMIXTURE
    CREATE_POP_FILE(PREPARE_PLINK.out.fam, pop_map_ch)

    // 8. Run supervised ADMIXTURE
    RUN_ADMIXTURE_SUPERVISED(
        PREPARE_PLINK.out.pruned,
        CREATE_POP_FILE.out,
        Channel.value(actual_k)
    )

    // 9. Summarize Q
    SUMMARIZE_Q(RUN_ADMIXTURE_SUPERVISED.out, pop_map_ch)

    // 10. Plot
    PLOT_ANCESTRY(SUMMARIZE_Q.out, pop_map_ch)

    // 11. PCA (all samples + references)
    RUN_PCA(PREPARE_PLINK.out.pruned)
    PLOT_PCA(RUN_PCA.out, pop_map_ch)

    // 12. Plot P matrix and annotate top SNPs
    PLOT_P_MATRIX(RUN_ADMIXTURE_SUPERVISED.out, PREPARE_PLINK.out.pruned)
}

// ── Processes ────────────────────────────────────────────────────────────────

process DOWNLOAD_CHAIN_FILE {
    storeDir "/datos/migccl/arriba_refs"

    output:
    path "hg19ToHg38.over.chain.gz", emit: chain

    script:
    """
    set -euo pipefail
    echo "Downloading hg19-to-hg38 liftover chain file from UCSC..." >&2
    wget -q -O hg19ToHg38.over.chain.gz \
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
    echo "Chain file downloaded successfully." >&2
    """
}

process DOWNLOAD_PICARD_JAR {
    storeDir "/datos/migccl/arriba_refs"

    output:
    path "picard-${params.picard_version}.jar", emit: picard

    script:
    """
    set -euo pipefail
    out="picard-${params.picard_version}.jar"
    url="https://github.com/broadinstitute/picard/releases/download/${params.picard_version}/picard.jar"

    echo "Downloading Picard ${params.picard_version} from GitHub releases..." >&2
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "\$url" -o "\$out"
    else
        wget -q -O "\$out" "\$url"
    fi

    # Quick sanity check: ensure the jar is executable and exposes LiftoverVcf.
    java -jar "\$out" LiftoverVcf --help >/dev/null 2>&1
    echo "Picard downloaded successfully: \$out" >&2
    """
}

/*
 * DETECT_AND_LIFTOVER
 *
 * For each reference VCF, detects the genome build by inspecting contig
 * lengths.  GRCh37 VCFs (chr1 ≈ 249 250 621 bp) are lifted over to GRCh38
 * using Picard LiftoverVcf.  GRCh38 VCFs pass through unchanged.
 *
 * Additionally, all-sites VCFs (ALT=".") are filtered to variant-only so
 * they don't inject false homozygous-reference genotypes into the merge.
 */
process DETECT_AND_LIFTOVER {
    tag { sample_id }
    publishDir "${params.outdir}/liftover", mode: 'copy', pattern: '*.liftover_report.txt'

    input:
    tuple val(sample_id), path(vcf)
    path chain
    path grch38_ref
    path grch38_fai
    path grch38_dict
    path picard_jar

    output:
    tuple val(sample_id), path("${sample_id}.harmonised.vcf.gz"), emit: vcf
    path "${sample_id}.liftover_report.txt",                      emit: report

    script:
    """
    set -euo pipefail

    REPORT="${sample_id}.liftover_report.txt"
    echo "=== Build detection report for ${sample_id} ===" > "\$REPORT"
    echo "Input VCF: ${vcf}" >> "\$REPORT"

    # ── 1. Detect genome build by chr1 contig length ──────────────────────
    # GRCh37/hg19 chr1 = 249250621 bp;  GRCh38/hg38 chr1 = 248956422 bp
    chr1_len=\$(bcftools view -h ${vcf} \
        | (grep -E '^##contig=<ID=(chr)?1,' || true) \
        | sed -E 's/.*length=([0-9]+).*/\\1/' \
        | head -1)

    if [ -z "\$chr1_len" ]; then
        echo "WARNING: No contig header found; assuming GRCh38 (pass-through)." >> "\$REPORT"
        chr1_len=248956422
    fi

    echo "Detected chr1 length: \$chr1_len" >> "\$REPORT"

    IS_GRCH37=0
    if [ "\$chr1_len" -eq 249250621 ]; then
        IS_GRCH37=1
        echo "" >> "\$REPORT"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >> "\$REPORT"
        echo "  BUILD: GRCh37 (hs37d5/hg19) detected" >> "\$REPORT"
        echo "  ACTION: Liftover to GRCh38 will be performed" >> "\$REPORT"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >> "\$REPORT"
        echo "" >> "\$REPORT"
    else
        echo "BUILD: GRCh38 detected  >>>  no liftover needed" >> "\$REPORT"
    fi

    # ── 2. Filter all-sites records (ALT=".") if present ──────────────────
    bcftools view -H ${vcf} | head -10000 > _sample_10k.tmp || true
    n_total=\$(wc -l < _sample_10k.tmp)
    n_refonly=\$(awk '\$5 == "."' _sample_10k.tmp | wc -l)
    rm -f _sample_10k.tmp
    allsites_pct=\$(( n_refonly * 100 / (n_total > 0 ? n_total : 1) ))

    if [ "\$allsites_pct" -gt 50 ]; then
        echo "ALL-SITES VCF detected (\${allsites_pct}% ref-only in first 10k records). Filtering to variant-only." >> "\$REPORT"
        FILTER_ALLSITES=1
    else
        FILTER_ALLSITES=0
    fi

    # ── 3. Prepare input VCF (filter all-sites if needed) ─────────────────
    if [ "\$FILTER_ALLSITES" -eq 1 ]; then
        bcftools view -i 'ALT!="."' ${vcf} -Oz -o filtered.vcf.gz
        bcftools index -t filtered.vcf.gz
        INPUT_VCF=filtered.vcf.gz
    else
        INPUT_VCF=${vcf}
        if [ ! -f "\${INPUT_VCF}.tbi" ]; then
            bcftools index -t "\$INPUT_VCF"
        fi
    fi

    n_variants=\$(bcftools view -H "\$INPUT_VCF" | wc -l)
    echo "Variants after all-sites filtering: \$n_variants" >> "\$REPORT"

    # ── 4. Ensure contigs have 'chr' prefix ──────────────────────────────
    first_chr=\$(bcftools query -f '%CHROM\\n' "\$INPUT_VCF" | head -1 || true)
    if [[ ! "\$first_chr" == chr* ]] && [[ -n "\$first_chr" ]]; then
        echo "    Renaming contigs: adding 'chr' prefix" >> "\$REPORT"
        cat > add_chr_map.tsv <<'CHRMAP'
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
CHRMAP
        bcftools annotate --rename-chrs add_chr_map.tsv "\$INPUT_VCF" -Oz -o renamed.vcf.gz
        bcftools index -t renamed.vcf.gz
        INPUT_VCF=renamed.vcf.gz
    fi

    # ── 5. Liftover if GRCh37 ────────────────────────────────────────────
    if [ "\$IS_GRCH37" -eq 1 ]; then

        echo ">>> PERFORMING LIFTOVER GRCh37 -> GRCh38 <<<" >> "\$REPORT"
        echo "    Chain file : ${chain}" >> "\$REPORT"
        echo "    Target ref : ${grch38_ref}" >> "\$REPORT"

        # Run Picard LiftoverVcf
        # Use a recent Picard jar to avoid htsjdk genotype-validation crashes seen in older releases.
        java -Xmx${task.memory.toGiga()}g -jar ${picard_jar} LiftoverVcf \
            I="\$INPUT_VCF" \
            O=lifted.vcf.gz \
            CHAIN=${chain} \
            REJECT=rejected.vcf.gz \
            R=${grch38_ref} \
            WRITE_ORIGINAL_POSITION=true \
            WARN_ON_MISSING_CONTIG=true \
            VALIDATION_STRINGENCY=LENIENT

        n_lifted=\$(bcftools view -H lifted.vcf.gz | wc -l)
        n_rejected=\$(bcftools view -H rejected.vcf.gz 2>/dev/null | wc -l || echo 0)
        pct_lifted=\$(( n_lifted * 100 / (n_lifted + n_rejected > 0 ? n_lifted + n_rejected : 1) ))
        echo "    Lifted variants   : \$n_lifted" >> "\$REPORT"
        echo "    Rejected variants : \$n_rejected" >> "\$REPORT"
        echo "    Liftover rate     : \${pct_lifted}%" >> "\$REPORT"

        # Sort the lifted VCF (liftover can scramble order)
        bcftools sort lifted.vcf.gz -Oz -o "${sample_id}.harmonised.vcf.gz"
        bcftools index -t "${sample_id}.harmonised.vcf.gz"

    else
        # ── GRCh38: pass through (just copy/link) ────────────────────────
        cp "\$INPUT_VCF" "${sample_id}.harmonised.vcf.gz"
        bcftools index -t "${sample_id}.harmonised.vcf.gz"
    fi

    echo "" >> "\$REPORT"
    echo "Output: ${sample_id}.harmonised.vcf.gz" >> "\$REPORT"
    cat "\$REPORT" >&2
    """
}

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

process HARMONIZE_CONTIGS {
    tag { sample_id }

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path study_sites

    output:
    tuple val(sample_id), path("${sample_id}.harmonized.vcf.gz"), path("${sample_id}.harmonized.vcf.gz.tbi"), emit: harmonized

    script:
    """
    set -euo pipefail

    if [ ! -s "${study_sites}" ]; then
        echo "ERROR: Study site list is empty." >&2
        exit 1
    fi

    if grep -q '^chr' "${study_sites}"; then
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

    first_chr=\$(bcftools query -f '%CHROM\\n' "${vcf}" | head -n 1 || true)
    needs_rename=0
    if [ "\$target_style" = "chr" ] && [[ ! "\$first_chr" == chr* ]]; then
        needs_rename=1
    fi
    if [ "\$target_style" = "nochr" ] && [[ "\$first_chr" == chr* ]]; then
        needs_rename=1
    fi

    if [ "\$needs_rename" -eq 1 ]; then
        bcftools annotate --threads ${task.cpus} --rename-chrs chr_map.tsv -Oz -o "${sample_id}.harmonized.vcf.gz" "${vcf}"
        bcftools index -t "${sample_id}.harmonized.vcf.gz"
    else
        cp "${vcf}" "${sample_id}.harmonized.vcf.gz"
        cp "${tbi}" "${sample_id}.harmonized.vcf.gz.tbi"
    fi
    """
}

process SPLIT_STUDY_SITES {
    input:
    path study_sites

    output:
    path "sites_*.tsv", emit: sites_files
    path "chr_list.txt", emit: chr_list

    script:
    """
    python3 - <<'PY'
import sys
from collections import defaultdict

study_sites = "${study_sites}"
sites_by_chr = defaultdict(list)

# Read study sites
with open(study_sites, 'r') as f:
    for line in f:
        idx = line.find('\t')
        if idx != -1:
            chrom = line[:idx]
            sites_by_chr[chrom].append(line)

# Write to per-chromosome site files and print chromosomes in order
with open("chr_list.txt", "w") as chr_list_f:
    for chrom, lines in sites_by_chr.items():
        with open(f"sites_{chrom}.tsv", "w") as out_f:
            out_f.writelines(lines)
        chr_list_f.write(chrom + "\\n")
PY
    """
}

process MERGE_CHROMOSOME {
    tag { chr }

    input:
    tuple val(chr), path(chr_sites)
    path all_harmonized_files

    output:
    tuple val(chr), path("merged_${chr}.vcf.gz"), path("merged_${chr}.vcf.gz.tbi"), emit: merged_chr

    script:
    """
    set -euo pipefail

    ls *.harmonized.vcf.gz | sort -V > all_files.txt

    : > file_list.txt
    declare -A seen_samples=()
    n_kept=0
    n_skipped=0

    while IFS= read -r vcf; do
        [ -n "\$vcf" ] || continue

        # Query samples in VCF
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
            continue
        fi

        printf '%s\\n' "\$vcf" >> file_list.txt
        n_kept=\$((n_kept + 1))
        for s in "\${samples[@]}"; do
            seen_samples["\$s"]=1
        done
    done < all_files.txt

    n_files=\$(wc -l < file_list.txt)
    if [ "\$n_files" -eq 0 ]; then
        echo "ERROR: No non-redundant VCF files available after sample deduplication." >&2
        exit 1
    fi

    if [ "\$n_files" -eq 1 ]; then
        bcftools view -T ${chr_sites} \$(cat file_list.txt) -Oz -o "merged_${chr}.vcf.gz"
        bcftools index -t "merged_${chr}.vcf.gz"
    else
        bcftools merge --threads ${task.cpus} -R ${chr_sites} -l file_list.txt -Oz -o "merged_${chr}.vcf.gz"
        bcftools index -t "merged_${chr}.vcf.gz"
    fi
    """
}

process CONCAT_VCFS {
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path all_merged_chrs
    path chr_list
    path study_sites

    output:
    tuple path("merged.vcf.gz"), path("merged.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    # Reconstruct concat list in correct chromosome order
    > concat_list.txt
    while IFS= read -r chr; do
        if [ -f "merged_\${chr}.vcf.gz" ]; then
            echo "merged_\${chr}.vcf.gz" >> concat_list.txt
        else
            echo "ERROR: File merged_\${chr}.vcf.gz not found" >&2
            exit 1
        fi
    done < ${chr_list}

    # Concatenate and index
    bcftools concat --threads ${task.cpus} -f concat_list.txt -Oz -o merged_raw.vcf.gz
    bcftools index -t merged_raw.vcf.gz

    # Filter/View to ensure correct final sites
    bcftools view --threads ${task.cpus} -T ${study_sites} merged_raw.vcf.gz -Oz -o merged.vcf.gz
    bcftools index -t merged.vcf.gz
    rm -f merged_raw.vcf.gz merged_raw.vcf.gz.tbi
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

    extract_cmd=""
    ${params.extract_snps ? """
    extract_file="${params.extract_snps}"
    if [[ "\$extract_file" != /* ]]; then
        extract_file="${workflow.launchDir}/\$extract_file"
    fi
    if [[ "\$extract_file" == *.csv ]]; then
        awk -F',' 'NR>1 {
            if (\$1 != "" && \$1 != ".") print \$1;
            if (\$2 != "" && \$3 != "") {
                print \$2":"\$3;
                print "chr"\$2":"\$3;
            }
        }' "\$extract_file" > extract_list.txt
    else
        cp "\$extract_file" extract_list.txt
    fi
    extract_cmd="--extract extract_list.txt"
    """ : ""}

    exclude_cmd=""
    ${params.exclude_snps ? """
    exclude_file="${params.exclude_snps}"
    if [[ "\$exclude_file" != /* ]]; then
        exclude_file="${workflow.launchDir}/\$exclude_file"
    fi
    if [[ "\$exclude_file" == *.csv ]]; then
        awk -F',' 'NR>1 {print \$1}' "\$exclude_file" > exclude_list.txt
    else
        cp "\$exclude_file" exclude_list.txt
    fi
    exclude_cmd="--exclude exclude_list.txt"
    """ : ""}

    ${params.exclude_regions ? """
    exclude_regions_file="${params.exclude_regions}"
    if [[ "\$exclude_regions_file" != /* ]]; then
        exclude_regions_file="${workflow.launchDir}/\$exclude_regions_file"
    fi
    cp "\$exclude_regions_file" exclude_regions.bed
    exclude_cmd="\$exclude_cmd --exclude bed0 exclude_regions.bed"
    """ : ""}

    # VCF -> PLINK2 pgen (autosomes only, sort variants)
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --vcf ${vcf} --max-alleles 2 --autosome --make-pgen --sort-vars --out sorted

    # pgen -> BED
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --pfile sorted --set-missing-var-ids @:# \$extract_cmd \$exclude_cmd --make-bed --out base

    # Keep biallelic SNPs only
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile base --snps-only just-acgt --max-alleles 2 --make-bed --out snps

    # Filter out ambiguous A/T and C/G SNPs
    awk '(\$5=="A" && \$6=="T") || (\$5=="T" && \$6=="A") || (\$5=="C" && \$6=="G") || (\$5=="G" && \$6=="C") {print \$2}' snps.bim > ambiguous_snps.txt
    \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile snps --exclude ambiguous_snps.txt --make-bed --out snps_filtered
    mv snps_filtered.bed snps.bed
    mv snps_filtered.bim snps.bim
    mv snps_filtered.fam snps.fam

    # QC (strict first): MAF, missingness, mind
    # Note: To avoid discarding WES study samples due to missingness in non-exonic WGS reference variants,
    # we filter variants by --geno first before applying the --mind filter.
    strict_ok=1
    if ! \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile snps --geno ${params.geno} --make-bed --out snps_geno_filtered || \\
       ! \$plink --threads ${task.cpus} --memory \$mem_mb \\
        --bfile snps_geno_filtered --maf ${params.maf} --mind ${params.mind} --make-bed --out clean_strict
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

process RUN_PCA {
    publishDir "${params.outdir}/pca", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("${bed.baseName}.eigenvec"), path("${bed.baseName}.eigenval")

    script:
    """
    set -euo pipefail
    mem_mb=${task.memory.toMega()}
    plink_bin=\${PLINK:-${params.plink}}

    # PCA from LD-pruned data generated upstream.
    \$plink_bin --threads ${task.cpus} --memory \$mem_mb \\
        --bfile ${bed.baseName} --pca header --out ${bed.baseName}
    """
}

process PLOT_PCA {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    tuple path(eigenvec), path(eigenval)
    path pop_map

    output:
    path "pca_all_samples.png"
    path "pca_all_samples.tsv"

    script:
    """
    python3 - <<'PY'
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

eigenvec_f = "${eigenvec}"
eigenval_f = "${eigenval}"
pop_map_f  = "${pop_map}"

# Read reference pop map (sample_id -> population)
pop_map = {}
with open(pop_map_f, "r") as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            parts = line.split("\\t")
            if len(parts) >= 2:
                pop_map[parts[0]] = parts[1]

# Read eigenvectors
df = pd.read_csv(eigenvec_f, sep="\\s+")
df.columns = [c.lstrip("#") for c in df.columns]

if "IID" not in df.columns:
    if len(df.columns) >= 2:
        df = df.rename(columns={df.columns[1]: "IID"})
    else:
        raise SystemExit("Malformed eigenvec: missing IID column")

pc_cols = [c for c in df.columns if c.upper().startswith("PC")]
if len(pc_cols) < 2:
    raise SystemExit("PCA output has fewer than 2 components")

for c in pc_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df["Population"] = df["IID"].map(pop_map).fillna("study")
df[["IID", "Population"] + pc_cols].to_csv("pca_all_samples.tsv", sep="\\t", index=False)

# Eigenvalue-based explained variance for axis labels (if available)
pc1, pc2 = pc_cols[0], pc_cols[1]
try:
    eig = pd.read_csv(eigenval_f, header=None).iloc[:, 0].astype(float).tolist()
    if len(eig) >= 2 and sum(eig) > 0:
        pc1 = f"{pc1} ({100.0 * eig[0] / sum(eig):.2f}%)"
        pc2 = f"{pc2} ({100.0 * eig[1] / sum(eig):.2f}%)"
except Exception:
    pass

fig, ax = plt.subplots(figsize=(11, 8))

ref_pops = sorted([p for p in df["Population"].unique() if p != "study"])
palette = plt.cm.tab10.colors
for i, pop in enumerate(ref_pops):
    sub = df[df["Population"] == pop]
    ax.scatter(sub[pc_cols[0]], sub[pc_cols[1]], s=28, alpha=0.8,
               color=palette[i % len(palette)], label=pop)

study = df[df["Population"] == "study"]
if not study.empty:
    ax.scatter(study[pc_cols[0]], study[pc_cols[1]], s=24, alpha=0.85,
               c="black", marker="x", label="study")

ax.set_title("PCA: Study Samples + Reference Populations")
ax.set_xlabel(pc1)
ax.set_ylabel(pc2)
ax.grid(True, alpha=0.25)
ax.legend(loc="best", frameon=True)
plt.tight_layout()
plt.savefig("pca_all_samples.png", dpi=300)
PY
    """
}

process PLOT_P_MATRIX {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    tuple path(q_matrix), path(p_matrix), path(fam_file), path(pop_file)
    tuple path(bed), path(bim), path("plink_fam")

    output:
    path "p_matrix_informative.png"
    path "p_matrix_top_markers.tsv"
    path "p_matrix_top50_gene_function.tsv"

    script:
    """
    # Copy the scripts to the working directory
    cp /datos/home/epaaso/ancestry/admix_whole/results/plots/make_p_matrix_plot.py .
    cp /datos/home/epaaso/ancestry/admix_whole/results/plots/annotate_top50_snps.py .
    
    # Ensure they are executable
    chmod +x make_p_matrix_plot.py annotate_top50_snps.py
    
    # We need to modify the script to use the local files instead of hardcoded paths
    sed -i 's|base = Path("/datos/home/epaaso/ancestry/admix_whole/results")|base = Path(".")|g' make_p_matrix_plot.py
    sed -i 's|admix_dir = base / "admixture"|admix_dir = base|g' make_p_matrix_plot.py
    sed -i 's|plink_dir = base / "plink"|plink_dir = base|g' make_p_matrix_plot.py
    sed -i 's|out_dir = base / "plots"|out_dir = base|g' make_p_matrix_plot.py
    sed -i 's|p_path = admix_dir / "merged_pruned.3.P"|p_path = Path("${p_matrix}")|g' make_p_matrix_plot.py
    sed -i 's|q_path = admix_dir / "merged_pruned.3.Q"|q_path = Path("${q_matrix}")|g' make_p_matrix_plot.py
    sed -i 's|pop_path = admix_dir / "merged_pruned.pop"|pop_path = Path("${pop_file}")|g' make_p_matrix_plot.py
    sed -i 's|fam_path = admix_dir / "merged_pruned.fam"|fam_path = Path("${fam_file}")|g' make_p_matrix_plot.py
    sed -i 's|bim_path = plink_dir / "merged_pruned.bim"|bim_path = Path("${bim}")|g' make_p_matrix_plot.py
    
    sed -i 's|BASE = Path("/datos/home/epaaso/ancestry/admix_whole/results/plots")|BASE = Path(".")|g' annotate_top50_snps.py
    
    # Run the scripts
    ./make_p_matrix_plot.py
    ./annotate_top50_snps.py
    
    # Run the plot script again to use the annotations
    ./make_p_matrix_plot.py
    """
}
