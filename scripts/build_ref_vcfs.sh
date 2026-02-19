#!/usr/bin/env bash
set -euo pipefail

# Build whole-genome per-sample VCFs from 1000 Genomes per-chromosome VCFs.
#
# For each sample list file provided, subsets the 1KG VCFs to those samples,
# concatenates all chromosomes, then splits into one VCF per sample.
# Output goes to a folder named after the sample list (without extension).
#
# Usage:
#   ./build_ref_vcfs.sh [OPTIONS] <sample_list1> [sample_list2] ...
#
# Options:
#   -r, --raw-vcfs DIR    Directory with chr*.1kg.vcf.gz files
#                          (default: nextflow_admix_structure/results/reference/raw_vcfs)
#   -o, --outdir DIR      Base output directory (default: current dir)
#   -t, --threads N       Threads for bcftools (default: 8)
#   -h, --help            Show this help
#
# Example:
#   ./build_ref_vcfs.sh -r /path/to/raw_vcfs \
#       meta/afr_40.samples meta/eur_40.samples meta/mxl_40.samples

RAW_VCFS="/datos/home/epaaso/ancestry/nextflow_admix_structure/results/reference/raw_vcfs"
OUTDIR="."
THREADS=8

usage() {
    sed -n '3,/^$/p' "$0" | sed 's/^# \?//'
    exit "${1:-0}"
}

# ── Parse arguments ──────────────────────────────────────────────────────────
SAMPLE_LISTS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--raw-vcfs) RAW_VCFS="$2"; shift 2 ;;
        -o|--outdir)   OUTDIR="$2";   shift 2 ;;
        -t|--threads)  THREADS="$2";  shift 2 ;;
        -h|--help)     usage 0 ;;
        -*)            echo "Unknown option: $1" >&2; usage 1 ;;
        *)             SAMPLE_LISTS+=("$1"); shift ;;
    esac
done

if [[ ${#SAMPLE_LISTS[@]} -eq 0 ]]; then
    echo "Error: provide at least one sample list file" >&2
    usage 1
fi

# ── Validate ─────────────────────────────────────────────────────────────────
for f in "${SAMPLE_LISTS[@]}"; do
    [[ -f "$f" ]] || { echo "Error: sample list not found: $f" >&2; exit 1; }
done

chrom_vcfs=()
for c in $(seq 1 22); do
    vcf="${RAW_VCFS}/chr${c}.1kg.vcf.gz"
    [[ -f "$vcf" ]] || { echo "Error: missing ${vcf}" >&2; exit 1; }
    chrom_vcfs+=("$vcf")
done

echo "Raw VCF dir : $RAW_VCFS"
echo "Output dir  : $OUTDIR"
echo "Threads     : $THREADS"
echo "Sample lists: ${SAMPLE_LISTS[*]}"
echo ""

# ── Process each sample list ─────────────────────────────────────────────────
for sample_list in "${SAMPLE_LISTS[@]}"; do
    list_name="$(basename "$sample_list" .samples)"
    list_name="$(basename "$list_name" .txt)"
    list_name="$(basename "$list_name" .ids)"
    dest="${OUTDIR}/${list_name}"
    mkdir -p "$dest"

    n_samples=$(grep -c . "$sample_list" || true)
    echo "=== ${list_name}: ${n_samples} samples → ${dest}/ ==="

    # 1. Subset + concat all chromosomes for this sample list
    merged="${dest}/${list_name}_all_chroms.vcf.gz"
    if [[ -f "$merged" ]]; then
        echo "  Merged VCF already exists, skipping subset+concat"
    else
        tmpdir=$(mktemp -d "${dest}/tmp.XXXXXX")
        trap "rm -rf '$tmpdir'" EXIT

        for vcf in "${chrom_vcfs[@]}"; do
            chrom=$(basename "$vcf" .1kg.vcf.gz)
            out_chrom="${tmpdir}/${chrom}.vcf.gz"
            echo "  Subsetting ${chrom}..."
            if bcftools view --threads "$THREADS" --force-samples \
                -S "$sample_list" "$vcf" -Oz -o "$out_chrom" 2>&1; then
                bcftools index -t "$out_chrom"
            else
                echo "  WARNING: skipping ${chrom} (file may be truncated)" >&2
                rm -f "$out_chrom"
            fi
        done

        echo "  Concatenating chromosomes..."
        ls "${tmpdir}"/chr*.vcf.gz | sort -V > "${tmpdir}/concat_list.txt"
        bcftools concat --threads "$THREADS" \
            -f "${tmpdir}/concat_list.txt" -Oz -o "$merged"
        bcftools index -t "$merged"

        rm -rf "$tmpdir"
        trap - EXIT
    fi

    # 2. Split into per-sample VCFs
    echo "  Splitting into per-sample VCFs..."
    while IFS= read -r sample_id; do
        [[ -z "$sample_id" ]] && continue
        sample_vcf="${dest}/${sample_id}.vcf.gz"
        if [[ -f "$sample_vcf" ]]; then
            continue
        fi
        bcftools view --threads 2 -s "$sample_id" \
            --min-ac 1 "$merged" -Oz -o "$sample_vcf"
        bcftools index -t "$sample_vcf"
    done < "$sample_list"

    n_done=$(ls "${dest}"/*.vcf.gz 2>/dev/null | grep -v all_chroms | wc -l)
    echo "  Done: ${n_done} per-sample VCFs in ${dest}/"
    echo ""
done

echo "All done."
