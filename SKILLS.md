# SKILLS.md

Reusable operating procedures for `admix_whole`.

## Current Handoff

For the paused leukemia exome ancestry run, read this file first:

```bash
sed -n '1,260p' AGY_CLI_HANDOFF.md
```

It documents the stopped run, the `MERGE_ALL` failure history, why a single giant
threaded `bcftools merge` is risky, and the recommended chromosome-sharded merge
continuation strategy.

## Preflight

Run from:

```bash
cd /datos/home/epaaso/ancestry/admix_whole
```

Check study VCF count:

```bash
find /datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names \
  -maxdepth 1 -name '*.vcf.gz' | wc -l
```

Expected: `50`.

Check required tools:

```bash
command -v nextflow
command -v bcftools
command -v plink2
command -v /home/epaaso/bin/plink
command -v admixture
command -v agy
```

Check required files:

```bash
test -s /datos/migccl/arriba_refs/GRCh38.primary_assembly.genome.fa
test -s /datos/migccl/arriba_refs/GRCh38.primary_assembly.genome.fa.fai
test -s /datos/migccl/arriba_refs/GRCh38.primary_assembly.genome.dict
test -s /datos/home/epaaso/ancestry/hm3_combined_ids_and_loci.txt
test -s /datos/home/epaaso/ancestry/leukemia_genes_hg38.bed
```

Check reference globs:

```bash
find /datos/migccl/ancestry_refs/vcfs/eur_40 -maxdepth 1 -name '*.vcf.gz' | wc -l
find /datos/migccl/ancestry_refs/vcfs/afr_40 -maxdepth 1 -name '*.vcf.gz' | wc -l
find /datos/migccl/ancestry_refs/vcfs/SimonsVCFs -maxdepth 1 -name '*.vcf.gz' | wc -l
find /datos/migccl/ancestry_refs/vcfs/maya_hgdp -maxdepth 1 -name '*.vcf.gz' | wc -l
```

## Index Study VCFs

The final VCF folder contains symlinks. Create adjacent indexes if missing:

```bash
for vcf in /datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz; do
  if [ ! -s "${vcf}.tbi" ] && [ ! -s "${vcf}.csi" ]; then
    bcftools index -t "$vcf"
  fi
done
```

## Primary HM3 Run

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

## Fallback No-HM3 Run

Use only if the primary run completes but marker retention is too poor for reliable ADMIXTURE/PCA.

```bash
nextflow run main.nf -resume \
  --input "/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz" \
  --outdir /datos/migccl/leukemia_exoma/outs-results/005-admix_whole_no_hm3 \
  --ref_pops "EUR:/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz;AFR:/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz" \
  --k 3 \
  --exclude_regions /datos/home/epaaso/ancestry/leukemia_genes_hg38.bed \
  --maf 0.05 --geno 0.05 --mind 0.8 --relaxed_geno 0.99 \
  --ld_window 50 --ld_step 5 --ld_r2 0.2 \
  --max_cpus 18 --max_memory "120 GB"
```

## Output Validation

For a candidate output directory:

```bash
test -s "$outdir/ancestry/ancestry.tsv"
test -s "$outdir/plots/ancestry_plot.png"
test -s "$outdir/plots/pca_all_samples.png"
wc -l "$outdir/ancestry/ancestry.tsv"
```

Check marker retention:

```bash
test -s "$outdir/plink/merged_pruned.bim"
wc -l "$outdir/plink/merged_pruned.bim"
```

The autowatcher treats fewer than 100 pruned markers as a poor-retention signal and tries the no-HM3 fallback.

## Autowatcher Status

The watcher status file is:

```bash
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch/last_status.md
```

Refresh it without launching or stopping any pipeline process:

```bash
cd /datos/home/epaaso/ancestry/admix_whole
bash scripts/run_admix_autowatch.sh --status-only
```

Read the live Nextflow log:

```bash
tail -f /datos/migccl/leukemia_exoma/outs-results/005-admix_whole_autowatch/logs/nextflow_primary_hm3.log
```

If the current process is `MERGE_ALL`, watch the active merge file from the work directory shown in `last_status.md`:

```bash
watch -n 60 'date; ls -lh /datos/home/epaaso/ancestry/admix_whole/work/b3/b7ea0b244889c9b7f0cfb043df49c9/merged_raw.vcf.gz'
```

## `agy` Recovery Prompt Pattern

Use the installed supported form:

```bash
agy --dangerously-skip-permissions --print --print-timeout 15m \
  --add-dir /datos/migccl/leukemia_exoma \
  --add-dir /datos/migccl/ancestry_refs \
  "PROMPT TEXT"
```

Prompt must include:

- the error fingerprint,
- relevant `.nextflow.log` tail,
- failed work directory path if available,
- instruction to read `AGENTS.md` and `SKILLS.md`,
- instruction to make targeted fixes and rerun with `-resume`,
- instruction not to delete input/reference/old output files.

## Notification

Preferred recipient:

```text
ernesto.paas@ciencias.unam.mx
```

If `mail` is missing, install `mailutils`:

```bash
sudo apt-get update
sudo apt-get install -y mailutils
```

Send blocked notification:

```bash
printf '%s\n' "message" | mail -s "admix_whole autowatch blocked" ernesto.paas@ciencias.unam.mx
```

If mail cannot be installed, use `/usr/bin/wall` as a local fallback and write the status log.
