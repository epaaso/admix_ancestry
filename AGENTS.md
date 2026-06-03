# AGENTS.md

Operational guide for agents working on `admix_whole`.

## Goal

Run supervised ADMIXTURE ancestry inference for the leukemia exome cohort using the 50 final VCF symlinks produced upstream.

Primary study input:

```text
/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz
```

Primary output:

```text
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole
```

Fallback no-HM3 output:

```text
/datos/migccl/leukemia_exoma/outs-results/005-admix_whole_no_hm3
```

## Pipeline Context

Run from:

```text
/datos/home/epaaso/ancestry/admix_whole
```

For the current paused leukemia exome run, read the continuation handoff before
making changes:

```text
AGY_CLI_HANDOFF.md
```

Main entrypoint:

```text
nextflow run main.nf -resume
```

The pipeline:

1. indexes study/reference VCFs,
2. collects study sites,
3. harmonizes/lifts over references to GRCh38,
4. merges study and reference VCFs at study sites,
5. converts to PLINK,
6. applies QC and LD pruning,
7. runs supervised ADMIXTURE,
8. writes ancestry and PCA outputs.

## Fixed Inputs

Study VCF glob:

```text
/datos/migccl/leukemia_exoma/outs-results/final_vcfs_new_names/*.vcf.gz
```

Reference populations:

```text
EUR:/datos/migccl/ancestry_refs/vcfs/eur_40/*.vcf.gz;AFR:/datos/migccl/ancestry_refs/vcfs/afr_40/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/SimonsVCFs/*.vcf.gz;MXL:/datos/migccl/ancestry_refs/vcfs/maya_hgdp/*.vcf.gz
```

Reference genome:

```text
/datos/migccl/arriba_refs/GRCh38.primary_assembly.genome.fa
```

Primary ancestry SNP list:

```text
/datos/home/epaaso/ancestry/hm3_combined_ids_and_loci.txt
```

Leukemia regions to exclude:

```text
/datos/home/epaaso/ancestry/leukemia_genes_hg38.bed
```

SSH entry:

```text
inmegen-gpu -> epaaso@192.168.112.101
```

## Default Parameters

Use these unless the user explicitly changes the analysis:

```text
--k 3
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

Primary run also includes:

```text
--extract_snps /datos/home/epaaso/ancestry/hm3_combined_ids_and_loci.txt
--exclude_regions /datos/home/epaaso/ancestry/leukemia_genes_hg38.bed
```

Fallback no-HM3 run omits only `--extract_snps`.

## Agent Rules

- Do not delete FASTQs, reference files, old workflow outputs, or final VCF symlinks.
- Do not change the study input glob unless the user explicitly requests it.
- Always inspect `.nextflow.log`, the watcher logs, and the failed work directory before editing code.
- If continuing the paused leukemia exome ancestry run, read `AGY_CLI_HANDOFF.md`
  first and preserve the documented context.
- Use `-resume` for every rerun.
- Make narrow, targeted fixes in `main.nf`, `nextflow.config`, or helper scripts only when the error requires it.
- Do not run `git reset --hard`, `git checkout --`, or destructive cleanup commands.
- If a fix changes behavior, record it in the watcher status file.
- If the same error fingerprint repeats 3 times, or if 5 total agent repair attempts are used, stop and notify the user.

## Success Criteria

The run is successful when either the primary HM3 run or the fallback no-HM3 run contains:

```text
ancestry/ancestry.tsv
plots/ancestry_plot.png
plots/pca_all_samples.png
```

The watcher should also record:

- final mode: `primary_hm3` or `fallback_no_hm3`,
- number of study VCFs,
- number of output rows in `ancestry.tsv`,
- final status and completion time.
