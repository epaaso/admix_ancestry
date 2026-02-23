# admix_whole

Pipeline Nextflow (DSL2) para inferencia de ancestría con ADMIXTURE supervisado a partir de VCFs completos.

## Qué hace

1. Carga VCFs de estudio y de referencia.
2. Indexa VCFs.
3. Extrae sitios del estudio (`CHROM`, `POS`).
4. Pre-armoniza referencias (opcional) y cachea resultados.
5. Hace `merge` restringido a sitios de estudio.
6. Convierte a PLINK, aplica QC, pruning LD.
7. Ejecuta ADMIXTURE supervisado.
8. Resume componentes y genera tabla/plot.
9. Calcula PCA sobre datos podados por LD y grafica estudio + referencias.

## Requisitos

- Nextflow 25.x
- `bcftools`
- `plink2`
- `plink` (1.9)
- `admixture`
- Python 3 con `pandas` y `matplotlib`

`nextflow.config` tiene `conda.enabled = true`, y `environment.yml` cubre dependencias de Python/BCFtools/PLINK para pasos auxiliares.

## Ejecución rápida

Desde `admix_whole/`:

```bash
nextflow run main.nf \
  --input "/ruta/estudio/*.vcf.gz" \
  --ref_pops "EUR:/ruta/eur/*.vcf.gz;AFR:/ruta/afr/*.vcf.gz;MXL:/ruta/mxl/*.vcf.gz" \
  --outdir results \
  --k 3
```

## Modos de referencia

Solo usa uno de estos modos:

- `--ref_pops "Pop1:glob;Pop2:glob;..."`
- `--ref_vcf_dir "/ruta/*.vcf.gz" --ref_pop_map pop_map.tsv`
- `--ref_vcf_dir "/ruta/*.vcf.gz" --ref_manifests "Pop1:manifest.tsv;Pop2:manifest.tsv"`

Importante: en `--ref_pops`/`--ref_manifests`, el pipeline infiere `sample_id` desde el nombre del archivo VCF. Ahora limpia sufijos técnicos comunes (por ejemplo `.hard-filtered`, `.annotated.nh`) para que coincida con el `IID` real del VCF.
Si tus archivos usan convenciones distintas, usa `--ref_pop_map` con IDs exactos.

## Parámetros importantes

- `--input`: glob de VCFs de estudio.
- `--outdir`: carpeta de salida (default `results`).
- `--k`: K esperado (si difiere del número real de poblaciones de referencia, el pipeline usa el real).
- `--maf`, `--geno`, `--mind`: filtros QC.
- `--relaxed_geno`: filtro alterno si QC estricto elimina poblaciones de referencia.
- `--max_cpus`, `--max_memory`: recursos por proceso (según `nextflow.config`).
- `--preharmonize_refs` (default `true`): pre-armoniza contigs en referencias antes del merge.
- `--ref_harmonized_cache` (default `${outdir}/ref_harmonized_cache`): caché persistente de referencias armonizadas.

## Salidas

Se publican en `${outdir}`:

- `merged/`: `merged.vcf.gz` final.
- `plink/`: archivos PLINK intermedios/finales.
- `admixture/`: archivos `.Q`, `.P`, `.pop`, log de admixture.
- `ancestry/`: `ancestry.tsv`.
- `pca/`: `merged_pruned.eigenvec`, `merged_pruned.eigenval`.
- `plots/`: `ancestry_plot.png`, `ancestry_study_only.tsv`, `pca_all_samples.png`, `pca_all_samples.tsv`.

## Rendimiento

- El `merge` de bcftools suele ser I/O-bound (CPU baja es normal).
- Para corridas repetidas, deja activado `--preharmonize_refs true` para reutilizar el caché.
- Si cambias el universo de sitios de estudio, el prearmonizado puede regenerarse porque depende de esos sitios.

## Notas

- El pipeline asume VCFs bgzip + index (`.tbi`), y si falta índice lo crea.
- Hay deduplicación de muestras antes del merge; si hay VCFs con traslape parcial de muestras, falla de forma explícita.
- En los modos de referencia por glob, el `sample_id` de referencia se infiere del nombre del archivo VCF.
