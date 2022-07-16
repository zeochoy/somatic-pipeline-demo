# Somatic variant calling snakemake pipeline demo
A somatic variant calling snakemake pipeline based on GATK4 Best Practice for human paired targeted sequencing.
**Intended for _DEMO_ only, not for production use.**

## Prerequisite
Assume you have your:
1. miniconda3 environment set up with snakemake and pandas installed;
2. fastq placed under data/fastq;
3. [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) initialised and set up;
4. clone the latest [mskcc-vcf2maf](https://github.com/mskcc/vcf2maf), [oncokb-annotator](https://github.com/oncokb/oncokb-annotator), [VEP](https://github.com/Ensembl/ensembl-vep);
5. set up data/config.yaml - incl. ONCOKB token, path to GATK4, ANNOVAR, VEP, ONCOKB etc.

## Dependencies
snakemake pandas

## Usage
```
run_pipeline.sh SAMPLE_ID YOUR_SNAKEMAKE_CONDA_ENV
```

### Expected output
1. aligned and duplicates marked BAM under `data/bam`.
2. VCF under `data/vcf`.
3. VEP, ANNOVAR, OncoKB annotated variant call results under `data/annotated`.

## Directory
```
.
|-- README.md
|-- data
|   |-- Snakefile
|   |-- annotated
|   |-- bam
|   |-- combine_filt_sort_anno.py
|   |-- config.yaml
|   |-- fastq
|   `-- vcf
|-- run_pipeline.sh
`-- tools

```
