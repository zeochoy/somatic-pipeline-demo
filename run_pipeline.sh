#!/usr/bin/env bash

set +eu

#---
## set up
#---
SAMPLE_ID=$1
CONDA_ENV_NAME=$2

HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"
CONDA_ACT_PATH="${HOME}/miniconda3/bin/activate"
CONDA_DEACT_PATH="${HOME}/miniconda3/bin/deactivate"

#---
## help script
#---
usage="
$(basename "$0") [SAMPLE_ID] [CONDA ENV NAME] [-h] - shell script to run the somatic variant calling pipeline

where:
    -h  show this help text
"

if [ "$1" == "-h" ]; then
  echo "$usage"
  exit 0
fi

#---
## main script
#---
cd ./data

### run snakemake
source ${CONDA_ACT_PATH} ${CONDA_ENV_NAME}

snakemake bam/${SAMPLE_ID}T_sorted.bam.bai
snakemake bam/${SAMPLE_ID}T_sorted_markdup.bam.bai
snakemake bam/${SAMPLE_ID}B_sorted.bam.bai
snakemake bam/${SAMPLE_ID}B_sorted_markdup.bam.bai
snakemake vcf/${SAMPLE_ID}_gatkfilt_vanfilt.vcf
snakemake annotated/${SAMPLE_ID}.hg19_multianno.txt
snakemake annotated/${SAMPLE_ID}_oncokb_maf.txt
snakemake annotated/${SAMPLE_ID}_annotated.txt

source ${CONDA_DEACT_PATH}
set -eu
