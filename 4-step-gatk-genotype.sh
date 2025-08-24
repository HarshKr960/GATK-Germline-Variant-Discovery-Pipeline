#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_Geno_pipeline.log) 2>&1


#Running: 4-step.sh 

gatk GenotypeGVCFs \
  --reference ../reference/Homo_sapiens_assembly38.fasta \
  -V NA12878.g.vcf.gz \
  -O NA12878.raw_variants.vcf.gz
