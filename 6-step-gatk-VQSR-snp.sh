#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_VQSR_SNP.log) 2>&1

### === Running: 6-step.sh === ###
   
gatk ApplyVQSR \
  -R ../reference/Homo_sapiens_assembly38.fasta \
  -V NA12878.raw_variants.vcf.gz \
  --recal-file NA12878.snp.recal \
  --tranches-file NA12878.snp.tranches \
  --truth-sensitivity-filter-level 99.5 \
  --mode SNP \
  -O NA12878.recalibrated_snps.vcf.gz
