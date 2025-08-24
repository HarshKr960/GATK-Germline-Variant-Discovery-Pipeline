#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_VQSR_SNP.log) 2>&1

gatk ApplyVQSR \
  -R ../reference/Homo_sapiens_assembly38.fasta \
  -V NA12878.recalibrated_snps.vcf.gz \
  --recal-file NA12878.indel.recal \
  --tranches-file NA12878.indel.tranches \
  --truth-sensitivity-filter-level 99.0 \
  --mode INDEL \
  -O NA12878.recalibrated_variants.vcf.gz
