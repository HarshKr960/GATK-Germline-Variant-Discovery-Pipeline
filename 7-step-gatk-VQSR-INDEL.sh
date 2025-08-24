#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_VQSR_Indel.log) 2>&1

gatk VariantRecalibrator \
  -R ../reference/Homo_sapiens_assembly38.fasta \
  -V NA12878.recalibrated_snps.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 \
      ../resources/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
      ../resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  -an QD \
  -an FS \
  -an SOR \
  -mode INDEL \
  --output NA12878.indel.recal \
  --tranches-file NA12878.indel.tranches \
  --rscript-file NA12878.indel.plots.R

