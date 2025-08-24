#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_VQSR_pipeline.log) 2>&1

gatk VariantRecalibrator \
  -R ../reference/Homo_sapiens_assembly38.fasta \
  -V NA12878.raw_variants.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
    ../resources/hg38_v0_hapmap_3.3.hg38.vcf.gz \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 \
    ../resources/hg38_v0_1000G_omni2.5.hg38.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
    ../resources/hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
    ../resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an MQ \
  -mode SNP \
  --output NA12878.snp.recal \
  --tranches-file NA12878.snp.tranches \
  --rscript-file NA12878.snp.plots.R
