#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_Haplo_pipeline.log) 2>&1

echo "=== Starting GATK Germline Variant Discovery Pipeline ==="
date

# Step 3: Generate GVCF using HaplotypeCaller
gatk --java-options "-Xmx8G -XX:ParallelGCThreads=16" HaplotypeCaller \
    --reference ../reference/Homo_sapiens_assembly38.fasta \
    --input NA12878_analysis_ready.bam \
    --output NA12878.g.vcf.gz \
    --emit-ref-confidence GVCF \
    --dbsnp ../resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --intervals ../reference/tutorial_intervals.list \
    --native-pair-hmm-threads 4 \
    --standard-min-confidence-threshold-for-calling 10.0 \
    --annotation QualByDepth \
    --annotation Coverage \
    --annotation FisherStrand \
    --annotation StrandOddsRatio \
    --annotation MappingQualityRankSumTest \
    --annotation ReadPosRankSumTest \
    --annotation RMSMappingQuality

echo "✅ HaplotypeCaller completed successfully!"

# Verify GVCF generation
echo "=== GVCF Statistics ==="
bcftools stats NA12878.g.vcf.gz | grep -E "(records|SNPs|indels)"

# Examine GVCF structure
echo "=== Sample GVCF Records ==="
bcftools view -H NA12878.g.vcf.gz | head -5

