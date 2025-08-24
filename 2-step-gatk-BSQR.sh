#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_BQSR_pipeline.log) 2>&1

echo "=== Starting GATK Germline Variant Discovery Pipeline ==="
date

### === Running: 2A-step-BQAR-gatk.sh === ###
#!/bin/bash

# Step 2A: Generate recalibration table
gatk --java-options "-Xmx8G" BaseRecalibrator \
    --input NA12878_marked_duplicates.bam \
    --reference ../reference/Homo_sapiens_assembly38.fasta \
    --known-sites ../resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ../resources/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites ../resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --output NA12878_recal_data.table \
    --intervals ../reference/tutorial_intervals.list
# Examine recalibration table
head -20 NA12878_recal_data.table

### === Running: 2B-step-BQSR-gatk.sh === ###
#!/bin/bash

# Step 2B: Apply base quality score recalibration
gatk --java-options "-Xmx8G" ApplyBQSR \
    --input NA12878_marked_duplicates.bam \
    --reference ../reference/Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file NA12878_recal_data.table \
    --output NA12878_analysis_ready.bam \
    --intervals ../reference/tutorial_intervals.list
# Verify the analysis-ready BAM
samtools flagstat NA12878_analysis_ready.bam

### === Running: 2C-step-BQSR-gatk.sh === ###
#!/bin/bash

# Generate post-recalibration table for comparison
gatk --java-options "-Xmx8G" BaseRecalibrator \
    --input NA12878_analysis_ready.bam \
    --reference ../reference/Homo_sapiens_assembly38.fasta \
    --known-sites ../resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ../resources/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites ../resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --output NA12878_post_recal_data.table \
    --intervals ../reference/tutorial_intervals.list
# Analyze BQSR effectiveness
gatk AnalyzeCovariates \
    --before-report-file NA12878_recal_data.table \
    --after-report-file NA12878_post_recal_data.table \
    --plots-report-file NA12878_BQSR_report.pdf

