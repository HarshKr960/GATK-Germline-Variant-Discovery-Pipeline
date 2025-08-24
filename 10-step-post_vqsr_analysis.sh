#!/bin/bash
set -euo pipefail
exec > >(tee -a gatk_post_VQSR_analysis.log) 2>&1

echo "=== Starting Post-VQSR Analysis Pipeline ==="
date

# Step 9A: Basic Variant Statistics
echo "=== Generating Basic Variant Statistics ==="
bcftools stats NA12878.recalibrated_variants.vcf.gz > NA12878_variant_stats.txt

# Count variants by type
echo "=== Counting Variants by Type ==="
echo "SNPs:"
bcftools view -v snps NA12878.recalibrated_variants.vcf.gz | bcftools stats | grep "number of SNPs"
echo "INDELs:"
bcftools view -v indels NA12878.recalibrated_variants.vcf.gz | bcftools stats | grep "number of indels"

# Check VQSR filter status
echo "=== VQSR Filter Summary ==="
echo "Total variants:"
bcftools view -H NA12878.recalibrated_variants.vcf.gz | wc -l
echo "PASS variants:"
bcftools view -f PASS NA12878.recalibrated_variants.vcf.gz -H | wc -l

# Step 9B: Extract High-Quality Variants
echo "=== Extracting High-Quality PASS Variants ==="
bcftools view -f PASS NA12878.recalibrated_variants.vcf.gz -O z -o NA12878.high_quality.vcf.gz
bcftools index NA12878.high_quality.vcf.gz

# Step 9C: Separate SNPs and INDELs for further analysis
echo "=== Separating SNPs and INDELs ==="
bcftools view -v snps NA12878.high_quality.vcf.gz -O z -o NA12878.high_quality_SNPs.vcf.gz
bcftools view -v indels NA12878.high_quality.vcf.gz -O z -o NA12878.high_quality_INDELs.vcf.gz

# Index the separated files
bcftools index NA12878.high_quality_SNPs.vcf.gz
bcftools index NA12878.high_quality_INDELs.vcf.gz

# Step 9D: Quality-based filtering
echo "=== Applying Additional Quality Filters ==="
bcftools view -i 'QUAL>30 && INFO/DP>10' NA12878.high_quality.vcf.gz -O z -o NA12878.filtered.vcf.gz
bcftools index NA12878.filtered.vcf.gz

# Step 9E: Generate Ti/Tv ratio (important QC metric)
echo "=== Calculating Ti/Tv Ratio ==="
bcftools stats NA12878.high_quality_SNPs.vcf.gz | grep "TSTV" | head -2

# Step 9F: Variant density analysis
echo "=== Generating Variant Density by Chromosome ==="
bcftools query -f '%CHROM\n' NA12878.high_quality.vcf.gz | sort | uniq -c > chromosome_variant_counts.txt

# Step 9G: Prepare for annotation (create summary for manual review)
echo "=== Creating Variant Summary for Review ==="
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/DP\n' \
    NA12878.filtered.vcf.gz | head -20 > variant_preview.tsv

# Step 9H: Generate MultiQC report if possible
if command -v multiqc &> /dev/null; then
    echo "=== Generating MultiQC Report ==="
    multiqc . --title "NA12878_Post_VQSR_Analysis" --filename NA12878_analysis_report.html
else
    echo "MultiQC not found - install with: pip install multiqc"
fi

# Final summary
echo "=== Analysis Complete ==="
echo "Files generated:"
echo "- NA12878.high_quality.vcf.gz (PASS variants only)"
echo "- NA12878.high_quality_SNPs.vcf.gz (PASS SNPs only)"  
echo "- NA12878.high_quality_INDELs.vcf.gz (PASS INDELs only)"
echo "- NA12878.filtered.vcf.gz (Additional quality filters applied)"
echo "- NA12878_variant_stats.txt (Comprehensive statistics)"
echo "- chromosome_variant_counts.txt (Variant density by chromosome)"

echo "Next steps:"
echo "1. Annotate variants with VEP or ANNOVAR"
echo "2. Filter by population frequency (if databases available)"
echo "3. Perform functional impact analysis"
echo "4. Generate final clinical/research report"

date
echo "=== Post-VQSR Analysis Pipeline Complete ==="
