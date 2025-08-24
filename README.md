GATK Germline Variant Discovery Pipeline (GRCh38)

This repository contains a stepwise GATK 4.x pipeline for germline variant discovery on Illumina whole-genome/exome sequencing data.
The workflow follows GATK Best Practices, from read alignment to variant quality score recalibration (VQSR).

📂 Repository Structure
```
alignbwa.sh → Align reads to GRCh38 using BWA-MEM, sort & index with SAMtools

1-step-gatk.sh → Mark duplicates with GATK MarkDuplicates

2-step-gatk-BSQR.sh → Base Quality Score Recalibration (BQSR) – BaseRecalibrator + ApplyBQSR + post-recalibration analysis

3-step-gatk-Haplo.sh → Variant calling with HaplotypeCaller (GVCF mode)

4-step-gatk-genotype.sh → Joint genotyping with GenotypeGVCFs

5-step-gatk-VQSR.sh → Build recalibration model for SNPs

6-step-gatk-VQSR-snp.sh → Apply VQSR to SNPs

7-step-gatk-VQSR-INDEL.sh → Build recalibration model for indels

9-step-gatk-QC-PASS-Variant.sh → Apply VQSR to indels, final PASS variants
```
⚙️ Workflow Steps

Alignment → Align reads to GRCh38 with BWA-MEM, sort & index BAM.

MarkDuplicates → Remove PCR/optical duplicates.

BQSR → Recalibrate base quality scores with known sites (dbSNP, Mills, 1000G indels).

HaplotypeCaller → Call variants per-sample in GVCF mode.

GenotypeGVCFs → Joint genotyping across samples.

VQSR (SNPs & Indels) → Build & apply variant recalibration models.

Final QC → Filtered, high-confidence VCF file of germline variants.

🛠️ Tools Used

BWA

SAMtools

GATK 4.x

bcftools

📊 Output Files
```
NA12878_aligned_sorted.bam → Aligned BAM + index

NA12878_marked_duplicates.bam → Deduplicated BAM

NA12878_analysis_ready.bam → BQSR-corrected BAM

NA12878.g.vcf.gz → GVCF file from HaplotypeCaller

NA12878.raw_variants.vcf.gz → Raw variants

NA12878.recalibrated_variants.vcf.gz → Final VQSR-filtered variants
```
▶️ Usage

Each step is run sequentially. Example:
```
# Step 1: Align reads
bash alignbwa.sh

# Step 2: Mark duplicates
bash 1-step-gatk.sh

# Step 3: Base Quality Score Recalibration
bash 2-step-gatk-BSQR.sh

# Step 4: HaplotypeCaller
bash 3-step-gatk-Haplo.sh

# Step 5: Genotyping
bash 4-step-gatk-genotype.sh

# Step 6: SNP Recalibration
bash 5-step-gatk-VQSR.sh
bash 6-step-gatk-VQSR-snp.sh

# Step 7: INDEL Recalibration
bash 7-step-gatk-VQSR-INDEL.sh

# Step 8: Apply final filters
bash 9-step-gatk-QC-PASS-Variant.sh
```
📜 References

GATK Best Practices (Broad Institute)

Van der Auwera GA, et al. (2013) From FastQ to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Nat Protoc.
