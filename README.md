GATK Germline Variant Discovery Pipeline (GRCh38)

This repository contains a stepwise GATK 4.x pipeline for germline variant discovery on Illumina whole-genome/exome sequencing data. The workflow follows GATK Best Practices, from read trimming and alignment through variant quality score recalibration (VQSR) and optional post-VQSR QC.

📂 Repository Structure 
.
├── Fastp-adapter.sh                  # Step 0: adapter/quality trimming + HTML/JSON report (fastp)
├── alignbwa.sh                       # Step 1: align to GRCh38 with BWA-MEM, sort & index
├── 1-step-gatk.sh                    # Step 2: MarkDuplicates (GATK)
├── 2-step-gatk-BSQR.sh               # Step 3: BQSR (BaseRecalibrator + ApplyBQSR + post-recalibration check)
├── 3-step-gatk-Haplo.sh              # Step 4: HaplotypeCaller (GVCF mode)
├── 4-step-gatk-genotype.sh           # Step 5: GenotypeGVCFs (joint genotyping)
├── 5-step-gatk-VQSR.sh               # Step 6: Build VQSR model for SNPs
├── 6-step-gatk-VQSR-snp.sh           # Step 7: Apply VQSR to SNPs
├── 7-step-gatk-VQSR-INDEL.sh         # Step 8: Build VQSR model for INDELs
├── 9-step-gatk-QC-PASS-Variant.sh    # Step 9: Apply VQSR to INDELs, emit final PASS variants
├── post_vqsr_analysis.sh             # Step 10 (optional): post-VQSR QC & reporting (bcftools/MultiQC)
└── tutorial_intervals.list           # Optional target regions for quick tutorial runs / chr20 demo


Fastp-adapter.sh runs paired-end adapter/quality trimming and emits fastp_report.html/json. 

post_vqsr_analysis.sh summarizes PASS variants, splits SNPs/INDELs, computes Ti/Tv, density by chromosome, and (if available) builds a MultiQC HTML report. 

tutorial_intervals.list currently contains a chr20 demo region (chr20:10000000-11000000) you can pass to GATK steps with -L. 

⚙️ Workflow Steps (updated)
```
0) Pre-processing (fastp) → Adapter/quality trimming + QC HTML/JSON. 

1) Alignment → Align reads to GRCh38 with BWA-MEM, sort & index BAM.
2) MarkDuplicates → Remove PCR/optical duplicates (GATK).
3) BQSR → Recalibrate base qualities using known sites (dbSNP, Mills, 1000G indels).
4) HaplotypeCaller → Per-sample GVCF mode.
5) GenotypeGVCFs → Joint genotyping across samples.
6–9) VQSR (SNPs & INDELs) → Build & apply recalibration models; produce final PASS VCF.
10) Post-VQSR QC (optional) → Stats, PASS extraction, SNP/INDEL split, Ti/Tv, density, MultiQC. 

Tip: For tutorial/quick runs, use -L tutorial_intervals.list with BQSR/HaplotypeCaller/GenotypeGVCFs to restrict processing to the demo region. 
```

🛠️ Tools Used

fastp (trimming + QC) 

BWA, SAMtools

GATK 4.x

bcftools (stats/filters/splitting) & MultiQC (optional aggregation) 

📊 Output Files (expanded)
```
NA12878_R1.trimmed.fastq.gz, NA12878_R2.trimmed.fastq.gz, fastp_report.html, fastp_report.json (pre-processing) 

NA12878_aligned_sorted.bam → aligned BAM + index

NA12878_marked_duplicates.bam → dedup BAM

NA12878_analysis_ready.bam → BQSR-corrected BAM

NA12878.g.vcf.gz → GVCF from HaplotypeCaller

NA12878.raw_variants.vcf.gz → raw variants

NA12878.recalibrated_variants.vcf.gz → final VQSR-filtered variants

(Optional post-VQSR) NA12878.high_quality*.vcf.gz, NA12878_variant_stats.txt, chromosome_variant_counts.txt, NA12878_analysis_report.html (if MultiQC available). 
```

▶️ Usage (updated)
```
# Make scripts executable once
chmod +x Fastp-adapter.sh alignbwa.sh 1-step-gatk.sh 2-step-gatk-BSQR.sh \
  3-step-gatk-Haplo.sh 4-step-gatk-genotype.sh 5-step-gatk-VQSR.sh \
  6-step-gatk-VQSR-snp.sh 7-step-gatk-VQSR-INDEL.sh 9-step-gatk-QC-PASS-Variant.sh \
  post_vqsr_analysis.sh

# Step 0: Pre-processing (fastp)
# (Adjust input FASTQs; emits trimmed FASTQs + HTML/JSON report)
./Fastp-adapter.sh    # :contentReference[oaicite:10]{index=10}

# Step 1: Align reads
bash alignbwa.sh

# Step 2: Mark duplicates
bash 1-step-gatk.sh

# Step 3: Base Quality Score Recalibration (optionally restrict to demo intervals)
bash 2-step-gatk-BSQR.sh            # add: -L tutorial_intervals.list   # :contentReference[oaicite:11]{index=11}

# Step 4: HaplotypeCaller (optionally restrict to demo intervals)
bash 3-step-gatk-Haplo.sh           # add: -L tutorial_intervals.list   # :contentReference[oaicite:12]{index=12}

# Step 5: Genotyping (optionally restrict to demo intervals)
bash 4-step-gatk-genotype.sh        # add: -L tutorial_intervals.list   # :contentReference[oaicite:13]{index=13}

# Step 6–9: VQSR
bash 5-step-gatk-VQSR.sh
bash 6-step-gatk-VQSR-snp.sh
bash 7-step-gatk-VQSR-INDEL.sh
bash 9-step-gatk-QC-PASS-Variant.sh

# Step 10 (optional): Post-VQSR QC & reporting
./post_vqsr_analysis.sh             # produces stats, PASS-only sets, Ti/Tv, density, MultiQC (if installed)  # :contentReference[oaicite:14]{index=14}
```
📜 References

GATK Best Practices (Broad Institute)

Van der Auwera GA, et al. 2013. From FastQ to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Nat Protoc.
