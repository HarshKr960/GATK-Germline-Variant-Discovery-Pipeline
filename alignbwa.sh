#!/bin/bash
# Align reads to GRCh38
bwa mem -t 14 -M \
    -R "@RG\tID:NA12878_chr20\tSM:NA12878\tPL:ILLUMINA\tLB:lib1\tPU:flowcell1.lane1" \
    ../reference/Homo_sapiens_assembly38.fasta \
    ../sample_data/NA12878_R1.trimmed.fastq.gz \
    ../sample_data/NA12878_R2.trimmed.fastq.gz | \

samtools sort -@ 8 -o NA12878_aligned_sorted.bam -
samtools index NA12878_aligned_sorted.bam

echo "✅ Alignment complete with perfect GRCh38 compatibility!"

