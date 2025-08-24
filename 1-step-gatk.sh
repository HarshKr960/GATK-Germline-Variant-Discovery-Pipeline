#!/bin/bash

# Step 1: Mark Duplicates with GATK 4.6
gatk --java-options "-Xmx8G -XX:ParallelGCThreads=8" MarkDuplicates \
    --INPUT NA12878_aligned_sorted.bam \
    --OUTPUT NA12878_marked_duplicates.bam \
    --METRICS_FILE NA12878_duplicate_metrics.txt \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER coordinate

# Check duplicate metrics
echo "=== Duplicate Metrics ==="
head -8 NA12878_duplicate_metrics.txt | tail -2
