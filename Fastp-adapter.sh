fastp -i NA12878_chr20_R1.fastq.gz -I NA12878_chr20_R2.fastq.gz \
  -o NA12878_R1.trimmed.fastq.gz -O NA12878_R2.trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --thread 8 \
  --html fastp_report.html \
  --json fastp_report.json
