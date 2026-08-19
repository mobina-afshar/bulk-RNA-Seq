###3)Alignment with STAR
#!/bin/bash
# -----------------------------
# Paths
# -----------------------------
GENOME_DIR="./GRCh38_STAR_INDEX"
FASTQ_DIR="./fastq"
OUT_DIR="./STAR_alignment"
CSV_FILE="$OUT_DIR/mapping_summary.csv"

# -----------------------------
# Check and create directories
# -----------------------------
if [ ! -d "$GENOME_DIR" ]; then
  echo "ERROR: Genome directory $GENOME_DIR does not exist!"
  exit 1
fi

if [ ! -d "$FASTQ_DIR" ]; then
  echo "ERROR: FASTQ directory $FASTQ_DIR does not exist!"
  exit 1
fi

mkdir -p $OUT_DIR

# -----------------------------
# Increase file descriptor limit
# -----------------------------
ulimit -n 10000

# Sample accessions are omitted from this public repository.
# Provide sample IDs in a local file before running the pipeline.

SAMPLE_FILE="sra_ids.txt"

while read -r sample
do
    # analysis commands
done < "$SAMPLE_FILE"

# -----------------------------
# Create CSV header
# -----------------------------
echo "Sample,Total_Reads,Uniquely_Mapped,%Uniquely_Mapped,%Multi_Mapped,Total_Reads_Flagstat,Mapped_Flagstat,%Mapped_Flagstat" > $CSV_FILE

# -----------------------------
# STAR alignment loop
# -----------------------------
for sample in "${SAMPLES[@]}"
do
  echo "Processing $sample ..."

  # Check if FASTQ files exist
  if [ ! -f "$FASTQ_DIR/${sample}_1.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${sample}_2.fastq.gz" ]; then
    echo "WARNING: FASTQ files for $sample not found, skipping."
    continue
  fi

  # Run STAR
  STAR --runThreadN 24 \
       --genomeDir $GENOME_DIR \
       --readFilesIn $FASTQ_DIR/${sample}_1.fastq.gz $FASTQ_DIR/${sample}_2.fastq.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix $OUT_DIR/${sample}_ \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMultimapNmax 20 \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNoverReadLmax 0.04 \
       --outFilterType BySJout \
       --outReadsUnmapped Fastx

  LOG_FILE="$OUT_DIR/${sample}_Log.final.out"
  BAM_FILE="$OUT_DIR/${sample}_Aligned.sortedByCoord.out.bam"

  if [ ! -f "$LOG_FILE" ]; then
    echo "ERROR: STAR failed for $sample, Log.final.out not found."
    continue
  fi

  # -----------------------------
  # Extract STAR mapping stats
  # -----------------------------
  TOTAL=$(grep "Number of input reads" $LOG_FILE | awk '{print $6}')
  UNIQUE=$(grep "Uniquely mapped reads number" $LOG_FILE | awk '{print $6}')
  PCT_UNIQUE=$(grep "Uniquely mapped reads %" $LOG_FILE | awk '{print $6}')
  MULTI=$(grep "Number of reads mapped to multiple loci" $LOG_FILE | awk '{print $8}')
  PCT_MULTI=$(grep "% of reads mapped to multiple loci" $LOG_FILE | awk '{print $1}')

  # -----------------------------
  # Run flagstat
  # -----------------------------
  FLAGSTAT_OUTPUT=$(samtools flagstat "$BAM_FILE")

  TOTAL_FLAG=$(echo "$FLAGSTAT_OUTPUT" | head -n 1 | awk '{print $1}')
  MAPPED_FLAG=$(echo "$FLAGSTAT_OUTPUT" | grep "mapped (" | head -n 1 | awk '{print $1}')
  PCT_MAPPED_FLAG=$(echo "$FLAGSTAT_OUTPUT" | grep "mapped (" | head -n 1 | awk '{print $5}' | tr -d '()%')

  # -----------------------------
  # Write to CSV
  # -----------------------------
  echo "$sample,$TOTAL,$UNIQUE,$PCT_UNIQUE,$PCT_MULTI,$TOTAL_FLAG,$MAPPED_FLAG,$PCT_MAPPED_FLAG" >> $CSV_FILE

done

echo "All samples processed. Summary saved to $CSV_FILE"