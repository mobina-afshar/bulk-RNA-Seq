#3)Alignment with STAR
#!/bin/bash

GENOME_DIR="./GRCh38_STAR_INDEX"
FASTQ_DIR="./fastq"
OUT_DIR="./STAR_alignment_single"
CSV_FILE="$OUT_DIR/mapping_summary.tsv"

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

# Sample accessions are omitted from this public repository.
# Provide sample IDs in a local file before running the pipeline.

SAMPLE_FILE="sra_ids.txt"

while read -r sample
do
    # analysis commands
done < "$SAMPLE_FILE"

# -----------------------------
# Create TSV header
# -----------------------------
echo -e "Sample\tTotal_Reads\tUniquely_Mapped\t%Uniquely_Mapped\t%Multi_Mapped\tFlagstat_Total\tFlagstat_Mapped\tFlagstat_Mapped%" > $CSV_FILE

# -----------------------------
# STAR + FLAGSTAT loop
# -----------------------------
for sample in "${SAMPLES[@]}"
do
  echo "Processing $sample ..."

  # Detect single-end fastq file
  if [ -f "$FASTQ_DIR/${sample}.fastq.gz" ]; then
    READ_FILE="$FASTQ_DIR/${sample}.fastq.gz"
  elif [ -f "$FASTQ_DIR/${sample}_1.fastq.gz" ]; then
    READ_FILE="$FASTQ_DIR/${sample}_1.fastq.gz"
  else
    echo "WARNING: FASTQ for $sample not found! Skipping."
    continue
  fi

  # Run STAR (SINGLE-END)
  STAR --runThreadN 24 \
       --genomeDir $GENOME_DIR \
       --readFilesIn $READ_FILE \
       --readFilesCommand zcat \
       --outFileNamePrefix $OUT_DIR/${sample}_ \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMultimapNmax 20 \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNoverReadLmax 0.04 \
       --outFilterType BySJout \
       --outReadsUnmapped Fastx

  LOG="$OUT_DIR/${sample}_Log.final.out"
  BAM="$OUT_DIR/${sample}_Aligned.sortedByCoord.out.bam"

  if [ ! -f "$LOG" ]; then
    echo "ERROR: STAR failed for $sample"
    continue
  fi

  # -----------------------------
  # Extract STAR mapping stats
  # -----------------------------
  TOTAL=$(grep "Number of input reads" "$LOG" | awk '{print $6}')
  UNIQUE=$(grep "Uniquely mapped reads number" "$LOG" | awk '{print $6}')
  PCT_UNIQUE=$(grep "Uniquely mapped reads %" "$LOG" | awk '{print $6}')
  PCT_MULTI=$(grep "% of reads mapped to multiple loci" "$LOG" | awk '{print $1}')

  # -----------------------------
  # samtools flagstat
  # -----------------------------
  FLAGSTAT=$(samtools flagstat "$BAM")

  FS_TOTAL=$(echo "$FLAGSTAT" | grep "in total" | awk '{print $1}')
  FS_MAPPED=$(echo "$FLAGSTAT" | grep "mapped (" | awk '{print $1}')
  FS_MAPPED_PCT=$(echo "$FLAGSTAT" | grep "mapped (" | awk '{print $5}' | tr -d '()%')

  # -----------------------------
  # Write to TSV
  # -----------------------------
  echo -e "${sample}\t${TOTAL}\t${UNIQUE}\t${PCT_UNIQUE}\t${PCT_MULTI}\t${FS_TOTAL}\t${FS_MAPPED}\t${FS_MAPPED_PCT}" >> $CSV_FILE

done

echo "DONE! Final summary saved to:"
echo "$CSV_FILE"
