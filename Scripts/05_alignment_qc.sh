#!/bin/bash
# -----------------------------
# Paths
# -----------------------------
OUT_DIR="./STAR_alignment"
CSV_FILE="$OUT_DIR/mapping_summary.csv"

# Sample accessions are omitted from this public repository.
# Provide sample IDs in a local file before running the pipeline.

SAMPLES="sra_ids.txt"

while read -r sample
do
    # analysis commands
done < "$SAMPLES"

# -----------------------------
# Write CSV header
# -----------------------------
echo "Sample,TOTAL,UNIQUE,PCT_UNIQUE,PCT_MULTI,TOTAL_FLAG,MAPPED_FLAG,PCT_MAPPED_FLAG" > "$CSV_FILE"

# -----------------------------
# Loop over samples
# -----------------------------
for sample in "${SAMPLES[@]}"
do
  # Note: Ensure these filenames match your actual output. 
  # Sometimes STAR outputs "Log.final.out" without the prefix depending on how you ran it.
  LOG_FILE="$OUT_DIR/${sample}_Log.final.out"
  BAM_FILE="$OUT_DIR/${sample}_Aligned.sortedByCoord.out.bam"

  if [ ! -f "$LOG_FILE" ]; then
    echo "WARNING: Log file for $sample not found, skipping."
    continue
  fi

  if [ ! -f "$BAM_FILE" ]; then
    echo "WARNING: BAM file for $sample not found, skipping."
    continue
  fi

  echo "Processing $sample..."

  # -----------------------------
  # Extract STAR mapping stats (FIXED)
  # -----------------------------
  # We use awk -F "|" to split by the pipe, getting the 2nd part, 
  # and tr -d ' \t' to remove whitespace.
  
  TOTAL=$(grep "Number of input reads" "$LOG_FILE" | awk -F "|" '{print $2}' | tr -d ' \t')
  UNIQUE=$(grep "Uniquely mapped reads number" "$LOG_FILE" | awk -F "|" '{print $2}' | tr -d ' \t')
  PCT_UNIQUE=$(grep "Uniquely mapped reads %" "$LOG_FILE" | awk -F "|" '{print $2}' | tr -d ' \t')
  MULTI=$(grep "Number of reads mapped to multiple loci" "$LOG_FILE" | head -n 1 | awk -F "|" '{print $2}' | tr -d ' \t')
  PCT_MULTI=$(grep "% of reads mapped to multiple loci" "$LOG_FILE" | head -n 1 | awk -F "|" '{print $2}' | tr -d ' \t')

  # -----------------------------
  # Run flagstat
  # -----------------------------
  FLAGSTAT_OUTPUT=$(samtools flagstat "$BAM_FILE")
  
  # Extracting values from flagstat
  TOTAL_FLAG=$(echo "$FLAGSTAT_OUTPUT" | head -n 1 | awk '{print $1}')
  MAPPED_FLAG=$(echo "$FLAGSTAT_OUTPUT" | grep "mapped (" | head -n 1 | awk '{print $1}')
  
  # Standard flagstat output format: "123 + 0 mapped (95.43% : N/A)"
  # awk $5 gets "(95.43%", then we delete ( ) and %
  PCT_MAPPED_FLAG=$(echo "$FLAGSTAT_OUTPUT" | grep "mapped (" | head -n 1 | awk '{print $5}' | tr -d '()%')

  # -----------------------------
  # Write to CSV
  # -----------------------------
  echo "$sample,$TOTAL,$UNIQUE,$PCT_UNIQUE,$PCT_MULTI,$TOTAL_FLAG,$MAPPED_FLAG,$PCT_MAPPED_FLAG" >> "$CSV_FILE"

done

echo "All samples processed. Summary saved to $CSV_FILE"
