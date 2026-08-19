###4)Featurecounts
#Install app
conda create -n rnaseq -y
conda activate rnaseq
conda install -y -c bioconda -c conda-forge subread
which featureCounts
featureCounts -v
 

 #Run
#!/bin/bash
# Number of threads
THREADS=32

# Path to annotation GTF file
GTF_FILE="./home/user/GRCh38_STAR_INDEX/Homo_sapiens.GRCh38.115.gtf"

# Path to genome FASTA file
GENOME_FASTA="./GRCh38_STAR_INDEX/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Path to the folder containing BAM files
BAM_DIR="./STAR_alignment"

# Output file path
OUT_FILE="all_counts.txt"

# Find all BAM files
BAM_FILES=("$BAM_DIR"/*.bam)

# Run featureCounts
featureCounts -T $THREADS -p --countReadPairs -C -F GTF -a "$GTF_FILE" -G "$GENOME_FASTA" -o "$OUT_FILE" "${BAM_FILES[@]}"
