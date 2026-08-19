#ّ##2)FastQC 
OUTPUT_DIR=./fastqc_results
THREADS=4

mkdir -p $OUTPUT_DIR

for fq in ./*_1.fastq.gz ./*_2.fastq.gz
do
    echo "Running FastQC on $fq ..."
    fastqc -o $OUTPUT_DIR -t $THREADS $fq
done

echo "All FastQC analyses are complete!"
