###1)Download SRA samples and convert to FASTQ 
cat sra_ids.txt | xargs -n 1 -P 4 -I{} parallel-fastq-dump --sra-id {} --threads 20 --split-files --gzip --outdir fastq


#(if it needs to download and convert one sample:)
parallel-fastq-dump --sra-id SRA ID --threads 20 --gzip --outdir fastq

