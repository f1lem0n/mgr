mkdir -p output/reads_filtered
mkdir -p output/QC_premap
for f in $(ls data/reads/*.fastq.gz); do
    fastp -i $f \
        -o output/reads_filtered/$(basename $f) \
        -h output/QC_premap/$(basename $f).html
done

mkdir -p output/genome_aln
for f in $(ls output/reads_filtered/*.fastq.gz); do
    bowtie2 -p 2 -q --local \
        -x output/index/S228C \
        -U $f \
        -S output/genome_aln/$(basename $f .fastq.gz).sam
    samtools view -h -S -b \
        -o output/genome_aln/$(basename $f .fastq.gz).bam \
        output/genome_aln/$(basename $f .fastq.gz).sam
    samtools sort -o output/genome_aln/$(basename $f .fastq.gz).sorted.bam \
        output/genome_aln/$(basename $f .fastq.gz).bam
    samtools index output/genome_aln/$(basename $f .fastq.gz).sorted.bam
    rm output/genome_aln/$(basename $f .fastq.gz).sam
done
