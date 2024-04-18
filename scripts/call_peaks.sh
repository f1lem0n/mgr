mkdir output/peaks -p
macs2 callpeak \
    -t output/genome_aln/SRR815629.bam \
    -c output/genome_aln/SRR815623.bam \
    -n macs2.out \
    --outdir output/peaks \
    -g 12359296 \
    --nomodel \
    --extsize 300
