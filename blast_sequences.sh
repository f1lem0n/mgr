mkdir -p output/blastn
for f in $(ls output/sequences/*.fasta); do
    blastn \
        -query $f \
        -db refseq_select_rna \
        -outfmt 5 \
        -evalue 1e-30 \
        -word_size 6 \
        -max_target_seqs 5 \
        -out output/blastn/$(basename $f .fasta).xml
done