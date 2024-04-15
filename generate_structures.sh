python scripts/split_query.py $1

mkdir -p output/structures_0
for f in $(ls output/sequences/*.fasta); do
    RNAfold < $f > output/structures_0/$(basename $f .fasta).out
done
mv *_ss.ps output/structures_0/.

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

python scripts/blast_xml2fasta.py

mkdir -p output/alignments
for f in $(ls output/blastn/*.fasta); do
    clustalw \
        -INFILE=$f \
        -OUTFILE="output/alignments/$(basename $f .fasta).aln"
done

mkdir -p output/structures_consensus
for f in $(ls output/alignments/*.aln); do
    RNAalifold < $f > output/structures_consensus/$(basename $f .fasta).out
done
mv *.ps output/structures_consensus/.