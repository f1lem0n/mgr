BLAST_DB="refseq_select_rna"

root=$(pwd)

for gene in $(ls output/); do
    echo "processing [ $gene ]"
    cd output/$gene
    mkdir -p invivo insilico consensus

    cd insilico
    echo "generating in silico structure..."
    RNAfold --noLP -p -d2 -i ../seq.fasta > MFEs.txt
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=ss.png *_ss.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=dp.png *_dp.ps > /dev/null

    cd ../invivo
    echo "generating in vivo structure..."
    RNAfold --noLP -p -d2 -C -i ../constrained.fasta > MFEs.txt
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=ss.png *_ss.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=dp.png *_dp.ps > /dev/null

    cd ../consensus
    echo "finding homologs for consensus structure..."
    blastn \
        -query ../seq.fasta \
        -db $BLAST_DB \
        -outfmt 5 \
        -evalue 1e-30 \
        -word_size 6 \
        -max_target_seqs 5 \
        -out blast_results.xml
    python ../../../scripts/blast_xml2fasta.py
    if [ ! -f homologs.fasta ]; then
        cd $root
        continue
    fi
    echo "aligning homologs..."
    clustalw \
        -INFILE="homologs.fasta" \
        -OUTFILE="homologs.aln"
    echo "generating consensus structure..."
    RNAalifold --noLP -p -d2 < homologs.aln > MFEs.txt
    cd $root
done

exit 0






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