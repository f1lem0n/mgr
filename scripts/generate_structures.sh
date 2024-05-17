BLAST_DB="refseq_select_rna"

root=$(pwd)

for gene in $(ls output/); do
    echo "processing [ $gene ]"
    cd output/$gene
    mkdir -p \
        vienna/invivo \
        vienna/insilico \
        vienna/consensus \
        vienna/CT \
        mfold/invivo \
        mfold/insilico \
        mfold/consensus \
        mfold/CT

    cd vienna/insilico
    echo "Generating in silico structure..."
    RNAfold --noLP -p -d2 -i ../seq.fasta > MFEs.txt
    b2ct < MFEs.txt > ../CT/insilico.ct
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=ss.png *_ss.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=dp.png *_dp.ps > /dev/null

    cd ../invivo
    echo "Generating in vivo structure..."
    RNAfold --noLP -p -d2 -C -i ../constrained.fasta > MFEs.txt
    b2ct < MFEs.txt > ../CT/invivo.ct
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=ss.png *_ss.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=dp.png *_dp.ps > /dev/null

    cd ..
    # echo "Searching for homologs..."
    # blastn \
    #     -query seq.fasta \
    #     -out blastx.xml \
    #     -db $BLAST_DB \
    #     -num_threads 12 \
    #     -outfmt "5" \
    #     -max_target_seqs 500 \
    #     -evalue 0.05
    # python $root/scripts/blast_xml2fasta.py blastx.xml > seqdump.txt
    if [ ! -s seqdump.txt ]; then
        echo "No homologs found! Skipping..."
        cd $root
        # rm -rf output/$gene
        continue
    fi
    echo "Aligning homologs..."
    cat seq.fasta seqdump.fasta > homologs.fasta
    clustalw \
        -INFILE="homologs.fasta" \
        -OUTFILE="homologs.aln"

    cd consensus
    echo "Generating consensus structure..."
    RNAalifold --noLP -p -d2 < ../homologs.aln > MFEs.txt
    b2ct < MFEs.txt > ../CT/consensus.ct
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alidot.png alidot.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alirna.png alirna.ps > /dev/null
    cd $root
done

exit 0
