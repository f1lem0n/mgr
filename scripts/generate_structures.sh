BLAST_DB="refseq_select_rna"

root=$(pwd)
for type in mRNA tRNA rRNA; do
    for gene in $(ls output/structures/$type/); do
        echo "processing [ $gene ]"
        cd output/$gene
        mkdir -p \
            vienna/invivo \
            vienna/insilico \
            vienna/consensus_nn \
            vienna/consensus_aa \
            vienna/CT \
            mfold/invivo \
            mfold/insilico \
            mfold/consensus_nn \
            mfold/consensus_aa \
            mfold/CT

        # if [ ! -s seqdump.txt ]; then
        #     echo "No homologs found! Skipping..."
        #     cd $root
        #     continue
        # fi

        echo "Aligning homologs..."
        cat nucleotide.fasta seqdump.txt > homologs_nn.fasta
        clustalw \
            -INFILE="homologs_nn.fasta" \
            -OUTFILE="homologs_nn.aln"

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

        cd ../consensus_nn
        echo "Generating consensus structure..."
        RNAalifold --noLP -p -d2 < ../homologs_nn.aln > MFEs.txt
        b2ct < MFEs.txt > ../CT/consensus_nn.ct
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alidot.png alidot.ps > /dev/null
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alirna.png alirna.ps > /dev/null
        cd $root
    done
done

for gene in $(ls output/structures/mRNA/); do
    echo "processing [ $gene ]"
    cd output/$gene

    echo "Aligning homologs..."
    cat nucleotide.fasta seqdump.txt > homologs_nn.fasta
    clustalw \
        -INFILE="homologs_nn.fasta" \
        -OUTFILE="homologs_nn.aln"

    cd ../consensus_nn
    echo "Generating consensus structure..."
    RNAalifold --noLP -p -d2 < ../homologs_nn.aln > MFEs.txt
    b2ct < MFEs.txt > ../CT/consensus_nn.ct
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alidot.png alidot.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alirna.png alirna.ps > /dev/null
    cd $root
    done
exit 0
