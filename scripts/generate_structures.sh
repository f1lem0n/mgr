BLAST_DB="refseq_select_rna"

root=$(pwd)

for type in mRNA rRNA tRNA; do
    for gene in $(ls output/structures/$type/); do
        echo "Processing [ $gene ]"
        cd output/structures/$type/$gene
        mkdir -p \
            vienna/invivo \
            vienna/insilico \
            vienna/consensus \
            vienna/consensus_guided \
            vienna/CT \
            mfold/invivo \
            mfold/insilico \
            mfold/consensus \
            mfold/consensus_guided \
            mfold/CT

        cd vienna/insilico
        echo "Generating in silico structure..."
        RNAfold --noLP -p -d2 -i ../../nucleotide.fasta > MFEs.txt
        b2ct < MFEs.txt > ../CT/insilico.ct
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=ss.png *_ss.ps > /dev/null
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=dp.png *_dp.ps > /dev/null

        cd ../invivo
        echo "Generating in vivo structure..."
        RNAfold --noLP -p -d2 -C -i ../../constrained.fasta > MFEs.txt
        b2ct < MFEs.txt > ../CT/invivo.ct
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=ss.png *_ss.ps > /dev/null
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=dp.png *_dp.ps > /dev/null

        cd ../..

        if [ ! -s seqdump.txt ]; then
            echo "No homologs found! Skipping..."
            cd $root
            continue
        fi

        echo "Aligning homologs..."
        cat nucleotide.fasta seqdump.txt > homologs_nn.fasta
        clustalw \
            -INFILE="homologs_nn.fasta" \
            -OUTFILE="nn.aln"

        cd vienna/consensus
        echo "Generating consensus structure..."
        RNAalifold --noLP -p -d2 < ../../nn.aln > MFEs.txt
        b2ct < MFEs.txt > ../CT/consensus.ct
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alidot.png alidot.ps > /dev/null
        gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alirna.png alirna.ps > /dev/null
        cd $root
    done
done

for gene in $(ls output/structures/mRNA/); do
    echo "Processing [ $gene ]"
    cd output/structures/mRNA/$gene

    if [ ! -s hits.csv ]; then
        echo "No homologs found! Skipping..."
        cd $root
        continue
    fi

    echo "Fetching sequences..."
    python $root/scripts/hits2fasta.py

    echo "Aligning homologs..."
    cat protein.fasta seqdump_aa.txt > homologs_aa.fasta
    cat nucleotide.fasta pal2nal_guides.txt > pal2nal_guides.fasta
    rm -f pal2nal_guides.txt
    clustalw \
        -INFILE="homologs_aa.fasta" \
        -OUTFILE="homologs_aa.aln"
    pal2nal.pl homologs_aa.aln pal2nal_guides.fasta > guided.aln

    cd vienna/consensus_guided
    echo "Generating consensus structure..."
    RNAalifold --noLP -p -d2 < ../../guided.aln > MFEs.txt
    b2ct < MFEs.txt > ../CT/consensus_guided.ct
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alidot.png alidot.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alirna.png alirna.ps > /dev/null
    cd $root
done

exit 0
