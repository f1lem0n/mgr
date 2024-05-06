BLAST_DB="refseq_select_rna"

root=$(pwd)

for gene in $(ls output/); do
    echo "processing [ $gene ]"
    cd output/$gene
    mkdir -p invivo insilico consensus CT

    cd insilico
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
    if [ ! -f seqdump.txt ]; then
        echo "No homologs found! Skipping..."
        cd $root
        continue
    fi
    echo "aligning homologs..."
    clustalw \
        -INFILE="seqdump.txt" \
        -OUTFILE="homologs.aln"
    echo "Generating consensus structure..."
    cd consensus
    RNAalifold --noLP -p -d2 < ../homologs.aln > MFEs.txt
    b2ct < MFEs.txt > ../CT/consensus.ct
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alidot.png alidot.ps > /dev/null
    gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -sOutputFile=alirna.png alirna.ps > /dev/null
    cd $root
done

exit 0
