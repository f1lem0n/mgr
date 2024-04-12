mkdir -p output/RNAfold_0
for f in $(ls output/sequences/*.fasta); do
    RNAfold < $f > output/RNAfold_0/$(basename $f .fasta).out
done
mv *_ss.ps output/RNAfold_0/.