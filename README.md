# Protokół analizy danych z DMS-seq

## Struktura projektu

```
tree -d
.
├── output
│   ├── index
│   ├── peaks
│   ├── QC_premap
│   ├── bowtie2_alignment
│   └── reads_trimmed
├── data
│   ├── reads
│   └── S288C_reference_genome_R62-1-1_20090218
└── scripts
```

(*wyświetlono tylko foldery*)

## Środowisko

UWAGA! Wszystkie poniższe skrypty należy uruchamiać z poziomu roota projektu
z aktywnym środowiskiem `mgr`.

Aby stworzyć i aktywować środowisko należy uruchomić:
```bash
conda env create -f environment.yml
conda activate mgr
```



## Dane

Odczyty pozyskano z serwera FTP portalu [ENA](https://www.ebi.ac.uk/ena/browser/home)
poprzez uruchomienie poniższego skryptu:

```bash
mkdir -p data/reads
cd data/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR815/SRR815623/SRR815623.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR815/SRR815629/SRR815629.fastq.gz
```

Należy pobrać przynajmniej dwa pliki z odczytami **(a)** próbka kontrolna (DMS-)
i **(b)** próbka eksperymentalna (DMS+). W tym przypadku pobrano pliki
o numerach SRR815623 (DMS-) i SRR815629 (DMS+).

Genom referencyjny *Saccharomyces cerevisiae* pobrano z serwera
[SGD](https://www.yeastgenome.org/):

```bash
cd data
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R62-1-1_20090218.tgz
tar -xf S288C_reference_genome_R62-1-1_20090218.tgz
rm S288C_reference_genome_R62-1-1_20090218.tgz
```


## Kontrola jakości odczytów

Do kontroli jakości i przycinania sekwencji adapterowych użyto narzędzia `fastp`:

```bash
mkdir -p output/reads_trimmed
mkdir -p output/QC_premap
for f in $(ls data/reads/*.fastq.gz); do
    fastp -i $f \
        -o output/reads_trimmed/$(basename $f) \
        -h output/QC_premap/$(basename $f).html
done
```

## Indeksowanie genomu referencyjnego

Do budowy indeksu i późniejszego mapowania odczytów do genomu referencyjnego
użyto narzędzia `STAR`. Najpierw należy zbudować indeks:

STAR:

```bash
mkdir -p output/STAR_index
STAR \
    --runMode genomeGenerate \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles data/S288C_reference_genome_R62-1-1_20090218/S288C_reference_sequence_R62-1-1_20090218.fsa \
    --genomeDir output/STAR_index
```


## Mapowanie odczytów do genomu referencyjnego

Do mapowania odczytów do genomu referencyjnego użyto narzędzia `STAR`.
Następnie przy pomocy `samtools view` przekonwertowano pliki
`sam` na format `bam`.

STAR:

```bash
mkdir -p output/STAR_alignment

for f in $(ls output/reads_trimmed/*.fastq.gz); do
    STAR \
        --genomeDir output/STAR_index \
        --readFilesCommand zcat \
        --readFilesIn output/reads_trimmed/$f \
        --runThreadN 12 \
        --outFileNamePrefix output/STAR_alignment/$(basename $f .fastq.gz)_ \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax 100 \
        --seedSearchStartLmax 15 \
        -outSAMattributes All \

    samtools view -h -S -b \
        -o output/STAR_alignment/$(basename $f .fastq.gz).bam \
        output/STAR_alignment/$(basename $f .fastq.gz)_Aligned.out.sam

    rm output/STAR_alignment/$(basename $f .fastq.gz)_Aligned.out.sam
done
```

## Wykrywanie sygnału od DMS

Aby wykryć sygnał od DMS, najpierw użyto `rf-count` z pakietu `RNAFramework`
do zliczenia liczby zatrzymań odwrotnej transkryptazy na każdej pozycji w genomie:

```bash
scripts/RNAFramework/rf-count \
    output/STAR_alignment/SRR81562*.bam \
    -f data/S288C_reference_genome_R62-1-1_20090218/S288C_reference_sequence_R62-1-1_20090218.fsa \
    -o output/RTS_counts/
```
