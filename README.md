# Protokół analizy danych z DMS-seq

## Struktura projektu

```
tree -d
.
├── output
│   ├── index
│   ├── peaks
│   ├── QC_premap
│   ├── genome_aln
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
użyto narzędzia `bowtie2`. Najpierw należy zbudować indeks:

```bash
bowtie2-build \
    data/S288C_reference_genome_R62-1-1_20090218/S288C_reference_sequence_R62-1-1_20090218.fsa \
    output/index/S288C
```

## Mapowanie odczytów do genomu referencyjnego

Do mapowania odczytów do genomu referencyjnego użyto narzędzia `bowtie2`.
Następnie przy pomocy `samtools view` przekonwertowano pliki
`sam` na format `bam`.

```bash
mkdir -p output/genome_aln
for f in $(ls output/reads_trimmed/*.fastq.gz); do
    bowtie2 \
	-p 4 --no-mixed --no-discordant --very-sensitive \
        -x output/index/S228C \
        -U $f \
        -S output/genome_aln/$(basename $f .fastq.gz).sam

    samtools view -h -S -b \
        -o output/genome_aln/$(basename $f .fastq.gz).bam \
        output/genome_aln/$(basename $f .fastq.gz).sam

    rm output/genome_aln/$(basename $f .fastq.gz).sam
done
```

## Wykrywanie sygnału od DMS

Do wykrywania sygnału od DMS użyto narzędzia `MACS2` (analogicznie
jak w przypadku protokołów ChIPseq). Trzeba mieć pliki `bam` dla
próbki kontrolnej i eksperymentalnej oraz oszacować wielkość genomu,
do którego mapowane były odczyty. W tym przypadku wielkość genomu
*Saccharomyces cerevisiae* wynosi 12359296 bp. Nie zastosowano modelu
fragmentacji (`--nomodel`), ponieważ `MACS2` nie znajduje odpowiedniej
liczby sparowanych sygnałów.

```bash
macs2 callpeak \
    -t output/genome_aln/SRR815629.bam \
    -c output/genome_aln/SRR815623.bam \
    --outdir output/peaks \
    -g 12359296 \
    --nomodel \
    -n DMS \
```
