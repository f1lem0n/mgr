from Bio import Entrez, SeqIO
from io import StringIO

with open("seqids.txt") as f:
    ids = f.read().splitlines()

aa_records = []
nn_records = []

for id in ids:
    with Entrez.efetch(db="protein", id=id, rettype="gb", retmode="text") as handle:
        gb = handle.read()
    aa_rec = SeqIO.read(StringIO(gb), "genbank").format("fasta")
    aa_rec = SeqIO.read(StringIO(aa_rec), "fasta")
    aa_records.append(aa_rec)
    cds_coords = gb.split("/coded_by=")[1].split("\n")[0]
    if "(" in cds_coords:
        cds_coords = cds_coords.split("(")[1].split(")")[0]
    seqid, coords = cds_coords.split(":")[0], cds_coords.split(":")[1]
    start, stop = coords.split("..")
    start, stop = int(start), int(stop)
    if start > stop:
        start, stop = stop, start
    start -= 1
    with Entrez.efetch(db="nucleotide", id=seqid, rettype="fasta", retmode="text") as handle:
        nn_rec = list(SeqIO.parse(handle, "fasta"))[0]
    nn_rec.seq = nn_rec.seq[start:stop]
    nn_records.append(nn_rec)

SeqIO.write(aa_records, "seqdump_aa.fasta", "fasta")
SeqIO.write(nn_records, "pal2nal_guides.fasta", "fasta")