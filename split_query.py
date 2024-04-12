from sys import argv
import pathlib

from Bio import SeqIO


def main():
    pathlib.Path("output/sequences").mkdir(exist_ok=True, parents=True)
    for seq in SeqIO.parse(argv[1], "fasta"):
        with open(f"output/sequences/{seq.id}.fasta", "w") as handle:
            SeqIO.write(seq, handle, "fasta")


if __name__ == "__main__":
    main()
