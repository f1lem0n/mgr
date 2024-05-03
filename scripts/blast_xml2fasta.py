import pathlib
import xml.etree.ElementTree as ET


def find_query(tree):
    qid = tree.find(".//BlastOutput_query-def").text
    qseq = tree.find(".//Hsp_qseq").text
    return qid, qseq


def find_hits(tree):
    hids = tree.findall(".//Hit_def")
    hseqs = tree.findall(".//Hsp_hseq")
    return (
        [hid.text for hid in hids],
        [hseq.text.replace("-", "") for hseq in hseqs],
    )


def write_fasta(ids, seqs, outfile):
    with open(f"{outfile}.fasta", "w") as f:
        for i, s in zip(ids, seqs):
            f.write(f">{i}\n{s}\n")


def main():
    try:
        tree = ET.parse("blast_results.xml")
        qid, qseq = find_query(tree)
        hids, hseqs = find_hits(tree)
        ids = [qid] + hids
        seqs = [qseq] + hseqs
        write_fasta(ids, seqs, outfile="homologs.fasta")
    except AttributeError:
        print(f"[!] No hits found. Skipping...")


if __name__ == "__main__":
    main()
