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
    with open(f"output/blastn/{outfile}.fasta", "w") as f:
        for i, s in zip(ids, seqs):
            f.write(f">{i}\n{s}\n")


def main():
    for path in pathlib.Path("output/blastn").iterdir():
        if path.suffix != ".xml":
            continue
        try:
            tree = ET.parse(path)
            qid, qseq = find_query(tree)
            hids, hseqs = find_hits(tree)
            ids = [qid] + hids
            seqs = [qseq] + hseqs
            write_fasta(ids, seqs, outfile=path.stem)
        except AttributeError:
            print(f"[{path.stem}] No hits found.")
            continue


if __name__ == "__main__":
    main()
