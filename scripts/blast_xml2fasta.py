from sys import argv
import xml.etree.ElementTree as ET


def find_query(tree):
    qid = tree.find(".//BlastOutput_query-def").text
    qseq = tree.find(".//Hsp_qseq").text
    return qid, qseq


def find_hits(tree):
    hids = tree.findall(".//Hit_def")
    hseqs = tree.findall(".//Hsp_hseq")
    hdefs = tree.findall(".//Hit_def")
    return (
        [hid.text for hid in hids],
        [hseq.text.replace("-", "") for hseq in hseqs],
        [hdef.text for hdef in hdefs]
    )


def filter_hits(ids, seqs, defs, num_hits=10):
    n = 0
    fasta = ""
    for i, s, d in zip(ids, seqs, defs):
        if n == num_hits:
            break
        if "query" not in d and "Saccharomyces cerevisiae" not in d:
            fasta += f">seq{n}\n{s}\n"
            n += 1
    return fasta


def main():
    try:
        tree = ET.parse(str(argv[1]))
        qid, qseq = find_query(tree)
        hids, hseqs, hdefs = find_hits(tree)
        ids = [qid] + hids
        seqs = [qseq] + hseqs
        defs = ["query"] + hdefs
        fasta = filter_hits(ids, seqs, defs)
        print(fasta)
    except ET.ParseError:
        pass


if __name__ == "__main__":
    main()
