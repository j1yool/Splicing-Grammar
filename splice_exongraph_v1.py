"""
splice_exongraph_v1.py  (Day 97, Splicing Grammar, Supporting A)
Data preparation only: GTF -> exon table, transcript exon-path table, exon-graph edge table.
No inference. Nothing is pre-registered (exploratory data prep, R6).

Usage:
  python splice_exongraph_v1.py --gtf PATH --genes FGFR1,FGFR2 [--outdir .] [--show FGFR2] [--overwrite]

Hard stops (nothing written): gtf missing; gene absent or ambiguous; zero exon rows;
mixed chromosome or strand inside a gene; overlapping exons inside a transcript;
a transcript assigned to two genes; reconstruction or edge-endpoint mismatch; output exists.
Exon IDs are unique (start,end) intervals per gene, numbered 5'->3' in transcript direction.
"""
import argparse, csv, gzip, hashlib, json, os, re, sys
from collections import defaultdict, Counter

ATTR = re.compile(r'(\S+) "([^"]*)"')


def stop(msg):
    print("STOP: " + msg)
    sys.exit(1)


def opener(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def order_key(strand):
    if strand == "+":
        return lambda e: (e[0], e[1])
    return lambda e: (-e[1], -e[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--show", default=None)
    ap.add_argument("--show-limit", type=int, default=40)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    if not os.path.isfile(a.gtf):
        stop("gtf not found: " + a.gtf)
    if not os.path.isdir(a.outdir):
        stop("outdir not found: " + a.outdir)
    genes = [g.strip() for g in a.genes.split(",") if g.strip()]
    if not genes or len(set(genes)) != len(genes):
        stop("--genes empty or has duplicates")
    if a.show and a.show not in genes:
        stop("--show gene not in --genes")
    needles = {g: 'gene_name "%s"' % g for g in genes}

    rows = []
    with opener(a.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            hit = None
            for g, nd in needles.items():
                if nd in line:
                    hit = g
                    break
            if hit is None:
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "exon":
                continue
            attrs = defaultdict(list)
            for k, v in ATTR.findall(cols[8]):
                attrs[k].append(v)
            if attrs.get("gene_name", [None])[0] != hit:
                continue
            tx = attrs.get("transcript_id", [None])[0]
            gid = attrs.get("gene_id", [None])[0]
            if tx is None or gid is None:
                stop("exon row without transcript_id/gene_id for " + hit)
            bt = (attrs.get("transcript_type") or attrs.get("transcript_biotype") or ["NA"])[0]
            rows.append({"gene": hit, "gene_id": gid, "tx": tx, "biotype": bt,
                         "tags": attrs.get("tag", []), "chrom": cols[0],
                         "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6]})
    if not rows:
        stop("zero exon rows matched any requested gene; nothing written")

    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    for g in genes:
        if g not in by_gene:
            stop("gene absent from GTF exon rows: " + g)
        ids = set(r["gene_id"] for r in by_gene[g])
        if len(ids) != 1:
            stop("gene_name %s maps to %d gene_ids: %s" % (g, len(ids), sorted(ids)))
        if len(set(r["chrom"] for r in by_gene[g])) != 1:
            stop("mixed chromosomes in " + g)
        if len(set(r["strand"] for r in by_gene[g])) != 1:
            stop("mixed strands in " + g)

    tx_gene = {}
    for r in rows:
        if tx_gene.setdefault(r["tx"], r["gene"]) != r["gene"]:
            stop("transcript %s assigned to two genes" % r["tx"])

    exon_rows, tx_out, edge_out, summary = [], [], [], {}
    for g in genes:
        rs = by_gene[g]
        strand = rs[0]["strand"]
        chrom = rs[0]["chrom"]
        key = order_key(strand)
        uniq = sorted(set((r["start"], r["end"]) for r in rs), key=key)
        w = max(2, len(str(len(uniq))))
        eid = {iv: "E%0*d" % (w, i + 1) for i, iv in enumerate(uniq)}

        txs = defaultdict(list)
        for r in rs:
            txs[r["tx"]].append(r)
        paths, pc_paths = {}, {}
        n_ex_rows = 0
        exon_tx = Counter()
        exon_pc = Counter()
        edge_all, edge_pc = Counter(), Counter()
        for t in sorted(txs):
            ex = sorted(txs[t], key=lambda r: key((r["start"], r["end"])))
            for i in range(len(ex) - 1):
                p, q = ex[i], ex[i + 1]
                if strand == "+" and not q["start"] > p["end"]:
                    stop("overlapping exons inside transcript " + t)
                if strand == "-" and not q["end"] < p["start"]:
                    stop("overlapping exons inside transcript " + t)
            if len(set((r["start"], r["end"]) for r in ex)) != len(ex):
                stop("duplicate exon interval inside transcript " + t)
            path = [eid[(r["start"], r["end"])] for r in ex]
            if path != sorted(path):
                stop("path not monotone in transcript " + t)
            n_ex_rows += len(path)
            bt = ex[0]["biotype"]
            ispc = (bt == "protein_coding")
            tags = sorted(set(x for r in ex for x in r["tags"]))
            paths[t] = path
            if ispc:
                pc_paths[t] = path
            for e in path:
                exon_tx[e] += 1
                if ispc:
                    exon_pc[e] += 1
            chain = ["START"] + path + ["END"]
            for i in range(len(chain) - 1):
                edge_all[(chain[i], chain[i + 1])] += 1
                if ispc:
                    edge_pc[(chain[i], chain[i + 1])] += 1
            tx_out.append([g, t, bt, ";".join(tags), len(path), path[0], path[-1], ">".join(path)])

        if n_ex_rows != len(rs):
            stop("reconstruction mismatch (exon rows) in " + g)
        if sum(exon_tx.values()) != len(rs):
            stop("reconstruction mismatch (exon counts) in " + g)
        if sum(c for (u, v), c in edge_all.items() if u == "START") != len(txs):
            stop("START edge count != transcript count in " + g)
        valid = set(eid.values()) | {"START", "END"}
        for (u, v) in edge_all:
            if u not in valid or v not in valid:
                stop("edge endpoint not an exon id in " + g)

        n_ov = 0
        for iv in uniq:
            ov = [eid[o] for o in uniq if o != iv and o[0] <= iv[1] and iv[0] <= o[1]]
            n_ov += len(ov)
            exon_rows.append([g, eid[iv], chrom, iv[0], iv[1], strand, iv[1] - iv[0] + 1,
                              exon_tx[eid[iv]], exon_pc[eid[iv]], ";".join(ov)])
        for (u, v) in sorted(edge_all, key=lambda e: (e[0] != "START", e[0], e[1])):
            edge_out.append([g, u, v, edge_all[(u, v)], edge_pc.get((u, v), 0)])

        pc_set = Counter(">".join(p) for p in pc_paths.values())
        all_set = Counter(">".join(p) for p in paths.values())
        n_pc = len(pc_paths)
        var_pc = sum(1 for e in eid.values() if 0 < exon_pc[e] < n_pc)
        summary[g] = {"gene_id": rs[0]["gene_id"], "chrom": chrom, "strand": strand,
                      "n_transcripts": len(txs), "n_protein_coding_transcripts": n_pc,
                      "n_distinct_paths": len(all_set), "n_distinct_paths_pc": len(pc_set),
                      "n_unique_exons": len(uniq), "n_edges": len(edge_all),
                      "n_overlapping_exon_pairs": n_ov // 2, "n_variable_exons_pc": var_pc,
                      "_pc_paths": pc_set}

    outs = {"exons": "splice_exons_v1.tsv", "transcripts": "splice_transcripts_v1.tsv",
            "edges": "splice_edges_v1.tsv", "summary": "splice_exongraph_v1_summary.json"}
    for fn in outs.values():
        if os.path.exists(os.path.join(a.outdir, fn)) and not a.overwrite:
            stop("output exists (use --overwrite): " + fn)

    def wtsv(fn, header, data):
        with open(os.path.join(a.outdir, fn), "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
            wr.writerow(header)
            wr.writerows(data)

    wtsv(outs["exons"], ["gene", "exon_id", "chrom", "start", "end", "strand", "length",
                         "n_transcripts", "n_protein_coding", "overlaps"], exon_rows)
    wtsv(outs["transcripts"], ["gene", "transcript_id", "biotype", "tags", "n_exons",
                               "first_exon", "last_exon", "path"], tx_out)
    wtsv(outs["edges"], ["gene", "from", "to", "n_transcripts", "n_protein_coding"], edge_out)

    me = open(__file__, "rb").read()
    js = {"gtf": os.path.abspath(a.gtf), "gtf_sha256": sha256_file(a.gtf),
          "script_sha256_raw": sha256_bytes(me),
          "script_sha256_lf": sha256_bytes(me.replace(b"\r\n", b"\n")),
          "genes": {g: {k: v for k, v in s.items() if not k.startswith("_")} for g, s in summary.items()}}
    with open(os.path.join(a.outdir, outs["summary"]), "w", encoding="utf-8") as fh:
        json.dump(js, fh, indent=2)

    print("gtf_sha256 %s" % js["gtf_sha256"])
    print("script_sha256_raw %s  script_sha256_lf %s" % (js["script_sha256_raw"], js["script_sha256_lf"]))
    print("%-8s %-3s %6s %6s %8s %8s %6s %6s %6s %6s" % ("gene", "str", "n_tx", "n_pc", "paths", "paths_pc",
                                                       "exons", "edges", "ovl", "var_pc"))
    for g in genes:
        s = summary[g]
        print("%-8s %-3s %6d %6d %8d %8d %6d %6d %6d %6d" % (
            g, s["strand"], s["n_transcripts"], s["n_protein_coding_transcripts"], s["n_distinct_paths"],
            s["n_distinct_paths_pc"], s["n_unique_exons"], s["n_edges"], s["n_overlapping_exon_pairs"],
            s["n_variable_exons_pc"]))
    for fn in outs.values():
        p = os.path.join(a.outdir, fn)
        print("WROTE %s  bytes=%d" % (p, os.path.getsize(p)))

    if a.show:
        print("\n== exon table: %s ==" % a.show)
        print("exon_id\tchrom\tstart\tend\tlen\tn_tx\tn_pc\toverlaps")
        for r in exon_rows:
            if r[0] == a.show:
                print("\t".join(str(x) for x in [r[1], r[2], r[3], r[4], r[6], r[7], r[8], r[9]]))
        print("\n== distinct protein-coding paths: %s (count, path) ==" % a.show)
        items = summary[a.show]["_pc_paths"].most_common()
        for pth, n in items[:a.show_limit]:
            print("%3d  %s" % (n, pth))
        if len(items) > a.show_limit:
            print("... %d more distinct paths not shown" % (len(items) - a.show_limit))


if __name__ == "__main__":
    main()