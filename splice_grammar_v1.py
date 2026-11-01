"""
splice_grammar_v1.py  (Day 98, Splicing Grammar, Supporting A)
Rebuilds a probabilistic grammar for ONE gene from protein-coding edge counts.
Descriptive reconstruction only: no inference, no variant data, nothing pre-registered (R6).

Grammar (right-linear PCFG over an acyclic exon graph):
  terminals     = exon ids (E03, E07, ...)
  nonterminals  = N[START] (start symbol) and N[exon] for every exon with out-edges
  productions   = N[u] -> v N[v]  with P = count(u,v) / out_count(u)     (v an exon)
                  N[u] -> EPS     with P = count(u,END) / out_count(u)
Inputs (schemas fixed by splice_exongraph_v1.py):
  splice_edges_v1.tsv        gene from to n_transcripts n_protein_coding
  splice_transcripts_v1.tsv  gene transcript_id biotype tags n_exons first_exon last_exon path
Usage:
  python splice_grammar_v1.py --selftest
  python splice_grammar_v1.py --edges splice_edges_v1.tsv --transcripts splice_transcripts_v1.tsv \
      --gene SRSF6 --expect-paths 18 --expect-observed 7 --expect-overgen 11 [--outdir .] [--overwrite]

Hard stops (nothing written):
  input missing / gene absent / bad schema or counts; no protein-coding edges or transcripts;
  edge counts do not equal the edge counts rebuilt from the observed protein-coding paths;
  flow not conserved at an exon; cycle in graph;
  H1 probabilities out of any state do not sum to exactly 1;
  H2 probabilities over all paths do not sum to 1 (exact) or float sum off by more than 1e-9;
  H3 an observed path does not derive from the grammar;
  H4/H5 |L(G)|, observed distinct paths, or overgeneration differ from the --expect values;
  an output file already exists (use --overwrite).
"""
import argparse, csv, hashlib, json, os, sys
from collections import Counter, defaultdict
from fractions import Fraction


class GrammarError(Exception):
    pass


def stop(msg):
    raise GrammarError(msg)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_edges(path, gene):
    if not os.path.isfile(path):
        stop("edges file not found: " + path)
    need = ["gene", "from", "to", "n_transcripts", "n_protein_coding"]
    edges = {}
    seen_gene = False
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames != need:
            stop("edges header mismatch: %s" % rd.fieldnames)
        for r in rd:
            if r["gene"] != gene:
                continue
            seen_gene = True
            try:
                npc = int(r["n_protein_coding"])
            except ValueError:
                stop("non-integer n_protein_coding for %s>%s" % (r["from"], r["to"]))
            if npc < 0:
                stop("negative count")
            key = (r["from"], r["to"])
            if key in edges:
                stop("duplicate edge row %s>%s" % key)
            if npc > 0:
                edges[key] = npc
    if not seen_gene:
        stop("gene absent from edges file: " + gene)
    if not edges:
        stop("no protein-coding edges for " + gene)
    return edges


def read_observed(path, gene):
    if not os.path.isfile(path):
        stop("transcripts file not found: " + path)
    need = ["gene", "transcript_id", "biotype", "tags", "n_exons", "first_exon", "last_exon", "path"]
    obs, ids = [], set()
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames != need:
            stop("transcripts header mismatch: %s" % rd.fieldnames)
        for r in rd:
            if r["gene"] != gene or r["biotype"] != "protein_coding":
                continue
            if r["transcript_id"] in ids:
                stop("duplicate transcript id " + r["transcript_id"])
            ids.add(r["transcript_id"])
            p = tuple(x for x in r["path"].split(">") if x)
            if not p or len(p) != int(r["n_exons"]):
                stop("path/n_exons mismatch in " + r["transcript_id"])
            obs.append(p)
    if not obs:
        stop("no protein-coding transcripts for " + gene)
    return obs


def build_grammar(edges, observed):
    """edges: {(u,v): count>0}; observed: list of exon-id tuples, one per protein-coding transcript."""
    n_tx = len(observed)
    rebuilt = Counter()
    for p in observed:
        chain = ("START",) + tuple(p) + ("END",)
        for i in range(len(chain) - 1):
            rebuilt[(chain[i], chain[i + 1])] += 1
    if dict(rebuilt) != dict(edges):
        stop("edge counts differ from edge counts rebuilt from observed paths")

    out_c, in_c = Counter(), Counter()
    for (u, v), c in edges.items():
        out_c[u] += c
        in_c[v] += c
    if out_c["START"] != n_tx or in_c["END"] != n_tx:
        stop("START out-count or END in-count != number of transcripts")
    if "END" in out_c or "START" in in_c:
        stop("edge out of END or into START")
    for n in set(out_c) | set(in_c):
        if n in ("START", "END"):
            continue
        if in_c[n] != out_c[n]:
            stop("flow not conserved at %s (in %d, out %d)" % (n, in_c[n], out_c[n]))

    succ = defaultdict(list)
    for (u, v), c in sorted(edges.items()):
        succ[u].append((v, c))
    prods = []
    for u in sorted(succ, key=lambda x: (x != "START", x)):
        tot = Fraction(0)
        for v, c in succ[u]:
            p = Fraction(c, out_c[u])
            tot += p
            prods.append({"lhs": "N[%s]" % u, "rhs": ("EPS" if v == "END" else "%s N[%s]" % (v, v)),
                          "u": u, "v": v, "count": c, "out_count": out_c[u], "prob": p})
        if tot != 1:
            stop("H1 failed: probabilities out of %s sum to %s" % (u, tot))

    paths = []
    onstack = set()

    def walk(u, acc, p):
        if u == "END":
            paths.append((tuple(acc), p))
            return
        if u in onstack:
            stop("cycle at " + u)
        if u not in succ:
            stop("dead end at " + u)
        onstack.add(u)
        for v, c in succ[u]:
            walk(v, acc if v == "END" else acc + [v], p * Fraction(c, out_c[u]))
        onstack.discard(u)

    walk("START", [], Fraction(1))
    tot = sum((p for _, p in paths), Fraction(0))
    if tot != 1 or abs(sum(float(p) for _, p in paths) - 1.0) > 1e-9:
        stop("H2 failed: path probabilities sum to %s" % tot)

    lang = set(p for p, _ in paths)
    if len(lang) != len(paths):
        stop("duplicate paths in L(G)")
    emap = {(pr["u"], pr["v"]) for pr in prods}
    for p in observed:
        chain = ("START",) + tuple(p) + ("END",)
        if not all((chain[i], chain[i + 1]) in emap for i in range(len(chain) - 1)) or tuple(p) not in lang:
            stop("H3 failed: observed path does not derive: " + ">".join(p))

    obs_cnt = Counter(observed)
    return {"prods": prods, "paths": paths, "n_tx": n_tx, "obs_cnt": obs_cnt,
            "n_L": len(paths), "n_obs_distinct": len(obs_cnt),
            "overgen": len(paths) - len(obs_cnt)}


def check_expect(g, e_paths, e_obs, e_over):
    for name, got, exp in (("|L(G)|", g["n_L"], e_paths), ("observed distinct", g["n_obs_distinct"], e_obs),
                           ("overgeneration", g["overgen"], e_over)):
        if exp is not None and got != exp:
            stop("%s = %d, expected %d" % (name, got, exp))


def selftest():
    def expect_fail(edges, obs, why):
        try:
            build_grammar(edges, obs)
        except GrammarError:
            return
        raise AssertionError("expected STOP: " + why)

    # Toy 1: unequal branch counts; probabilities known by hand.
    obs1 = [("E1", "E3"), ("E1", "E3"), ("E1", "E4"), ("E2", "E3")]
    ed1 = {("START", "E1"): 3, ("START", "E2"): 1, ("E1", "E3"): 2, ("E1", "E4"): 1, ("E2", "E3"): 1,
           ("E3", "END"): 3, ("E4", "END"): 1}
    g = build_grammar(ed1, obs1)
    pm = {p: pr for p, pr in g["paths"]}
    assert g["n_L"] == 3 and g["n_obs_distinct"] == 3 and g["overgen"] == 0
    assert pm[("E1", "E3")] == Fraction(1, 2) and pm[("E1", "E4")] == Fraction(1, 4) and pm[("E2", "E3")] == Fraction(1, 4)
    # Toy 2: shared bottleneck M (same shape as SRSF6 E14>E15>E16): 2 observed, 4 generated.
    obs2 = [("A", "M", "C"), ("B", "M", "D")]
    ed2 = {("START", "A"): 1, ("START", "B"): 1, ("A", "M"): 1, ("B", "M"): 1, ("M", "C"): 1, ("M", "D"): 1,
           ("C", "END"): 1, ("D", "END"): 1}
    g = build_grammar(ed2, obs2)
    assert g["n_L"] == 4 and g["n_obs_distinct"] == 2 and g["overgen"] == 2
    assert all(pr == Fraction(1, 4) for _, pr in g["paths"])
    check_expect(g, 4, 2, 2)
    # Negative tests: each hard stop must fire.
    expect_fail(ed2, [("A", "M", "C")], "edge counts != rebuilt counts")
    bad = dict(ed2); bad[("M", "C")] = 2
    expect_fail(bad, obs2, "flow / rebuilt mismatch")
    cyc_obs = [("A", "B", "A")]
    expect_fail({("START", "A"): 1, ("A", "B"): 1, ("B", "A"): 1, ("A", "END"): 1}, cyc_obs, "cycle or flow")
    try:
        check_expect(g, 5, None, None)
    except GrammarError:
        pass
    else:
        raise AssertionError("expected STOP on wrong |L|")
    print("SELFTEST OK: toy1 (3 paths, probs 1/2,1/4,1/4), toy2 (4 paths, overgen 2), 4 negative stops fired")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--edges")
    ap.add_argument("--transcripts")
    ap.add_argument("--gene", default="SRSF6")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--expect-paths", type=int)
    ap.add_argument("--expect-observed", type=int)
    ap.add_argument("--expect-overgen", type=int)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        selftest()
        return
    if not a.edges or not a.transcripts:
        stop("--edges and --transcripts required (or --selftest)")
    if not os.path.isdir(a.outdir):
        stop("outdir not found: " + a.outdir)
    outs = {"prods": "splice_grammar_v1_productions.tsv", "paths": "splice_grammar_v1_paths.tsv",
            "json": "splice_grammar_v1_results.json"}
    for fn in outs.values():
        if os.path.exists(os.path.join(a.outdir, fn)) and not a.overwrite:
            stop("output exists (use --overwrite): " + fn)

    edges = read_edges(a.edges, a.gene)
    observed = read_observed(a.transcripts, a.gene)
    g = build_grammar(edges, observed)
    check_expect(g, a.expect_paths, a.expect_observed, a.expect_overgen)

    def f12(x):
        return "%.12g" % float(x)

    with open(os.path.join(a.outdir, outs["prods"]), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "lhs", "rhs", "count", "out_count", "probability", "fraction"])
        for pr in g["prods"]:
            w.writerow([a.gene, pr["lhs"], pr["rhs"], pr["count"], pr["out_count"], f12(pr["prob"]), str(pr["prob"])])
    ranked = sorted(g["paths"], key=lambda t: (-t[1], t[0]))
    with open(os.path.join(a.outdir, outs["paths"]), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "path_id", "path", "n_exons", "probability", "fraction", "observed",
                    "n_observed_transcripts", "expected_transcripts"])
        for i, (p, pr) in enumerate(ranked, 1):
            w.writerow([a.gene, "P%02d" % i, ">".join(p), len(p), f12(pr), str(pr),
                        1 if p in g["obs_cnt"] else 0, g["obs_cnt"].get(p, 0), "%.4f" % (float(pr) * g["n_tx"])])

    me = open(__file__, "rb").read()
    res = {"gene": a.gene, "edges_file": os.path.abspath(a.edges), "edges_sha256": sha256_file(a.edges),
           "transcripts_file": os.path.abspath(a.transcripts), "transcripts_sha256": sha256_file(a.transcripts),
           "script_sha256_raw": sha256_bytes(me), "script_sha256_lf": sha256_bytes(me.replace(b"\r\n", b"\n")),
           "n_protein_coding_transcripts": g["n_tx"], "n_productions": len(g["prods"]),
           "n_paths_L": g["n_L"], "n_observed_distinct": g["n_obs_distinct"], "overgeneration": g["overgen"],
           "path_prob_sum_float": sum(float(p) for _, p in g["paths"]),
           "hard_stops_passed": ["flow", "H1", "H2", "H3", "H4", "H5"] if None not in (
               a.expect_paths, a.expect_observed, a.expect_overgen) else ["flow", "H1", "H2", "H3"],
           "expect": {"paths": a.expect_paths, "observed": a.expect_observed, "overgen": a.expect_overgen}}
    with open(os.path.join(a.outdir, outs["json"]), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)

    print("script_sha256_raw %s  script_sha256_lf %s" % (res["script_sha256_raw"], res["script_sha256_lf"]))
    print("gene %s | pc transcripts %d | productions %d | |L(G)| %d | observed distinct %d | overgeneration %d"
          % (a.gene, g["n_tx"], len(g["prods"]), g["n_L"], g["n_obs_distinct"], g["overgen"]))
    print("path prob sum (float) %.12f | stops passed: %s" % (res["path_prob_sum_float"], ",".join(res["hard_stops_passed"])))
    print("\n%-4s %-10s %-9s %s" % ("id", "prob", "observed", "path"))
    for i, (p, pr) in enumerate(ranked, 1):
        print("P%02d  %-10s %-9s %s" % (i, f12(pr), g["obs_cnt"].get(p, 0), ">".join(p)))
    for fn in outs.values():
        pth = os.path.join(a.outdir, fn)
        print("WROTE %s  bytes=%d" % (pth, os.path.getsize(pth)))


if __name__ == "__main__":
    try:
        main()
    except GrammarError as e:
        print("STOP: " + str(e))
        sys.exit(1)