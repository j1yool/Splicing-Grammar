"""
dmd_splice_prereg_scorer_v1.py  (Day 99, Splicing Grammar, DMD validation)
PRE-REGISTRATION HEADER AND SCORER. NO CLINVAR VARIANT HAS BEEN READ. The first commit containing this text
precedes any ClinVar variant row being opened (R6: this result is headed to a manuscript).

QUESTION. For single-nucleotide DMD splice-region variants, does a grammar-based violation score separate
ClinVar pathogenic from benign variants better than the reading-frame rule?

PLAN IMPACT. Brings forward the Days 241-270 step "test grammar against known pathogenic splicing mutations".
Not a pivot (no R7). FGFR2 remains the fallback gene.

DESIGN (D-SG1..D-SG6 approved Day 99; D-SG7 added to make D-SG4 computable)
D-SG1 Variants: ClinVar variant_summary, GRCh38, gene DMD, single-nucleotide. In scope if the position is
      intronic +1..+8 after an exon's donor, intronic -1..-20 before an exon's acceptor, or within 3 nt of
      an exon edge on the exonic side, for an INTERNAL exon (see D-SG4). Strand is handled per exon row.
D-SG2 Labels: pathogenic + likely pathogenic = 1; benign + likely benign = 0; VUS and conflicting excluded;
      review status of at least one star. ClinVar parsing is written AFTER this commit, after reading only
      the file's header line, and applies D-SG1/D-SG2 unchanged.
D-SG3 Modeled consequence: every in-scope variant is modeled as skipping the single exon whose window holds
      it. A position inside the windows of more than one internal exon is excluded and counted. A position
      in the window of an internal exon that is not on the reference path (D-SG7) is excluded and counted.
D-SG4 Grammar and score. Right-linear PCFG as in splice_grammar_v1.py (P(u,v) = count(u,v)/out_count(u)),
      built from the protein-coding transcripts of the gene with each path's first and last exon collapsed
      to START and END (internal exons only). Floor for an unseen edge = (smallest observed edge
      probability) / 10, fixed by the grammar before any variant is read. Violation score
      V = -ln P(post-variant path) with the floor applied to unseen edges. Higher V = more violation,
      predicted pathogenic. Ties get midranks in every AUC.
D-SG5 Baseline: reading-frame rule (Monaco 1988). Skipping an exon whose length mod 3 is 0 is predicted
      benign (score 0); any other length is predicted pathogenic (score 1).
D-SG6 Outcome: AUC(V) minus AUC(frame), pathogenic = positive class. 95% percentile interval from 10,000
      bootstrap replicates that resample EXONS with replacement (all variants of a drawn exon come along),
      seed 99; replicates in which either class is empty are discarded and counted. PASS only if the
      lower bound is above 0 and the point difference is at least 0.01 (the Option C rule). If either class
      has fewer than 10 variants or fewer than 5 exons, the result is DESCRIPTIVE and gets no verdict.
      No threshold is tuned. A failed rule is reported as computed.
D-SG7 Reference path: one reference path, the protein-coding path with the most exons (ties: smallest exon-id
      tuple). Only internal exons on it are scored; the post-variant path is the reference path with that
      exon removed, so every V differs from the reference value by local edge terms only.

ASSUMPTIONS AND LIMITS (page numbers are from Edward's reading notes; verify against the PDFs before use)
A1 "Variant -> exon skipped" is a default, not a demonstrated outcome. Robberson 1990 (MCB 10:84-94)
   shows no skipping in its own data; Wang & Burge 2008 (RNA 14:802-813) cite Robberson for skipping, so
   the two are not independent support.
A2 Cryptic-site use and intron retention are alternative outcomes (Robberson p.92; Wang & Burge pp.804-805)
   and are not modeled.
A3 Monaco 1988 (Genomics 2:90-95) states the frame rule for genomic deletions: 3 DMD and 5 BMD patients
   (the abstract says three BMD; the body gives five), 0 of 8 exceptions, retrospective, with 2 patients
   not sequence-confirmed. Extending it to splice variants is speculation there (p.94).
A4 Effects of a sequence change depend on cell type and context (Wang & Burge pp.806-807, 809); the
   measured outcome in one assay may not transfer to muscle.
A5 Collapsing first and last exons to START/END merges alternative-promoter and alternative-3'-end isoforms
   of DMD; internal-exon structure only is modeled.
A6 The exon table is assumed to be GRCh38; the build step prints the gene's lowest start and highest end for
   a check against Ensembl before any ClinVar position is mapped.

Inputs: splice_transcripts_v1.tsv, splice_exons_v1.tsv (schemas fixed by splice_exongraph_v1.py);
        splice_grammar_v1.py in the same folder (read_observed is reused).
Usage:
  python dmd_splice_prereg_scorer_v1.py --selftest
  python dmd_splice_prereg_scorer_v1.py --build --transcripts splice_transcripts_v1.tsv --exons splice_exons_v1.tsv
      --gene DMD --expect-pc-transcripts N --expect-collapsed-edges N --expect-ref-exons N
      --expect-ref-internal N --expect-internal-set N [--outdir .] [--overwrite]
"""
import argparse, csv, hashlib, json, math, os, random, sys
from collections import Counter, defaultdict

DONOR_INTRON, ACCEPTOR_INTRON, EXON_EDGE = 8, 20, 3
FLOOR_DIV = 10.0
N_BOOT, SEED = 10000, 99
MIN_PER_CLASS, MIN_EXONS, MIN_DIFF = 10, 5, 0.01


class ScorerError(Exception):
    pass


def stop(msg):
    raise ScorerError(msg)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- grammar
def build_chains(observed):
    if not observed:
        stop("no observed paths")
    return [("START",) + tuple(p[1:-1]) + ("END",) for p in observed]


def edge_model(chains):
    cnt, out = Counter(), Counter()
    for c in chains:
        for i in range(len(c) - 1):
            cnt[(c[i], c[i + 1])] += 1
            out[c[i]] += 1
    probs = {e: n / out[e[0]] for e, n in cnt.items()}
    return cnt, out, probs, min(probs.values()) / FLOOR_DIV


def reference_path(observed):
    return min(observed, key=lambda p: (-len(p), p))


def chain_v(chain, probs, floor):
    return -sum(math.log(probs.get((chain[i], chain[i + 1]), floor)) for i in range(len(chain) - 1))


def skip_scores(ref, probs, floor):
    inner = list(ref[1:-1])
    out = []
    for i, e in enumerate(inner):
        post = ["START"] + inner[:i] + inner[i + 1:] + ["END"]
        out.append((e, chain_v(post, probs, floor)))
    return out


# ---------------------------------------------------------------- exon windows and mapping
def exon_windows(ex):
    s, e, st = ex["start"], ex["end"], ex["strand"]
    if st == "-":
        return [(e + 1, e + ACCEPTOR_INTRON, "acceptor_intron"), (e - EXON_EDGE + 1, e, "acceptor_exon"),
                (s, s + EXON_EDGE - 1, "donor_exon"), (s - DONOR_INTRON, s - 1, "donor_intron")]
    if st == "+":
        return [(s - ACCEPTOR_INTRON, s - 1, "acceptor_intron"), (s, s + EXON_EDGE - 1, "acceptor_exon"),
                (e - EXON_EDGE + 1, e, "donor_exon"), (e + 1, e + DONOR_INTRON, "donor_intron")]
    stop("strand must be + or -: " + str(st))


def map_position(chrom, pos, exons):
    """exons: list of dict(id, chrom, start, end, strand, on_ref); internal exons only.
    Returns (status, exon_id, window)."""
    hit = {}
    for ex in exons:
        if ex["chrom"] != chrom:
            continue
        for lo, hi, kind in exon_windows(ex):
            if lo <= pos <= hi:
                hit.setdefault(ex["id"], kind)
    if not hit:
        return ("out_of_scope", None, None)
    if len(hit) > 1:
        return ("ambiguous", None, None)
    eid = next(iter(hit))
    on_ref = [x for x in exons if x["id"] == eid][0]["on_ref"]
    return ("scored" if on_ref else "off_reference", eid, hit[eid])


# ---------------------------------------------------------------- statistics
def auc(scores, labels):
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return None
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    rsum = sum(ranks[k] for k in range(len(scores)) if labels[k] == 1)
    return (rsum - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def evaluate(records, n_boot=N_BOOT, seed=SEED):
    """records: dicts with exon, label (1 pathogenic, 0 benign), v, frame."""
    labels = [r["label"] for r in records]
    a_g = auc([r["v"] for r in records], labels)
    a_f = auc([r["frame"] for r in records], labels)
    n_pos, n_neg = sum(labels), len(labels) - sum(labels)
    n_ex = len({r["exon"] for r in records})
    res = {"n_variants": len(records), "n_pathogenic": n_pos, "n_benign": n_neg, "n_exons": n_ex,
           "auc_grammar": a_g, "auc_frame": a_f}
    if a_g is None or n_pos < MIN_PER_CLASS or n_neg < MIN_PER_CLASS or n_ex < MIN_EXONS:
        res.update({"verdict": "DESCRIPTIVE", "diff": None if a_g is None else a_g - a_f})
        return res
    by_exon = defaultdict(list)
    for r in records:
        by_exon[r["exon"]].append(r)
    keys = sorted(by_exon)
    rng = random.Random(seed)
    diffs, dropped = [], 0
    for _ in range(n_boot):
        samp = [r for k in (rng.choice(keys) for _ in keys) for r in by_exon[k]]
        lab = [r["label"] for r in samp]
        g = auc([r["v"] for r in samp], lab)
        if g is None:
            dropped += 1
            continue
        diffs.append(g - auc([r["frame"] for r in samp], lab))
    diffs.sort()
    lo = diffs[int(0.025 * (len(diffs) - 1))]
    hi = diffs[int(round(0.975 * (len(diffs) - 1)))]
    diff = a_g - a_f
    res.update({"diff": diff, "ci_lo": lo, "ci_hi": hi, "boot_used": len(diffs), "boot_dropped": dropped,
                "verdict": "PASS" if (lo > 0 and diff >= MIN_DIFF) else "FAIL"})
    return res


# ---------------------------------------------------------------- inputs
def read_exons(path, gene):
    if not os.path.isfile(path):
        stop("exons file not found: " + path)
    need = ["gene", "exon_id", "chrom", "start", "end", "strand", "length", "n_transcripts",
            "n_protein_coding", "overlaps"]
    ex = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames != need:
            stop("exons header mismatch: %s" % rd.fieldnames)
        for r in rd:
            if r["gene"] != gene:
                continue
            d = {"id": r["exon_id"], "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                 "strand": r["strand"], "length": int(r["length"])}
            if d["end"] - d["start"] + 1 != d["length"]:
                stop("length != end-start+1 for " + d["id"])
            if d["strand"] not in "+-" or len(d["strand"]) != 1:
                stop("bad strand for " + d["id"])
            if d["id"] in ex:
                stop("duplicate exon id " + d["id"])
            ex[d["id"]] = d
    if not ex:
        stop("gene absent from exons file: " + gene)
    if len({(d["chrom"], d["strand"]) for d in ex.values()}) != 1:
        stop("exons of one gene span more than one chromosome or strand")
    return ex


def build_scope(observed, ex):
    ref = reference_path(observed)
    chains = build_chains(observed)
    cnt, out, probs, floor = edge_model(chains)
    internal_set = {e for p in observed for e in p[1:-1]}
    for e in internal_set | set(ref):
        if e not in ex:
            stop("exon %s in a path but not in the exon table" % e)
    rows = []
    for e, v in skip_scores(ref, probs, floor):
        d = ex[e]
        rows.append({"exon_id": e, "chrom": d["chrom"], "start": d["start"], "end": d["end"],
                     "strand": d["strand"], "length": d["length"], "mod3": d["length"] % 3,
                     "frame_pred": 0 if d["length"] % 3 == 0 else 1, "v_skip": v})
    ref_ids = {r["exon_id"] for r in rows}
    exons = [dict(ex[e], on_ref=(e in ref_ids)) for e in sorted(internal_set)]
    return {"ref": ref, "cnt": cnt, "probs": probs, "floor": floor, "rows": rows, "exons": exons,
            "internal_set": internal_set, "n_edges": len(cnt)}


# ---------------------------------------------------------------- selftest
def selftest():
    res = []

    def check(name, ok):
        res.append(ok)
        print("%s: %s" % ("PASS" if ok else "FAIL", name))

    near = lambda a, b: abs(a - b) < 1e-9
    obs = [("F1", "A", "B", "C", "L1"), ("F2", "A", "B", "C", "L2"), ("F3", "A", "C", "L3"), ("F4", "B", "C", "L4")]
    ex = {"A": dict(id="A", chrom="X", start=1000, end=1029, strand="-", length=30),
          "B": dict(id="B", chrom="X", start=2000, end=2030, strand="-", length=31),
          "C": dict(id="C", chrom="X", start=3000, end=3031, strand="-", length=32)}
    for k in ("F1", "F2", "F3", "F4", "L1", "L2", "L3", "L4"):
        ex[k] = dict(id=k, chrom="X", start=9000, end=9009, strand="-", length=10)
    sc = build_scope(obs, ex)
    check("reference path = longest, smallest tuple", sc["ref"] == ("F1", "A", "B", "C", "L1"))
    check("edge counts and floor (0.25/10)", sc["n_edges"] == 6 and near(sc["floor"], 0.025))
    check("P(START,A) = 3/4 and P(A,C) = 1/3", near(sc["probs"][("START", "A")], 0.75) and near(sc["probs"][("A", "C")], 1 / 3))
    v = {r["exon_id"]: r["v_skip"] for r in sc["rows"]}
    check("V(skip A) = -ln 0.25", near(v["A"], -math.log(0.25)))
    check("V(skip B) = -ln 0.25 (tie with A)", near(v["B"], -math.log(0.25)))
    check("V(skip C) = -ln(0.75*2/3*0.025), floor used", near(v["C"], -math.log(0.75 * (2 / 3) * 0.025)))
    fp = {r["exon_id"]: r["frame_pred"] for r in sc["rows"]}
    check("frame: len 30 -> 0, 31 -> 1, 32 -> 1", (fp["A"], fp["B"], fp["C"]) == (0, 1, 1))
    check("one-exon and two-exon transcripts collapse to START>END", build_chains([("Q",), ("Q", "R")]) == [("START", "END")] * 2)
    check("AUC midranks: 0.875 by hand", near(auc([0.9, 0.5, 0.5, 0.1], [1, 1, 0, 0]), 0.875))
    check("AUC perfect 1.0, reversed 0.0, all ties 0.5",
          auc([3, 2, 1, 0], [1, 1, 0, 0]) == 1.0 and auc([0, 1, 2, 3], [1, 1, 0, 0]) == 0.0 and auc([1, 1, 1, 1], [1, 1, 0, 0]) == 0.5)
    check("AUC with an empty class is None", auc([1, 2], [1, 1]) is None)

    mi = dict(id="X1", chrom="X", start=1000, end=1099, strand="-", on_ref=True)
    pl = dict(id="Y1", chrom="X", start=5000, end=5099, strand="+", on_ref=True)
    zz = dict(id="Z1", chrom="X", start=5110, end=5200, strand="+", on_ref=True)
    of = dict(id="W1", chrom="X", start=8000, end=8100, strand="+", on_ref=False)
    exs = [mi, pl, zz, of]
    st = lambda pos: map_position("X", pos, exs)
    check("minus: donor +1 = 999, +8 = 992, +9 = 991 out",
          st(999) == ("scored", "X1", "donor_intron") and st(992)[0] == "scored" and st(991)[0] == "out_of_scope")
    check("minus: acceptor -1 = 1100, -20 = 1119, -21 = 1120 out",
          st(1100) == ("scored", "X1", "acceptor_intron") and st(1119)[0] == "scored" and st(1120)[0] == "out_of_scope")
    check("minus: exonic edges 1000..1002 and 1097..1099 in, 1003 and 1096 out",
          st(1000)[2] == "donor_exon" and st(1002)[0] == "scored" and st(1099)[2] == "acceptor_exon"
          and st(1097)[0] == "scored" and st(1003)[0] == "out_of_scope" and st(1096)[0] == "out_of_scope")
    check("plus: acceptor -1 = 4999, -20 = 4980, -21 out; donor +1 = 5100 sits in an overlap",
          st(4999) == ("scored", "Y1", "acceptor_intron") and st(4980)[0] == "scored" and st(4979)[0] == "out_of_scope" and st(5100)[0] == "ambiguous")
    check("plus: exonic 5000..5002 and 5097..5099 in; 5003 out",
          st(5000)[2] == "acceptor_exon" and st(5097)[0] != "out_of_scope" and st(5003)[0] == "out_of_scope")
    check("two windows of different exons -> ambiguous (5105)", st(5105)[0] == "ambiguous")
    check("exon not on the reference path -> off_reference; other chromosome -> out_of_scope",
          st(8000) == ("off_reference", "W1", "acceptor_exon") and map_position("Y", 999, exs)[0] == "out_of_scope")

    recs = [{"exon": "e%d" % (i // 2), "label": 1 if i < 12 else 0, "v": 10.0 - i, "frame": 0.5} for i in range(24)]
    r1, r2 = evaluate(recs, 500, SEED), evaluate(recs, 500, SEED)
    check("perfect grammar vs constant baseline: diff 0.5, PASS, deterministic",
          near(r1["diff"], 0.5) and r1["verdict"] == "PASS" and r1["ci_lo"] == r2["ci_lo"] and r1["ci_hi"] == r2["ci_hi"])
    recs2 = [{"exon": "e%d" % (i // 2), "label": 1 if i < 12 else 0, "v": float(i % 2), "frame": float(i % 2)} for i in range(24)]
    r3 = evaluate(recs2, 500, SEED)
    check("identical scores: diff 0 and FAIL", near(r3["diff"], 0.0) and r3["verdict"] == "FAIL")
    check("fewer than 10 per class or 5 exons -> DESCRIPTIVE",
          evaluate(recs[:6] + recs[-6:], 500, SEED)["verdict"] == "DESCRIPTIVE")

    def must_stop(fn, why):
        try:
            fn()
        except ScorerError:
            check("hard stop fires: " + why, True)
            return
        check("hard stop fires: " + why, False)

    must_stop(lambda: build_chains([]), "no observed paths")
    must_stop(lambda: exon_windows(dict(start=1, end=2, strand="?")), "bad strand")
    must_stop(lambda: build_scope(obs, {k: v for k, v in ex.items() if k != "B"}), "exon missing from the table")
    n = sum(res)
    print("SELFTEST %s %d/%d" % ("PASS" if n == len(res) else "FAIL", n, len(res)))
    if n != len(res):
        sys.exit(1)


# ---------------------------------------------------------------- build
def build(a):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import splice_grammar_v1 as sg
    outs = {"scope": "dmd_scorer_v1_scope.tsv", "json": "dmd_scorer_v1_build.json"}
    for fn in outs.values():
        if os.path.exists(os.path.join(a.outdir, fn)) and not a.overwrite:
            stop("output exists (use --overwrite): " + fn)
    observed = sg.read_observed(a.transcripts, a.gene)
    ex = read_exons(a.exons, a.gene)
    sc = build_scope(observed, ex)
    got = {"pc_transcripts": len(observed), "collapsed_edges": sc["n_edges"], "ref_exons": len(sc["ref"]),
           "ref_internal": len(sc["rows"]), "internal_set": len(sc["internal_set"])}
    exp = {"pc_transcripts": a.expect_pc_transcripts, "collapsed_edges": a.expect_collapsed_edges,
           "ref_exons": a.expect_ref_exons, "ref_internal": a.expect_ref_internal, "internal_set": a.expect_internal_set}
    for k in got:
        if exp[k] is None or got[k] != exp[k]:
            stop("%s = %s, expected %s" % (k, got[k], exp[k]))
    with open(os.path.join(a.outdir, outs["scope"]), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "exon_id", "chrom", "start", "end", "strand", "length", "mod3", "frame_pred", "v_skip"])
        for r in sc["rows"]:
            w.writerow([a.gene, r["exon_id"], r["chrom"], r["start"], r["end"], r["strand"], r["length"],
                        r["mod3"], r["frame_pred"], "%.6f" % r["v_skip"]])
    me = open(__file__, "rb").read()
    res = {"gene": a.gene, "transcripts_sha256": sha256_file(a.transcripts), "exons_sha256": sha256_file(a.exons),
           "grammar_script_sha256": sha256_file(sg.__file__), "script_sha256_raw": hashlib.sha256(me).hexdigest(),
           "script_sha256_lf": hashlib.sha256(me.replace(b"\r\n", b"\n")).hexdigest(), "counts": got,
           "floor_probability": sc["floor"], "reference_path": list(sc["ref"]),
           "gene_min_start": min(d["start"] for d in ex.values()), "gene_max_end": max(d["end"] for d in ex.values()),
           "constants": {"donor_intron": DONOR_INTRON, "acceptor_intron": ACCEPTOR_INTRON, "exon_edge": EXON_EDGE,
                         "floor_div": FLOOR_DIV, "n_boot": N_BOOT, "seed": SEED, "min_per_class": MIN_PER_CLASS,
                         "min_exons": MIN_EXONS, "min_diff": MIN_DIFF}}
    with open(os.path.join(a.outdir, outs["json"]), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    vs = [r["v_skip"] for r in sc["rows"]]
    print("script_sha256_raw %s  script_sha256_lf %s" % (res["script_sha256_raw"], res["script_sha256_lf"]))
    print("gene %s | pc transcripts %d | collapsed edges %d | reference path exons %d | scored internal exons %d | internal set %d"
          % (a.gene, got["pc_transcripts"], got["collapsed_edges"], got["ref_exons"], got["ref_internal"], got["internal_set"]))
    print("floor probability %.6g | distinct V values %d | V min %.4f max %.4f" % (sc["floor"], len(set(round(x, 9) for x in vs)), min(vs), max(vs)))
    print("frame_pred 0 (in-frame): %d | frame_pred 1: %d" % (sum(1 for r in sc["rows"] if r["frame_pred"] == 0), sum(1 for r in sc["rows"] if r["frame_pred"] == 1)))
    print("gene lowest start %d | highest end %d (check against Ensembl GRCh38 before mapping any ClinVar position)" % (res["gene_min_start"], res["gene_max_end"]))
    for fn in outs.values():
        p = os.path.join(a.outdir, fn)
        print("WROTE %s  bytes=%d" % (p, os.path.getsize(p)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--transcripts")
    ap.add_argument("--exons")
    ap.add_argument("--gene", default="DMD")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--overwrite", action="store_true")
    for k in ("pc-transcripts", "collapsed-edges", "ref-exons", "ref-internal", "internal-set"):
        ap.add_argument("--expect-" + k, type=int)
    a = ap.parse_args()
    if a.selftest:
        selftest()
    elif a.build:
        if not a.transcripts or not a.exons or not os.path.isdir(a.outdir):
            stop("--transcripts, --exons and an existing --outdir are required")
        build(a)
    else:
        stop("use --selftest or --build")


if __name__ == "__main__":
    try:
        main()
    except ScorerError as e:
        print("STOP: " + str(e))
        sys.exit(1)