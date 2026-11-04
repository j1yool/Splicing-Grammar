"""dmd_splice_eval_v1.py (Day 100, Splicing Grammar, DMD validation)
EVALUATION DRIVER for dmd_splice_prereg_scorer_v1.py (commit 66f4c97). It imports map_position, read_exons,
build_scope and evaluate from that file UNCHANGED and adds nothing to the design: D-SG1..D-SG7 are applied
exactly as registered. The verdict is whatever evaluate() computes (PASS / FAIL / DESCRIPTIVE), reported as computed.

Order of operations (hard stops on every pin):
 1. scorer file LF sha256 starts with 3452a563 (the committed file)
 2. extractor LF sha256 equals the pinned value; extract table sha256 equals ce00c917... (extractor run output)
 3. transcripts/exons/grammar-script sha256 equal those recorded in dmd_scorer_v1_build.json
 4. the rebuilt scope reproduces the recorded counts, floor and reference path
 5. map every extracted variant (D-SG1/D-SG3), score scored exons (D-SG4/D-SG5), evaluate (D-SG6)

Usage:
  python dmd_splice_eval_v1.py --selftest
  python dmd_splice_eval_v1.py --transcripts splice_transcripts_v1.tsv --exons splice_exons_v1.tsv
      [--extract dmd_clinvar_extract_v1.tsv] [--build-json dmd_scorer_v1_build.json] [--outdir .] [--overwrite]
"""
import argparse, csv, hashlib, json, os, sys, tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

SCORER_LF_PREFIX = "3452a563"
EXTRACT_TSV_SHA256 = "ce00c9175e96d4f81d438ea4565034fdfb6f58f039a16936c5cb7f85c1d2454f"
EXTRACTOR_LF_SHA256 = "ba5e5c4455b90031d8ecd158ee1427e9b02f59bea4558f25d753a6b461643768"
EXTRACT_COLS = ["VariationID", "AlleleID", "Name", "GeneSymbol", "Chromosome", "Start", "Stop", "PositionVCF",
                "ReferenceAlleleVCF", "AlternateAlleleVCF", "ClinicalSignificance", "ReviewStatus", "label"]


class EvalError(Exception):
    pass


def stop(msg):
    raise EvalError(msg)


def lf_bytes(path, strip_trailing=False):
    b = open(path, "rb").read().replace(b"\r\n", b"\n")
    return b.rstrip(b"\n") if strip_trailing else b


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_extract(path):
    out = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames != EXTRACT_COLS:
            stop("extract header mismatch: " + str(rd.fieldnames))
        for r in rd:
            lab = int(r["label"])
            if lab not in (0, 1):
                stop("bad label for VariationID " + r["VariationID"])
            out.append({"id": r["VariationID"], "chrom": r["Chromosome"], "pos": int(r["Start"]), "label": lab,
                        "name": r["Name"]})
    if not out:
        stop("extract table is empty")
    if len({v["id"] for v in out}) != len(out):
        stop("duplicate VariationID in extract table")
    return out


def run_eval(sc, observed, ex, variants):
    """Pure function: scope -> mapping -> records -> sc.evaluate. Returns a dict."""
    scope = sc.build_scope(observed, ex)
    gene_chrom = scope["exons"][0]["chrom"]
    if gene_chrom.replace("chr", "") != "X":
        stop("exon table chromosome is " + gene_chrom + ", expected X for DMD")
    by_exon = {r["exon_id"]: r for r in scope["rows"]}
    status, kinds, mapped, records = Counter(), Counter(), [], []
    for v in variants:
        if v["chrom"].replace("chr", "") != gene_chrom.replace("chr", ""):
            stop("variant on another chromosome: " + v["id"])
        st, eid, kind = sc.map_position(gene_chrom, v["pos"], scope["exons"])
        status[st] += 1
        mapped.append((v, st, eid, kind))
        if st == "scored":
            row = by_exon[eid]
            kinds[(kind, v["label"])] += 1
            records.append({"exon": eid, "label": v["label"], "v": row["v_skip"], "frame": row["frame_pred"],
                            "id": v["id"], "kind": kind, "pos": v["pos"]})
    res = sc.evaluate(records)
    return {"scope": scope, "status": dict(status), "kinds": kinds, "mapped": mapped, "records": records,
            "res": res}


def pins(scorer_path, extractor_path, extract_tsv):
    h = hashlib.sha256(lf_bytes(scorer_path)).hexdigest()
    if not h.startswith(SCORER_LF_PREFIX):
        stop("scorer LF sha256 " + h[:16] + " does not start with " + SCORER_LF_PREFIX)
    e = hashlib.sha256(lf_bytes(extractor_path, True)).hexdigest()
    if e != EXTRACTOR_LF_SHA256:
        stop("extractor LF sha256 " + e + " != pinned " + EXTRACTOR_LF_SHA256)
    t = sha256_file(extract_tsv)
    if t != EXTRACT_TSV_SHA256:
        stop("extract table sha256 " + t + " != recorded " + EXTRACT_TSV_SHA256)
    return {"scorer_lf_sha256": h, "extractor_lf_sha256": e, "extract_tsv_sha256": t}


# ---------------------------------------------------------------- selftest
def _synthetic():
    def ex(i, s, e):
        return {"id": "E%d" % i, "chrom": "X", "start": s, "end": e, "strand": "+", "length": e - s + 1}
    ex_all = {d["id"]: d for d in [ex(1, 1000, 1100), ex(2, 2000, 2100), ex(3, 2110, 2200), ex(4, 3000, 3050),
                                  ex(5, 4000, 4100), ex(6, 5000, 5100)]}
    observed = [("E1", "E2", "E4", "E5", "E6"), ("E1", "E2", "E3", "E4", "E5", "E6"),
                ("E1", "E2", "E3", "E4", "E5", "E6")]
    return ex_all, observed


def selftest():
    import dmd_splice_prereg_scorer_v1 as sc
    ok = [0, 0]

    def check(name, cond):
        ok[0 if cond else 1] += 1
        print(("PASS " if cond else "FAIL ") + name)

    ex_all, observed = _synthetic()
    V = lambda i, pos, lab: {"id": str(i), "chrom": "X", "pos": pos, "label": lab, "name": "v" + str(i)}
    variants = [V(1, 2105, 1),     # E2 donor intron +5 AND E3 acceptor intron -5: two windows -> ambiguous
                V(2, 3051, 0),     # E4 donor intron +1 only -> scored
                V(3, 3999, 1),     # E5 acceptor intron -1 only -> scored
                V(4, 3500, 0),     # nowhere -> out_of_scope
                V(5, 2050, 1)]     # interior of E2, outside every window -> out_of_scope
    out = run_eval(sc, observed, ex_all, variants)
    s = out["status"]
    check("ambiguous variant counted", s.get("ambiguous") == 1)
    check("two scored variants", s.get("scored") == 2)
    check("two out_of_scope", s.get("out_of_scope") == 2)
    ids = {r["id"]: r for r in out["records"]}
    check("E4 donor-intron variant scored on E4", ids["2"]["exon"] == "E4" and ids["2"]["kind"] == "donor_intron")
    check("E5 acceptor-intron variant scored on E5", ids["3"]["exon"] == "E5" and ids["3"]["kind"] == "acceptor_intron")
    row4 = [r for r in out["scope"]["rows"] if r["exon_id"] == "E4"][0]
    check("frame prediction from scorer (E4 length 51 -> mod3 0 -> 0)", ids["2"]["frame"] == 0 and row4["mod3"] == 0)
    check("tiny set is DESCRIPTIVE", out["res"]["verdict"] == "DESCRIPTIVE")
    # off-reference exon: E3 is in the internal set but not on the reference path when ref = the longer path;
    # here the reference path is the longest (6 exons, includes E3), so make a case with an extra exon
    ex2 = dict(ex_all)
    ex2["E7"] = {"id": "E7", "chrom": "X", "start": 6000, "end": 6100, "strand": "+", "length": 101}
    obs2 = observed + [("E1", "E7", "E6")]
    out2 = run_eval(sc, obs2, ex2, [V(9, 5999, 1)])
    check("variant in window of off-reference exon is off_reference", out2["status"].get("off_reference") == 1)
    try:
        run_eval(sc, observed, ex_all, [dict(V(10, 3051, 0), chrom="1")])
        check("variant on wrong chromosome stops", False)
    except EvalError:
        check("variant on wrong chromosome stops", True)
    with tempfile.TemporaryDirectory() as td:
        a, b = os.path.join(td, "a.py"), os.path.join(td, "b.tsv")
        open(a, "wb").write(b"x = 1\r\n")
        open(b, "wb").write(b"t\n")
        try:
            pins(a, a, b)
            check("scorer pin stops on a wrong file", False)
        except EvalError:
            check("scorer pin stops on a wrong file", True)
    print("SELFTEST: " + str(ok[0]) + " passed, " + str(ok[1]) + " failed")
    return ok[1] == 0


# ---------------------------------------------------------------- real run
def fmt(x):
    return "None" if x is None else ("%.4f" % x if isinstance(x, float) else str(x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--transcripts")
    ap.add_argument("--exons")
    ap.add_argument("--extract", default="dmd_clinvar_extract_v1.tsv")
    ap.add_argument("--build-json", default="dmd_scorer_v1_build.json")
    ap.add_argument("--gene", default="DMD")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    try:
        if a.selftest:
            sys.exit(0 if selftest() else 1)
        for p in (a.transcripts, a.exons, a.extract, a.build_json):
            if not p or not os.path.isfile(p):
                stop("missing input file: " + str(p))
        if not os.path.isdir(a.outdir):
            stop("outdir is not a folder")
        outs = [os.path.join(a.outdir, f) for f in ("dmd_eval_v1_variants.tsv", "dmd_eval_v1_results.json")]
        for fn in outs:
            if os.path.exists(fn) and not a.overwrite:
                stop("output exists (use --overwrite): " + fn)
        if not selftest():
            stop("selftest failed; real inputs not read")
        import dmd_splice_prereg_scorer_v1 as sc
        import splice_grammar_v1 as sg
        pinned = pins(os.path.join(HERE, "dmd_splice_prereg_scorer_v1.py"),
                      os.path.join(HERE, "dmd_clinvar_extract_v1.py"), a.extract)
        bj = json.load(open(a.build_json, encoding="utf-8"))
        for key, path in (("transcripts_sha256", a.transcripts), ("exons_sha256", a.exons),
                          ("grammar_script_sha256", sg.__file__)):
            if sha256_file(path) != bj[key]:
                stop(key + " differs from the value recorded in " + a.build_json)
        observed = sg.read_observed(a.transcripts, a.gene)
        ex = sc.read_exons(a.exons, a.gene)
        variants = read_extract(a.extract)
        out = run_eval(sc, observed, ex, variants)
        scope = out["scope"]
        got = {"pc_transcripts": len(observed), "collapsed_edges": scope["n_edges"], "ref_exons": len(scope["ref"]),
               "ref_internal": len(scope["rows"]), "internal_set": len(scope["internal_set"])}
        if got != bj["counts"]:
            stop("rebuilt scope " + str(got) + " != recorded " + str(bj["counts"]))
        if abs(scope["floor"] - bj["floor_probability"]) > 1e-15 or list(scope["ref"]) != bj["reference_path"]:
            stop("rebuilt floor or reference path differs from the recorded build")
        res = out["res"]
        print("pins OK:", {k: v[:16] for k, v in pinned.items()})
        print("scope reproduced:", got)
        print("variants read:", len(variants), "| mapping status:", dict(sorted(out["status"].items())))
        print("D-SG3 exclusions: ambiguous =", out["status"].get("ambiguous", 0), "| off_reference =",
              out["status"].get("off_reference", 0), "| out_of_scope =", out["status"].get("out_of_scope", 0))
        print("SCORED: n variants %d (pathogenic %d, benign %d) over %d exons" %
              (res["n_variants"], res["n_pathogenic"], res["n_benign"], res["n_exons"]))
        print("AUC(V) %s | AUC(frame) %s | diff %s" % (fmt(res["auc_grammar"]), fmt(res["auc_frame"]), fmt(res.get("diff"))))
        if "ci_lo" in res:
            print("95%% interval [%s, %s] | boot used %d, dropped %d" %
                  (fmt(res["ci_lo"]), fmt(res["ci_hi"]), res["boot_used"], res["boot_dropped"]))
        print("VERDICT (as computed): " + res["verdict"])
        print("scored variants by window and label (descriptive, not a test):")
        for (kind, lab), n in sorted(out["kinds"].items()):
            print("  %-16s label %d : %d" % (kind, lab, n))
        with open(outs[0], "w", encoding="utf-8", newline="\n") as f:
            f.write("\t".join(["VariationID", "Start", "label", "status", "exon_id", "window"]) + "\n")
            for v, st, eid, kind in out["mapped"]:
                f.write("\t".join([v["id"], str(v["pos"]), str(v["label"]), st, eid or "", kind or ""]) + "\n")
        summary = {"results": res, "mapping_status": out["status"], "pins": pinned, "scope_counts": got,
                   "kinds": {k[0] + "|" + str(k[1]): n for k, n in out["kinds"].items()}}
        with open(outs[1], "w", encoding="utf-8", newline="\n") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
            f.write("\n")
        for fn in outs:
            print("WROTE", fn, "bytes=" + str(os.path.getsize(fn)))
    except EvalError as e:
        print("STOP: " + str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()