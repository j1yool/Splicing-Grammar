"""dmd_clinvar_extract_v1.py (Day 100, Splicing Grammar, DMD validation)
Applies D-SG1 and D-SG2 of dmd_splice_prereg_scorer_v1.py (commit 66f4c97) unchanged.
Written after reading ONLY the header line of variant_summary.txt.gz (clinvar_vs_fetch_header_v1.py).
No exon mapping, no scoring, no outcome. Output = labelled variant table plus counts.

D-SG1 here: Assembly == GRCh38, GeneSymbol (split on ';') contains DMD, Type == 'single nucleotide variant',
Start inside the Ensembl-verified DMD span (31097677-33339609, checked Day 100). The exon-window test of
D-SG1 and the modeled consequence of D-SG3 belong to the evaluation driver, not to this file.
D-SG2 here: ClinicalSignificance EXACTLY one of Pathogenic / Likely pathogenic / Pathogenic/Likely pathogenic = 1;
Benign / Likely benign / Benign/Likely benign = 0. Every other string (uncertain significance, conflicting,
compound strings such as 'Pathogenic, risk factor') is excluded and counted by string. ReviewStatus must be
in ONE_PLUS (at least one star). A ReviewStatus string that is in none of ONE_PLUS, ZERO or CONFLICT, on a row
with a labelled significance, stops the run: nothing is guessed.
FIXED BEFORE ANY ROW IS READ: exact-match handling of compound significance strings. Changing it later is a
dated amendment in the manuscript.

Usage:
  python dmd_clinvar_extract_v1.py --selftest
  python dmd_clinvar_extract_v1.py --input variant_summary.txt.gz [--outdir .] [--overwrite]
"""
import argparse, csv, gzip, hashlib, json, os, sys, tempfile
from collections import Counter

EXPECT_SHA256 = "7c39f2d875f3cd1ccd2d568662f06a7df647dd91fc0f51a072a758b08004c72e"
EXPECT_NCOLS = 43
SPAN = (31097677, 33339609)
SNV = "single nucleotide variant"
POS = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
NEG = {"Benign", "Likely benign", "Benign/Likely benign"}
ONE_PLUS = {"criteria provided, single submitter",
            "criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel",
            "practice guideline"}
ZERO = {"no assertion criteria provided", "no classification provided", "no assertion provided",
        "no classification for the single variant"}
CONFLICT = {"criteria provided, conflicting classifications", "criteria provided, conflicting interpretations"}
NEEDED = ["#AlleleID", "Type", "Name", "GeneSymbol", "ClinicalSignificance", "Assembly", "Chromosome", "Start",
          "Stop", "ReviewStatus", "VariationID", "ReferenceAlleleVCF", "AlternateAlleleVCF", "PositionVCF"]
OUT_COLS = ["VariationID", "AlleleID", "Name", "GeneSymbol", "Chromosome", "Start", "Stop", "PositionVCF",
            "ReferenceAlleleVCF", "AlternateAlleleVCF", "ClinicalSignificance", "ReviewStatus", "label"]


class ExtractError(Exception):
    pass


def stop(msg):
    raise ExtractError(msg)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def classify(r, counts):
    """r: dict keyed by NEEDED names. Returns an output record dict or None (reason counted)."""
    if r["Assembly"] != "GRCh38":
        counts["drop_assembly_not_GRCh38"] += 1
        return None
    if "DMD" not in r["GeneSymbol"].split(";"):
        counts["drop_gene_not_exactly_DMD"] += 1
        return None
    if r["GeneSymbol"] != "DMD":
        counts["note_multigene_symbol_kept_in_funnel"] += 1
    if r["Chromosome"] != "X":
        stop("DMD row on chromosome " + r["Chromosome"] + " (VariationID " + r["VariationID"] + ")")
    if r["Type"] != SNV:
        counts["drop_type_not_SNV"] += 1
        return None
    try:
        start, end = int(r["Start"]), int(r["Stop"])
    except ValueError:
        stop("non-integer Start/Stop at VariationID " + r["VariationID"])
    if start != end:
        counts["drop_SNV_start_ne_stop"] += 1
        return None
    if not (SPAN[0] <= start <= SPAN[1]):
        counts["drop_outside_DMD_span"] += 1
        return None
    sig = r["ClinicalSignificance"]
    if sig in POS:
        label = 1
    elif sig in NEG:
        label = 0
    else:
        counts["excl_significance: " + sig] += 1
        return None
    rs = r["ReviewStatus"]
    if rs in CONFLICT:
        counts["excl_review_conflicting"] += 1
        return None
    if rs in ZERO:
        counts["excl_review_zero_star"] += 1
        return None
    if rs not in ONE_PLUS:
        stop("UNKNOWN ReviewStatus string on a labelled row: " + repr(rs) + " (VariationID " + r["VariationID"] + ")")
    return {"VariationID": r["VariationID"], "AlleleID": r["#AlleleID"], "Name": r["Name"],
            "GeneSymbol": r["GeneSymbol"], "Chromosome": r["Chromosome"], "Start": start, "Stop": end,
            "PositionVCF": r["PositionVCF"], "ReferenceAlleleVCF": r["ReferenceAlleleVCF"],
            "AlternateAlleleVCF": r["AlternateAlleleVCF"], "ClinicalSignificance": sig, "ReviewStatus": rs,
            "label": label}


def run(path, outdir, expect_sha, overwrite, ncols=EXPECT_NCOLS):
    out_tsv = os.path.join(outdir, "dmd_clinvar_extract_v1.tsv")
    out_json = os.path.join(outdir, "dmd_clinvar_extract_v1_counts.json")
    for fn in (out_tsv, out_json):
        if os.path.exists(fn) and not overwrite:
            stop("output exists (use --overwrite): " + fn)
    got = sha256_file(path)
    if got != expect_sha:
        stop("input sha256 " + got + " != recorded " + expect_sha)
    counts, recs, n_rows, n_prefilter = Counter(), [], 0, 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as g:
        rd = csv.reader(g, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(rd)
        if len(header) != ncols or header[0] != "#AlleleID":
            stop("header is not the recorded " + str(ncols) + "-column header")
        miss = [c for c in NEEDED if c not in header]
        if miss:
            stop("needed columns missing: " + str(miss))
        ix = {c: header.index(c) for c in NEEDED}
        gi = ix["GeneSymbol"]
        for ln, row in enumerate(rd, start=2):
            n_rows += 1
            if len(row) != ncols:
                stop("line " + str(ln) + " has " + str(len(row)) + " fields, expected " + str(ncols))
            if "DMD" not in row[gi]:
                continue
            n_prefilter += 1
            rec = classify({c: row[ix[c]] for c in NEEDED}, counts)
            if rec:
                recs.append(rec)
    seen, kept = {}, []
    for rec in recs:
        k = rec["VariationID"]
        if k in seen:
            counts["dup_VariationID_dropped"] += 1
            if seen[k] != rec["label"]:
                stop("same VariationID with different labels: " + k)
            continue
        seen[k] = rec["label"]
        kept.append(rec)
    with open(out_tsv, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(OUT_COLS) + "\n")
        for rec in kept:
            f.write("\t".join(str(rec[c]) for c in OUT_COLS) + "\n")
    summ = {"input": os.path.basename(path), "input_sha256": got, "n_rows_read": n_rows,
            "n_rows_with_DMD_substring": n_prefilter, "n_pathogenic_label1": sum(r["label"] == 1 for r in kept),
            "n_benign_label0": sum(r["label"] == 0 for r in kept),
            "kept_review_status": dict(Counter(r["ReviewStatus"] for r in kept)),
            "kept_significance": dict(Counter(r["ClinicalSignificance"] for r in kept)),
            "counts": dict(sorted(counts.items())), "span": list(SPAN),
            "output_tsv_sha256": sha256_file(out_tsv)}
    with open(out_json, "w", encoding="utf-8", newline="\n") as f:
        json.dump(summ, f, indent=2, sort_keys=True)
        f.write("\n")
    return summ


# ---------------------------------------------------------------- selftest
def _row(**kw):
    base = {"#AlleleID": "1", "Type": SNV, "Name": "NM_x(DMD):c.1A>G", "GeneSymbol": "DMD",
            "ClinicalSignificance": "Pathogenic", "Assembly": "GRCh38", "Chromosome": "X", "Start": "32000000",
            "Stop": "32000000", "ReviewStatus": "criteria provided, single submitter", "VariationID": "100",
            "ReferenceAlleleVCF": "A", "AlternateAlleleVCF": "G", "PositionVCF": "32000000"}
    base.update(kw)
    return base


def selftest():
    ok = [0, 0]

    def check(name, cond):
        ok[0 if cond else 1] += 1
        print(("PASS " if cond else "FAIL ") + name)

    def one(**kw):
        c = Counter()
        return classify(_row(**kw), c), c

    r, _ = one()
    check("P + one star kept, label 1", r is not None and r["label"] == 1)
    r, _ = one(ClinicalSignificance="Likely pathogenic")
    check("LP kept, label 1", r is not None and r["label"] == 1)
    r, _ = one(ClinicalSignificance="Pathogenic/Likely pathogenic")
    check("P/LP kept, label 1", r is not None and r["label"] == 1)
    r, _ = one(ClinicalSignificance="Benign", ReviewStatus="reviewed by expert panel")
    check("B kept, label 0", r is not None and r["label"] == 0)
    r, _ = one(ClinicalSignificance="Likely benign")
    check("LB kept, label 0", r is not None and r["label"] == 0)
    r, c = one(ClinicalSignificance="Uncertain significance")
    check("VUS excluded and counted", r is None and c["excl_significance: Uncertain significance"] == 1)
    r, c = one(ClinicalSignificance="Conflicting classifications of pathogenicity")
    check("conflicting excluded", r is None and sum(c.values()) == 1)
    r, c = one(ClinicalSignificance="Pathogenic, risk factor")
    check("compound string excluded (exact match)", r is None and c["excl_significance: Pathogenic, risk factor"] == 1)
    r, c = one(ReviewStatus="no assertion criteria provided")
    check("zero-star excluded", r is None and c["excl_review_zero_star"] == 1)
    r, c = one(Assembly="GRCh37")
    check("GRCh37 dropped", r is None and c["drop_assembly_not_GRCh38"] == 1)
    r, c = one(Type="Deletion")
    check("non-SNV dropped", r is None and c["drop_type_not_SNV"] == 1)
    r, c = one(GeneSymbol="DMDX")
    check("symbol DMDX is not DMD", r is None and c["drop_gene_not_exactly_DMD"] == 1)
    r, c = one(GeneSymbol="DMD;OTHER")
    check("DMD;OTHER kept and noted", r is not None and c["note_multigene_symbol_kept_in_funnel"] == 1)
    r, c = one(Start="31097676", Stop="31097676")
    check("one base below span dropped", r is None and c["drop_outside_DMD_span"] == 1)
    r, c = one(Start="33339609", Stop="33339609")
    check("span end inclusive", r is not None)
    r, c = one(Start="32000000", Stop="32000001")
    check("start != stop dropped", r is None and c["drop_SNV_start_ne_stop"] == 1)
    try:
        one(ReviewStatus="some new string")
        check("unknown ReviewStatus stops", False)
    except ExtractError:
        check("unknown ReviewStatus stops", True)
    try:
        one(Chromosome="1")
        check("DMD on wrong chromosome stops", False)
    except ExtractError:
        check("DMD on wrong chromosome stops", True)
    # end to end on a synthetic 5-column-count file
    ncols = EXPECT_NCOLS
    header = NEEDED + ["pad" + str(i) for i in range(ncols - len(NEEDED))]

    def line(**kw):
        d = _row(**kw)
        return "\t".join([d[c] for c in NEEDED] + [""] * (ncols - len(NEEDED)))
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "fake.txt.gz")
        rows = [line(VariationID="1"), line(VariationID="2", ClinicalSignificance="Benign"),
                line(VariationID="3", ClinicalSignificance="Uncertain significance"),
                line(VariationID="4", GeneSymbol="TTN"), line(VariationID="1")]
        with gzip.open(p, "wt", encoding="utf-8", newline="") as f:
            f.write("\t".join(header) + "\n")
            for x in rows:
                f.write(x + "\n")
        try:
            run(p, td, "0" * 64, False)
            check("wrong sha256 stops", False)
        except ExtractError:
            check("wrong sha256 stops", True)
        s = run(p, td, sha256_file(p), False)
        check("end-to-end: 1 P, 1 B, duplicate dropped", s["n_pathogenic_label1"] == 1 and s["n_benign_label0"] == 1
              and s["counts"].get("dup_VariationID_dropped") == 1)
        try:
            run(p, td, sha256_file(p), False)
            check("existing output stops without --overwrite", False)
        except ExtractError:
            check("existing output stops without --overwrite", True)
    print("SELFTEST: " + str(ok[0]) + " passed, " + str(ok[1]) + " failed")
    return ok[1] == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--input")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    try:
        if a.selftest:
            sys.exit(0 if selftest() else 1)
        if not a.input or not os.path.isfile(a.input) or not os.path.isdir(a.outdir):
            stop("--input must be an existing file and --outdir an existing folder")
        if not selftest():
            stop("selftest failed; real file not opened")
        s = run(a.input, a.outdir, EXPECT_SHA256, a.overwrite)
        print("rows read:", s["n_rows_read"], "| rows with DMD substring:", s["n_rows_with_DMD_substring"])
        print("LABELLED: pathogenic (1) =", s["n_pathogenic_label1"], "| benign (0) =", s["n_benign_label0"])
        print("kept ReviewStatus:", s["kept_review_status"])
        print("kept ClinicalSignificance:", s["kept_significance"])
        for k, v in s["counts"].items():
            print("  " + k + ": " + str(v))
        print("output sha256:", s["output_tsv_sha256"])
    except ExtractError as e:
        print("STOP: " + str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()