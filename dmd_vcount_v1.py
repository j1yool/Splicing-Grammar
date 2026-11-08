import csv
rows = list(csv.DictReader(open("dmd_scorer_v1_scope.tsv", encoding="utf-8"), delimiter="\t"))
v = [float(r["v_skip"]) for r in rows]
print("rows", len(rows))
print("frame_pred 0 (in-frame):", sum(1 for r in rows if r["frame_pred"] == "0"), "| frame_pred 1:", sum(1 for r in rows if r["frame_pred"] == "1"))
print("distinct V (round 9, as the scorer prints):", len(set(round(x, 9) for x in v)))
print("distinct V (strings as stored):", len(set(r["v_skip"] for r in rows)))
print("V min %.4f max %.4f" % (min(v), max(v)))