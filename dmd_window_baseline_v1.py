# dmd_window_baseline_v1.py  EXPLORATORY (R6): no prereg, no self-test, descriptive only.
import hashlib, json
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

SRC = "dmd_eval_v1_variants.tsv"
sha = hashlib.sha256(open(SRC, "rb").read()).hexdigest()
d = pd.read_csv(SRC, sep="\t")
s = d[d.status == "scored"].reset_index(drop=True)
assert len(s) == 1084 and int(s.label.sum()) == 386 and s.exon_id.nunique() == 77, "scored set differs from dmd_eval_v1_results.json"

rate = s.groupby("window").label.mean()
auc_insample = roc_auc_score(s.label, s.window.map(rate))

pred = np.empty(len(s))
for e in s.exon_id.unique():
    m = (s.exon_id == e).values
    r = s[~m].groupby("window").label.mean()
    pred[m] = s.loc[m, "window"].map(r).fillna(s[~m].label.mean()).values
auc_loeo = roc_auc_score(s.label, pred)

rng = np.random.default_rng(101)
exons = s.exon_id.unique()
groups = {e: np.where(s.exon_id.values == e)[0] for e in exons}
aucs = []
for _ in range(2000):
    pick = rng.choice(exons, len(exons), replace=True)
    idx = np.concatenate([groups[e] for e in pick])
    y = s.label.values[idx]
    if y.min() == y.max():
        continue
    aucs.append(roc_auc_score(y, pred[idx]))
lo, hi = np.percentile(aucs, [2.5, 97.5])

out = {"source": SRC, "source_sha256": sha, "n": len(s),
       "window_rates": {k: float(v) for k, v in rate.items()},
       "auc_insample": float(auc_insample), "auc_loeo": float(auc_loeo),
       "loeo_ci95": [float(lo), float(hi)], "boots_used": len(aucs), "seed": 101,
       "label": "EXPLORATORY"}
json.dump(out, open("dmd_window_baseline_v1_results.json", "w"), indent=2)
print(json.dumps(out, indent=2))