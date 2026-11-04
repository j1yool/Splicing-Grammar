"""clinvar_vs_fetch_header_v1.py (Day 100, Splicing Grammar)
Downloads variant_summary.txt.gz, records date/size/sha256, prints ONLY the header line,
and checks the DMD gene span against Ensembl REST. No variant row is read or printed."""
import gzip, hashlib, json, os, sys, urllib.request

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
OUT = "variant_summary.txt.gz"
ENS = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/DMD?content-type=application/json"
EXPECT = (31097677, 33339609)  # from dmd_scorer_v1_build.json (exon-table span)

if os.path.exists(OUT):
    print("STOP: " + OUT + " exists; move it first.")
    sys.exit(1)

req = urllib.request.Request(URL, headers={"User-Agent": "Splicing-Grammar-fetch"})
with urllib.request.urlopen(req) as r, open(OUT, "wb") as f:
    last_mod = r.headers.get("Last-Modified")
    h = hashlib.sha256()
    n = 0
    while True:
        chunk = r.read(1 << 20)
        if not chunk:
            break
        f.write(chunk)
        h.update(chunk)
        n += len(chunk)
print("file:", OUT)
print("server Last-Modified:", last_mod)
print("bytes:", n)
print("sha256:", h.hexdigest())

with gzip.open(OUT, "rt", encoding="utf-8", newline="") as g:
    header = g.readline().rstrip("\r\n")
cols = header.split("\t")
print("HEADER (" + str(len(cols)) + " columns):")
for i, c in enumerate(cols):
    print(i, repr(c))

with urllib.request.urlopen(urllib.request.Request(ENS, headers={"User-Agent": "Splicing-Grammar-fetch"})) as r:
    d = json.loads(r.read().decode("utf-8"))
print("ENSEMBL:", {k: d.get(k) for k in ("id", "seq_region_name", "start", "end", "strand", "assembly_name")})
ok = (d.get("start"), d.get("end")) == EXPECT
print("DMD span matches the exon table span:", ok, "| exon table:", EXPECT)