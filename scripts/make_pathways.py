# #!/usr/bin/env python
# import argparse, json
# import gseapy as gp
# p = argparse.ArgumentParser()
# p.add_argument("--out", default="data/pathways.json")
# p.add_argument("--library", default="Hallmark")
# p.add_argument("--organism", default="Human")
# args = p.parse_args()
# gmt = gp.get_library(args.library, organism=args.organism)
# with open(args.out,"w") as f: json.dump(gmt,f,indent=2)
# print(f"Wrote {args.out} with {len(gmt)} pathways")

#!/usr/bin/env python
# Robust MSigDB/Enrichr library fetcher: auto-detects a Hallmark library for the given organism.
import argparse, json, re, sys
try:
    import gseapy as gp
except Exception as e:
    print("Please install gseapy (pip install gseapy)", file=sys.stderr)
    raise

p = argparse.ArgumentParser()
p.add_argument("--out", default="data/pathways.json")
p.add_argument("--organism", default="Human")
p.add_argument("--prefer", default="hallmark", help="substring to select library (case-insensitive)")
args = p.parse_args()

# List available libraries for this organism
libs = gp.get_library_name(organism=args.organism)

# Prefer Hallmark; pick the newest year if multiple
cands = [l for l in libs if args.prefer.lower() in l.lower()]
def year_key(l):
    m = re.search(r"(20\d{2})", l)
    return int(m.group(1)) if m else -1

if cands:
    lib = sorted(cands, key=year_key)[-1]
else:
    # Fallback to Reactome, then KEGG
    fallback = next((l for l in libs if "reactome" in l.lower()), None) \
           or next((l for l in libs if "kegg" in l.lower()), None)
    if not fallback:
        # Last-resort message with a peek at options
        sample = ", ".join(libs[:10])
        raise SystemExit(f"No library containing '{args.prefer}' found for organism {args.organism}. "
                         f"Example libraries: {sample} ...")
    lib = fallback

# Fetch and write
gmt = gp.get_library(lib, organism=args.organism)
with open(args.out, "w") as f:
    json.dump(gmt, f, indent=2)
print(f"Wrote {args.out} using library '{lib}' with {len(gmt)} pathways")
