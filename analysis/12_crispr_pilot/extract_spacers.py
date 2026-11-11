#!/usr/bin/env python3
"""
De novo CRISPR array extraction around the two type I-E cas loci.

No external dependencies: minced is not available, and the arrays are easy
enough to find directly. Within a window either side of the cas operon we
look for the most frequent k-mer (k = 28-30) occurring at least three times
with regular spacing; in a CRISPR array that k-mer is the direct repeat.
The sequences between consecutive repeats are the spacers.

Writes one FASTA of spacers per locus plus a per-genome summary.
"""
import csv, json, os, sys
from collections import Counter, defaultdict

ROOT   = os.path.expanduser("~/mnt/kpne_ags")
FASTA  = os.path.join(ROOT, "input_data/kpne_260_chr_fasta")
ANNO   = os.path.join(ROOT, "output/data/all_pirate_anno_full.csv")
OUT    = os.path.join(ROOT, "analysis/12_crispr_pilot/out")

WIN        = 25000      # bp searched either side of the cas operon
K_RANGE    = (28, 31)   # candidate direct-repeat lengths
MIN_COPIES = 3          # a repeat must occur at least this often
SPACER_LEN = (25, 50)   # plausible spacer length

LOCI = {
 "A": {"g002579","g003424","g004669","g004793","g005291","g005762","g005961","g006512"},
 "B": {"g003179","g004638","g004752","g005162","g005653","g005811","g006892","g014756"},
}

def read_fasta(p):
    seqs, name, buf = {}, None, []
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def find_array(seq):
    """Most frequent regularly spaced k-mer = the direct repeat."""
    best = None
    for k in range(*K_RANGE):
        pos = defaultdict(list)
        for i in range(len(seq) - k + 1):
            pos[seq[i:i+k]].append(i)
        for rep, ps in pos.items():
            if len(ps) < MIN_COPIES or len(set(rep)) < 3: continue
            gaps = [b - a - k for a, b in zip(ps, ps[1:])]
            good = [g for g in gaps if SPACER_LEN[0] <= g <= SPACER_LEN[1]]
            if len(good) < MIN_COPIES - 1: continue
            score = len(good)
            if best is None or score > best[0]:
                best = (score, rep, k, ps)
    return best

def main():
    os.makedirs(OUT, exist_ok=True)
    want = set().union(*LOCI.values())

    # cas gene coordinates per genome
    span = defaultdict(lambda: defaultdict(lambda: [None, None, None]))
    csv.field_size_limit(10**9)
    with open(ANNO, newline="") as f:
        rd = csv.reader(f); h = next(rd)
        I = {c: h.index(c) for c in ("geno_id","gene_family","seqnames","start","end")}
        for row in rd:
            gf = row[I["gene_family"]]
            if gf not in want: continue
            L = "A" if gf in LOCI["A"] else "B"
            g = row[I["geno_id"]]
            s, e = int(row[I["start"]]), int(row[I["end"]])
            rec = span[L][g]
            rec[0] = row[I["seqnames"]]
            rec[1] = s if rec[1] is None else min(rec[1], s)
            rec[2] = e if rec[2] is None else max(rec[2], e)

    summary = []
    for L in ("A", "B"):
        out_fa = open(os.path.join(OUT, f"spacers_locus_{L}.fasta"), "w")
        n_sp = 0
        for g, (ctg, lo, hi) in sorted(span[L].items()):
            p = os.path.join(FASTA, g + ".fasta")
            if not os.path.exists(p): continue
            seqs = read_fasta(p)
            s = seqs.get(ctg) or max(seqs.values(), key=len)
            w = s[max(0, lo - WIN): hi + WIN].upper()
            hit = find_array(w)
            if not hit:
                summary.append((L, g, 0, "", 0)); continue
            _, rep, k, ps = hit
            sp = []
            for a, b in zip(ps, ps[1:]):
                frag = w[a + k: b]
                if SPACER_LEN[0] <= len(frag) <= SPACER_LEN[1]: sp.append(frag)
            for i, x in enumerate(sp, 1):
                out_fa.write(f">{L}|{g}|sp{i}\n{x}\n"); n_sp += 1
            summary.append((L, g, len(sp), rep, k))
        out_fa.close()
        print(f"locus {L}: {n_sp} spacers from {sum(1 for r in summary if r[0]==L and r[2])} genomes")

    with open(os.path.join(OUT, "array_summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["locus","geno_id","n_spacers","direct_repeat","repeat_len"])
        w.writerows(summary)
    print("wrote", OUT)

if __name__ == "__main__":
    main()
