#!/usr/bin/env python3
"""
Subtree tip count (k) and summed subtree branch length, per internal node,
under an exponentially growing population.

Direct analogue of the ape::rcoal script: the root is dropped so k < n always,
and singletons (k = 1) are excluded.

Units
-----
initial_size = 1.0 with ploidy = 1 puts time in coalescent units, matching
rcoal exactly (pairwise coalescence rate = 1/(ploidy * N) = 1 per unit time).
So `growth_rate` here is alpha in COALESCENT units, i.e. the dimensionless
product alpha_per_gen * N0 for a haploid population. Set growth_rate = 0.0 to
reproduce rcoal.
"""

import argparse
import csv
import os
import sys
import time

import msprime


def subtree_stats(ts, n_samples):
    """Return (node, k, subtree_branch_length) for internal non-root nodes.

    Identical algorithm to the R post-order version: an edge's length is
    credited to its PARENT's subtree, so the stem branch above a clade is
    never counted in that clade's total.
    """
    tree = ts.first()
    num_nodes = ts.num_nodes

    k = [0] * num_nodes
    bl = [0.0] * num_nodes
    for u in tree.samples():
        k[u] = 1

    out = []
    root = tree.root
    for u in tree.nodes(order="postorder"):
        ku = k[u]
        blu = bl[u]
        for c in tree.children(u):
            ku += k[c]
            blu += bl[c] + tree.branch_length(c)
        k[u] = ku
        bl[u] = blu
        # keep internal, non-root nodes only -> 2 <= k < n_samples
        if ku >= 2 and u != root:
            out.append((u, ku, blu))
    return out


def build_demography(initial_size, growth_rate):
    dem = msprime.Demography()
    dem.add_population(name="pop", initial_size=initial_size,
                       growth_rate=growth_rate)
    return dem


def simulate(n_genomes, n_trees, out_file, growth_rate, initial_size=1.0,
             chunk_size=10000, seed=42):

    dem = build_demography(initial_size, growth_rate)

    with open(out_file, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["internal_node", "freq", "subtree_br_len", "tree_id"])

        replicates = msprime.sim_ancestry(
            samples=n_genomes,
            demography=dem,
            ploidy=1,               # haploid: one sample node per genome
            sequence_length=1,
            recombination_rate=0,   # exactly one tree per replicate
            num_replicates=n_trees,
            random_seed=seed,
        )

        buf = []
        t0 = time.time()
        for tree_id, ts in enumerate(replicates, start=1):
            for node, kk, blen in subtree_stats(ts, n_genomes):
                buf.append((node, kk, blen, tree_id))

            if tree_id % chunk_size == 0:
                writer.writerows(buf)
                buf.clear()
                fh.flush()
                print(f"Trees 1 - {tree_id} written "
                      f"({time.time() - t0:.0f}s)", flush=True)

        if buf:
            writer.writerows(buf)

    return out_file


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-genomes", type=int, default=260)
    ap.add_argument("--n-trees", type=int, default=500000)
    ap.add_argument("--growth-rate", type=float, required=True,
                    help="alpha in coalescent units; 0 reproduces rcoal")
    ap.add_argument("--initial-size", type=float, default=1.0)
    ap.add_argument("--chunk-size", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default=os.environ.get("OUTDIR_DAT", "."))
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fname = (f"expanding_sim_{args.n_genomes}_alpha{args.growth_rate:g}"
             f"_ntree_{args.n_trees}.csv")
    out_file = os.path.join(args.outdir, fname)

    t0 = time.time()
    simulate(args.n_genomes, args.n_trees, out_file,
             growth_rate=args.growth_rate,
             initial_size=args.initial_size,
             chunk_size=args.chunk_size,
             seed=args.seed)
    print(f"Done: {out_file} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
