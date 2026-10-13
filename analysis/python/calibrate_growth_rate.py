#!/usr/bin/env python3
"""
Calibrate the exponential growth rate alpha against an observed pi / theta_W.

Both estimators are functionals of the branch-length site frequency spectrum
L_i = total length of branches subtending exactly i tips (i = 1 .. n-1):

    theta_W  proportional to  sum_i L_i / a_n,     a_n = sum_{i=1}^{n-1} 1/i
    pi       proportional to  sum_i w_i L_i,       w_i = 2 i (n-i) / (n (n-1))

so the ratio

    R = a_n * sum_i w_i L_i / sum_i L_i

has expectation 1 under a constant-sized population and decreases under
expansion. Every scale factor (mutation rate, N0, ploidy) cancels, so R is a
function of alpha and n only.

NOTE: this uses BRANCH lengths keyed by descendant count, which is NOT the
subtree total used in simulate_genealogies_msprime.py. Singletons (i = 1) are
included here and must be -- that is where the expansion signal sits.

Usage:
    python calibrate_growth_rate.py --target 0.65 --n 260 --reps 2000
"""

import argparse

import msprime
import numpy as np


def branch_sfs(ts, n):
    """L[i] = total branch length subtending exactly i tips, i = 1 .. n-1."""
    tree = ts.first()
    L = np.zeros(n, dtype=np.float64)
    k = [0] * ts.num_nodes
    for u in tree.samples():
        k[u] = 1
    root = tree.root
    for u in tree.nodes(order="postorder"):
        ku = k[u]
        for c in tree.children(u):
            ku += k[c]
        k[u] = ku
        if u != root:
            L[ku] += tree.branch_length(u)
    return L


def ratio_pi_thetaW(L, n):
    i = np.arange(1, n)
    w = 2.0 * i * (n - i) / (n * (n - 1.0))
    a_n = np.sum(1.0 / i)
    tot = L[1:n].sum()
    if tot <= 0:
        return np.nan
    return a_n * np.sum(w * L[1:n]) / tot


def mean_ratio(alpha, n, reps, seed):
    dem = msprime.Demography()
    dem.add_population(name="pop", initial_size=1.0, growth_rate=alpha)
    reps_iter = msprime.sim_ancestry(
        samples=n, demography=dem, ploidy=1, sequence_length=1,
        recombination_rate=0, num_replicates=reps, random_seed=seed,
    )
    vals = [ratio_pi_thetaW(branch_sfs(ts, n), n) for ts in reps_iter]
    vals = np.asarray(vals)
    return np.nanmean(vals), np.nanstd(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=260)
    ap.add_argument("--target", type=float, required=True,
                    help="observed pi / theta_W, e.g. 0.013/0.020 = 0.65")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--alphas", type=float, nargs="+",
                    default=[0, 1, 2, 5, 10, 20, 50, 100, 200, 500])
    args = ap.parse_args()

    print(f"n = {args.n}, target R = {args.target:.3f}, "
          f"{args.reps} replicates per alpha\n")
    print(f"{'alpha':>8}  {'mean R':>8}  {'sd':>7}")
    print("-" * 27)

    results = []
    for a in args.alphas:
        m, s = mean_ratio(a, args.n, args.reps, args.seed)
        results.append((a, m))
        flag = "  <-- closest so far" if (
            not results[:-1] or
            abs(m - args.target) < min(abs(x[1] - args.target)
                                       for x in results[:-1])
        ) else ""
        print(f"{a:>8.4g}  {m:>8.4f}  {s:>7.4f}{flag}")

    # linear interpolation between the two bracketing alphas
    arr = sorted(results, key=lambda x: x[0])
    for (a0, r0), (a1, r1) in zip(arr, arr[1:]):
        if (r0 - args.target) * (r1 - args.target) <= 0 and r0 != r1:
            a_hat = a0 + (args.target - r0) * (a1 - a0) / (r1 - r0)
            print(f"\nTarget bracketed between alpha = {a0:g} and {a1:g}; "
                  f"interpolated alpha ~ {a_hat:.3g}")
            print("Re-run with a finer grid around this value to refine.")
            break
    else:
        print("\nTarget not bracketed -- widen --alphas.")


if __name__ == "__main__":
    main()
