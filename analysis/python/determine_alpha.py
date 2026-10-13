#!/usr/bin/env python3
"""
Calibrate exponential growth rate alpha from core-gene diversity.

R = pi / theta_W is free of theta and depends only on alpha and n.
    Simulate a grid of alpha, interpolate R(alpha) = R_obs in log10(1+alpha).

Scaling conventions, which the AG simulation script must match:
    ploidy = 1, initial_size = 1, so one generation = one coalescent unit,
    alpha is per coalescent unit.

Usage:
    python determine_alpha_and_theta.py --n 260 --pi-obs 0.012 --thetaw-obs 0.020 \
        --reps 5000 --out alpha_theta_curve.csv --alpha-out alpha_hat.csv
"""

import argparse

import msprime
import numpy as np


def sim_alpha(alpha, n, reps, seed):
    """Return (R, se_R, pi_branch) at a single alpha.

    R is the mean of per-replicate ratios
    """
    dem = msprime.Demography()
    dem.add_population(name="pop", initial_size=1.0, growth_rate=alpha)
    reps_iter = msprime.sim_ancestry(
        samples=n, demography=dem, ploidy=1, sequence_length=1,
        recombination_rate=0, num_replicates=reps, random_seed=seed,
    )

    i = np.arange(1, n)
    w = 2.0 * i * (n - i) / (n * (n - 1.0))
    a_n = np.sum(1.0 / i)

    num = np.empty(reps)
    den = np.empty(reps)
    for j, ts in enumerate(reps_iter):
        afs = ts.allele_frequency_spectrum(mode="branch", polarised=True, span_normalise=True)
        num[j] = np.sum(w * afs[1:n])
        den[j] = afs[1:n].sum()

    with np.errstate(divide="ignore", invalid="ignore"):
        vals = np.where(den > 0, a_n * num / den, np.nan)

    R = float(np.nanmean(vals))
    se_R = float(np.nanstd(vals, ddof=1) / np.sqrt(np.sum(~np.isnan(vals))))
    pi_branch = float(num.mean())
    return R, se_R, pi_branch


def interp_alpha(alphas, Rs, target):
    """Interpolate R(alpha) = target in log10(1 + alpha). Returns nan if
    the target is not bracketed by the grid."""
    alphas = np.asarray(alphas, dtype=float)
    Rs = np.asarray(Rs, dtype=float)
    if target > Rs.max() or target < Rs.min():
        return float("nan")
    x = np.log10(1.0 + alphas)
    o = np.argsort(Rs)              # np.interp requires increasing x argument
    return float(10.0 ** np.interp(target, Rs[o], x[o]) - 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=260)
    ap.add_argument("--pi-obs", type=float, required=True,
                    help="observed core-gene pi_S, e.g. 0.012")
    ap.add_argument("--thetaw-obs", type=float, required=True,
                    help="observed core-gene theta_W, e.g. 0.020")
    ap.add_argument("--reps", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--alphas", type=float, nargs="+",
                    default=[0, 1, 2, 2.5, 3, 3.5, 4, 4.5,
                             5, 10, 15, 20, 30, 50, 75,
                             100, 150, 200, 350, 500])
    ap.add_argument("--out", type=str, default="alpha_theta_curve.csv")
    ap.add_argument("--alpha-out", type=str, default="alpha_hat.csv")
    args = ap.parse_args()

    target = args.pi_obs / args.thetaw_obs

    print(f"n = {args.n}, pi_obs = {args.pi_obs:g}, "
          f"theta_W_obs = {args.thetaw_obs:g}, R_obs = {target:.4f}")
    print(f"{args.reps} replicates per alpha, seed {args.seed}\n")
    print(f"{'alpha':>8}  {'R':>8}  {'se_R':>8}  {'ci_lo':>8}  "
          f"{'ci_hi':>8}  {'pi_branch':>10}")
    print("-" * 60)

    rows = []
    for a in args.alphas:
        R, se, pib = sim_alpha(a, args.n, args.reps, args.seed)
        rows.append((a, R, se, R - 1.96 * se, R + 1.96 * se, pib))
        print(f"{a:>8.4g}  {R:>8.4f}  {se:>8.5f}  {R - 1.96*se:>8.4f}  "
              f"{R + 1.96*se:>8.4f}  {pib:>10.4f}")

    # ---- validation at alpha = 0 -------------------------------------------
    if args.alphas and args.alphas[0] == 0:
        R0, _, pib0 = rows[0][1], rows[0][2], rows[0][5]
        print(f"\nCheck at alpha = 0: R = {R0:.4f} (expect 1.00), "
              f"pi_branch = {pib0:.4f} (expect 2.00)")
        if abs(pib0 - 2.0) > 0.05:
            print("WARNING: constant-size check failed, "
                  "the scaling convention is wrong. Stop here.")

    # ---- step 1, alpha -----------------------------------------------------
    alphas = [r[0] for r in rows]
    Rs = [r[1] for r in rows]
    a_hat = interp_alpha(alphas, Rs, target)

    if not np.isfinite(a_hat):
        print(f"\nR_obs = {target:.4f} is not bracketed by the grid "
              f"(simulated range {min(Rs):.4f} to {max(Rs):.4f}). "
              f"Widen --alphas.")
    else:
        print(f"\nalpha_hat = {a_hat:.4g}")
        with open(args.alpha_out, "w") as fh:
            fh.write("pi_obs,thetaw_obs,alpha_hat\n")
            fh.write(f"{args.pi_obs:g},{args.thetaw_obs:g},{a_hat:.6g}\n")
        print(f"Wrote {args.alpha_out}")        

    # ---- write the curve ---------------------------------------------------
    header = "alpha,log10_1p_alpha,R,se_R,ci_lo,ci_hi,pi_branch"
    with open(args.out, "w") as fh:
        fh.write(header + "\n")
        for a, R, se, lo, hi, pib in rows:
            fh.write(f"{a:g},{np.log10(1.0 + a):.6f},{R:.6f},"
                     f"{se:.8f},{lo:.6f},{hi:.6f},{pib:.6f}\n")
    print(f"\nWrote {args.out}")
    print(f"R_obs = {target:.6f}  (add as a horizontal line when plotting)")


if __name__ == "__main__":
    main()
