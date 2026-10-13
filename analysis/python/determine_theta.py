#!/usr/bin/env python3
"""
Calibrate theta_hat at a fixed growth rate alpha from core-gene diversity.

Step 2 of the two-step calibration (step 1 is determine_alpha.py).
Alpha is fixed, so the genealogy shape is set and theta only scales the
number of mutations. Simulate a grid of theta, record mean pi and theta_W,
then interpolate E[pi](theta) = pi_obs (and E[theta_W](theta) = theta_W_obs
as a check).

Scaling conventions, which must match determine_alpha.py and the AG null:
    ploidy = 1, initial_size = 1, so one generation = one coalescent unit.
    alpha is per coalescent unit (msprime growth_rate = alpha).
    theta = 2 * N0 * mu per site, so msprime mutation rate = theta / 2.
    Infinite sites (discrete_genome = False).

Under infinite sites E[pi] = (theta / 2) * pi_branch(alpha) exactly, so the
curve should be a straight line through the origin. The analytic value
2 * pi_obs / pi_branch is reported alongside the interpolated values.

New trees are simulated for every theta (independent seeds per theta).

Usage:
    python determine_theta.py --alpha 2.88 --pi-obs 0.01226 --thetaw-obs 0.019726 \
        --reps 10000 --out theta_pi_curve.csv --theta-out theta_hat.csv
"""

import argparse

import msprime
import numpy as np


def to_msprime_seed(x):
    """Map a uint32 to the valid msprime seed range 1 .. 2^31 - 1."""
    return int(x) % (2**31 - 2) + 1


def sim_theta(theta, alpha, n, reps, seq_len, anc_seed, mut_seed):
    """Return per-replicate arrays (pi, theta_W, pi_branch, tbl) at one theta.

    pi and theta_W are per site, from simulated mutations.
    pi_branch (mean pairwise branch length = 2 E[T_pair]) and tbl (total
    branch length) come from the trees and feed the analytic check.
    """
    dem = msprime.Demography()
    dem.add_population(name="pop", initial_size=1.0, growth_rate=alpha)
    reps_iter = msprime.sim_ancestry(
        samples=n, demography=dem, ploidy=1, sequence_length=seq_len,
        recombination_rate=0, num_replicates=reps, random_seed=anc_seed,
    )

    a_n = np.sum(1.0 / np.arange(1, n))
    rng = np.random.default_rng(mut_seed)
    mut_seeds = rng.integers(1, 2**31 - 1, size=reps)

    pi = np.zeros(reps)
    tw = np.zeros(reps)
    pib = np.empty(reps)
    tbl = np.empty(reps)
    for j, ts in enumerate(reps_iter):
        # one tree spanning the sequence, so span-normalised = tree value
        pib[j] = ts.diversity(mode="branch")
        tbl[j] = ts.segregating_sites(mode="branch")
        if theta > 0:
            mts = msprime.sim_mutations(
                ts, rate=theta / 2.0, discrete_genome=False,
                random_seed=int(mut_seeds[j]),
            )
            pi[j] = mts.diversity(mode="site")                  # per site
            tw[j] = mts.segregating_sites(mode="site") / a_n    # per site
    return pi, tw, pib, tbl


def interp_linear(thetas, means, target):
    """Linear interpolation of mean(theta) = target between the two
    bracketing grid points. Returns nan if the target is not bracketed."""
    thetas = np.asarray(thetas, dtype=float)
    means = np.asarray(means, dtype=float)
    if target > means.max() or target < means.min():
        return float("nan")
    o = np.argsort(means)           # np.interp requires increasing x argument
    return float(np.interp(target, means[o], thetas[o]))


def fit_quad0(thetas, means, ses, target):
    """Weighted least-squares fit of mean = b1*theta + b2*theta^2 (through
    the origin), weights 1/se^2, theta = 0 excluded. Solve for target.

    Returns (theta_hat, b1, b2). Uses the numerically stable form of the
    positive root, which reduces to target / b1 when b2 = 0.
    """
    t = np.asarray(thetas, dtype=float)
    y = np.asarray(means, dtype=float)
    s = np.asarray(ses, dtype=float)
    keep = (t > 0) & (s > 0)
    w = 1.0 / s[keep]
    X = np.column_stack([t[keep], t[keep] ** 2]) * w[:, None]
    b1, b2 = np.linalg.lstsq(X, y[keep] * w, rcond=None)[0]

    disc = b1 ** 2 + 4.0 * b2 * target
    if b1 <= 0 or disc < 0:
        return float("nan"), float(b1), float(b2)
    root = 2.0 * target / (b1 + np.sqrt(disc))
    if root <= 0 or root > t.max():
        return float("nan"), float(b1), float(b2)
    return float(root), float(b1), float(b2)


def bootstrap(thetas, pi_list, tw_list, se_pi, se_tw, pi_obs, tw_obs, B, rng):
    """Resample replicates within each theta (pi and theta_W paired) and
    re-solve all four estimates. Fit weights are held at the original se.
    Columns: linear_pi, quad_pi, linear_thetaw, quad_thetaw."""
    K = len(thetas)
    out = np.full((B, 4), np.nan)
    for b in range(B):
        mp = np.empty(K)
        mt = np.empty(K)
        for k in range(K):
            p, w = pi_list[k], tw_list[k]
            idx = rng.integers(0, p.size, p.size)
            mp[k] = p[idx].mean()
            mt[k] = w[idx].mean()
        out[b, 0] = interp_linear(thetas, mp, pi_obs)
        out[b, 1] = fit_quad0(thetas, mp, se_pi, pi_obs)[0]
        out[b, 2] = interp_linear(thetas, mt, tw_obs)
        out[b, 3] = fit_quad0(thetas, mt, se_tw, tw_obs)[0]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=260)
    ap.add_argument("--alpha", type=float, required=True,
                    help="alpha_hat from determine_alpha.py, e.g. 2.88")
    ap.add_argument("--pi-obs", type=float, required=True,
                    help="observed core-gene pi_S, e.g. 0.01226")
    ap.add_argument("--thetaw-obs", type=float, required=True,
                    help="observed core-gene theta_W, e.g. 0.019726")
    ap.add_argument("--reps", type=int, default=10000,
                    help="new genealogies per theta")
    ap.add_argument("--seq-len", type=float, default=1000.0,
                    help="sites per replicate; sets mutational noise only, "
                         "not the expectation")
    ap.add_argument("--boot", type=int, default=1000,
                    help="bootstrap resamples for confidence intervals")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--thetas", type=float, nargs="+",
                    default=list(np.round(np.arange(0.0, 0.1 + 1e-12, 0.005), 6)))
    ap.add_argument("--out", type=str, default="theta_pi_curve.csv")
    ap.add_argument("--theta-out", type=str, default="theta_hat.csv")
    ap.add_argument("--reps-out", type=str, default=None,
                    help="optional per-replicate output file")
    args = ap.parse_args()

    thetas = np.asarray(sorted(args.thetas), dtype=float)
    a_n = np.sum(1.0 / np.arange(1, args.n))

    # independent seeds per theta, plus one stream for the bootstrap
    children = np.random.SeedSequence(args.seed).spawn(len(thetas) + 1)
    boot_rng = np.random.default_rng(children[-1])

    print(f"n = {args.n}, alpha = {args.alpha:g}, pi_obs = {args.pi_obs:g}, "
          f"theta_W_obs = {args.thetaw_obs:g}")
    print(f"{args.reps} new genealogies per theta, L = {args.seq_len:g}, "
          f"seed {args.seed}\n")
    print(f"{'theta':>7}  {'mean_pi':>9}  {'se_pi':>9}  {'mean_tW':>9}  "
          f"{'se_tW':>9}  {'R_rom':>7}  {'R_mor':>7}  {'pi/E[pi]':>8}")
    print("-" * 80)

    rows, pi_list, tw_list, pib_all, tbl_all = [], [], [], [], []
    rep_fh = open(args.reps_out, "w") if args.reps_out else None
    if rep_fh:
        rep_fh.write("theta,rep,pi,thetaw,pi_branch,total_branch_length\n")

    for k, th in enumerate(thetas):
        anc_seed, mut_seed = (to_msprime_seed(x)
                              for x in children[k].generate_state(2))
        pi, tw, pib, tbl = sim_theta(th, args.alpha, args.n, args.reps,
                                     args.seq_len, anc_seed, mut_seed)
        pi_list.append(pi)
        tw_list.append(tw)
        pib_all.append(pib)
        tbl_all.append(tbl)

        m_pi, m_tw = pi.mean(), tw.mean()
        se_pi = pi.std(ddof=1) / np.sqrt(pi.size)
        se_tw = tw.std(ddof=1) / np.sqrt(tw.size)
        r_rom = m_pi / m_tw if m_tw > 0 else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            r_mor = float(np.nanmean(np.where(tw > 0, pi / tw, np.nan))) \
                if th > 0 else np.nan
        exp_pi = (th / 2.0) * pib.mean()
        check = m_pi / exp_pi if exp_pi > 0 else np.nan

        rows.append((th, m_pi, se_pi, m_tw, se_tw, r_rom, r_mor,
                     pib.mean(), check))
        print(f"{th:>7.4g}  {m_pi:>9.6f}  {se_pi:>9.6f}  {m_tw:>9.6f}  "
              f"{se_tw:>9.6f}  {r_rom:>7.4f}  {r_mor:>7.4f}  {check:>8.4f}")

        if rep_fh:
            for j in range(pi.size):
                rep_fh.write(f"{th:g},{j + 1},{pi[j]:.8g},{tw[j]:.8g},"
                             f"{pib[j]:.8g},{tbl[j]:.8g}\n")

    if rep_fh:
        rep_fh.close()
        print(f"\nWrote {args.reps_out}")

    # ---- check at theta = 0 ------------------------------------------------
    if thetas[0] == 0 and rows[0][1] != 0:
        print("WARNING: non-zero pi at theta = 0. Stop here.")

    # ---- estimates ---------------------------------------------------------
    m_pi = np.array([r[1] for r in rows])
    se_pi = np.array([r[2] for r in rows])
    m_tw = np.array([r[3] for r in rows])
    se_tw = np.array([r[4] for r in rows])

    lin_pi = interp_linear(thetas, m_pi, args.pi_obs)
    quad_pi, b1_pi, b2_pi = fit_quad0(thetas, m_pi, se_pi, args.pi_obs)
    lin_tw = interp_linear(thetas, m_tw, args.thetaw_obs)
    quad_tw, b1_tw, b2_tw = fit_quad0(thetas, m_tw, se_tw, args.thetaw_obs)

    pib_mean = np.concatenate(pib_all).mean()
    tbl_mean = np.concatenate(tbl_all).mean()
    an_pi = 2.0 * args.pi_obs / pib_mean
    an_tw = 2.0 * a_n * args.thetaw_obs / tbl_mean

    print(f"\nRunning {args.boot} bootstrap resamples ...")
    boot = bootstrap(thetas, pi_list, tw_list, se_pi, se_tw,
                     args.pi_obs, args.thetaw_obs, args.boot, boot_rng)
    lo = np.nanpercentile(boot, 2.5, axis=0)
    hi = np.nanpercentile(boot, 97.5, axis=0)

    print(f"\npi_branch(alpha) = {pib_mean:.4f}   "
          f"(pooled over {len(thetas) * args.reps} trees)")
    print(f"Quadratic fit, pi:      b1 = {b1_pi:.5f}  b2 = {b2_pi:.5f}  "
          f"(expect b1 = pi_branch/2 = {pib_mean / 2:.5f}, b2 = 0)")
    print(f"Quadratic fit, theta_W: b1 = {b1_tw:.5f}  b2 = {b2_tw:.5f}\n")

    results = [
        ("linear_pi",       args.pi_obs,     lin_pi,  lo[0], hi[0]),
        ("quadratic_pi",    args.pi_obs,     quad_pi, lo[1], hi[1]),
        ("analytic_pi",     args.pi_obs,     an_pi,   np.nan, np.nan),
        ("linear_thetaw",   args.thetaw_obs, lin_tw,  lo[2], hi[2]),
        ("quadratic_thetaw", args.thetaw_obs, quad_tw, lo[3], hi[3]),
        ("analytic_thetaw", args.thetaw_obs, an_tw,   np.nan, np.nan),
    ]
    print(f"{'method':>17}  {'target':>9}  {'theta_hat':>9}  "
          f"{'ci_lo':>9}  {'ci_hi':>9}")
    print("-" * 62)
    for name, tgt, est, l, h in results:
        print(f"{name:>17}  {tgt:>9.6f}  {est:>9.6f}  {l:>9.6f}  {h:>9.6f}")

    with open(args.theta_out, "w") as fh:
        fh.write("method,alpha,target,theta_hat,ci_lo,ci_hi\n")
        for name, tgt, est, l, h in results:
            fh.write(f"{name},{args.alpha:g},{tgt:g},{est:.6g},"
                     f"{l:.6g},{h:.6g}\n")
    print(f"\nWrote {args.theta_out}")

    # ---- write the curve ---------------------------------------------------
    header = ("theta,mean_pi,se_pi,ci_lo_pi,ci_hi_pi,mean_thetaw,se_thetaw,"
              "R_ratio_of_means,R_mean_of_ratios,mean_pi_branch,pi_over_expected")
    with open(args.out, "w") as fh:
        fh.write(header + "\n")
        for th, mp, sp, mt, st, rr, rm, pb, ck in rows:
            fh.write(f"{th:g},{mp:.8f},{sp:.8f},{mp - 1.96 * sp:.8f},"
                     f"{mp + 1.96 * sp:.8f},{mt:.8f},{st:.8f},{rr:.6f},"
                     f"{rm:.6f},{pb:.6f},{ck:.6f}\n")
    print(f"Wrote {args.out}")
    print(f"pi_obs = {args.pi_obs:.6f}  (add as a horizontal line when plotting)")


if __name__ == "__main__":
    main()
