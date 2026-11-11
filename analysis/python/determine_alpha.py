#!/usr/bin/env python3
"""
Calibrate exponential growth rate alpha from core-gene diversity.

R = pi / theta_W is free of theta and depends only on alpha and n.
    Simulate a grid of alpha, interpolate R(alpha) = R_obs in log10(1+alpha).

Scaling conventions, which the AG simulation script must match:
    ploidy = 1, initial_size = 1, so one generation = one coalescent unit,
    alpha is per coalescent unit.

Uncertainty on alpha_hat is obtained by resampling the per-replicate R values
within each alpha and re-interpolating (same scheme as determine_theta.py).
Simulation noise alone is propagated by default. Pass --r-obs-se to also
propagate uncertainty in the observed ratio, which is usually the larger term.

Usage:
    python determine_alpha.py --n 260 --pi-obs 0.012 --thetaw-obs 0.020 \
        --reps 5000 --boot 1000 --out alpha_theta_curve.csv --alpha-out alpha_hat.csv
"""

import argparse

import msprime
import numpy as np


def to_msprime_seed(x):
    """Map a uint32 to the valid msprime seed range 1 .. 2^31 - 1."""
    return int(x) % (2**31 - 2) + 1


def sim_alpha(alpha, n, reps, seed):
    """Return (vals, pi_branch) at a single alpha.

    vals is the per-replicate ratio a_n * pi_branch / total_branch_length,
    i.e. one R per genealogy. Replicates with no branches are nan.
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

    pi_branch = float(num.mean())
    return vals, pi_branch


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


def bootstrap(alphas, val_list, target, B, rng, r_obs_se=None):
    """Resample replicates within each alpha and re-solve alpha_hat.

    R_obs is held fixed unless r_obs_se is given, in which case it is drawn
    from N(target, r_obs_se) on each resample so that uncertainty in the
    observed ratio is propagated as well.
    """
    K = len(alphas)
    out = np.full(B, np.nan)
    for b in range(B):
        mR = np.empty(K)
        for k in range(K):
            v = val_list[k]
            idx = rng.integers(0, v.size, v.size)
            mR[k] = np.nanmean(v[idx])
        tgt = target if r_obs_se is None else rng.normal(target, r_obs_se)
        out[b] = interp_alpha(alphas, mR, tgt)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=260)
    ap.add_argument("--pi-obs", type=float, required=True,
                    help="observed core-gene pi_S, e.g. 0.012")
    ap.add_argument("--thetaw-obs", type=float, required=True,
                    help="observed core-gene theta_W, e.g. 0.020")
    ap.add_argument("--reps", type=int, default=5000)
    ap.add_argument("--boot", type=int, default=1000,
                    help="bootstrap resamples for the confidence interval")
    ap.add_argument("--r-obs-se", type=float, default=None,
                    help="standard error of the observed pi_S / theta_W ratio; "
                         "if given, propagated into the interval")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--alphas", type=float, nargs="+",
                    default=[0, 1, 2, 2.5, 3, 3.5, 4, 4.5,
                             5, 10, 15, 20, 30, 50, 75,
                             100, 150, 200, 350, 500])
    ap.add_argument("--out", type=str, default="alpha_theta_curve.csv")
    ap.add_argument("--alpha-out", type=str, default="alpha_hat.csv")
    ap.add_argument("--reps-out", type=str, default=None,
                    help="optional per-replicate output file")
    args = ap.parse_args()

    target = args.pi_obs / args.thetaw_obs
    alphas = np.asarray(args.alphas, dtype=float)

    # independent seeds per alpha, plus one stream for the bootstrap
    children = np.random.SeedSequence(args.seed).spawn(len(alphas) + 1)
    boot_rng = np.random.default_rng(children[-1])

    print(f"n = {args.n}, pi_obs = {args.pi_obs:g}, "
          f"theta_W_obs = {args.thetaw_obs:g}, R_obs = {target:.4f}")
    if args.r_obs_se is not None:
        print(f"se(R_obs) = {args.r_obs_se:g}, propagated into the interval")
    print(f"{args.reps} replicates per alpha, seed {args.seed}\n")
    print(f"{'alpha':>8}  {'R':>8}  {'se_R':>8}  {'ci_lo':>8}  "
          f"{'ci_hi':>8}  {'pi_branch':>10}")
    print("-" * 60)

    rows, val_list = [], []
    rep_fh = open(args.reps_out, "w") if args.reps_out else None
    if rep_fh:
        rep_fh.write("alpha,rep,R\n")

    for k, a in enumerate(alphas):
        anc_seed = to_msprime_seed(children[k].generate_state(1)[0])
        vals, pib = sim_alpha(a, args.n, args.reps, anc_seed)
        val_list.append(vals)

        R = float(np.nanmean(vals))
        se = float(np.nanstd(vals, ddof=1) / np.sqrt(np.sum(~np.isnan(vals))))
        rows.append((a, R, se, R - 1.96 * se, R + 1.96 * se, pib))
        print(f"{a:>8.4g}  {R:>8.4f}  {se:>8.5f}  {R - 1.96*se:>8.4f}  "
              f"{R + 1.96*se:>8.4f}  {pib:>10.4f}")

        if rep_fh:
            for j in range(vals.size):
                rep_fh.write(f"{a:g},{j + 1},{vals[j]:.8g}\n")

    if rep_fh:
        rep_fh.close()
        print(f"\nWrote {args.reps_out}")

    # ---- validation at alpha = 0 -------------------------------------------
    if alphas[0] == 0:
        R0, pib0 = rows[0][1], rows[0][5]
        print(f"\nCheck at alpha = 0: R = {R0:.4f} (expect 1.00), "
              f"pi_branch = {pib0:.4f} (expect 2.00)")
        if abs(pib0 - 2.0) > 0.05:
            print("WARNING: constant-size check failed, "
                  "the scaling convention is wrong. Stop here.")

    # ---- step 1, alpha -----------------------------------------------------
    Rs = [r[1] for r in rows]
    a_hat = interp_alpha(alphas, Rs, target)

    if not np.isfinite(a_hat):
        print(f"\nR_obs = {target:.4f} is not bracketed by the grid "
              f"(simulated range {min(Rs):.4f} to {max(Rs):.4f}). "
              f"Widen --alphas.")
        ci_lo = ci_hi = float("nan")
        n_ok = 0
    else:
        print(f"\nRunning {args.boot} bootstrap resamples ...")
        boot = bootstrap(alphas, val_list, target, args.boot, boot_rng,
                         r_obs_se=args.r_obs_se)
        n_ok = int(np.sum(np.isfinite(boot)))
        ci_lo = float(np.nanpercentile(boot, 2.5))
        ci_hi = float(np.nanpercentile(boot, 97.5))
        if n_ok < args.boot:
            print(f"NOTE: {args.boot - n_ok} of {args.boot} resamples fell "
                  f"outside the alpha grid and were dropped. Widen --alphas "
                  f"if this is a large fraction.")
        print(f"\nalpha_hat = {a_hat:.4g}  "
              f"[{ci_lo:.4g}, {ci_hi:.4g}]  (95% bootstrap)")

    with open(args.alpha_out, "w") as fh:
        fh.write("method,pi_obs,thetaw_obs,R_obs,alpha_hat,ci_lo,ci_hi,"
                 "n_boot,n_boot_ok\n")
        fh.write(f"linear_R,{args.pi_obs:g},{args.thetaw_obs:g},{target:g},"
                 f"{a_hat:.6g},{ci_lo:.6g},{ci_hi:.6g},"
                 f"{args.boot:d},{n_ok:d}\n")
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
