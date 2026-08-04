"""
detection_index.py -- does the DETECTION INDEX improve by the predicted factor T?

The lens law of "Matching the Lens to the Zero" is a statement about the exponential
RATE:  log|q_a(rho)| is maximal at d = |rho - 1/2|, and exceeds its value at the
classical lens d = 1/2 by a factor T.  What a search actually cares about is the
INDEX n at which the coefficient first goes negative, and the two are not the same
quantity: a single earlier experiment gave a gain of 30.89 where the rate predicts
30.01, a 3% residue the paper attributes but does not explain.

This script measures the index directly, at several heights, so the residue can be
seen as a function of T rather than at one point.

    S_n(a) = sum_rho [ 1 - q_a(rho)^n ],   q_a(rho) = (rho - a)/(rho - (1-a))

Zeros on the critical line give q = e^{i theta}, theta = 2 arctan(gamma/d) - pi, and
contribute 1 - cos(n theta) >= 0.  A planted quartet at 1/2 +- beta +- iT contributes
a term growing like R^n with R = |q_a| > 1, which eventually drives the sum negative.

Usage
    python detection_index.py                          # default sweep
    python detection_index.py --T 10 30 100 --beta 0.1
    ZETA_ZEROS=/path/to/zeros.txt python detection_index.py

Without a zero file it falls back to mpmath.zetazero, which is slow; supply a list.
"""
from __future__ import annotations
import argparse
import os
import numpy as np

TWO_PI = 2.0 * np.pi


# ----------------------------------------------------------------- the zeros
def load_zeros(n_wanted: int) -> np.ndarray:
    path = os.environ.get("ZETA_ZEROS")
    if path and os.path.exists(path):
        g = np.loadtxt(path)
        return g[:n_wanted]
    import mpmath as mp
    mp.mp.dps = 15
    return np.array([float(mp.im(mp.zetazero(k))) for k in range(1, n_wanted + 1)])


# ------------------------------------------------------------------ the sum
def theta_online(gam: np.ndarray, d: float) -> np.ndarray:
    """q_a = e^{i theta} for a zero at 1/2 + i gamma; d = 1/2 - a."""
    return 2.0 * np.arctan(gam / d) - np.pi


def quartet_terms(beta: float, T: float, d: float):
    """the four planted zeros 1/2 +- beta +- iT, as (log|q|, arg q) pairs."""
    out = []
    for sb in (+1.0, -1.0):
        for sT in (+1.0, -1.0):
            rho = 0.5 + sb * beta + 1j * sT * T
            q = (rho - (0.5 - d)) / (rho - (0.5 + d))
            out.append((np.log(abs(q)), np.angle(q)))
    return out


def S_of_n(ns: np.ndarray, th: np.ndarray, quart, chunk: int = 4000) -> np.ndarray:
    """Re S_n for each n in ns.  Online zeros contribute 1 - cos(n theta) each
    (both signs of gamma), the quartet contributes 1 - |q|^n cos(n arg q)."""
    out = np.empty(len(ns), dtype=float)
    for i0 in range(0, len(ns), chunk):
        nn = ns[i0:i0 + chunk].astype(float)[:, None]           # (k,1)
        # online, doubled for +-gamma
        out[i0:i0 + chunk] = 2.0 * np.sum(1.0 - np.cos(nn * th[None, :]), axis=1)
        for lg, ph in quart:
            out[i0:i0 + chunk] += 1.0 - np.exp(nn[:, 0] * lg) * np.cos(nn[:, 0] * ph)
    return out


def first_negative(th, quart, n_max: int = 400_000) -> int | None:
    """First n with Re S_n < 0.  Brackets geometrically, then scans the decade."""
    lo = 8
    while lo < n_max:
        hi = min(lo * 2, n_max)
        ns = np.arange(lo, hi, max(1, (hi - lo) // 400))
        if np.any(S_of_n(ns, th, quart) < 0):
            fine = np.arange(lo, hi)
            vals = S_of_n(fine, th, quart)
            k = np.argmax(vals < 0)
            return int(fine[k]) if vals[k] < 0 else None
        lo = hi
    return None


# ----------------------------------------------------------------- the sweep
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, nargs="*", default=[10.0, 30.0, 100.0])
    ap.add_argument("--beta", type=float, default=0.1)
    ap.add_argument("--nzeros", type=int, default=2000)
    ap.add_argument("--nmax", type=int, default=400_000)
    args = ap.parse_args()

    gam = load_zeros(args.nzeros)
    print(f"detection index vs the predicted factor-T gain")
    print(f"  {len(gam)} critical-line zeros, gamma up to {gam[-1]:.1f}; "
          f"planted quartet at 1/2 +- {args.beta} +- iT\n")
    print(f"  {'T':>6} {'n (classical)':>14} {'n (matched)':>12} {'measured gain':>14} "
          f"{'predicted T':>12} {'ratio':>8}")

    for T in args.T:
        R = np.hypot(args.beta, T)
        rows = {}
        for label, d in (("classical", 0.5), ("matched", R)):
            th = theta_online(gam, d)
            rows[label] = first_negative(th, quartet_terms(args.beta, T, d), args.nmax)
        nc, nm = rows["classical"], rows["matched"]
        if nc is None or nm is None:
            print(f"  {T:>6.0f} {str(nc):>14} {str(nm):>12} "
                  f"{'(no crossing below n_max)':>36}")
            continue
        gain = nc / nm
        print(f"  {T:>6.0f} {nc:>14d} {nm:>12d} {gain:>14.2f} {T:>12.0f} "
              f"{gain/T:>8.3f}")

    print("\n  the last column is what the paper does not predict: the rate law gives")
    print("  a gain of exactly T, and any departure is the background's own dependence")
    print("  on the lens, which section 3.4 identifies but does not compute.")
    print("\n  [C] 'first negative' is a crossing of an oscillating quantity; the index")
    print("  is therefore sensitive to where the oscillation happens to sit, and two")
    print("  runs differing only in the zero count can differ by a few percent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
