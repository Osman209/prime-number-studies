"""
edge_precision.py -- the edge limit in extended precision.

The overlap R(t) is O((L-t)^3) as t -> L, obtained as a difference of two O(L-t)
quantities. In float64 that cancellation destroys the answer before delta gets
small: the residual of the leading-order expansion stops falling at delta ~ 1e-4
and RISES again by 1e-5, and the parity-block ratio sigma2/sigma1 turns around in
the same place. Neither turn is mathematics; both are the arithmetic.

This module removes the problem twice over. First it replaces the quadrature-style
closed form of builder_sine.overlap by an exact reduction with no trigonometric
cancellation left in its structure, and then it evaluates that reduction in mpmath
at a precision set by the target delta.

    THE REDUCTION.  With a_j = j pi / L, delta = L - t and S_j = sin(j pi delta / L),

        R_jk = ((-1)^k / pi) [ (S_j - S_k)/(j - k)  -  (S_j + S_k)/(j + k) ]

    and on the diagonal the first bracket is its limit (pi delta / L) cos(j pi delta / L).
    Symmetrising and expanding S_j = x_j - x_j^3/6 + ... with x_j = j pi delta / L,
    the O(delta) parts cancel identically and

        R_jk = (pi^2 delta^3 / 6 L^3) j k ( (-1)^{j+1} + (-1)^{k+1} ) + O(delta^5),

    which is the expansion of A. Groskin (see the attribution in
    papers/prime_edge_two_paths.md), here re-derived from the closed form rather
    than assumed.

Usage:
    python edge_precision.py                 # all four checks
    python edge_precision.py --m 20 40 80    # several basis sizes
    python edge_precision.py --dps 16        # see the float64 failure reproduced

That last one is the point of the module in one command. Pinning the working
precision at 16 digits instead of setting it from delta makes check 3 print
ratios 100.07, 100.015, 6.06, 0.0097 -- the same turnaround float64 shows, at the
same place. Raise the precision and it goes away. The turnaround is arithmetic.

Requires mpmath and numpy. The eigenproblem runs in software arithmetic but the
matrices are small: about 6 s for m = 20 and 40 together, 10 s for m = 80.
"""
from __future__ import annotations
import argparse
import math
import mpmath as mp
import numpy as np

from conventions import L_DEFAULT


# ----------------------------------------------------------------- the overlap
def overlap(delta, L, m):
    """Symmetrised overlap at t = L - delta, in mpmath at the current precision."""
    d = mp.mpf(delta)
    Ls = mp.mpf(L)
    x = [j * mp.pi * d / Ls for j in range(m + 1)]          # x[j] = j pi delta / L
    S = [mp.sin(v) for v in x]
    C = [mp.cos(v) for v in x]
    R = mp.zeros(m, m)
    for j in range(1, m + 1):
        for k in range(1, m + 1):
            if j == k:
                first = (mp.pi * d / Ls) * C[j]
            else:
                first = (S[j] - S[k]) / (j - k)
            val = ((-1) ** k / mp.pi) * (first - (S[j] + S[k]) / (j + k))
            R[j - 1, k - 1] = val
    for j in range(m):
        for k in range(j + 1, m):
            avg = (R[j, k] + R[k, j]) / 2
            R[j, k] = avg
            R[k, j] = avg
    return R


def leading(delta, L, m):
    """The O(delta^3) leading matrix."""
    d = mp.mpf(delta)
    Ls = mp.mpf(L)
    pre = mp.pi ** 2 * d ** 3 / (6 * Ls ** 3)
    A = mp.zeros(m, m)
    for j in range(1, m + 1):
        for k in range(1, m + 1):
            A[j - 1, k - 1] = pre * j * k * ((-1) ** (j + 1) + (-1) ** (k + 1))
    return A


def abs_eigs(M):
    """|eigenvalues| of a real symmetric mpmath matrix, descending."""
    E, _ = mp.eigsy(M)
    return sorted((abs(E[i]) for i in range(M.rows)), reverse=True)


def fro(M):
    return mp.sqrt(sum(M[i, j] ** 2 for i in range(M.rows) for j in range(M.cols)))


def dps_for(delta, m, headroom=25, order=2):
    """Precision needed at a given delta.

    The matrix is O(delta^3) assembled from O(delta) pieces, so about
    2*log10(1/delta) digits vanish to cancellation before anything is computed.
    Quantities that are themselves small relative to the matrix need more: the
    parity ratio sigma2/sigma1 is O(delta^4), so it is requested with order=4."""
    return int(order * math.log10(1.0 / delta) + math.log10(max(m, 2)) + headroom)


# ------------------------------------------------------------- float64 baseline
def overlap_f64(delta, L, m):
    d = float(delta)
    j = np.arange(1, m + 1)
    S = np.sin(j * np.pi * d / L)
    C = np.cos(j * np.pi * d / L)
    J = j[:, None] * np.ones(m)[None, :]
    K = np.ones(m)[:, None] * j[None, :]
    Sj = S[:, None] * np.ones(m)[None, :]
    Sk = np.ones(m)[:, None] * S[None, :]
    first = np.where(J == K,
                     (np.pi * d / L) * C[:, None] * np.ones(m)[None, :],
                     (Sj - Sk) / np.where(J == K, 1.0, J - K))
    R = ((-1.0) ** K / np.pi) * (first - (Sj + Sk) / (J + K))
    return 0.5 * (R + R.T)


# --------------------------------------------------------------------- checks
def check_agreement(L, m):
    print("\n=== 1. the reduction agrees with builder_sine.overlap where float64 still works ===")
    try:
        from builder_sine import overlap as builder_overlap
    except Exception as e:                                   # pragma: no cover
        print(f"   (builder_sine not importable: {e})")
        return
    print(f"   {'delta':>8} {'max abs diff':>14} {'matrix scale':>14} {'relative':>11}")
    for d in (1e-1, 1e-2, 1e-3):
        A = builder_overlap(L - d, L, m)
        B = overlap_f64(d, L, m)
        sc = np.abs(A).max()
        print(f"   {d:8.0e} {np.abs(A - B).max():14.3e} {sc:14.3e} "
              f"{np.abs(A - B).max() / sc:11.3e}")
    print("   Same object. What follows changes only the arithmetic it is evaluated in.")


def check_law(L, ms, epsilons):
    print("\n=== 2. the (m-1)/(m+2) law, as eps -> 0 ===")
    print(f"   {'m':>4} " + " ".join(f"{'eps=' + f'{e:.0e}':>14}" for e in epsilons) + f" {'(m-1)/(m+2)':>14}")
    for m in ms:
        row = []
        for eps in epsilons:
            delta = L * eps
            mp.mp.dps = dps_for(delta, m)
            ev = abs_eigs(overlap(delta, L, m))
            row.append(ev[1] / ev[0])
        pred = mp.mpf(m - 1) / (m + 2)
        print(f"   {m:>4} " + " ".join(f"{mp.nstr(r, 10):>14}" for r in row) +
              f" {mp.nstr(pred, 10):>14}")
        print(f"        diff " + " ".join(f"{mp.nstr(abs(r - pred), 3):>14}" for r in row))
    print("   In float64 the discrepancy bottoms out near 1e-6 and then grows; here it")
    print("   falls monotonically, which is what a delta -> 0 limit is supposed to do.")


def check_expansion(L, m, deltas):
    print("\n=== 3. the O(delta^5) remainder really is O(delta^5) ===")
    print(f"   {'delta':>9} {'||R-lead||/||R||':>20} {'ratio':>10}   (delta^2 scaling predicts 100)")
    prev = None
    for d in deltas:
        mp.mp.dps = dps_for(d, m)
        R = overlap(d, L, m)
        A = leading(d, L, m)
        res = fro(R - A) / fro(R)
        line = f"   {d:9.0e} {mp.nstr(res, 8):>20}"
        if prev is not None:
            line += f" {mp.nstr(prev / res, 6):>10}"
        print(line)
        prev = res
    print("   float64 gives 100.3 then 6.4 and then turns around; the ratio here stays")
    print("   at 100 for as long as the precision is raised to match delta.")


def check_parity(L, m, epsilons):
    print("\n=== 4. each parity block approaches rank one, without turning around ===")
    odd = [j for j in range(1, m + 1) if j % 2 == 1]
    even = [j for j in range(1, m + 1) if j % 2 == 0]
    print(f"   {'eps':>9} {'odd s2/s1':>18} {'even s2/s1':>18} {'odd sign':>9} {'even sign':>10}")
    for eps in epsilons:
        delta = L * eps
        mp.mp.dps = dps_for(delta, m, order=6)
        R = overlap(delta, L, m)
        out = []
        signs = []
        for idx in (odd, even):
            B = mp.matrix(len(idx), len(idx))
            for a, j in enumerate(idx):
                for b, k in enumerate(idx):
                    B[a, b] = R[j - 1, k - 1]
            E, _ = mp.eigsy(B)
            vals = sorted(((abs(E[i]), E[i]) for i in range(B.rows)), reverse=True)
            out.append(vals[1][0] / vals[0][0])
            signs.append("+" if vals[0][1] > 0 else "-")
        print(f"   {eps:9.0e} {mp.nstr(out[0], 10):>18} {mp.nstr(out[1], 10):>18} "
              f"{signs[0]:>9} {signs[1]:>10}")
    print("   float64 gives 3.6e-7, 9.8e-9, then 2.6e-5 and 1.4e-2 -- the rise was the")
    print("   arithmetic. The dominant eigenvalue is positive on the odd block and")
    print("   negative on the even one, so the full matrix is their difference.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=L_DEFAULT)
    ap.add_argument("--m", type=int, nargs="*", default=[20, 40])
    ap.add_argument("--dps", type=int, default=None,
                    help="override the automatic precision (not recommended)")
    args = ap.parse_args()
    if args.dps:
        global dps_for
        dps_for = lambda *a, **k: args.dps          # noqa: E731

    L = args.L
    print(f"edge_precision.py   L = {L}   m = {args.m}")
    print("precision is set per delta from dps_for(); float64 is shown only for contrast.")
    check_agreement(L, min(args.m))
    check_law(L, args.m, (1e-3, 1e-4, 1e-5, 1e-6))
    check_expansion(L, min(args.m), (1e-2, 1e-3, 1e-4, 1e-5, 1e-6))
    check_parity(L, min(args.m), (1e-3, 1e-4, 1e-5, 1e-6, 1e-8))
    print("\ndone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
