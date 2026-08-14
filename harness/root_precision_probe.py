#!/usr/bin/env python3
"""
root_precision_probe.py — is a branch's deviation from gamma REAL, or is it the
precision of the scan?

THE QUESTION.  At c = 100, N = 150, T = 400 the branch with lambda = -3.32e-248
places its first root 5.76e-244 from gamma_1.  Three explanations were on the table
and sizing them from the outside did not separate them:

  (a) the eigenvector from the dps = 500 diagonalisation.  The spectral estimate is
      ||Q|| * 10^-500 / gap = 7.8641 * 1e-500 / 5.77e-240 ~ 1.4e-260.  BUT a root
      position is not a vector error: |dr| ~ |dF| / |F'(r)|, so a condition number
      of ~1e16 turns 1.4e-260 into 1e-244.  **Not excluded without measuring F'.**
  (b) the scan at ref_dps.  Cancellation destroys RELATIVE precision, not absolute:
      summing terms of total size ~10 at 300 digits gives ~1e-298, not 1e-244.  To
      blame the scan you need |F'(r)| ~ 1e-55.  **Also a measurement, not a guess.**
  (c) the deviation is REAL — first-order perturbation, r - gamma ~ lambda * D/Z',
      which at |lambda| = 3.32e-248 needs an amplification of only 1.7e4.

The bisection tolerance at ref_dps = 300 is 1e-280, thirty-six orders BELOW the
observed 1e-244, so the refinement floor is not the ceiling either.

THE DECISIVE TEST, and it needs no rebuild.  Recompute the SAME root from the SAME
cached eigenvector at several ref_dps and compare the roots TO EACH OTHER rather
than to zeta:

  * if r(260) = r(300) = r(350) = r(400) agree past digit 244, the deviation is real
    and (c) is the answer;
  * if the root itself moves near digit 244 as ref_dps changes, the scan is the limit.

Alongside it the script measures |F(r)|, |F'(r)| and the condition number

    kappa(r) = sum_k |c_k * term_k(r)| / |F'(r)|,

which turns the scan bound into an actual number: |dr|_scan <~ kappa(r) * 10^-ref_dps.

    python3 root_precision_probe.py Q_c100_N150_T400_dps500_even_eig.json --c 100 --N 150
    python3 root_precision_probe.py <eig.json> --c 100 --N 150 --dps-list 260,300,350,400

NOTE ON SCOPE.  This concerns the ACCURACY of a root's position, not the COUNT of
roots.  Matching an ordinate needs only 1e-3 and the spurious roots are 1e-2 to 1
apart, so 300 digits is far more than enough to see their signs.  The count is
threatened by the sweep step h, not by ref_dps, and that is tested separately.
"""
from __future__ import annotations
import argparse, json, sys, time
import mpmath as mp


def make_terms(coef, N, L):
    """Return (F, terms) where terms(tau) gives every summand, for kappa."""
    def summands(tau, dps):
        old = mp.mp.dps
        mp.mp.dps = dps
        try:
            t, Lm = mp.mpf(tau), mp.mpf(L)
            e = mp.exp(-1j * t * Lm)
            out = []
            for k in range(-N, N + 1):
                ck = coef[k + N]
                if ck == 0:
                    continue
                den = 2 * mp.pi * k / Lm - t
                out.append(ck * (mp.mpc(Lm, 0)
                                 if abs(den) < mp.mpf(10) ** (-(dps - 10))
                                 else (e - 1) / (1j * den)))
            return out, mp.exp(1j * t * Lm / 2) / mp.sqrt(Lm)
        finally:
            mp.mp.dps = old

    def F(tau, dps):
        s, pre = summands(tau, dps)
        return mp.re(pre * sum(s))

    def kappa(tau, dps, Fp):
        s, pre = summands(tau, dps)
        scale = sum(abs(pre * x) for x in s)
        return scale / abs(Fp) if Fp != 0 else mp.inf
    return F, kappa


def refine(F, lo, hi, dps, tol):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    fl = F(lo, dps)
    n = int(mp.ceil(mp.log((hi - lo) / tol, 2))) + 10
    for _ in range(n):
        mid = (lo + hi) / 2
        fm = F(mid, dps)
        if fm == 0:
            return mid
        if mp.sign(fm) != mp.sign(fl):
            hi = mid
        else:
            lo, fl = mid, fm
        if hi - lo < tol:
            break
    return (lo + hi) / 2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("eig")
    ap.add_argument("--c", type=int, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--branch", default=None,
                    help="key in the eig cache; default is the one nearest zero")
    ap.add_argument("--zero", type=int, default=1, help="which gamma_k to examine")
    ap.add_argument("--dps-list", default="260,300,350,400")
    a = ap.parse_args()

    dpss = [int(x) for x in a.dps_list.split(",")]
    top = max(dpss) + 60
    ec = json.load(open(a.eig, encoding="utf-8"))
    mp.mp.dps = top
    E = [mp.mpf(x) for x in ec["eigenvalues"]]
    key = a.branch or min(ec["vectors"], key=lambda k: abs(E[int(k)]))
    coef = [mp.mpf(x) for x in ec["vectors"][key]]
    lam = E[int(key)]
    L = mp.log(a.c)
    print(f"branch {key}: lambda = {mp.nstr(lam, 10)}  (|lambda| = {mp.nstr(abs(lam),6)})")

    g = mp.im(mp.zetazero(a.zero))
    print(f"gamma_{a.zero} = {mp.nstr(g, 25)}\n")

    F, kappa = make_terms(coef, a.N, L)
    # a bracket around gamma that does not use gamma to more than 3 digits
    lo, hi = g - mp.mpf("0.01"), g + mp.mpf("0.01")

    roots = {}
    hdr = ("ref_dps", "|r - gamma|", "|F(r)|", "|Fprime(r)|", "kappa", "kappa*10^-dps")
    print("%8s %16s %14s %14s %12s %15s" % hdr)
    for d in dpss:
        t0 = time.time()
        tol = mp.mpf(10) ** (-(d - 20))
        r = refine(F, lo, hi, d, tol)
        # DIFFERENTIATE AT THE SAME PRECISION F IS EVALUATED AT. mp.diff picks its
        # step from the CURRENT dps; at `top` that step is far below what F(., d)
        # can resolve, and the quotient collapses to 0. Set dps to d and give the
        # step explicitly.
        mp.mp.dps = d
        h = mp.mpf(10) ** (-(d // 3))
        Fp = (F(r + h, d) - F(r - h, d)) / (2 * h)
        mp.mp.dps = top
        kap = kappa(r, d, Fp)
        roots[d] = r
        bound = kap * mp.mpf(10) ** (-d)
        print(f"{d:>8} {mp.nstr(abs(r-g),6):>16} {mp.nstr(abs(F(r,d)),4):>14} "
              f"{mp.nstr(abs(Fp),4):>14} {mp.nstr(kap,4):>12} {mp.nstr(bound,4):>15}"
              f"   ({time.time()-t0:.0f}s)")

    print("\nTHE ROOTS AGAINST EACH OTHER — this is the decisive comparison:")
    ds = sorted(roots)
    for i in range(len(ds) - 1):
        d = abs(roots[ds[i + 1]] - roots[ds[i]])
        print(f"  |r({ds[i]}) - r({ds[i+1]})| = {mp.nstr(d, 6)}")
    # THE STABILITY MEASURE IS THE LAST PAIR, not the widest spread. A low ref_dps
    # in the list is expected to be wrong; including it would report the failure of
    # the worst setting as though it were the uncertainty of the best.
    spread = abs(roots[ds[-1]] - roots[ds[-2]])
    dev = abs(roots[ds[-1]] - g)
    print(f"\n  movement between the top two settings ({ds[-2]} -> {ds[-1]}) : {mp.nstr(spread,6)}")
    print(f"  deviation of the root from gamma_{a.zero}                  : {mp.nstr(dev,6)}")
    if spread == 0 or dev / spread > mp.mpf(100):
        print("\n  ==> the root is STABLE under ref_dps while sitting far from gamma.")
        print("      The deviation is REAL, not a scan artefact.")
    elif dev / spread < mp.mpf(10):
        print("\n  ==> the root MOVES with ref_dps by an amount comparable to its")
        print("      distance from gamma. The scan is the limit; raise ref_dps.")
    else:
        print("\n  ==> inconclusive: the movement is within two orders of the deviation.")
        print("      Extend --dps-list upward before reading anything into it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
