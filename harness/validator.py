"""
validator.py -- non-circular check of the assembled form against the zero side.

The zeta ordinates come from mpmath.zetazero and enter no part of the
construction, so agreement certifies the assembly rather than restating it.

    W(f)  ?=  sum over ALL nontrivial zeros of |fhat(gamma)|^2
           =  2 * sum_{gamma > 0} |fhat(gamma)|^2          (f real)

Only valid on the saturated diagonal p_c = floor(exp(L)); see conventions.py §5.
Off the diagonal prime terms are missing by construction and the identity is not
expected to hold, so the validator refuses to run there unless forced.

Usage:
    python validator.py                 # default diagonal case, prints a table
    from validator import validate      # returns a dict for the harness
"""
from __future__ import annotations
import math, os, sys
import numpy as np
import mpmath as mp

from conventions import GAMMA_MAX_DEFAULT
from builder_sine import build, phihat

ZEROS_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zeros.npy")


def zeta_ordinates(n: int = 700) -> np.ndarray:
    """First n positive ordinates, cached to disk. Computed at 12 digits, which
    is far more than the double-precision comparison can use."""
    if os.path.exists(ZEROS_CACHE):
        g = np.load(ZEROS_CACHE)
        if len(g) >= n:
            return g[:n]
    mp.mp.dps = 12
    g = np.array([float(mp.im(mp.zetazero(k))) for k in range(1, n + 1)])
    np.save(ZEROS_CACHE, g)
    return g


def tail_bound(c: np.ndarray, L: float, m: int, gamma_max: float) -> float:
    """Bound on the omitted zero-sum tail beyond gamma_max.

    |phihat_j(r)| <= sqrt(2/L) * a_j * 2 / (r^2 - a_j^2), so for r > 2 a_max
    |fhat(r)|^2 <= K / r^4 with K = (2/L) * (2 a_max)^2 * 4 * ||c||_1^2.
    Zero density is dN = log(gamma / 2pi) / (2pi) dgamma. Integrating
    K log(g/2pi) / (2pi g^4) from gamma_max to infinity gives the bound below."""
    a_max = m * math.pi / L
    if gamma_max <= 2.0 * a_max:
        return float("inf")
    K = (2.0 / L) * (2.0 * a_max) ** 2 * 4.0 * float(np.sum(np.abs(c))) ** 2
    G = gamma_max
    return 2.0 * K * (math.log(G / (2 * math.pi)) / 3.0 + 1.0 / 9.0) / (2 * math.pi * G ** 3)


def validate(L: float, m: int, n_zeros: int = 700, ranks=None, force: bool = False) -> dict:
    p_c = math.floor(math.exp(L))
    B = build(L, p_c, m)
    if not B["saturated"] and not force:
        raise ValueError("validator is only meaningful on the saturated diagonal")
    W = B["W"]
    gam = zeta_ordinates(n_zeros)
    gmax = float(gam[-1])

    wv, V = np.linalg.eigh(W)
    if ranks is None:
        ranks = [int(0.6 * m), int(0.75 * m), int(0.9 * m), m - 1]

    rows = []
    for k in ranks:
        c = V[:, k]
        fh = (phihat(gam, L, m) @ c)
        zs = 2.0 * float(np.sum(np.abs(fh) ** 2))
        tb = tail_bound(c, L, m, gmax)
        lam = float(wv[k])
        rows.append({
            "rank": k, "lambda": lam, "zero_side": zs,
            "abs_residual": lam - zs, "tail_bound": tb,
            "rel_residual": (lam - zs) / zs if zs != 0 else float("nan"),
            "within_tail": abs(lam - zs) <= max(tb, 1e-12),
        })
    return {"L": L, "p_c": p_c, "m": m, "n_zeros": n_zeros,
            "gamma_max": gmax, "arch_tail_bound": B["arch_tail_bound"], "rows": rows}


if __name__ == "__main__":
    L = float(sys.argv[1]) if len(sys.argv) > 1 else 4.5
    m = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    out = validate(L, m)
    print(f"L = {out['L']}  p_c = {out['p_c']} (saturated)  m = {out['m']}")
    print(f"zeros used: {out['n_zeros']}, gamma_max = {out['gamma_max']:.1f}")
    print(f"archimedean truncation tail bound: {out['arch_tail_bound']:.2e}\n")
    print(f"{'rank':>6} {'lambda':>16} {'zero side':>16} {'abs resid':>12} {'rel resid':>12} {'tail bd':>11}")
    for r in out["rows"]:
        print(f"{r['rank']:>6} {r['lambda']:>16.9f} {r['zero_side']:>16.9f} "
              f"{r['abs_residual']:>12.2e} {r['rel_residual']:>12.2e} {r['tail_bound']:>11.2e}")
