"""
builder_sine.py -- fixed-support (Dirichlet sine) constructor for the truncated
Weil quadratic form, with every block exposed separately.

Conventions: see conventions.py. Nothing here chooses a convention on its own.

    build(L, p_c, m) -> dict with keys
        'arch', 'pole', 'prime', 'W'          (m x m float arrays)
        'L', 'p_c', 'm', 'saturated'          (parameters and diagonal flag)
        'prime_powers'                        (list of (q, Lambda(q)) actually used)

The two parameters L and p_c are independent. Passing p_c = floor(exp(L)) puts you
on the saturated diagonal, where the prime sum is complete up to the support edge
and the zero-side identity of conventions.py section 5 is expected to hold.
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from sympy import primerange

from conventions import R_MAX, NR


# ------------------------------------------------------------------ basis
def modes(L: float, m: int) -> np.ndarray:
    """a_j = j pi / L for j = 1..m."""
    return np.arange(1, m + 1) * math.pi / L


def phihat(r: np.ndarray, L: float, m: int) -> np.ndarray:
    """(len(r), m) array of phihat_j(r). Removable poles at r = a_j are not
    evaluated -- callers must supply an r-grid avoiding them (see arch_block)."""
    a = modes(L, m)
    sgn = (-1.0) ** np.arange(1, m + 1)
    E = np.exp(-1j * r * L)
    den = a[None, :] ** 2 - r[:, None] ** 2
    return (math.sqrt(2.0 / L) * a[None, :] * (1.0 - sgn[None, :] * E[:, None]) / den
            * np.exp(1j * r * L / 2.0)[:, None])


# ------------------------------------------------------------------ blocks
def arch_block(L: float, m: int, r_max: float = R_MAX, nr: int = NR):
    """Archimedean block, plus a bound on the truncation tail beyond r_max.

    Returns (A, tail_bound). Psi(r) = Re psi(1/4 + i r/2) - log pi grows like
    log(r/2); |phihat_j|^2 = O(1/r^4), so the omitted tail is bounded by
    integrating C log(r) / r^4 from r_max to infinity."""
    r = np.linspace(-r_max, r_max, nr)
    a = modes(L, m)
    keep = np.min(np.abs(r[:, None] - a[None, :]), axis=1) > 1e-9
    r = r[keep]
    Psi = np.array([float(mp.re(mp.digamma(mp.mpc(0.25, x / 2.0)))) for x in r]) - math.log(math.pi)
    Phi = phihat(r, L, m)
    w = Psi * (2.0 * r_max / (nr - 1))
    A = np.real((np.conj(Phi).T * w) @ Phi) / (2.0 * math.pi)

    # tail bound: |phihat_j(r)| <= sqrt(2/L) * a_j * 2 / (r^2 - a_j^2) for r > a_max
    a_max = a[-1]
    if r_max > 2.0 * a_max:
        C = (2.0 / L) * (2.0 * a_max) ** 2 * 4.0
        tail = C * (math.log(r_max / 2.0) + 1.0 / 3.0) / (3.0 * r_max ** 3) / (2.0 * math.pi)
    else:
        tail = float("inf")   # r_max too small for the bound to be meaningful
    return A, tail


def pole_block(L: float, m: int) -> np.ndarray:
    """Rank-two pole block: m+ (m-)^T + m- (m+)^T."""
    a = modes(L, m)
    j = np.arange(1, m + 1)
    lo = -L / 2.0
    out = {}
    for s in (0.5, -0.5):
        num = math.exp(s * lo) * (math.exp(s * L) * (-a * (-1.0) ** j) + a)
        out[s] = math.sqrt(2.0 / L) * num / (s ** 2 + a ** 2)
    return np.outer(out[0.5], out[-0.5]) + np.outer(out[-0.5], out[0.5])


def _I(w, ph, S):
    w = np.asarray(w, float)
    out = np.empty_like(w)
    small = np.abs(w) < 1e-13
    out[small] = S * np.cos(ph[small] if np.ndim(ph) else ph)
    ws = w[~small]
    phs = ph[~small] if np.ndim(ph) else ph
    out[~small] = (np.sin(ws * S + phs) - np.sin(phs)) / ws
    return out


def overlap(t: float, L: float, m: int) -> np.ndarray:
    """Symmetrised Int_0^{L-t} phi_j(y) phi_k(y+t) dy, in closed form.

    NOTE: cancellation-limited in float64 as t -> L. For the delta -> 0 regime
    use edge_precision.py instead of this routine."""
    if t >= L:
        return np.zeros((m, m))
    a = modes(L, m)
    A = a[:, None] * np.ones(m)[None, :]
    B = np.ones(m)[:, None] * a[None, :]
    S = L - t
    R = (2.0 / L) * 0.5 * (_I(A - B, -B * t, S) - _I(A + B, B * t, S))
    return 0.5 * (R + R.T)


def prime_powers(p_c: float, L: float):
    """(q, Lambda(q)) for prime powers q <= p_c with log q < L."""
    out = []
    for p in primerange(2, int(p_c) + 2):
        k, q = 1, p
        while q <= p_c:
            if math.log(q) < L:
                out.append((q, math.log(p)))
            k += 1
            q = p ** k
    return sorted(out)


def prime_block(L: float, p_c: float, m: int):
    pps = prime_powers(p_c, L)
    P = np.zeros((m, m))
    for q, Lam in pps:
        P += 2.0 * Lam / math.sqrt(q) * overlap(math.log(q), L, m)
    return P, pps


# ------------------------------------------------------------------ assembly
def build(L: float, p_c: float, m: int, r_max: float = R_MAX, nr: int = NR) -> dict:
    A, tail = arch_block(L, m, r_max, nr)
    Po = pole_block(L, m)
    Pr, pps = prime_block(L, p_c, m)
    W = A + Po - Pr
    W = 0.5 * (W + W.T)
    return {
        "arch": A, "pole": Po, "prime": Pr, "W": W,
        "arch_tail_bound": tail,
        "L": L, "p_c": p_c, "m": m,
        "saturated": int(p_c) == int(math.floor(math.exp(L))),
        "prime_powers": pps,
        "r_max": r_max, "nr": nr,
    }


def inertia(W: np.ndarray, tau: float) -> dict:
    """Full inertia triple at an explicitly stated absolute tolerance."""
    w = np.linalg.eigvalsh(0.5 * (W + W.T))
    return {
        "lambda_min": float(w[0]),
        "lambda_max": float(w[-1]),
        "n_minus": int((w < -tau).sum()),
        "n_zero": int(((w >= -tau) & (w <= tau)).sum()),
        "n_plus": int((w > tau).sum()),
        "tau": tau,
    }
