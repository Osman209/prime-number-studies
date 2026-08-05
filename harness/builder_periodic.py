"""
builder_periodic.py -- the CvS (periodic Fourier) constructor behind the same
interface as builder_sine.build, so both families run under one conventions.py.

This file adds no mathematics of its own. It calls
connes_cvs.operator.build_galerkin_matrix (A. Groskin, PyPI `connes-cvs`), which
implements CvS Proposition 4.1, and re-presents the result in the shape the rest
of this harness expects.

--------------------------------------------------------------------------------
THE NORMALISATION MAP, read off the constructor's own source rather than assumed.

`build_galerkin_matrix(c, N, T, dps)` sets L = mp.log(c) internally and sums prime
powers n <= c with weight Lambda(n)/sqrt(n).  For the basis, take the periodic
side from what it actually reconstructs: connes_cvs.extract_zeros integrates
against exp(2 pi i k x / L) on [0, L], so the periodic family is the FULL Fourier
basis of period L, frequencies 2 pi k / L.  conventions.py fixes the support edge
L, sums prime powers with log q < L, and works in

    phi_j(x) = sqrt(2/L) sin(j pi (x + L/2)/L),   j = 1 .. m,

frequencies j pi / L.  So the periodic sine half (e_k - e_{-k}) ~ sin(2 pi k x/L)
is exactly phi_{2k} -- this family's EVEN-INDEXED modes -- while the periodic
cosine half and this family's odd-indexed modes have no counterpart on the other
side at all.  Hence the parameters correspond as

    L = log c,      p_c = c   (the saturated diagonal),      m = 2N.

The prime q = c is included by the periodic side (n <= c) and excluded by this
one (log q < L, strict).  Its overlap integrates over an interval of length
L - log c = 0, so it contributes nothing and the two prime sums agree.

--------------------------------------------------------------------------------
WHAT THE TWO FAMILIES SHARE, AND WHAT THEY DO NOT.  Both forms commute with
parity, so each splits into an even and an odd sector, and the two sectors behave
differently:

  * The periodic odd combinations (e_k - e_{-k}) ~ sin(2 pi k x / L) vanish at the
    centre AND at both endpoints.  They ARE the Dirichlet sine modes, and the two
    constructors compute the SAME operator there -- measured difference falls
    like the archimedean truncation, 1.7e-4 -> 2.4e-5 -> 3.2e-6 at T = 100, 200,
    400 (c = 13, N = 16, dps = 80).

  * The two "even" sectors share no vector at FINITE size -- one is spanned by
    cosines, the other by the half-integer-frequency sines -- and the measured gap
    is flat in T: 5.0e-2, 4.7e-2, 4.9e-2 over the same sweep.  They are not
    disjoint in the limit: the odd-indexed sine modes carry the same parity and are
    complete in the same subspace, with an L2 residual falling like 1/m.  So the
    difference should converge in m rather than in T -- an inference from
    completeness, since the m-dependence of the gap was not swept.

  * The shared sector carries the zeros. At c = 13, N = 100, T = 400, dps = 80
    the ODD-sector near-null vector reconstructs the first three ordinates to
    8.70e-53, 6.67e-50, 4.28e-48 against the even sector's 1.45e-55, 2.69e-52,
    2.49e-50.  The endpoint condition costs about three decimal places; it is not
    what makes the reconstruction work.

  * The near-null counts nevertheless agree in BOTH sectors (7/7 and 7/7 at
    N = 20; 5/5 and 5/5 at N = 12).

So "one set of conventions" is possible, but only if the report is made sector by
sector.  A single matrix-to-matrix comparison of the two families mixes a sector
where they agree with a sector where they must not.

--------------------------------------------------------------------------------
PRECISION.  Two traps, both measured:

  * dps must be raised with T.  At dps = 40 and T = 800 the comparison degrades
    to 5.5e-1, worse than T = 200 gives; at dps = 80 the T-trend is monotone.
  * `extract_zeros` must receive L as an mpmath number.  Passing a Python float
    caps the reconstruction at ~1e-17 instead of ~1e-55 -- 38 orders of
    magnitude, silently.  See upstream issue #3.  This module always passes
    mp.log(c) and never a float.

Requires: pip install connes-cvs   (python-flint strongly recommended; without it
the build is roughly a hundred times slower).
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp

try:
    from connes_cvs.operator import build_galerkin_matrix
except ImportError as exc:                                    # pragma: no cover
    raise ImportError(
        "builder_periodic requires the `connes-cvs` package: pip install connes-cvs"
    ) from exc


# ------------------------------------------------------------------ parameter map
def N_for(m: int) -> int:
    """Basis half-size N of the periodic family matching m sine modes."""
    if m % 2:
        raise ValueError(f"m must be even for the m = 2N correspondence, got {m}")
    return m // 2


def c_for(L: float) -> int:
    """Prime cutoff c matching a support edge L on the saturated diagonal.

    floor(exp(L)) is WRONG here: in float64, exp(log(19)) rounds just below 19 and
    the floor returns 18 -- likewise 36 for 37 and 66 for 67, silently building at
    the wrong cutoff. Round instead, and refuse the call unless L really is log of
    the integer returned.
    """
    c = int(round(math.exp(L)))
    if c < 2 or abs(L - math.log(c)) > 1e-9:
        raise ValueError(
            f"L = {L!r} is not log of an integer cutoff to 1e-9; nearest is {c} "
            f"with log {math.log(max(c,2))!r}. The saturated diagonal is only "
            f"defined at L = log c for integer c."
        )
    return c


# ------------------------------------------------------------------ parity split
def parity_projectors(N: int):
    """(V_even, V_odd) mapping the (2N+1) trigonometric basis onto its two sectors.

    V_even spans e_0 and (e_k + e_{-k})/sqrt2; V_odd spans (e_k - e_{-k})/sqrt2,
    the sine modes.  Returned as real arrays, since Q is real symmetric."""
    D = 2 * N + 1
    r = 1.0 / math.sqrt(2.0)
    Ve = np.zeros((D, N + 1))
    Vo = np.zeros((D, N))
    Ve[N, 0] = 1.0
    for k in range(1, N + 1):
        Ve[N + k, k] = Ve[N - k, k] = r
        Vo[N + k, k - 1] = r
        Vo[N - k, k - 1] = -r
    return Ve, Vo


def sine_parity_index(m: int):
    """(even_idx, odd_idx) into the sine basis.

    phi_j(-x) = (-1)^{j+1} phi_j(x), so j odd gives functions even about the
    centre and j even gives odd ones.  Returned as 0-based column indices."""
    return np.arange(0, m, 2), np.arange(1, m, 2)


# ------------------------------------------------------------------ assembly
def build(L: float, p_c: float, m: int, T: int = 400, dps: int = 80) -> dict:
    """Periodic-family counterpart of builder_sine.build, same return keys.

    Only the saturated diagonal is defined for this family: the periodic
    constructor derives L from c, so L and p_c are not independent here.  Off the
    diagonal there is no periodic constructor to call, and the caller is told so
    rather than being handed a matrix built at the wrong support.

    'arch', 'pole' and 'prime' are returned as None: the upstream constructor
    assembles Q from psi(n) and psi'(n) and does not expose the three blocks
    separately.  Everything else has the same meaning as in builder_sine.
    """
    c = c_for(L)
    if int(p_c) != c:
        raise ValueError(
            f"the periodic family is defined only on the saturated diagonal: "
            f"p_c must be floor(exp(L)) = {c}, got {int(p_c)}"
        )
    N = N_for(m)
    Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    D = Q.rows
    W = np.array([[float(Q[i, j]) for j in range(D)] for i in range(D)])
    W = 0.5 * (W + W.T)
    return {
        "arch": None, "pole": None, "prime": None, "W": W,
        "arch_tail_bound": None,
        "L": float(mp.log(c)), "p_c": c, "m": m, "N": N, "dim": D,
        "saturated": True,
        "prime_powers": None,
        "T": T, "dps": dps,
        "family": "periodic",
    }


def sectors(W: np.ndarray, N: int) -> dict:
    """Eigenvalues of the periodic matrix restricted to each parity sector."""
    Ve, Vo = parity_projectors(N)
    return {
        "even": np.sort(np.linalg.eigvalsh(Ve.T @ W @ Ve)),
        "odd": np.sort(np.linalg.eigvalsh(Vo.T @ W @ Vo)),
    }


def sine_sectors(W: np.ndarray, m: int) -> dict:
    """The same split for a sine-family matrix, for a like-for-like comparison."""
    ie, io = sine_parity_index(m)
    return {
        "even": np.sort(np.linalg.eigvalsh(W[np.ix_(ie, ie)])),
        "odd": np.sort(np.linalg.eigvalsh(W[np.ix_(io, io)])),
    }


def inertia(W: np.ndarray, tau: float) -> dict:
    """Full inertia triple, same contract as builder_sine.inertia."""
    w = np.linalg.eigvalsh(0.5 * (W + W.T))
    return {
        "lambda_min": float(w[0]),
        "lambda_max": float(w[-1]),
        "n_minus": int((w < -tau).sum()),
        "n_zero": int(((w >= -tau) & (w <= tau)).sum()),
        "n_plus": int((w > tau).sum()),
        "tau": tau,
    }
