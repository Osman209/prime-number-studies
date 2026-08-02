#!/usr/bin/env python3
"""
verify_redheffer.py -- regenerates every number printed in section 1 of
"A Numerical Study of the Division Table" concerning the 0-1 shadow of the
sheet: the empty spectrum of the raw divisibility matrix, the Redheffer
modification whose determinant is the Mertens function, the weighted version
whose determinant is a truncated Dirichlet series for 1/zeta(s), and the count
of non-trivial eigenvalues.

Exits 0 if every check passes, 1 otherwise.  Run:  python3 verify_redheffer.py
Optional:  --fast   (skips the slowest exact determinants)
"""
import sys, math
from fractions import Fraction

FAIL = []
FAST = "--fast" in sys.argv


def note(ok, label, detail=""):
    print(("  ok   " if ok else "  FAIL ") + label + (("   " + detail) if detail else ""))
    if not ok:
        FAIL.append(label)


# ---------------------------------------------------------------- utilities
def divisibility_matrix(n):
    """K(d,m) = 1 iff d | m, rows and columns 1..n."""
    return [[1 if m % d == 0 else 0 for m in range(1, n + 1)] for d in range(1, n + 1)]


def redheffer(n, weights=None):
    """K with the first column overwritten: entry (d,1) = weights[d] (default 1)."""
    A = [[Fraction(x) for x in row] for row in divisibility_matrix(n)]
    for d in range(1, n + 1):
        A[d - 1][0] = Fraction(1) if weights is None else weights[d]
    return A


def det_exact(A):
    """Fraction-free-enough Gaussian elimination over Q.  Exact."""
    A = [row[:] for row in A]
    n = len(A)
    det = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det *= A[c][c]
        inv = Fraction(1) / A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] * inv
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return det


def mobius_sieve(N):
    mu = [1] * (N + 1)
    primes, comp = [], [False] * (N + 1)
    mu[0] = 0
    for i in range(2, N + 1):
        if not comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > N:
                break
            comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu


MU = mobius_sieve(200000)
def M(n):
    return sum(MU[1:n + 1])


# ---------------------------------------------------------------- check 1
print("\n[1] the raw table is triangular: spectrum {1}, det 1, and not diagonalisable")
try:
    import numpy as np
except ImportError:
    np = None
    print("  (numpy absent -- eigenvalue checks skipped)")

for n in (20, 60, 121):
    K = divisibility_matrix(n)
    upper = all(K[d][m] == 0 for d in range(n) for m in range(d))
    diag = all(K[d][d] == 1 for d in range(n))
    d1 = det_exact([[Fraction(x) for x in row] for row in K])
    note(upper and diag and d1 == 1,
         f"n={n:3d}: upper triangular, unit diagonal, det = {d1}")
    if np is not None:
        A = np.array(K, dtype=float)
        ker = n - np.linalg.matrix_rank(A - np.eye(n))
        odd = (n + 1) // 2                      # number of odd integers <= n
        note(ker == odd,
             f"n={n:3d}: dim ker(K-I) = {ker}, number of odd integers <= n = {odd}")

# ---------------------------------------------------------------- check 2
print("\n[2] one column filled in: det R_n = M(n)   (Redheffer 1977)")
printed = [1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2]     # the paper's list, n = 1..12
got = [int(det_exact(redheffer(n))) for n in range(1, 13)]
note(got == printed, "n = 1..12 matches the list printed in the paper", str(got))

ns = [15, 20, 30, 50] + ([] if FAST else [65, 80])
for n in ns:
    d = int(det_exact(redheffer(n)))
    note(d == M(n), f"n={n:3d}: det R_n = {d:4d} = M(n)")

# ---------------------------------------------------------------- check 3
print("\n[3] weighted column w_k = k^{-s}: det = truncated Dirichlet series for 1/zeta(s)")
for s in (2.0, 3.0):
    for n in (12, 25, 40):
        w = {k: Fraction(1) for k in range(1, n + 1)}          # placeholder, filled below
        # use floats here: the weights are irrational for non-integer s
        A = [[float(x) for x in row] for row in divisibility_matrix(n)]
        for d in range(1, n + 1):
            A[d - 1][0] = d ** (-s)
        # float Gaussian elimination is fine for the comparison, the matrix is small
        det, m = 1.0, [row[:] for row in A]
        for c in range(n):
            p = max(range(c, n), key=lambda r: abs(m[r][c]))
            if abs(m[p][c]) < 1e-14:
                det = 0.0
                break
            if p != c:
                m[c], m[p] = m[p], m[c]
                det = -det
            det *= m[c][c]
            for r in range(c + 1, n):
                f = m[r][c] / m[c][c]
                for k in range(c, n):
                    m[r][k] -= f * m[c][k]
        trunc = sum(MU[k] * k ** (-s) for k in range(1, n + 1))
        note(abs(det - trunc) < 1e-9,
             f"s={s}, n={n:2d}: det = {det:.12f}, sum mu(k)k^-s = {trunc:.12f}")

# ---------------------------------------------------------------- check 4
print("\n[4] only floor(log2 n)+1 eigenvalues differ from 1; two are large")
if np is None:
    print("  (skipped: numpy absent)")
else:
    for n in (64, 200):
        A = np.array([[float(x) for x in row] for row in redheffer(n)])
        ev = np.linalg.eigvals(A)
        nontriv = sorted([e for e in ev if abs(e - 1) > 0.05], key=lambda z: -abs(z))
        pred = int(math.log2(n)) + 1
        note(len(nontriv) == pred,
             f"n={n:3d}: {len(nontriv)} eigenvalues off 1, floor(log2 n)+1 = {pred}")
        bj = math.sqrt(n) + math.log(math.sqrt(n)) + 0.5772156649015329 - 0.5
        note(abs(nontriv[0].real - bj) / bj < 0.03,
             f"n={n:3d}: largest = {nontriv[0].real:+.3f}, Barrett-Jarvis = {bj:+.3f}")
    print("  NOTE: the 0.05 tolerance undercounts for large n -- Vaughan proved there are")
    print("        non-trivial eigenvalues arbitrarily close to 1.  This is not a bug.")

# ---------------------------------------------------------------- verdict
print()
if FAIL:
    print("FAILED:", len(FAIL), "check(s):")
    for f in FAIL:
        print("   -", f)
    sys.exit(1)
print("all checks passed" + ("  (--fast: largest exact determinants skipped)" if FAST else ""))
sys.exit(0)
