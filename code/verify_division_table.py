#!/usr/bin/env python3
"""
verify_division_table.py -- regenerates the numbers printed in
"A Numerical Study of the Division Table" and fails loudly if any of them moves.

Sections covered: 1 (table, harmonic sums, Redheffer determinant), 2-3 (arches and
the five laws), 4 (Fermat step counts, arch spectrum), 5 (the 2041 table, striping,
Dickman, the three residue families), 6 (rho^ord and its summatory function, error
statistics), 7 (wobbles, gamma, the boundary table, Voronoi), 8 (soft edge), 9 (pi/4).
Sections 10 and 11 have their own scripts: verify_knowledge_clock.py and
verify_row_inheritance.py.

    python3 verify_division_table.py            full run, sieve limit 10^7  (~10 min)
    python3 verify_division_table.py --fast     sieve limit 10^6            (~1 min)

Exits 0 only if every check passes.
"""
import math, sys
from math import isqrt, gcd
from fractions import Fraction
import numpy as np

FAST = "--fast" in sys.argv
XLIM = 10**6 if FAST else 10**7
G = np.euler_gamma
FAILS = []


def chk(label, got, want, tol=0.0, rel=False):
    if isinstance(got, (list, tuple)):
        good = list(got) == list(want)
    elif isinstance(got, float) or isinstance(want, float):
        d = abs(got - want)
        good = d <= (tol * abs(want) if rel else tol)
    else:
        good = got == want
    print(("  ok   " if good else "  FAIL ") + f"{label}: {got}" +
          ("" if good else f"   (paper: {want})"))
    if not good:
        FAILS.append(label)


def tau_sieve(n):
    t = np.zeros(n + 1, dtype=np.int32)
    for d in range(1, n + 1):
        t[d::d] += 1
    return t


def D_hyp(n):
    r = isqrt(n)
    return 2 * sum(n // d for d in range(1, r + 1)) - r * r


# ----------------------------------------------------------------- section 1
print("\n=== 1. the table itself ===")
t100 = tau_sieve(100)
chk("D(100)", int(t100[1:].sum()), 482)
chk("whole entries D(100)-100", int(t100[1:].sum()) - 100, 382)
Hn = lambda n: sum(1 / d for d in range(1, n + 1))
chk("R(5) = 5(H_5-1)", round(5 * (Hn(5) - 1), 5), round(6 + 5 / 12, 5), 1e-9)
chk("R(100)", round(100 * (Hn(100) - 1), 4), 418.7378, 1e-4)

def divis(n):
    return [[1 if m % d == 0 else 0 for m in range(1, n + 1)] for d in range(1, n + 1)]

def det_exact(A):
    A = [[Fraction(x) for x in row] for row in A]
    n, det = len(A), Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            if A[r][c]:
                f = A[r][c] * inv
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return det

def redheffer(n, w=None):
    A = divis(n)
    for d in range(1, n + 1):
        A[d - 1][0] = 1 if w is None else w(d)
    return A

chk("spectrum of the raw table is {1} (det = 1, n=40)", int(det_exact(divis(40))), 1)
K = np.array(divis(60), dtype=float)
chk("dim ker(K-I) at n=60 = #odd <= 60", 60 - np.linalg.matrix_rank(K - np.eye(60)), 30)
chk("det R_n for n=1..12", [int(det_exact(redheffer(n))) for n in range(1, 13)],
    [1, 0, -1, -1, -2, -1, -2, -2, -2, -1, -2, -2])
mu = np.zeros(20001, dtype=np.int8); mu[1] = 1
comp = np.zeros(20001, dtype=bool); pr = []
for i in range(2, 20001):
    if not comp[i]:
        pr.append(i); mu[i] = -1
    for p in pr:
        if i * p > 20000: break
        comp[i * p] = True
        if i % p == 0:
            mu[i * p] = 0; break
        mu[i * p] = -mu[i]
for n in (30, 50, 80):
    chk(f"det R_{n} = M({n})", int(det_exact(redheffer(n))), int(mu[1:n + 1].sum()))
for s in (2.0, 3.0):
    A = np.array([[float(x) for x in r] for r in divis(30)])
    for d in range(1, 31):
        A[d - 1][0] = d ** (-s)
    chk(f"det with w_k=k^-{s:.0f} = truncated 1/zeta", round(float(np.linalg.det(A)), 9),
        round(float(sum(int(mu[k]) * k ** (-s) for k in range(1, 31))), 9), 1e-8)

# ----------------------------------------------------------------- sections 2-3
print("\n=== 2-3. the arches and the five laws ===")
chk("arch S=12", [a * (12 - a) for a in range(2, 11)], [20, 27, 32, 35, 36, 35, 32, 27, 20])
chk("arch S=11", [a * (11 - a) for a in range(2, 10)], [18, 24, 28, 30, 30, 28, 24, 18])
Nf = lambda s, h: (s * s - h * h) // 4
chk("(a) N(14,H)-N(12,H)", [Nf(14, h) - Nf(12, h) for h in (0, 2, 4)], [13, 13, 13])
chk("(a) 49,48,45", [Nf(14, 0), Nf(14, 2), Nf(14, 4)], [49, 48, 45])
chk("(b) 35+45 = 48+32", (Nf(12, 2) + Nf(14, 4), Nf(14, 2) + Nf(12, 4)), (80, 80))
chk("(c) mean of the four neighbours of (12,2)",
    (Nf(10, 2) + Nf(12, 0) + Nf(12, 4) + Nf(14, 2)) // 4, 35)
chk("(d) S=20 row", [a * (20 - a) for a in range(10, 1, -1)],
    [100, 99, 96, 91, 84, 75, 64, 51, 36])
chk("(d) depths are squares", [100 - a * (20 - a) for a in range(10, 1, -1)],
    [0, 1, 4, 9, 16, 25, 36, 49, 64])

# ----------------------------------------------------------------- section 4
print("\n=== 4. Fermat ===")
for p, q, Nc, st in [(999983, 1000003, 999985999949, 0),
                     (100003, 10000019, 1000031900057, 4049995),
                     (1009, 1000000007, 1009000007063, 498996018),
                     (101, 10000000019, 1010000001919, 4998995072)]:
    chk(f"{p} x {q}", p * q, Nc)
    r = isqrt(Nc); cs = r if r * r == Nc else r + 1
    chk(f"  steps", (p + q) // 2 - cs, st)
t6 = tau_sieve(10**6)
for k in (3, 4, 5, 6):
    lo, hi = 10 ** (k - 1), 10 ** k
    avg = float(np.mean([(int(t6[n]) + 1) // 2 for n in range(lo, hi + 1)]))
    chk(f"avg arch-spectrum size on 10^{k-1}..10^{k}", round(avg, 2),
        {3: 3.67, 4: 4.81, 5: 5.96, 6: 7.11}[k], 0.02)

# ----------------------------------------------------------------- section 5
print("\n=== 5. the 2041 table and the repair ===")
def fac(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
gaps = {x: x * x - 2041 for x in range(46, 55)}
chk("gaps 46..54", list(gaps.values()), [75, 168, 263, 360, 459, 560, 663, 768, 875])
par = lambda n: "".join(str(fac(n).get(p, 0) % 2) for p in (2, 3, 5, 7))
chk("parity vectors 46,47,49,51,53,54",
    [par(gaps[x]) for x in (46, 47, 49, 51, 53, 54)],
    ["0100", "1101", "1010", "0011", "0100", "0011"])
chk("75*768 = 240^2", (75 * 768, isqrt(75 * 768) ** 2), (57600, 57600))
chk("gcd(2438-240, 2041)", gcd(46 * 53 - 240, 2041), 157)
chk("the 51,54 pair gives 13", gcd(51 * 54 - isqrt(560 * 875), 2041), 13)
x = 46
while (x * x - 2041) != isqrt(x * x - 2041) ** 2:
    x += 1
chk("Fermat needs 39 steps", x - 46, 39)

Nb = 2578692190013
chk("N_b = 5519 * 467239027", 5519 * 467239027, Nb)
x0 = isqrt(Nb) + 1
chk("x0 = ceil(sqrt N_b)", x0, 1605831)
chk("roots mod 7, indexed from x0", sorted({i for i in range(7) if ((x0 + i) ** 2 - Nb) % 7 == 0}), [0, 1])
import sympy
fb = [p for p in sympy.primerange(3, 3000) if pow(Nb, (p - 1) // 2, p) == 1]
chk("odd factor-base primes below 3000", len(fb), 222)
adds = sum(len({i for i in range(p) if ((x0 + i) ** 2 - Nb) % p == 0}) * (300000 // p) for p in fb)
chk("striped additions", adds, 484562)
chk("trial division / sieve ratio", round(300000 * 222 / adds), 137)

def dickman(u, h=1e-4):
    us = np.arange(0.0, max(u, 1.0) + 1 + h, h)
    rho = np.ones_like(us); integ = 0.0
    for i in range(1, len(us)):
        uu = us[i]
        if uu <= 1.0: continue
        u0 = us[i - 1]
        f1 = np.interp(max(u0 - 1, 0), us[:i], rho[:i]) / max(u0, 1e-12)
        f2 = np.interp(uu - 1, us[:i], rho[:i]) / uu
        integ += 0.5 * (f1 + f2) * h
        rho[i] = 1.0 - integ
    return float(np.interp(u, us, rho))
for u, want in [(1.5, 0.5945), (2.0, 0.3069), (2.5, 0.1303), (3.0, 0.0486), (3.5, 0.0162), (4.0, 0.0049)]:
    chk(f"rho({u})", round(dickman(u), 4), want, 5e-5)
for u, want in [(3.19, 0.0324), (2.61, 0.1062), (1.73, 0.4519)]:
    chk(f"rho({u}) in the measured table", round(dickman(u), 4), want, 5e-5)
gm = math.exp(sum(math.log(abs(x * x - Nb)) for x in range(x0 - 100000, x0 + 100000)) / 200000)
chk("QS geometric-mean gap over x0 +- 1e5", round(gm / 1e11, 1), 1.2, 0.05)
chk("  its u = log(gm)/log 3000", round(math.log(gm) / math.log(3000), 2), 3.19, 0.01)
chk("CFRAC / QS smooth-rate ratio", round(0.388 / 0.0284, 1), 13.7, 0.05)

# ----------------------------------------------------------------- section 6
print("\n=== 6. rho^ord and the error term ===")
tX = tau_sieve(XLIM)
idx = np.arange(1, XLIM + 1)
rho = np.zeros(XLIM + 1, dtype=np.int32)
odd = idx % 2 == 1
rho[1:][odd] = tX[1:][odd]
q4 = idx[idx % 4 == 0]
rho[q4] = tX[q4 // 4]
brute = lambda n: sum(1 for d in range(1, n + 1) if n % d == 0 and (n // d - d) % 2 == 0)
chk("rho^ord formula vs direct search, n < 3000",
    all(int(rho[n]) == brute(n) for n in range(1, 3000)), True)

def odd_le(y): return (y + 1) // 2
def A_oddodd(X):
    r = isqrt(X); s = 0; d = 1
    while d <= r:
        s += odd_le(X // d); d += 2
    return 2 * s - odd_le(r) ** 2
for e, want in ((6, 81), (9, 285), (12, 1504)):
    if FAST and e > 9: continue
    X = 10 ** e
    err = A_oddodd(X) + D_hyp(X // 4) - (0.5 * X * math.log(X) + (G - 0.5) * X)
    chk(f"sum rho^ord - main term at 1e{e}", round(err), want, 1)

n = np.arange(1, XLIM + 1, dtype=np.float64)
Dc = np.cumsum(tX[1:].astype(np.int64))
Delta = Dc - (n * np.log(n) + (2 * G - 1) * n)
Rc = np.cumsum(rho[1:].astype(np.int64))
Earch = Rc - (0.5 * n * np.log(n) + (G - 0.5) * n)
r2 = np.zeros(XLIM + 1, dtype=np.int32)
for d in range(1, XLIM + 1, 2):
    r2[d::d] += 4 if d % 4 == 1 else -4
Ecirc = np.cumsum(r2[1:].astype(np.int64)) + 1 - math.pi * n
chk("sanity: lattice points in the disc of area 100pi", int(np.cumsum(r2[1:101])[-1] + 1), 317)

def decade_max(E, ks):
    return [float(np.abs(E[10 ** k - 1:10 ** (k + 1)]).max()) for k in ks]
ks = [3, 4, 5] if FAST else [3, 4, 5, 6]
tops = [10 ** (k + 1) for k in ks]
for name, E, want, wexp in (("circle", Ecirc, [50.0, 104.9, 209.6, 419.9], 0.308),
                            ("divisor", Delta, [35.5, 85.2, 169.2, 296.6], 0.307),
                            ("arch", Earch, [27.17, 56.87, 119.12, 234.34], 0.312)):
    got = decade_max(E, ks)
    chk(f"{name}: raw decade maxima", [round(v, 1) for v in got],
        [round(v, 1) for v in want[:len(ks)]], 0)
    chk(f"{name}: normalised by the decade top",
        [round(v / t ** 0.25, 2) for v, t in zip(got, tops)],
        [round(v / t ** 0.25, 2) for v, t in zip(want[:len(ks)], tops)], 0)
    if not FAST:
        chk(f"{name}: decade-max exponent",
            round(math.log(got[-1] / got[0]) / math.log(1000), 3), wexp, 0.001)
    lo = 10 ** 5
    chk(f"{name}: std of E/x^(1/4) over 1e5..{XLIM:.0e}",
        round(float(np.std(E[lo:] / n[lo:] ** 0.25)), 2),
        {"circle": 1.59, "divisor": 0.98, "arch": 0.79}[name], 0.02)
chk("skew decomposition sums to", round(-0.004 + 0.048 + 0.005 + 0.056, 3), 0.105, 5e-4)
chk("third moment / sigma^3 = the skew", round(0.105 / 0.79 ** 3, 2), 0.21, 0.02)
chk("mean-square decomposition", round(0.1725 + 0.3232 - 0.0691, 4), 0.4266, 1e-4)

# ----------------------------------------------------------------- section 7
print("\n=== 7. wobbles, gamma, boundaries, Voronoi ===")
chk("count under ab <= 100", int(t100[1:].sum()), 482)
chk("area at x=100", round(100 * math.log(100) + (2 * G - 1) * 100, 1), 476.0, 0.05)
run = 0.0; tbl = {}
for a in range(1, 101):
    run += 100 / a - math.floor(100 / a) - 0.5
    tbl[a] = run
for a, want in ((1, -0.500), (3, -1.167), (6, -2.000), (8, -2.214), (11, -3.512)):
    chk(f"running wobble sum at a={a}", round(tbl[a], 3), want, 1e-3)
chk("running sum at x=1e4",
    round(sum(10**4 / a - (10**4 // a) - 0.5 for a in range(1, 101)), 2), -10.22, 0.01)
for e, want in ((3, 0.583529), (5, 0.577354), (7, 0.577225)):
    if FAST and e > 5: continue
    m = 10 ** e
    if e <= 5:
        W = sum(m / d - (m // d) for d in range(1, m + 1))
    else:
        W = m * (math.log(m) + G + 1 / (2 * m) - 1 / (12 * m * m)) - D_hyp(m)
    chk(f"1 - W(n)/n at 1e{e}", round(1 - W / m, 6), want, 1e-6)
a = np.arange(1, XLIM + 1, dtype=np.float64)
for slope, c, name, want, wexp in ((1.5, 0.25, "slope 3/2", [0.25, 0.25, 0.25, 0.25], 0.000),
                                   ((1 + 5 ** 0.5) / 2, 0.5, "slope phi", [1.01, 1.23, 1.51, 1.71], 0.076),
                                   (math.pi, 0.5, "slope pi", [32.64, 38.56, 39.17, 40.69], 0.032)):
    E = np.cumsum(np.floor(slope * a)) - (slope * a * (a + 1) / 2 - c * a)
    got = decade_max(E, ks)
    chk(f"{name}: decade maxima", [round(v, 2) for v in got],
        [round(v, 2) for v in want[:len(ks)]], 0)
    if not FAST and got[0] > 0:
        chk(f"{name}: exponent", round(math.log(got[-1] / got[0]) / math.log(1000), 3), wexp, 0.002)
xv = 5000.3
tv = tau_sieve(400000).astype(np.float64)
tr = []
for Kt in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 4 * 10 ** 5):
    m = np.arange(1, Kt + 1)
    tr.append(xv ** 0.25 / (math.pi * math.sqrt(2)) *
              float(np.sum(tv[1:Kt + 1] * m ** -0.75 * np.cos(4 * math.pi * np.sqrt(m * xv) - math.pi / 4))))
chk("Voronoi truncations at x=5000.3", [round(v, 2) for v in tr], [6.18, 9.26, 16.33, 15.05, 15.00], 0)
chk("measured Delta(5000.3)",
    round(sum(int(xv) // d for d in range(1, int(xv) + 1)) - (xv * math.log(xv) + (2 * G - 1) * xv), 3),
    14.976, 1e-3)

# ----------------------------------------------------------------- section 8
print("\n=== 8. the soft edge ===")
for e, want in ((3, -6.9e-6), (4, -6.9e-7), (5, -6.9e-8)):
    if FAST and e > 4: continue
    X = 10 ** e; M = 40 * X
    tt = tau_sieve(M).astype(np.float64)
    s = float(np.sum(tt[1:] * np.exp(-np.arange(1, M + 1) / X)))
    err = s - (X * (math.log(X) + G) + 0.25)
    chk(f"soft-edge error at X=1e{e}", float(f"{err:.1e}"), want, 5e-2, rel=True)
    chk(f"  equals -1/(144X) at X=1e{e}", round(err * 144 * X, 3), -1.0, 5e-3)

# ----------------------------------------------------------------- section 9
print("\n=== 9. circles ===")
for r, want in ((0.15, 0.070686), (0.25, 0.196350), (0.35, 0.384845), (0.50, 0.785398)):
    chk(f"pi r^2 at r={r}", round(math.pi * r * r, 6), want, 1e-6)

print()
if FAILS:
    print(f"FAILED: {len(FAILS)} check(s)")
    for f in FAILS:
        print("   -", f)
    sys.exit(1)
print("all checks passed" + ("   (--fast: sieve limit 10^6, the 10^7 and 10^12 rows are skipped)" if FAST else ""))
sys.exit(0)
