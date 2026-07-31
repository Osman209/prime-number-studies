"""
verify_angle_lock.py -- regenerates every table and number in
"The Rational-Angle Lock", the note that closes Appendix A item 2 of
"The Odd Division Table as a Single Column".

    python verify_angle_lock.py            # default cutoff 2e7
    python verify_angle_lock.py --fast     # cutoff 4e6, same conclusions, looser digits

Checks, in order:
   1  the Ramanujan sum identity c_q(1) = mu(q), which is where mu enters
   2  the residue itself: Psi_{a,q}(x)/x -> mu(q)/phi(q), section 3.2
   3  the transfer lemma isolated prime-free on a_n = 1, section 3.3
   4  the edge law at 2 pi a/q: fitted exponent beta = alpha
   5  the edge law's amplitude ratio C_q / C_2 = -mu(q)/phi(q)
   6  checks 4 and 5 again, in the odd convention of section 4 of the paper
   7  the universal prime-free cusp: rescaling by phi(q)/mu(q) collapses every q
   8  which rational angle hosts the global minimum, by enumeration
   9  the crossing alpha* where pi hands over to 2 pi/5
  10  the paper's own two numbers, reproduced in both conventions
  11  cutoff stability of alpha*, and of a fitted exponent

Every printed number carries its summation cutoff, since none of them is
cutoff-free -- see the third method rule in the repository README.

Exits nonzero if any check fails.

Python 3, numpy.
"""
from __future__ import annotations
import sys
from math import gcd, pi, gamma, cos
import numpy as np

FAIL: list[str] = []
FAST = "--fast" in sys.argv
CUTOFF = 4_000_000 if FAST else 20_000_000
TOL_BETA = 0.05 if FAST else 0.03
TOL_AMP = 0.03 if FAST else 0.01


# ------------------------------------------------------------------ helpers
def mu(n: int) -> int:
    r, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    return -r if m > 1 else r


def phi(n: int) -> int:
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    return r - r // m if m > 1 else r


def von_mangoldt_support(N: int):
    """(n, Lambda(n)) over the prime powers n <= N, sorted."""
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p:: p] = False
    ns, lams = [], []
    for p in np.nonzero(sieve)[0]:
        q, lp = int(p), float(np.log(p))
        while q <= N:
            ns.append(q)
            lams.append(lp)
            q *= int(p)
    ns = np.asarray(ns, dtype=np.float64)
    lams = np.asarray(lams, dtype=np.float64)
    o = np.argsort(ns)
    return ns[o], lams[o]


def weights(ns, lams, alpha):
    return 2.0 * lams * ns ** (-(1.0 + alpha))


def symbol(ns, w, theta):
    return float(np.cos(ns * theta) @ w)


def cusp(ns, w, theta0, eps):
    """symmetric second difference -- kills the smooth linear background."""
    return 0.5 * (symbol(ns, w, theta0 + eps) + symbol(ns, w, theta0 - eps)) \
        - symbol(ns, w, theta0)


def check(name, ok, detail=""):
    print(f"    {'ok  ' if ok else 'FAIL'}  {name}{('  ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(name)


# ------------------------------------------------------------------ setup
print(__doc__.strip().splitlines()[1])
print(f"\nsummation cutoff N = {CUTOFF:,}")
N_ALL, L_ALL = von_mangoldt_support(CUTOFF)
ODD = N_ALL % 2 == 1
N_ODD, L_ODD = N_ALL[ODD], L_ALL[ODD]
print(f"prime powers: {len(N_ALL):,}   odd only: {len(N_ODD):,}")

EPS = (1e-3, 3e-5)
QS = (2, 3, 5, 6, 7, 10, 11, 15, 30, 4, 9, 12)


# ------------------------------------------------------------------ 1
print("\n[1] Ramanujan sum: tau(chi_0 mod q) = c_q(1) = mu(q)")
bad = []
for q in range(2, 41):
    c = sum(np.exp(2j * pi * b / q) for b in range(1, q + 1) if gcd(b, q) == 1)
    if abs(c.real - mu(q)) > 1e-9 or abs(c.imag) > 1e-9:
        bad.append(q)
check("c_q(1) = mu(q) for q = 2..40", not bad, f"max |imag| checked, offenders {bad}")


# ------------------------------------------------------------------ 2
print("\n[2] the residue: Psi_{a,q}(x)/x -> mu(q)/phi(q)   (section 3.2)")
print(f"    {'q':>4} " + " ".join(f"{x:>11.0e}" for x in (1e5, 1e6, 3e6)) + f" {'mu/phi':>9}")
ok2 = True
for q in QS:
    row = []
    for x in (1e5, 1e6, 3e6):
        m = N_ALL <= x
        row.append(float(L_ALL[m] @ np.cos(2 * pi * N_ALL[m] / q)) / x)
    print(f"    {q:>4} " + " ".join(f"{v:>11.4f}" for v in row) + f" {mu(q) / phi(q):>9.4f}")
    ok2 &= abs(row[-1] - mu(q) / phi(q)) < 0.02
check("Psi_{a,q}(x)/x lands on mu(q)/phi(q), and on 0 when mu(q) = 0", ok2,
      "the pole of -L'/L at s = 1 through chi_0, weighted by tau(chi_0) = mu(q)")


# ------------------------------------------------------------------ 3
print("\n[3] the transfer lemma, prime-free on a_n = 1   (section 3.3)")
A_T = 0.5
tgt = gamma(-A_T) * cos(pi * A_T / 2)
vals = []
for e in (1e-3, 3e-4, 1e-4, 3e-5):
    s = 0.0
    for lo in range(1, CUTOFF + 1, 2_000_000):
        n = np.arange(lo, min(lo + 2_000_000, CUTOFF + 1), dtype=np.float64)
        s += float(np.sum(n ** (-(1.0 + A_T)) * (np.cos(n * e) - 1.0)))
    vals.append(s / e ** A_T)
print(f"    alpha = {A_T}, terms to {CUTOFF:,}, eps = 1e-3 .. 3e-5")
print("    measured / eps^alpha:  " + ",  ".join(f"{v:.3f}" for v in vals)
      + f"    against Gamma(-a)cos(pi a/2) = {tgt:.4f}")
check("approaches the prime-free constant, from above, slowly",
      all(tgt < v < 0.90 * tgt for v in vals) and vals[0] < vals[-1],
      "O(1) relative remainder on a vanishing leading term -- see section 3.3")


# ------------------------------------------------------------------ 4,5
def edge_table(ns, lams, alpha, label):
    w = weights(ns, lams, alpha)
    rows = {}
    for q in QS:
        t0 = 2 * pi / q
        d = [cusp(ns, w, t0, e) for e in EPS]
        b = float(np.log(abs(d[1] / d[0])) / np.log(EPS[1] / EPS[0])) if d[0] and d[1] else float("nan")
        rows[q] = (d[1] / EPS[1] ** alpha, b)
    c2 = rows[2][0]
    print(f"\n    {label}, alpha = {alpha}, cutoff {CUTOFF:,}")
    print(f"    {'a/q':>6} {'mu/phi':>8} {'beta':>7} {'C_q/C_2':>9} {'-mu/phi':>9}")
    okb = oka = True
    for q, (c, b) in rows.items():
        pred = -mu(q) / phi(q)
        # a fitted exponent is meaningless where mu(q) = 0: there is no leading
        # edge, so b is a fit to rounding noise.  Print a dash, as the note does.
        bs = f"{b:>7.3f}" if mu(q) else f"{'--':>7}"
        print(f"    {'1/' + str(q):>6} {mu(q) / phi(q):>8.4f} {bs} {c / c2:>9.4f} {pred:>9.4f}")
        if mu(q) != 0:
            okb &= abs(b - alpha) < TOL_BETA
            oka &= abs(c / c2 - pred) < TOL_AMP
        else:
            oka &= abs(c / c2) < 1e-2
    return okb, oka


print("\n[4,5] the edge law  f(2 pi a/q + eps) - f(2 pi a/q) ~ C_q |eps|^beta")
okb, oka = edge_table(N_ALL, L_ALL, 0.75, "all prime powers")
check("beta = alpha (squarefree q)", okb, f"tol {TOL_BETA}")
check("C_q/C_2 = -mu(q)/phi(q), and 0 when mu(q) = 0", oka, f"tol {TOL_AMP}")


# ------------------------------------------------------------------ 6
print("\n[6] the same law in the odd convention of section 4")
okb, oka = edge_table(N_ODD, L_ODD, 0.75, "odd prime powers")
check("odd convention: beta = alpha", okb)
check("odd convention: C_q/C_2 = -mu(q)/phi(q)", oka,
      "dropping the 2-tower subtracts a term holomorphic at s = 1")


# ------------------------------------------------------------------ 7
print("\n[7] one universal cusp, scaled by one arithmetic number")
alpha = 0.75
w = weights(N_ALL, L_ALL, alpha)
uni = 2.0 * gamma(-alpha) * cos(pi * alpha / 2)   # the 2 is the one in front of the sum
ratios = []
for q in (2, 3, 5, 6, 7, 11, 30):
    t0 = 2 * pi / q
    for e in (1e-3, 3e-4, 1e-4, 3e-5):
        c = cusp(N_ALL, w, t0, e) / e ** alpha
        ratios.append(c * phi(q) / mu(q) / uni)
r = np.array(ratios)
print(f"    rescaled measured / Gamma(-a)cos(pi a/2) form:"
      f" mean {r.mean():.4f}, sd {r.std():.4f}, range {r.min():.3f}..{r.max():.3f}")
check("phi(q)/mu(q) rescaling collapses every q onto one cusp",
      abs(r.mean() - 1) < 0.10 and r.std() < 0.12,
      "the shape carries no prime information; all of it is in mu(q)/phi(q)")


# ------------------------------------------------------------------ 8
print("\n[8] which rational angle hosts the global minimum")
QMAX = 40
cands = [(a, q, 2 * pi * a / q)
         for q in range(2, QMAX + 1) if mu(q) == -1
         for a in range(1, q // 2 + 1) if gcd(a, q) == 1]
COS = np.empty((len(cands), len(N_ALL)))
for i, (_, _, t) in enumerate(cands):
    COS[i] = np.cos(N_ALL * t)
print(f"    {len(cands)} downward-cusp candidates (mu(q) = -1, q <= {QMAX})")
print(f"    {'alpha':>6} {'winner':>9} {'f_min':>11}   runners-up")
wins = {}
for a_ in (0.50, 0.60, 0.70, 0.735, 0.745, 0.75, 0.90, 1.00):
    v = COS @ weights(N_ALL, L_ALL, a_)
    o = np.argsort(v)
    wins[a_] = cands[o[0]]
    rest = ", ".join(f"{cands[j][0]}/{cands[j][1]}:{v[j]:.4f}" for j in o[1:4])
    print(f"    {a_:>6.3f} {str(cands[o[0]][0]) + '/' + str(cands[o[0]][1]):>9}"
          f" {v[o[0]]:>11.6f}   {rest}")
check("theta = pi wins below the crossing", all(wins[a][:2] == (1, 2) for a in (0.50, 0.60, 0.70, 0.735)))
check("theta = 2 pi/5 wins above it", all(wins[a][:2] == (1, 5) for a in (0.745, 0.75, 0.90, 1.00)))

# and in the odd convention there is no handover at all
COS_O = np.empty((len(cands), len(N_ODD)))
for i, (_, _, tt) in enumerate(cands):
    COS_O[i] = np.cos(N_ODD * tt)
ratios, okodd = [], True
for a_ in (0.30, 0.50, 0.70, 0.75, 0.90, 0.99):
    v = COS_O @ weights(N_ODD, L_ODD, a_)
    o = np.argsort(v)
    okodd &= cands[o[0]][:2] == (1, 2)
    ratios.append(v[o[1]] / v[o[0]])
print(f"    odd convention, alpha = 0.30 .. 0.99: winner is theta = pi throughout,"
      f" runner-up depth {min(ratios):.2f}-{max(ratios):.2f} of the winner's")
check("odd convention: no handover in the tested range", okodd and max(ratios) < 0.7,
      "so the fourth paper's 0.71 is not the odd-convention crossing either")


# ------------------------------------------------------------------ 9,11
print("\n[9] the crossing alpha*")
c2v, c5v = np.cos(N_ALL * pi), np.cos(N_ALL * 2 * pi / 5)


def gap(a_, ns=N_ALL, lams=L_ALL, u=None, v=None):
    w_ = weights(ns, lams, a_)
    return float((u if u is not None else c2v) @ w_) - float((v if v is not None else c5v) @ w_)


def bisect(f, lo, hi, n=50):
    for _ in range(n):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if f(m) < 0 else (lo, m)
    return 0.5 * (lo + hi)


astar = bisect(gap, 0.70, 0.80)
print(f"    alpha* = {astar:.8f}   (cutoff {CUTOFF:,})")
print(f"    at alpha = 0.71 the gap is still {gap(0.71):+.4e}, i.e. pi is winning comfortably")
check("alpha* is near 0.7400, not 0.71", abs(astar - 0.74005) < 2e-3,
      "the paper's 'about alpha = 0.71' is corrected here")

# ------------------------------------------------------------------ 10
print("\n[10] the paper's own numbers, reproduced")
w75 = weights(N_ALL, L_ALL, 0.75)
v_all = symbol(N_ALL, w75, 2 * pi / 5)
v_odd = symbol(N_ODD, weights(N_ODD, L_ODD, 0.75), pi)
print(f"    all prime powers, alpha = 0.75:  f(2 pi/5) = {v_all:.6f}   paper: -0.58633")
print(f"    odd prime powers, alpha = 0.75:  f(pi)     = {v_odd:.6f}   paper: -1.158989")
check("all-prime-power minimum matches the paper", abs(v_all + 0.58633) < 5e-4)
check("odd-convention minimum matches the paper", abs(v_odd + 1.158989) < 5e-4)


print("\n[11] cutoff stability")
half = CUTOFF // 5
ns2, ls2 = von_mangoldt_support(half)
c2b, c5b = np.cos(ns2 * pi), np.cos(ns2 * 2 * pi / 5)
astar2 = bisect(lambda a_: gap(a_, ns2, ls2, c2b, c5b), 0.70, 0.80)
print(f"    alpha* = {astar:.6f} at cutoff {CUTOFF:,},"
      f"  {astar2:.6f} at cutoff {half:,}   (shift {abs(astar - astar2):.2e})")
check("alpha* stable under a fivefold change of cutoff", abs(astar - astar2) < 5e-4)
w_a, w_b = weights(N_ALL, L_ALL, 0.75), weights(ns2, ls2, 0.75)
b_a = float(np.log(abs(cusp(N_ALL, w_a, pi, EPS[1]) / cusp(N_ALL, w_a, pi, EPS[0])))
            / np.log(EPS[1] / EPS[0]))
b_b = float(np.log(abs(cusp(ns2, w_b, pi, EPS[1]) / cusp(ns2, w_b, pi, EPS[0])))
            / np.log(EPS[1] / EPS[0]))
print(f"    beta at theta = pi, alpha = 0.75:  {b_a:.4f} at cutoff {CUTOFF:,},"
      f"  {b_b:.4f} at cutoff {half:,}   (target alpha = 0.75)")
check("the fitted exponent approaches alpha as the cutoff grows",
      abs(b_a - 0.75) <= abs(b_b - 0.75),
      "the fourth paper saw the same drift and read it correctly")


# ------------------------------------------------------------------ verdict
print("\n" + "-" * 68)
if FAIL:
    print("FAILURES:", "; ".join(FAIL))
    sys.exit(1)
print("all checks passed")
