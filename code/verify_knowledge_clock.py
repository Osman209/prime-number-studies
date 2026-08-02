#!/usr/bin/env python3
"""
verify_knowledge_clock.py
=========================
Regenerates every number quoted in section 10, "The sieve as a clock of
knowledge", of *A Numerical Study of the Division Table*.

The section's claim, in one line: cutting the division table off at row z
turns "prime" into a time-dependent verdict, and the clock is

        u = log x / log z,

with Buchstab's function omega(u) governing the density of survivors and
1/(u*omega(u)) governing the share of survivors that are genuinely prime.

Checks
------
  1  omega solver against its three closed-form values
  2  omega(u) -> e^{-gamma}
  3  frontier identity  Phi(x,z) = pi(x) - pi(z) + 1  for z >= sqrt(x)
  4  set identity: at u = 2 the candidate set IS the set of primes in (z,x]
  5  first false candidate after columns <= p equals p_next^2
  6  small-z regime: Phi(x,z)/x equals the wheel density prod(1-1/p)
  7  Buchstab table: Phi(x,z) / (x prod_{p<=z}(1-1/p))  vs  e^gamma omega(u)
  8  the deviation for u < 2 is exactly the prime-counting secondary term
  9  prime share vs 1/(u omega(u)): ordering, sign of the finite-size excess
 10  share at fixed u = 3 drifts toward the limit as x grows

Usage
-----
    python verify_knowledge_clock.py            # full,  x = 10^8   (~3-6 min, ~1.5 GB peak)
    python verify_knowledge_clock.py --fast     # quick, x = 10^7   (~20 s)

Exits 0 only if every check passes; exits 1 and names the failing checks
otherwise.  Every printed table carries the x it was computed at.
"""

import sys
import math
import numpy as np

FAST = "--fast" in sys.argv
XMAX = 10**7 if FAST else 10**8
GAMMA = 0.5772156649015329
EGAMMA = math.exp(GAMMA)

FAIL = []


def fail(tag, msg):
    FAIL.append(tag)
    print(f"    FAIL [{tag}] {msg}")


# ----------------------------------------------------------------------
# base sieve
# ----------------------------------------------------------------------
print(f"building sieve to {XMAX:,} ...", flush=True)
is_p = np.ones(XMAX + 1, bool)
is_p[:2] = False
for p in range(2, math.isqrt(XMAX) + 1):
    if is_p[p]:
        is_p[p * p :: p] = False
PRIMES = np.nonzero(is_p)[0]
PI = np.cumsum(is_p)  # PI[n] = pi(n)
print(f"  pi({XMAX:,}) = {int(PI[-1]):,}\n")


def Phi(x, z):
    """#{ 2 <= n <= x : n has no prime factor <= z }"""
    m = np.ones(x + 1, bool)
    m[0] = m[1] = False
    for p in PRIMES[PRIMES <= z]:
        m[p::p] = False
    return m


def wheel_density(z):
    ps = PRIMES[PRIMES <= z].astype(float)
    return float(np.prod(1.0 - 1.0 / ps))


# ----------------------------------------------------------------------
# Buchstab omega:   (u omega(u))' = omega(u-1),   omega(u) = 1/u on [1,2]
# ----------------------------------------------------------------------
def omega_grid(umax=14.0, h=1e-4):
    n = int(umax / h) + 1
    u = np.arange(n) * h
    w = np.zeros(n)
    sel = (u >= 1.0) & (u <= 2.0)
    w[sel] = 1.0 / u[sel]
    i2 = int(2.0 / h)
    acc = 1.0  # = u*omega(u) at u = 2
    for i in range(i2 + 1, n):
        f1 = np.interp(u[i - 1] - 1.0, u, w)
        f2 = np.interp(u[i] - 1.0, u, w)
        acc += 0.5 * (f1 + f2) * h
        w[i] = acc / u[i]
    return u, w


UG, WG = omega_grid()
omega = lambda t: float(np.interp(t, UG, WG))

print("CHECK 1 - omega against closed forms")
c1 = [("omega(2)", omega(2.0), 0.5),
      ("omega(3)", omega(3.0), (1 + math.log(2)) / 3),
      ("omega(2.5)", omega(2.5), (1 + math.log(1.5)) / 2.5)]
for name, got, want in c1:
    print(f"    {name:>10} = {got:.8f}   exact {want:.8f}")
    if abs(got - want) > 1e-6:
        fail("1", f"{name} off by {got-want:.2e}")

print("\nCHECK 2 - omega(u) -> e^-gamma")
for u in (6.0, 8.0, 10.0, 12.0):
    print(f"    omega({u:4.1f}) = {omega(u):.9f}")
d = abs(omega(12.0) - math.exp(-GAMMA))
print(f"    e^-gamma      = {math.exp(-GAMMA):.9f}   |omega(12) - e^-gamma| = {d:.2e}")
if d > 1e-7:
    fail("2", f"omega(12) has not settled: {d:.2e}")

# ----------------------------------------------------------------------
print("\nCHECK 3 - frontier identity  Phi(x,z) = pi(x) - pi(z) + 1  for z >= sqrt(x)")
xs3 = [10**5, 10**6, 10**7] + ([] if FAST else [10**8])
for x in xs3:
    z = math.isqrt(x - 1) + 1
    m = Phi(x, z)
    got = int(m.sum())
    want = int(PI[x] - PI[z])
    u = math.log(x) / math.log(z)
    print(f"    x={x:>11,}  z={z:>6,}  u={u:.4f}  Phi={got:>9,}  pi(x)-pi(z)={want:>9,}")
    if got != want:
        fail("3", f"x={x}: {got} != {want}")
    if x == 10**8 and got != 5760226:
        fail("3", f"x=1e8: Phi = {got}, paper 5,760,226")

print("\nCHECK 4 - at u = 2 the candidate SET is exactly the primes in (z, x]")
x = 10**6
z = math.isqrt(x - 1) + 1
cand = np.nonzero(Phi(x, z))[0]
if not np.array_equal(cand, PRIMES[(PRIMES > z) & (PRIMES <= x)]):
    fail("4", "candidate set differs from the prime set")
else:
    print(f"    x={x:,}  z={z:,}: {len(cand):,} candidates, set-identical to the primes")

# ----------------------------------------------------------------------
print("\nCHECK 5 - first false candidate after columns <= p is p_next^2")
for i in range(0, 15):
    p, q = int(PRIMES[i]), int(PRIMES[i + 1])
    m = Phi(min(q * q + 10, XMAX), p)
    surv = np.nonzero(m)[0]
    surv = surv[surv > p]
    first_false = next((int(n) for n in surv if not is_p[n]), None)
    ok = first_false == q * q
    if i < 8:
        print(f"    columns <= {p:>3}: first false candidate = {first_false:>5}  = {q}^2  {'ok' if ok else 'MISMATCH'}")
    if not ok:
        fail("5", f"p={p}: got {first_false}, expected {q*q}")

# ----------------------------------------------------------------------
print("\nCHECK 6 - small z: the clock reading is exactly the wheel density")
x = XMAX
for z in (3, 7, 13, 31):
    frac = int(Phi(x, z).sum()) / x
    w = wheel_density(z)
    print(f"    z={z:>3}  Phi/x = {frac:.7f}   prod(1-1/p) = {w:.7f}   diff = {frac-w:+.2e}")
    if abs(frac - w) > 2e-5:
        fail("6", f"z={z}: wheel density off by {frac-w:.2e}")

# ----------------------------------------------------------------------
print(f"\nCHECK 7 - Buchstab table at x = {XMAX:,}")
print(f"    {'z':>10} {'u':>7} {'measured':>10} {'e^g*omega':>10} {'rel err':>9}")
zs = [7, 31, 101, 331, 1009, 3163, 10007]
if not FAST:
    zs += [31623, 100003, 316231]
rows7 = []
for z in zs:
    if z * z > x * 100:
        continue
    f = int(Phi(x, z).sum())
    u = math.log(x) / math.log(z)
    meas = f / (x * wheel_density(z))
    pred = EGAMMA * omega(u)
    rel = meas / pred - 1
    rows7.append((z, u, meas, pred, rel))
    print(f"    {z:>10,} {u:>7.3f} {meas:>10.5f} {pred:>10.5f} {100*rel:>8.2f}%")
PAPER7 = {7:(9.466,1.00000,1.00000), 31:(5.364,1.00000,1.00000), 101:(3.991,0.99977,0.99999),
          331:(3.175,1.00360,1.00070), 1009:(2.663,0.99870,1.00901), 3163:(2.286,0.95387,0.97504),
          10007:(2.000,0.94618,0.89060), 31623:(1.778,1.06302,1.00185)}
if not FAST:
    for z, u, meas, pred, rel in rows7:
        pu, pm, pp = PAPER7[z]
        if abs(u-pu) > 5e-4 or abs(meas-pm) > 5e-6 or abs(pred-pp) > 5e-6:
            fail("7", f"z={z}: ({u:.3f},{meas:.5f},{pred:.5f}) vs paper ({pu},{pm},{pp})")
for z, u, meas, pred, rel in rows7:
    if u >= 2.2 and abs(rel) > 0.03:
        fail("7", f"z={z} (u={u:.3f}) deviates {100*rel:.2f}% from e^gamma*omega(u)")
print("    (gate applies to u >= 2.2; rows with u < 2.2 are governed by check 8)")

print("\nCHECK 8 - for u < 2 the whole deviation splits into two named factors")
print("    ratio = A * B   with   A = (pi(x)-pi(z)) log x / x   (PNT secondary term)")
print("                          B = 1 / (e^gamma * prod(1-1/p) * log z)   (Mertens remainder)")
print(f"    {'z':>10} {'u':>7} {'ratio':>9} {'A':>9} {'B':>9} {'A*B - ratio':>13}")
zs8 = [int(x ** (1 / 1.9)), int(x ** (1 / 1.7)), int(x ** (1 / 1.4))]
for z in zs8:
    f = int(Phi(x, z).sum())
    u = math.log(x) / math.log(z)
    meas = f / (x * wheel_density(z))
    pred = EGAMMA * omega(u)
    ratio = meas / pred
    A = int(PI[x] - PI[z]) * math.log(x) / x
    B = 1.0 / (EGAMMA * wheel_density(z) * math.log(z))
    print(f"    {z:>10,} {u:>7.3f} {ratio:>9.5f} {A:>9.5f} {B:>9.5f} {A*B-ratio:>+13.2e}")
    PAPER8 = {16237:(1.06208,1.06095), 50802:(1.06068,1.06034), 517947:(1.05352,1.05340)}
    if not FAST and z in PAPER8:
        pr_, pa_ = PAPER8[z]
        if abs(ratio-pr_) > 5e-6 or abs(A-pa_) > 5e-6:
            fail("8", f"z={z}: (ratio {ratio:.5f}, A {A:.5f}) vs paper ({pr_}, {pa_})")
    if abs(A * B - ratio) > 1e-9:
        fail("8", f"z={z}: decomposition does not close ({A*B-ratio:+.2e})")
    if abs(B - 1.0) > 0.02:
        fail("8", f"z={z}: Mertens factor B = {B:.5f} is further from 1 than expected")

# ----------------------------------------------------------------------
print(f"\nCHECK 9 - prime share as the clock runs, x = {XMAX:,}")
print(f"    {'z':>10} {'u':>7} {'candidates':>12} {'primes':>11} {'share':>9} {'1/(u*omega)':>12}")
zs9 = [10007, 1009, 101, 31, 13, 7, 5, 3]
rows9 = []
rows9full = []
for z in zs9:
    f = int(Phi(x, z).sum())
    u = math.log(x) / math.log(z)
    pr = int(PI[x] - PI[z])
    share = pr / f
    lim = 1.0 / (u * omega(u))
    rows9.append((u, share, lim))
    rows9full.append((z, round(u,3), f, pr, share))
    print(f"    {z:>10,} {u:>7.3f} {f:>12,} {pr:>11,} {share:>9.5f} {lim:>12.5f}")
PAPER9 = {2.000:(5760225,5760225,1.00000), 2.663:(8078009,5761286,0.71321),
          3.991:(11909899,5761429,0.48375), 5.364:(15285209,5761444,0.37693),
          7.182:(19180819,5761449,0.30038), 9.466:(22857141,5761451,0.25206),
          16.767:(33333332,5761453,0.17284)}
if not FAST:
    for z, u, f_, pr_, share in rows9full:
        key = round(u, 3)
        if key in PAPER9:
            pf, pp_, ps_ = PAPER9[key]
            if f_ != pf or pr_ != pp_ or abs(share-ps_) > 5e-6:
                fail("9", f"z={z}: ({f_},{pr_},{share:.5f}) vs paper ({pf},{pp_},{ps_})")
if any(rows9[i][1] <= rows9[i + 1][1] for i in range(len(rows9) - 1)):
    fail("9", "prime share is not strictly decreasing in u")
if any(s < lim for _, s, lim in rows9):
    fail("9", "a measured share fell below the asymptotic limit (excess should be positive)")

# ----------------------------------------------------------------------
print("\nCHECK 10 - at fixed u = 3 the share drifts toward 1/(3*omega(3)) = "
      f"{1/(3*omega(3.0)):.5f}")
print(f"    {'x':>12} {'z':>7} {'u':>6} {'share':>9} {'excess':>9}")
xs10 = [10**5, 10**6, 10**7] + ([] if FAST else [10**8])
sh = []
for xx in xs10:
    z = int(round(xx ** (1 / 3)))
    f = int(Phi(xx, z).sum())
    u = math.log(xx) / math.log(z)
    pr = int(PI[xx] - PI[z])
    share = pr / f
    lim = 1.0 / (u * omega(u))
    sh.append(share - lim)
    PAPER10 = {100000:(0.08498,0.98), 1000000:(0.05922,0.82),
               10000000:(0.04893,0.79), 100000000:(0.04178,0.77)}
    pe, pp_ = PAPER10[xx]
    prod = (share - lim) * math.log(xx)
    if abs((share - lim) - pe) > 5e-5:
        fail("10", f"x={xx}: excess {share-lim:+.5f} vs paper {pe:+.5f}")
    if abs(round(prod, 2) - pp_) > 1e-9:
        fail("10", f"x={xx}: excess*log x = {prod:.4f}, paper {pp_}")
    print(f"    {xx:>12,} {z:>7,} {u:>6.3f} {share:>9.5f} {share-lim:>+9.5f}"
          f"   excess x log x = {prod:.4f} (paper {pp_})")
if not all(sh[i] > sh[i + 1] for i in range(len(sh) - 1)):
    fail("10", f"excess did not decrease monotonically: {['%+.5f' % e for e in sh]}")
if sh[-1] <= 0:
    fail("10", "excess changed sign")

# ----------------------------------------------------------------------
print("\n" + "=" * 62)
if FAIL:
    print(f"FAILED checks: {sorted(set(FAIL))}   (run at x = {XMAX:,})")
    sys.exit(1)
print(f"all checks passed   (run at x = {XMAX:,}"
      f"{', --fast' if FAST else ''})")
sys.exit(0)
