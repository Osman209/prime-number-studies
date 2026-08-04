"""
verify_phase_masking.py -- regenerates every table and number in
"When the Lens Goes Blind".

    ZETA_ZEROS=/path/to/zeros.txt python verify_phase_masking.py
    ZETA_ZEROS=/path/to/zeros.txt python verify_phase_masking.py --fast

Checks, in order:
   1  the closed rate form, tanh r(d) = 2 beta d / (d^2 + beta^2 + T^2)
   2  the exact symmetry d -> R^2/d of the rate
   3  the sech form on the logarithmic axis
   4  the phase at the matched lens, phi(R) = -pi/2 exactly
   5  the quartet identity C_n = 4 - 4 cosh(n r) cos(n phi)
   6  the n mod 4 cycle at d = R
   7  the same cycle in the FULL coefficient, with real ordinates
   8  the large-d asymptotics, and why masking works
   9  N(d) linear in d, and the U shape with its minimum at R
  10  N(d) r(d) = log M + O(1)
  11  the symmetry in N(d): approximate, and its sign at large M
  12  NEGATIVE: phase-tuned lenses do not beat d = R on the first index
  13  NEGATIVE: the far detection pockets are swamped by the background

Every number the note prints comes from here.  Exits nonzero if any check that
could be run fails.

THE ORDINATE LIST.  Read from ZETA_ZEROS: one ordinate per line, increasing.
Without one, checks 7 and 9-13 are skipped and the rest still run.  The note's
tables use the first 100,000 ordinates (gamma = 14.1347 to 74920.8275); such a
list can be taken from Odlyzko's tables (www.dtc.umn.edu/~odlyzko/zeta_tables)
or the LMFDB, or generated with mpmath.zetazero, which is far slower.

A SHORT LIST DOES NOT SILENTLY BECOME A LONG ONE.  numpy slicing is quiet:
gam[:100000] on a 2000-element list returns 2000.  Every slice here is guarded,
and a table needing more ordinates than the list holds is SKIPPED and listed
under "REDUCED RUN" at the end -- never computed on fewer and labelled with the
count it asked for.  A short list is not a failure and does not set the exit
code; it is reported so that no printed number carries a label it did not earn.

Python 3, numpy, sympy.
"""
from __future__ import annotations
import argparse
import os
import sys
import numpy as np
import sympy as sp

FAIL: list[str] = []
SHORT: list[str] = []          # tables skipped for want of ordinates, not failures
BETA, HEIGHT = 0.1, 30.0


def need(gam, M, tag):
    """True if the list really holds M ordinates.

    numpy slicing is silent: gam[:100000] on a 2000-element list returns 2000
    and every table downstream would then be labelled with a count it never
    used.  Nothing in this script slices without asking here first."""
    if len(gam) >= M:
        return True
    SHORT.append((tag, M, len(gam)))
    return False


# ------------------------------------------------------------------ helpers
def q_of(beta: float, T: float, d: float) -> complex:
    """q_a(rho) for rho = 1/2 + beta + iT and lens half-width d = 1/2 - a."""
    rho = 0.5 + beta + 1j * T
    return (rho - (0.5 - d)) / (rho - (0.5 + d))


def quartet_C(n: int, beta: float, T: float, d: float) -> float:
    """the four zeros 1/2 +- beta +- iT, summed as sum (1 - q^n), real part."""
    tot = 0.0 + 0j
    for sb in (1, -1):
        for sT in (1, -1):
            rho = 0.5 + sb * beta + 1j * sT * T
            tot += 1 - ((rho - (0.5 - d)) / (rho - (0.5 + d))) ** n
    return tot.real


def online_theta(gam: np.ndarray, d: float) -> np.ndarray:
    """q = e^{i theta} for a zero at 1/2 + i gamma."""
    return 2.0 * np.arctan(gam / d) - np.pi


def quartet_polar(beta: float, T: float, d: float):
    out = []
    for sb in (1, -1):
        for sT in (1, -1):
            rho = 0.5 + sb * beta + 1j * sT * T
            q = (rho - (0.5 - d)) / (rho - (0.5 + d))
            out.append((np.log(abs(q)), np.angle(q)))
    return out


def S_of(ns: np.ndarray, th: np.ndarray, quart, chunk: int = 400) -> np.ndarray:
    """Re S_n for each n: on-line zeros (both signs) plus the quartet."""
    out = np.empty(len(ns), dtype=float)
    for i in range(0, len(ns), chunk):
        nn = ns[i:i + chunk].astype(float)
        out[i:i + chunk] = 2.0 * np.sum(1.0 - np.cos(nn[:, None] * th[None, :]), axis=1)
        for lg, ph in quart:
            out[i:i + chunk] += 1.0 - np.exp(nn * lg) * np.cos(nn * ph)
    return out


def first_negative(gam: np.ndarray, beta: float, T: float, d: float,
                   n_max: int = 4_000_000) -> int | None:
    """min{n : Re S_n < 0}.  Bracketed by doubling, bisected, then scanned."""
    th, quart = online_theta(gam, d), quartet_polar(beta, T, d)
    lo = 8
    while lo < n_max:
        hi = min(lo * 2, n_max)
        probe = np.unique(np.linspace(lo, hi, 300).astype(np.int64))
        if np.any(S_of(probe, th, quart) < 0):
            break
        lo = hi
    else:
        return None
    a, b = lo, hi
    while b - a > 3000:
        m = (a + b) // 2
        probe = np.unique(np.linspace(a, m, 300).astype(np.int64))
        if np.any(S_of(probe, th, quart) < 0):
            b = m
        else:
            a = m
    fine = np.arange(a, b + 1)
    v = S_of(fine, th, quart)
    k = int(np.argmax(v < 0))
    return int(fine[k]) if v[k] < 0 else None


def head(n: int, title: str) -> None:
    print("\n" + "=" * 76)
    print(f"{n:>2}  {title}")
    print("=" * 76)


# ------------------------------------------------------------------- checks
def check1_rate_form():
    head(1, "the closed rate form:  tanh r(d) = 2 beta d / (d^2 + beta^2 + T^2)")
    beta, d, T = sp.symbols('beta d T', positive=True)
    N = (beta + d) ** 2 + T ** 2
    D = (beta - d) ** 2 + T ** 2
    resid = sp.simplify(sp.tanh(sp.log(N / D) / 2).rewrite(sp.exp)
                        - 2 * beta * d / (d ** 2 + beta ** 2 + T ** 2))
    print(f"   symbolic residual: {resid}")
    if resid != 0:
        FAIL.append("check 1: the closed form is not an identity")
    print(f"\n   numerically, against the approximate form 2 beta d/(d^2+T^2):")
    print(f"   {'d':>8} {'exact tanh r':>16} {'d^2+beta^2+T^2':>16} {'d^2+T^2':>16}")
    for dv in (0.5, 5.0, 30.0, 100.0):
        ex = np.tanh(np.log(abs(q_of(BETA, HEIGHT, dv))))
        f1 = 2 * BETA * dv / (dv ** 2 + BETA ** 2 + HEIGHT ** 2)
        f2 = 2 * BETA * dv / (dv ** 2 + HEIGHT ** 2)
        print(f"   {dv:>8.1f} {ex:>16.10f} {f1:>16.10f} {f2:>16.10f}")
        if abs(ex - f1) > 1e-14:
            FAIL.append(f"check 1: mismatch at d={dv}")
    print("   the last column is the beta << T approximation and is NOT the identity.")


def check2_symmetry():
    head(2, "the rate is invariant under d -> R^2/d, and d = R is its only fixed point")
    beta, d, T = sp.symbols('beta d T', positive=True)
    R2 = beta ** 2 + T ** 2
    f = 2 * beta * d / (d ** 2 + R2)
    resid = sp.simplify(f.subs(d, R2 / d) - f)
    print(f"   symbolic residual of f(R^2/d) - f(d): {resid}")
    if resid != 0:
        FAIL.append("check 2: the symmetry is not exact")
    fixed = sp.solve(sp.Eq(d, R2 / d), d)
    print(f"   fixed points of d -> R^2/d in d > 0: {fixed}   (= R)")
    R = np.hypot(BETA, HEIGHT)
    print(f"\n   {'k':>8} {'r(kR)':>16} {'r(R/k)':>16} {'difference':>13}")
    for k in (0.01, 0.1, 0.5, 2.0):
        r1 = np.log(abs(q_of(BETA, HEIGHT, k * R)))
        r2 = np.log(abs(q_of(BETA, HEIGHT, R / k)))
        print(f"   {k:>8} {r1:>16.12f} {r2:>16.12f} {abs(r1-r2):>13.2e}")
        if abs(r1 - r2) > 1e-13:
            FAIL.append(f"check 2: rates differ at k={k}")


def check3_sech():
    head(3, "on the logarithmic axis x = log(d/R) the rate is a sech dome")
    R = np.hypot(BETA, HEIGHT)
    print(f"   {'x':>8} {'2 beta d/(d^2+R^2)':>20} {'(beta/R) sech x':>18}")
    worst = 0.0
    for x in (-3.0, -1.0, -0.2, 0.0, 0.5, 2.0, 4.0):
        d = R * np.exp(x)
        lhs = 2 * BETA * d / (d ** 2 + R ** 2)
        rhs = (BETA / R) / np.cosh(x)
        worst = max(worst, abs(lhs - rhs))
        print(f"   {x:>8.1f} {lhs:>20.12e} {rhs:>18.12e}")
    print(f"\n   worst difference: {worst:.2e}   (and the profile is even in x)")
    if worst > 1e-15:
        FAIL.append("check 3: the sech form fails")

    # the half-maximum window quoted in section 2
    xh = np.arccosh(2.0)
    lo, hi = np.exp(-xh), np.exp(xh)
    ds = np.linspace(0.01 * R, 20 * R, 400001)
    rate = np.arctanh(2 * BETA * ds / (ds ** 2 + R ** 2))
    keep = ds[rate >= 0.5 * np.arctanh(BETA / R)]
    print(f"\n   half-maximum window, from sech x = 1/2 i.e. |x| <= arcosh 2 = {xh:.4f}:")
    print(f"     predicted   {lo:.4f} R  to  {hi:.4f} R")
    print(f"     measured    {keep.min()/R:.4f} R  to  {keep.max()/R:.4f} R")
    if abs(keep.min() / R - lo) > 2e-3 or abs(keep.max() / R - hi) > 2e-2:
        FAIL.append("check 3: the half-maximum window does not match exp(-+ arcosh 2)")


def check4_phase():
    head(4, "at the matched lens the phase is exactly -pi/2, for every beta and T")
    beta, d, T = sp.symbols('beta d T', positive=True)
    R = sp.sqrt(beta ** 2 + T ** 2)
    print(f"   Re numerator of q_a after rationalising is R^2 - d^2;"
          f" at d = R that is {sp.simplify((R**2 - d**2).subs(d, R))}")
    print(f"   the imaginary part is -2dT < 0, so q = -i|q| and phi = -pi/2.\n")
    print(f"   {'beta':>7} {'T':>8} {'phi(R)':>18} {'-pi/2':>18}")
    for bv, Tv in ((0.1, 30.0), (0.01, 500.0), (0.4, 10.0), (0.001, 3.0)):
        Rv = np.hypot(bv, Tv)
        ph = np.angle(q_of(bv, Tv, Rv))
        print(f"   {bv:>7} {Tv:>8.0f} {ph:>18.14f} {-np.pi/2:>18.14f}")
        if abs(ph + np.pi / 2) > 1e-12:
            FAIL.append(f"check 4: phase is not -pi/2 at beta={bv}, T={Tv}")


def check5_quartet_identity():
    head(5, "the quartet contributes C_n(d) = 4 - 4 cosh(n r) cos(n phi)")
    worst = 0.0
    R = np.hypot(BETA, HEIGHT)
    for d in (0.5, 10.0, R, 100.0):
        q = q_of(BETA, HEIGHT, d)
        r, ph = np.log(abs(q)), np.angle(q)
        for n in (3, 50, 101, 402):
            direct = quartet_C(n, BETA, HEIGHT, d)
            model = 4 - 4 * np.cosh(n * r) * np.cos(n * ph)
            worst = max(worst, abs(direct - model))
    print(f"   worst |direct sum - identity| over d and n: {worst:.2e}")
    if worst > 1e-10:
        FAIL.append("check 5: the quartet identity fails")


def check6_mod4():
    head(6, "at d = R the quartet runs on a strict n mod 4 cycle")
    R = np.hypot(BETA, HEIGHT)
    r = np.log(abs(q_of(BETA, HEIGHT, R)))
    print(f"   {'n':>6} {'n mod 4':>8} {'C_n(R)':>16} {'meaning':>24}")
    for n in (100, 101, 102, 103, 200, 201, 202, 203):
        c = quartet_C(n, BETA, HEIGHT, R)
        m = n % 4
        meaning = {0: "detection", 2: "masking", 1: "no exponential term",
                   3: "no exponential term"}[m]
        print(f"   {n:>6} {m:>8} {c:>16.6f} {meaning:>24}")
        if m == 0 and not c < 0:
            FAIL.append(f"check 6: n={n} (0 mod 4) is not negative")
        if m == 2 and not c > 4:
            FAIL.append(f"check 6: n={n} (2 mod 4) is not strongly positive")
        if m % 2 == 1 and abs(c - 4) > 1e-6:
            FAIL.append(f"check 6: odd n={n} is not exactly 4")


def check7_mod4_full(gam, fast):
    head(7, "the same cycle in the FULL coefficient, with real ordinates")
    if gam is None:
        print("   skipped: no ordinate list (set ZETA_ZEROS)")
        return
    heights = (10.0, 30.0, 50.0) if fast else (10.0, 20.0, 30.0, 50.0, 80.0, 120.0)
    betas = (0.1,) if fast else (0.1, 0.05)
    counts = (2000,) if fast else (2000, 8000)
    res, rows = [], []
    for beta in betas:
        for T in heights:
            for M in counts:
                if not need(gam, M, "check 7"):
                    continue
                R = np.hypot(beta, T)
                n = first_negative(gam[:M], beta, T, R)
                if n is not None:
                    res.append(n % 4)
                    rows.append((beta, T, M, n))
    if not rows:
        print("   skipped: the list is shorter than 2000 ordinates")
        return
    print(f"   {'beta':>6} {'T':>6} {'zeros':>7} {'first n < 0':>12} {'n mod 4':>9}")
    for beta, T, M, n in rows:
        print(f"   {beta:>6} {T:>6.0f} {M:>7} {n:>12} {n % 4:>9}")
    zeros = sum(1 for x in res if x == 0)
    print(f"\n   {zeros} of {len(res)} first-negative indices are 0 mod 4")
    if zeros != len(res):
        FAIL.append("check 7: not every first-negative index is 0 mod 4")


def check8_large_d():
    head(8, "large d: q -> -1, and every C_n turns positive -- the masking mechanism")
    print(f"   {'d':>10} {'q_a':>34} {'-exp(2z/d)':>34}")
    for d in (200.0, 1000.0, 5000.0):
        q = q_of(BETA, HEIGHT, d)
        approx = -np.exp(2 * (BETA + 1j * HEIGHT) / d)
        print(f"   {d:>10.0f} {str(np.round(q, 8)):>34} {str(np.round(approx, 8)):>34}")
    print(f"\n   C_n -> 4(1 - (-1)^n):  0 for even n, 8 for odd n")
    print(f"   {'d':>10} {'C_4':>14} {'C_5':>14} {'C_10':>14} {'C_11':>14}")
    for d in (1e3, 1e5, 1e7):
        print(f"   {d:>10.0e} " + " ".join(
            f"{quartet_C(n, BETA, HEIGHT, d):>14.6e}" for n in (4, 5, 10, 11)))
    print(f"\n   even n:  C_n = 8 n^2 (T^2 - beta^2)/d^2   [positive when T > beta]")
    print(f"   {'n':>5} {'d':>10} {'exact':>15} {'model':>15} {'ratio':>8}")
    for n in (2, 4, 10, 50):
        for d in (2e4, 2e5):
            ex = quartet_C(n, BETA, HEIGHT, d)
            mo = 8 * n * n * (HEIGHT ** 2 - BETA ** 2) / d ** 2
            print(f"   {n:>5} {d:>10.0e} {ex:>15.6e} {mo:>15.6e} {ex/mo:>8.4f}")
            if abs(ex / mo - 1) > 0.01:
                FAIL.append(f"check 8: even-n asymptotic off at n={n}, d={d}")


def check9_U_and_linear(gam, fast):
    head(9, "N(d) is U-shaped with its minimum at d = R, and grows linearly in d")
    if gam is None:
        print("   skipped: no ordinate list")
        return
    R = np.hypot(BETA, HEIGHT)
    M = 2000
    if not need(gam, M, "check 9"):
        print(f"   skipped: needs {M} ordinates, the list has {len(gam)}")
        return
    ks = (0.1, 1.0, 10.0) if fast else (0.005, 0.1, 0.6, 1.0, 2.0, 5.0, 20.0)
    print(f"   the U shape, M = {M}")
    print(f"   {'d/R':>10} {'N(d)':>12}")
    vals = {}
    for k in ks:
        n = first_negative(gam[:M], BETA, HEIGHT, k * R)
        vals[k] = n
        print(f"   {k:>10.3f} {str(n):>12}")
    if 1.0 in vals and vals[1.0] is not None:
        if any(v is not None and v < vals[1.0] for k, v in vals.items() if k != 1.0):
            FAIL.append("check 9: the minimum is not at d = R")
    if fast:
        return
    print(f"\n   linear growth at large d, M = {M}")
    print(f"   {'d/R':>8} {'N(d)':>12} {'N(d)/(d/R)':>12}")
    for k in (1, 2, 10, 50, 200):
        n = first_negative(gam[:M], BETA, HEIGHT, k * R)
        print(f"   {k:>8} {str(n):>12} {(n/k if n else float('nan')):>12.1f}")


def check10_logM(gam, fast):
    head(10, "N_M(d) r(d) = log M + O(1), at the matched lens")
    if gam is None:
        print("   skipped: no ordinate list")
        return
    R = np.hypot(BETA, HEIGHT)
    r = np.log(abs(q_of(BETA, HEIGHT, R)))
    Ms = (250, 2000, 8000) if fast else (250, 500, 1000, 2000, 4000, 8000,
                                         16000, 32000, 64000, 100000)
    Ms = [m for m in Ms if need(gam, m, "check 10")]
    if not Ms:
        print(f"   skipped: needs at least {250} ordinates, the list has {len(gam)}")
        return
    print(f"   {'M':>8} {'N_M(R)':>9} {'N r':>10} {'log M':>10} {'E_M':>9}")
    E = []
    for M in Ms:
        n = first_negative(gam[:M], BETA, HEIGHT, R)
        if n is None:
            continue
        e = n * r - np.log(M)
        E.append(e)
        print(f"   {M:>8} {n:>9} {n*r:>10.3f} {np.log(M):>10.3f} {e:>+9.3f}")
    if E:
        print(f"\n   residual range over {len(E)} values of M: "
              f"[{min(E):+.3f}, {max(E):+.3f}]")
        if max(abs(x) for x in E) > 1.0:
            FAIL.append("check 10: the residual is not O(1)")
    if len(E) >= 5:
        print("   bounded, but NOT settling to zero: the law is")
        print("   log M + O(1), not log M + o(1).")
    else:
        print("   too few values of M to say whether the residual settles;")
        print("   the note's O(1)-not-o(1) reading needs the full ten-point run.")


def check11_symmetry_in_N(gam, fast):
    head(11, "does the exact rate symmetry survive in N(d)?  Approximately, with a sign")
    if gam is None or fast:
        print("   skipped" + ("" if gam is None else " in fast mode"))
        return
    R = np.hypot(BETA, HEIGHT)
    rows = []
    for M in (2000, 4000, 8000, 16000, 32000, 48000, 64000, 100000):
        if not need(gam, M, "check 11"):
            continue
        a = first_negative(gam[:M], BETA, HEIGHT, 0.1 * R)
        b = first_negative(gam[:M], BETA, HEIGHT, 10.0 * R)
        if a is None or b is None:
            continue
        rows.append((M, a, b, a / b))
    if not rows:
        print("   skipped: needs at least 2000 ordinates, "
              f"the list has {len(gam)}")
        return
    print(f"   {'M':>8} {'N(0.1R)':>9} {'N(10R)':>9} {'ratio':>8} {'sign':>6}")
    for M, a, b, q in rows:
        print(f"   {M:>8} {a:>9} {b:>9} {q:>8.4f} {'+' if q > 1 else '-':>6}")

    # the claim is that the symmetry survives APPROXIMATELY -- that is a gate.
    worst = max(abs(q - 1) for _, _, _, q in rows)
    print(f"\n   worst departure from 1: {worst:.4f}")
    if worst > 0.10:
        FAIL.append(f"check 11: the symmetry does not survive in N(d) "
                    f"(worst departure {worst:.4f} > 0.10)")

    # and the turn point is read off the data, not assumed.
    neg = [M for M, _, _, q in rows if q < 1]
    if neg and max(neg) < rows[-1][0]:
        turn = min(M for M, _, _, q in rows if q > 1 and M > max(neg))
        tail = [q for M, _, _, q in rows if M >= turn]
        print(f"   the sign is negative up to M = {max(neg)} and positive at "
              f"every M >= {turn}")
        print(f"   ({len(tail)} values, departures "
              f"{min(tail) - 1:+.4f} to {max(tail) - 1:+.4f})")
        if any(q <= 1 for q in tail):
            FAIL.append("check 11: the late sign is not uniformly positive")
    else:
        print("   no sign change within the values of M available here")


def check12_phase_tuned_lenses(gam, fast):
    head(12, "NEGATIVE: a phase-tuned lens does not beat d = R on the first index")
    if gam is None or fast:
        print("   skipped" + ("" if gam is None else " in fast mode"))
        return
    if not need(gam, 2000, "check 12"):
        print(f"   skipped: needs 2000 ordinates, the list has {len(gam)}")
        return
    R = np.hypot(BETA, HEIGHT)
    base = first_negative(gam[:2000], BETA, HEIGHT, R)
    print(f"   at d = R exactly: N = {base}   (M = 2000)")
    print(f"   {'d':>10} {'d - R':>9} {'N(d)':>8}")
    best = (R, base)
    for dd in R + np.array([-3.0, -1.5, -0.75, -0.25, 0.25, 0.75, 1.5, 3.0]):
        n = first_negative(gam[:2000], BETA, HEIGHT, dd)
        if n is not None:
            print(f"   {dd:>10.4f} {dd-R:>+9.4f} {n:>8}")
            if n < best[1]:
                best = (dd, n)
    ties = sum(1 for dd in R + np.array([-3.0, -1.5, -0.75, -0.25,
                                          0.25, 0.75, 1.5, 3.0])
               if first_negative(gam[:2000], BETA, HEIGHT, dd) == base)
    print(f"\n   best over the scan: d = {best[0]:.4f}, N = {best[1]}"
          + (f"   ({ties} displaced lens(es) tie with it, none below)"
             if ties else "   (no displaced lens reaches it)"))
    if best[1] < base:
        FAIL.append("check 12: some displaced lens beat d = R")
    print("   optimising C_n at a FIXED n is not the same problem as minimising")
    print("   min{n : S_n < 0}, and the crossing is a first passage over all n.")


def check13_far_pockets(gam, fast):
    head(13, "NEGATIVE: the far detection pockets are swamped by the background")
    if gam is None or fast:
        print("   skipped" + ("" if gam is None else " in fast mode"))
        return
    if not need(gam, 2000, "check 13"):
        print(f"   skipped: needs 2000 ordinates, the list has {len(gam)}")
        return
    g = gam[:2000]
    print("   pocket centres d = nT/(pi k) for even n, relative full width 2 beta/T")
    print(f"   {'n':>5} {'k':>3} {'centre':>11} {'C_n':>13} {'B_n':>12} {'S_n':>12} {'<0?':>5}")
    any_neg = False
    for n in (40, 100, 400):
        for k in (1, 2):
            d0 = n * HEIGHT / (np.pi * k)
            ds = np.linspace(d0 * 0.995, d0 * 1.005, 2001)
            cv = np.array([quartet_C(n, BETA, HEIGHT, d) for d in ds])
            j = int(np.argmin(cv))
            d, cn = ds[j], cv[j]
            th = online_theta(g, d)
            bn = 2.0 * np.sum(1.0 - np.cos(n * th))
            sn = bn + cn
            any_neg |= sn < 0
            print(f"   {n:>5} {k:>3} {d:>11.3f} {cn:>13.3e} {bn:>12.4f} "
                  f"{sn:>12.4f} {str(sn < 0):>5}")
    print(f"\n   along a pocket line C_n ~ -8 (pi k beta/T)^2, independent of n, so the")
    print("   pockets never deepen however far out or however large the order.")
    if any_neg:
        FAIL.append("check 13: a far pocket produced a negative coefficient")

    # the note claims the measured centre approaches n T/(pi k); gate it.
    print("\n   centre against the prediction n T/(pi k)")
    print(f"   {'n':>5} {'k':>3} {'measured':>11} {'predicted':>11} {'rel':>10}")
    worst_far = 0.0
    for n in (40, 100, 400):
        for k in (1, 2):
            d0 = n * HEIGHT / (np.pi * k)
            ds = np.linspace(d0 * 0.995, d0 * 1.005, 2001)
            cv = np.array([quartet_C(n, BETA, HEIGHT, d) for d in ds])
            dm = ds[int(np.argmin(cv))]
            rel = abs(dm - d0) / d0
            worst_far = max(worst_far, rel) if n == 400 else worst_far
            print(f"   {n:>5} {k:>3} {dm:>11.3f} {d0:>11.3f} {rel:>10.1e}")
    print("   the departure shrinks with n: the minimum sits where cosh and cos")
    print("   balance, not exactly where the phase realigns.")
    if worst_far > 1e-4:
        FAIL.append(f"check 13: at n = 400 the pocket centre is off by {worst_far:.1e}, "
                    "more than the 1e-4 the note claims")

    print(f"\n   relative full width, measured against 2 beta/T = {2*BETA/HEIGHT:.5f}")
    print(f"   {'n':>6} {'k':>3} {'width':>10}")
    for n, k in ((100, 1), (400, 1), (400, 2)):
        d0 = n * HEIGHT / (np.pi * k)
        ds = np.linspace(0.98 * d0, 1.02 * d0, 40001)
        neg = ds[np.array([quartet_C(n, BETA, HEIGHT, d) for d in ds]) < 0]
        if len(neg) < 2:
            FAIL.append(f"check 13: no pocket found at n={n}, k={k}")
            continue
        w = (neg.max() - neg.min()) / ((neg.max() + neg.min()) / 2)
        print(f"   {n:>6} {k:>3} {w:>10.5f}")
        if abs(w - 2 * BETA / HEIGHT) > 3e-4:
            FAIL.append(f"check 13: pocket width {w:.5f} is not 2 beta/T")


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="skip the expensive sweeps; the symbolic checks still run")
    args = ap.parse_args()

    path = os.environ.get("ZETA_ZEROS")
    gam = None
    if path and os.path.exists(path):
        gam = np.loadtxt(path)
        print(f"ordinates: {len(gam)} from {path}, "
              f"gamma from {gam[0]:.4f} to {gam[-1]:.4f}")
    else:
        print("no ZETA_ZEROS set -- checks 7 and 9-13 will be skipped")
    print(f"planted quartet: rho = 1/2 +- {BETA} +- {HEIGHT}i")

    check1_rate_form()
    check2_symmetry()
    check3_sech()
    check4_phase()
    check5_quartet_identity()
    check6_mod4()
    check7_mod4_full(gam, args.fast)
    check8_large_d()
    check9_U_and_linear(gam, args.fast)
    check10_logM(gam, args.fast)
    check11_symmetry_in_N(gam, args.fast)
    check12_phase_tuned_lenses(gam, args.fast)
    check13_far_pockets(gam, args.fast)

    print("\n" + "=" * 76)
    if SHORT:
        print("REDUCED RUN -- the ordinate list is shorter than these tables need,")
        print("so they were skipped rather than computed on fewer ordinates and")
        print("labelled with the count they asked for:")
        worst = {}
        for tag, M, have in SHORT:
            worst[tag] = (max(M, worst.get(tag, (0, have))[0]), have)
        for tag in sorted(worst, key=lambda t: int(t.split()[1])):
            M, have = worst[tag]
            print(f"  - {tag}: needs up to {M} ordinates, the list has {have}")
        print("The note's tables use 100,000 ordinates; see its reproducibility")
        print("section for where to obtain a list.")
        print("=" * 76)
    if FAIL:
        print("FAILURES:")
        for f in FAIL:
            print("  -", f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
