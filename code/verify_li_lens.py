"""
verify_li_lens.py -- regenerates every table and number in
"Matching the Lens to the Zero".

    python verify_li_lens.py

Checks, in order:
  1  Proposition 1, S_n(a) = (1-2a) L_n(a), against the classical lambda_1
  2  the differentiation-step caveat of section 2.2, shown rather than asserted
  3  Proposition 2, the optimal lens width d = T
  4  Corollary 3, the factor-T gain over the classical lens
  5  section 3.3, the 30.9x measurement explained as T = 30
  6  section 4, the truncation constants as (1 + 1/b)

Exits nonzero if any check fails.

Python 3, numpy, sympy, mpmath.
"""
from __future__ import annotations
import sys
import numpy as np
import sympy as sp
import mpmath as mp

FAIL: list[str] = []


# ------------------------------------------------------------------ helpers
def logxi(s):
    """log xi(s), with the pole of zeta at s = 1 cancelled against (s-1)."""
    s = mp.mpmathify(s)
    z = mp.mpf(1) if abs(s - 1) < mp.mpf('1e-40') else (s - 1) * mp.zeta(s)
    return (mp.log(mp.mpf(0.5)) + mp.log(s) + mp.log(z)
            - (s / 2) * mp.log(mp.pi) + mp.loggamma(s / 2))


def L_n(n, a, h=mp.mpf('0.02')):
    """(n-1)!^-1 d^n/ds^n [(s-a)^(n-1) log xi(s)] at s = 1-a.

    The explicit step h is not optional -- see check 2."""
    f = lambda s: (s - a) ** (n - 1) * logxi(s)
    return mp.diff(f, mp.mpf(1) - a, n, h=h) / mp.factorial(n - 1)


def S_n(n, a, gammas):
    """sum over zeros of [1 - q_a(rho)^n], summed symmetrically over +-gamma."""
    tot = mp.mpf(0)
    for g in gammas:
        for sgn in (1, -1):
            rho = mp.mpc(mp.mpf(1) / 2, sgn * g)
            tot += 1 - ((rho - a) / (rho - (1 - a))) ** n
    return mp.re(tot)


def rate(beta, T, d):
    """log|q_a(rho)| for rho = 1/2 + beta + iT and lens width d = 1/2 - a."""
    return 0.5 * np.log(((beta + d) ** 2 + T ** 2) / ((beta - d) ** 2 + T ** 2))


# ------------------------------------------------------------------ checks
def check0_map_geometry():
    """section 2.1: the unit circle of q_a IS the critical line, and q_a(1-rho) = 1/q_a(rho).
    Elementary, but the paper asserts them, so they are checked."""
    print("=" * 74)
    print("0  section 2.1: the geometry of the map")
    print("=" * 74)
    rng = np.random.default_rng(0)
    worst_on, worst_sym, worst_rate, n_eq = 0.0, 0.0, 0.0, 0
    for a in (0.0, -0.5, -1.5, -3.5, 2.0, 7.25):
        d = 0.5 - a
        for _ in range(2000):
            T = rng.uniform(-500, 500)
            q_on = (0.5 + 1j * T - a) / (0.5 + 1j * T - (1 - a))
            worst_on = max(worst_on, abs(abs(q_on) - 1))
            beta = rng.uniform(0.01, 0.4) * rng.choice([-1, 1])
            r = 0.5 + beta + 1j * T
            q = (r - a) / (r - (1 - a))
            if abs(abs(q) - 1) == 0.0:
                n_eq += 1
            # the deviation is not merely nonzero: it must equal the closed form of Prop 2
            closed = 0.5 * np.log(((beta + d) ** 2 + T ** 2) / ((beta - d) ** 2 + T ** 2))
            worst_rate = max(worst_rate, abs(np.log(abs(q)) - closed))
            lhs = ((1 - r) - a) / ((1 - r) - (1 - a))
            worst_sym = max(worst_sym, abs(lhs - 1.0 / q))
    print(f"   on the critical line, worst | |q_a| - 1 |          = {worst_on:.3e}")
    print(f"   off it, draws with |q_a| exactly 1 (of 12000)      = {n_eq}")
    print(f"   worst | log|q_a| - closed form of Prop 2 |         = {worst_rate:.3e}")
    print(f"   worst | q_a(1-rho) - 1/q_a(rho) |                  = {worst_sym:.3e}")
    print("   (off the line the deviation can be tiny -- beta/T^2 at the classical lens --")
    print("    so what is checked is the identity, not the size.)")
    if worst_on > 1e-12:
        FAIL.append("check 0: |q_a| != 1 on the critical line")
    if n_eq:
        FAIL.append("check 0: |q_a| = 1 found off the critical line")
    if worst_rate > 1e-12:
        FAIL.append("check 0: log|q_a| does not match the closed form")
    if worst_sym > 1e-10:
        FAIL.append("check 0: the functional-equation symmetry fails")


def check1_normalisation(nzeros=400):
    print("=" * 74)
    print("1  Proposition 1:  S_n(a) = (1-2a) L_n(a)")
    print("=" * 74)
    mp.mp.dps = 30
    lam1 = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
    got = L_n(1, 0)
    print(f"   classical lambda_1 = 1 + gamma/2 - log(4 pi)/2 = {mp.nstr(lam1, 12)}")
    print(f"   L_1(0) from the derivative object              = {mp.nstr(got, 12)}")
    print(f"   difference                                     = {mp.nstr(abs(got-lam1), 3)}"
          f"   (finite-difference step h = 2e-2, error O(h^2))")
    if abs(got - lam1) > mp.mpf('1e-6'):
        FAIL.append("check 1: L_1(0) does not reproduce lambda_1")
    print()
    gam = [mp.im(mp.zetazero(k)) for k in range(1, nzeros + 1)]
    print(f"   zero sum against the derivative side, {nzeros} zeros available")
    print(f"   {'a':>6} {'n':>3} {'(1-2a) L_n':>16} {'S_n(100)':>13} "
          f"{'S_n(200)':>13} {'S_n(400)':>13} {'gap':>9}")
    for a in (0, -0.5, -1.5):
        for n in (1, 2):
            tgt = (1 - 2 * a) * L_n(n, a)
            vals = [S_n(n, a, gam[:k]) for k in (100, 200, nzeros)]
            gaps = [abs(tgt - v) for v in vals]
            closing = gaps[0] > gaps[1] > gaps[2]
            print(f"   {a:>6} {n:>3} {mp.nstr(mp.re(tgt), 10):>16} "
                  + " ".join(f"{mp.nstr(v, 8):>13}" for v in vals)
                  + f" {'closing' if closing else 'NOT':>9}")
            if not closing:
                FAIL.append(f"check 1: gap not closing at a={a}, n={n}")
    print("   (the zero sum converges slowly; what is checked is that the gap closes)")


def check2_step_caveat():
    print("\n" + "=" * 74)
    print("2  section 2.2: the differentiation step is not optional")
    print("=" * 74)
    mp.mp.dps = 40
    lam1 = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
    default = mp.diff(logxi, mp.mpf(1), 1)
    explicit = mp.diff(logxi, mp.mpf(1), 1, h=mp.mpf('0.01'))
    print(f"   lambda_1, closed form        {mp.nstr(lam1, 12)}")
    print(f"   mp.diff, default step        {mp.nstr(default, 12)}   <- wrong")
    print(f"   mp.diff, h = 1e-2            {mp.nstr(explicit, 12)}   <- right")
    if abs(explicit - lam1) > mp.mpf('1e-8'):
        FAIL.append("check 2: explicit step does not recover lambda_1")
    if abs(default - lam1) < mp.mpf('1e-3'):
        FAIL.append("check 2: the default step no longer fails; caveat may be stale")


def check3a_proof_steps():
    """section 3.1: every step of the proof of Proposition 2, not just its conclusion."""
    print("\n" + "=" * 74)
    print("3a  section 3.1: the proof of Proposition 2, step by step")
    print("=" * 74)
    beta, d, T, u, v, R = sp.symbols('beta d T u v R', positive=True)
    a = sp.Rational(1, 2) - d
    rho = sp.Rational(1, 2) + beta + sp.I * T
    steps = []

    # 1  the substitution a = 1/2 - d puts the two factors in the stated form
    steps.append(("1  rho - a          = (beta+d) + iT",
                  sp.simplify(sp.expand(rho - a) - ((beta + d) + sp.I * T))))
    steps.append(("1  rho - (1-a)      = (beta-d) + iT",
                  sp.simplify(sp.expand(rho - (1 - a)) - ((beta - d) + sp.I * T))))

    # 2  |q|^2 = N/D
    N = (beta + d) ** 2 + T ** 2
    D = (beta - d) ** 2 + T ** 2
    q2 = sp.simplify(sp.Abs((rho - a) / (rho - (1 - a))) ** 2)
    steps.append(("2  |q_a|^2          = N/D",
                  sp.simplify(sp.radsimp(q2 - N / D))))

    # 4  the derivative of N/D has numerator 2[(beta+d)D + (beta-d)N]
    dNdD = sp.together(sp.diff(N / D, d))
    steps.append(("4  d/dd (N/D) numerator = 2[(beta+d)D + (beta-d)N]",
                  sp.simplify(sp.numer(dNdD) - 2 * ((beta + d) * D + (beta - d) * N))))

    # 5  the bracket factors as (u+v)(uv + T^2)
    steps.append(("5  u(v^2+T^2) + v(u^2+T^2) = (u+v)(uv + T^2)",
                  sp.simplify(u * (v ** 2 + T ** 2) + v * (u ** 2 + T ** 2) - (u + v) * (u * v + T ** 2))))

    # 6  with u = beta+d, v = beta-d that is 2 beta (beta^2 - d^2 + T^2)
    steps.append(("6  (u+v)(uv+T^2)|_{u=beta+d, v=beta-d} = 2 beta(beta^2 - d^2 + T^2)",
                  sp.simplify(((u + v) * (u * v + T ** 2)).subs({u: beta + d, v: beta - d})
                              - 2 * beta * (beta ** 2 - d ** 2 + T ** 2))))

    # 8  at d = R the two forms collapse
    Rv0 = sp.sqrt(beta ** 2 + T ** 2)
    steps.append(("8  N|_{d=R}         = 2R(R+beta)",
                  sp.simplify(N.subs(d, Rv0) - 2 * Rv0 * (Rv0 + beta))))
    steps.append(("8  D|_{d=R}         = 2R(R-beta)",
                  sp.simplify(D.subs(d, Rv0) - 2 * Rv0 * (Rv0 - beta))))

    # 9  the artanh step, as pure algebra plus a stated domain condition.
    #    artanh(z) = (1/2) log((1+z)/(1-z)) is standard on |z| < 1, so what has to be
    #    checked is that N/D at d = R is (R+beta)/(R-beta), that (x-1)/(x+1) = beta/R
    #    for that x, and that beta/R < 1.
    steps.append(("9a N/D at d=R      = (R+beta)/(R-beta)",
                  sp.simplify((N / D).subs(d, Rv0) - (Rv0 + beta) / (Rv0 - beta))))
    x = (R + beta) / (R - beta)
    steps.append(("9b (x-1)/(x+1)     = beta/R   for that x",
                  sp.simplify((x - 1) / (x + 1) - beta / R)))
    steps.append(("9c R^2 - beta^2    = T^2 > 0, so beta/R < 1 and artanh applies",
                  sp.simplify(Rv0 ** 2 - beta ** 2 - T ** 2)))

    for label, residual in steps:
        ok = residual == 0
        print(f"   {'OK ' if ok else 'FAIL'}  {label}")
        if not ok:
            print(f"          residual: {residual}")
            FAIL.append(f"check 3a: step failed -- {label}")

    # 7  the root is unique in d > 0 and the sign changes + -> -
    root = sp.solve(sp.Eq(2 * beta * (beta ** 2 - d ** 2 + T ** 2), 0), d)
    print(f"\n   7  roots in d > 0: {root}")
    if root != [sp.sqrt(T ** 2 + beta ** 2)]:
        FAIL.append("check 3a: the stationary point is not unique / not R")
    f = sp.lambdify((beta, T, d), 2 * beta * (beta ** 2 - d ** 2 + T ** 2), 'numpy')
    lo, hi = f(0.1, 30.0, 30.00017 * 0.5), f(0.1, 30.0, 30.00017 * 1.5)
    print(f"   7  sign of the bracket below and above R: {np.sign(lo):+.0f}, {np.sign(hi):+.0f}"
          f"   (must be +1 then -1, so R is a maximum)")
    if not (lo > 0 > hi):
        FAIL.append("check 3a: the sign change at R is not + -> -")


def check3_optimal_width():
    print("\n" + "=" * 74)
    print("3  Proposition 2: the optimum is exactly d = R = |rho - 1/2|,")
    print("   with value artanh(beta/R)")
    print("=" * 74)
    b_, d_, T_ = sp.symbols('beta d T', positive=True)
    N = (b_ + d_) ** 2 + T_ ** 2
    D = (b_ - d_) ** 2 + T_ ** 2
    crit = sp.factor(sp.numer(sp.together(sp.diff(N / D, d_))))
    sol = sp.solve(sp.Eq(crit, 0), d_)
    print(f"   d/dd numerator factors as   {crit}")
    print(f"   stationary points           {sol}")
    if sol != [sp.sqrt(T_ ** 2 + b_ ** 2)]:
        FAIL.append("check 3: the stationary point is not sqrt(beta^2 + T^2)")
    print(f"\n   {'T':>6} {'beta':>6} {'argmax_d':>13} {'R':>13} "
          f"{'max rate':>14} {'artanh(beta/R)':>15}")
    for T, beta in ((10.0, 0.1), (30.0, 0.1), (100.0, 0.1), (30.0, 0.01), (500.0, 0.01)):
        R = np.hypot(beta, T)
        ds = np.linspace(0.05, 5 * T, 800_000)
        f = rate(beta, T, ds)
        dbest, fmax = ds[np.argmax(f)], f.max()
        print(f"   {T:>6.0f} {beta:>6} {dbest:>13.5f} {R:>13.5f} "
              f"{fmax:>14.6e} {np.arctanh(beta/R):>15.6e}")
        if abs(dbest / R - 1) > 1e-4:
            FAIL.append(f"check 3: argmax != R at T={T}, beta={beta}")
        if abs(fmax / np.arctanh(beta / R) - 1) > 1e-6:
            FAIL.append(f"check 3: max rate != artanh(beta/R) at T={T}, beta={beta}")


def check3b_online_and_far(zerofile="2_million.txt"):
    """section 3.4: on-line terms are 1 - cos(n theta), and the far zeros contribute
    a positive n^2 d^2 term matching the density model."""
    import os
    print("\n" + "=" * 74)
    print("3b  section 3.4: on-line terms, and the far-zero contribution")
    print("=" * 74)
    worst = 0.0
    for d in (0.5, 30.0):
        for g in (50.0, 500.0, 5000.0):
            for n in (3, 100):
                rho = 0.5 + 1j * g
                q = (rho - (0.5 - d)) / (rho - (0.5 + d))
                theta = np.angle(q)
                worst = max(worst, abs((1 - q ** n).real - (1 - np.cos(n * theta))),
                            abs(abs(q) - 1))
                if abs(theta - (2 * np.arctan(g / d) - np.pi)) > 1e-12:
                    FAIL.append(f"check 3b: theta != 2 atan(g/d) - pi at d={d}, g={g}")
    print(f"   worst |Re(1-q^n) - (1-cos n theta)| and | |q|-1 | on the line: {worst:.2e}")
    if worst > 1e-12:
        FAIL.append("check 3b: the on-line identity fails")

    path = os.environ.get("ZETA_ZEROS", zerofile)
    if not os.path.exists(path):
        print(f"   (far-zero table skipped: set ZETA_ZEROS to a zero list; looked for {path!r})")
        return
    g = np.loadtxt(path)
    print(f"\n   far-zero sum against the model 2 n^2 d^2/gamma^2, {len(g)} ordinates")
    print(f"   {'d':>6} {'n':>5} {'T':>9} {'exact':>13} {'model':>13} {'ratio':>7}")
    for d, n in ((0.5, 20), (0.5, 50), (30.0, 50)):
        for T in (5e3, 2e4, 1e5):
            if T < 20 * n * d:
                continue
            sel = g[g > T]
            th = 2 * np.arctan(sel / d) - np.pi
            exact = 2 * np.sum(1 - np.cos(n * th))
            model = 2 * (2 * n * n * d * d) * (np.log(T / (2 * np.pi)) + 1) / (2 * np.pi * T)
            print(f"   {d:>6} {n:>5} {T:>9.0e} {exact:>13.4e} {model:>13.4e} {exact/model:>7.3f}")
            if not (0.85 < exact / model < 1.05):
                FAIL.append(f"check 3b: far-zero model off at d={d}, n={n}, T={T}")
    print("   (the shortfall at large T is the finite zero list, not the model)")


def check4_gain():
    print("\n" + "=" * 74)
    print("4  Corollary 3: the gain over the classical lens d = 1/2 is a factor T")
    print("=" * 74)
    print(f"   {'T':>7} {'beta':>7} {'rate at d=R':>14} {'rate at d=1/2':>15} {'gain':>10}")
    for T in (10.0, 30.0, 100.0, 500.0):
        for beta in (0.1, 0.01):
            R = np.hypot(beta, T)
            g = rate(beta, T, R) / rate(beta, T, 0.5)
            print(f"   {T:>7.0f} {beta:>7} {rate(beta,T,R):>14.3e} "
                  f"{rate(beta,T,0.5):>15.3e} {g:>10.3f}")
            if abs(g / T - 1) > 0.01:
                FAIL.append(f"check 4: gain is not T at T={T}, beta={beta}")


def check5_the_measured_case():
    print("\n" + "=" * 74)
    print("5  section 3.3: the 30.9x improvement is T = 30")
    print("=" * 74)
    r_classical = rate(0.1, 30.0, 0.5)
    r_matched = rate(0.1, 30.0, 30.0)
    pred = r_matched / r_classical
    measured = 41646 / 1348
    print(f"   planted quartet rho = 1/2 +- 0.1 +- 30i")
    print(f"   rate at a = 0      (d = 1/2)   {r_classical:.6e}")
    print(f"   rate at a = -29.5  (d = 30)    {r_matched:.6e}")
    print(f"   predicted threshold ratio      {pred:.2f}")
    print(f"   measured 41646 / 1348          {measured:.2f}")
    print(f"   residue                        {abs(pred-measured)/measured*100:.1f}%"
          f"   (finite-n background, ~ (n/2) log n)")
    if abs(pred - 30.0) > 0.1:
        FAIL.append("check 5: predicted ratio is not T = 30")


def check6_truncation_constants():
    print("\n" + "=" * 74)
    print("6  section 4: the truncation constants are (1 + 1/b), b = -a")
    print("=" * 74)
    a, b, d = sp.symbols('a b d')
    expr = ((d + sp.Rational(1, 2)) / (d - sp.Rational(1, 2))).subs(d, sp.Rational(1, 2) - a)
    simplified = sp.simplify(expr)
    in_b = sp.simplify(simplified.subs(a, -b))
    print(f"   (d+1/2)/(d-1/2) with d = 1/2 - a   ->   {simplified}")
    print(f"   with b = -a                        ->   {in_b}")
    if sp.simplify(in_b - (b + 1) / b) != 0:
        FAIL.append("check 6: the rate does not reduce to (1+1/b)")
    print(f"\n   {'a':>7} {'d':>6} {'(d+1/2)/(d-1/2)':>17} {'1 + 1/|a|':>11}")
    for av in (-0.5, -1.5, -3.5):
        dv = 0.5 - av
        lhs, rhs = (dv + 0.5) / (dv - 0.5), 1 + 1 / abs(av)
        print(f"   {av:>7} {dv:>6} {lhs:>17.4f} {rhs:>11.4f}")
        if abs(lhs - rhs) > 1e-12:
            FAIL.append(f"check 6: mismatch at a={av}")
    print("\n   this is the term (1 + 1/b)^n of Sekatskii, arXiv:1404.7276 eq. (7),")
    print("   whose source is the formal zero of xi at s = 1 and which is cancelled")
    print("   by the pole of zeta once the prime sum is regrouped before truncation.")


def check8_prime_side(M=(500, 2000, 8000)):
    """section 2.3: the Laguerre prime formula, checked rather than quoted."""
    print("\n" + "=" * 74)
    print("8  section 2.3: the prime-side Laguerre formula")
    print("=" * 74)
    mp.mp.dps = 30

    def lhs(n, a):
        f = lambda s: (s - a) ** (n - 1) * mp.log(mp.zeta(s))
        return mp.diff(f, mp.mpf(1) - a, n, h=mp.mpf('0.02')) / mp.factorial(n - 1)

    def rhs(n, a, cutoff, over_logm=False):
        tot = mp.mpf(0)
        primes = [q for q in range(2, cutoff + 1)
                  if all(q % r for r in range(2, int(q ** 0.5) + 1))]
        for p in primes:
            k, q = 1, p
            while q <= cutoff:
                coef = mp.log(p) / mp.log(q) if over_logm else mp.log(p)
                tot += coef * mp.mpf(q) ** (a - 1) * mp.laguerre(n - 1, 1, (1 - 2 * a) * mp.log(q))
                k += 1
                q = p ** k
        return tot

    print("   residue |LHS + RHS(M)| / |LHS|, which must fall with the prime cutoff M")
    print(f"   {'a':>6} {'n':>3} " + " ".join(f"{'M=' + str(m):>12}" for m in M))
    for a in (-0.5, -1.5):
        for n in (1, 2):
            L = lhs(n, a)
            row = [abs(L + rhs(n, a, m)) / abs(L) for m in M]
            print(f"   {a:>6} {n:>3} " + " ".join(f"{float(x):>12.3e}" for x in row))
            if not (row[0] > row[-1]):
                FAIL.append(f"check 8: residue not falling with M at a={a}, n={n}")
            if row[-1] > 1e-1:
                FAIL.append(f"check 8: residue too large at a={a}, n={n}")
    print("\n   the same residue for the variant carrying an extra 1/log m,")
    print("   which must NOT converge to zero as the cutoff grows")
    print(f"   {'a':>6} {'n':>3} " + " ".join(f"{'M=' + str(m):>12}" for m in M))
    for a in (-0.5, -1.5):
        for n in (1, 2):
            L = lhs(n, a)
            row = [abs(L + rhs(n, a, m, over_logm=True)) / abs(L) for m in M]
            print(f"   {a:>6} {n:>3} " + " ".join(f"{float(x):>12.3e}" for x in row))
            correct_last = abs(L + rhs(n, a, M[-1])) / abs(L)
            if row[-1] <= correct_last:
                FAIL.append(f"check 8: the 1/log m variant is not worse at a={a}, n={n}")
    print("\n   the correct coefficient converges with M; the variant stalls or drifts.")


def check7_weil_ceiling():
    """section 6: the Paley-Wiener ceiling, and that it has no interior maximum."""
    print("\n" + "=" * 74)
    print("7  section 6: the Weil-side ceiling is monotone in L -- no optimum")
    print("=" * 74)

    def phihat(r, L, j):
        a = j * np.pi / L
        sg = (-1.0) ** j
        return (np.sqrt(2.0 / L) * a * (1 - sg * np.exp(-1j * r * L))
                / (a ** 2 - r ** 2) * np.exp(1j * r * L / 2))

    def l1(L, j):
        """||phi_j||_1 in closed form.  With u = j pi y / L the integral is
        (L/(j pi)) * int_0^{j pi} |sin u| du = (L/(j pi)) * 2j = 2L/pi,
        so the j cancels and the answer does not depend on the mode."""
        return np.sqrt(2.0 / L) * 2 * L / np.pi

    print(f"   {'L':>6} {'T':>5} {'beta':>6} {'ratio to ceiling':>18} {'ceiling e^(bL/2)':>18}")
    prev = {}
    for L in (4.5, 9.0, 18.0, 36.0):
        for T, beta in ((30.0, 0.1), (30.0, 0.5)):
            j = max(1, int(round(T * L / np.pi)))
            j = j if abs(j * np.pi / L - T) > 1e-9 else j + 1
            ratio = abs(phihat(T - 1j * beta, L, j)) / (l1(L, j) * np.exp(beta * L / 2))
            ceil = np.exp(beta * L / 2)
            print(f"   {L:>6.1f} {T:>5.0f} {beta:>6} {ratio:>18.6f} {ceil:>18.4f}")
            if ratio > 1.0:
                FAIL.append(f"check 7: Paley-Wiener ceiling violated at L={L}, beta={beta}")
            if beta in prev and ceil <= prev[beta]:
                FAIL.append(f"check 7: ceiling not increasing in L at beta={beta}")
            prev[beta] = ceil
    print("   every ratio <= 1, and the ceiling beta*L/2 is strictly increasing:")
    print("   the support length has no interior optimum, unlike the Li lens width.")


def main() -> int:
    check0_map_geometry()
    check1_normalisation()
    check2_step_caveat()
    check3a_proof_steps()
    check3_optimal_width()
    check3b_online_and_far()
    check4_gain()
    check5_the_measured_case()
    check6_truncation_constants()
    check7_weil_ceiling()
    check8_prime_side()
    print("\n" + "=" * 74)
    if FAIL:
        print("FAILURES:")
        for f in FAIL:
            print("  -", f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
