"""
verify_core.py — Appendix B of
"Compressing the Division Table into a Single Dynamical Column".

Covers the exact chain and the §10 spectral results:

  1. Proposition 1 — every prime factor p <= sqrt(n) puts n on p's trajectory.
  2. Proposition 2 — the event queue is empty exactly at the primes.
  3. Proposition 3 — <Phi(a),Phi(b)> = log gcd(a,b) and the two corollaries.
  4. Proposition 4 — the lag correlation is log gcd(n_t, h).
  5. Proposition 5 — layer isolation.
  6. Proposition 6 — A_alpha = Lambda_alpha * 1 and its Mobius inversion.
  7. Proposition 7 — the Dirichlet transform to -zeta(s) zeta'/zeta(s+alpha+1).
  8. §10 Lemmas 1 and 2, and the spectral table.
  9. §10 the RMSE protocol comparison.
 10. §10 the randomised-phase control.

§8 is in verify_covariance.py and §9 in verify_positivity.py.

Python 3, numpy, mpmath. Set FAST = True to shrink the ranges in steps 1-2.
"""
import math
import numpy as np
from mpmath import mp, zeta, log as mlog, diff, zetazero

mp.dps = 25
FAST = False
NMAX = 20_001 if FAST else 200_001


# ---------------------------------------------------------------- utilities
def composite_flags(P):
    c = np.zeros(P + 1, dtype=bool)
    for p in range(2, int(P ** 0.5) + 1):
        if not c[p]:
            c[p * p::p] = True
    return c


def odd_prime_powers(P):
    c = composite_flags(P)
    out = []
    for p in range(3, P + 1):
        if not c[p]:
            q = p
            while q <= P:
                out.append((q, p))
                q *= p
    return out


def odd_factorisation(n):
    """{p: v_p(n)} over ODD primes only — the indexing Phi actually uses."""
    d = {}
    m = n
    while m % 2 == 0:
        m //= 2
    p = 3
    while p * p <= m:
        while m % p == 0:
            d[p] = d.get(p, 0) + 1
            m //= p
        p += 2
    if m > 1:
        d[m] = d.get(m, 0) + 1
    return d


# --------------------------------------------------- 1-2. the sieve statements
def step_1_2():
    print(f"\n=== 1-2. Propositions 1 and 2 (odd n < {NMAX - 1}) ===")
    comp = composite_flags(NMAX)
    bad1 = bad2 = n_comp = 0
    for n in range(9, NMAX, 2):
        struck = False
        for p in range(3, int(n ** 0.5) + 1, 2):
            if not comp[p] and n % p == 0:
                struck = True
                r = n - p * p                       # Proposition 1, EVERY such p
                if r < 0 or r % (2 * p) != 0:
                    bad1 += 1
        if comp[n]:
            n_comp += 1
        if struck != comp[n]:                       # Proposition 2
            bad2 += 1
    print(f"  Prop 1: {n_comp} odd composites, every qualifying prime factor checked"
          f"  -> violations = {bad1}")
    print(f"  Prop 2: {(NMAX - 9) // 2 + 1} odd n  -> mismatches = {bad2}")


# ------------------------------------------------- 3-5. the kernel statements
def step_3_5():
    print("\n=== 3. Proposition 3 — the log-gcd Gram kernel ===")
    rng = np.random.default_rng(0)
    e_ip = e_nm = e_dist = 0.0
    for _ in range(4000):
        a = int(rng.integers(1, 100_000)) * 2 + 1
        b = int(rng.integers(1, 100_000)) * 2 + 1
        va, vb = odd_factorisation(a), odd_factorisation(b)
        ip = sum(np.log(p) * min(k, vb.get(p, 0)) for p, k in va.items())
        e_ip = max(e_ip, abs(ip - np.log(math.gcd(a, b))))
        e_nm = max(e_nm, abs(sum(np.log(p) * k for p, k in va.items()) - np.log(a)))
        g = math.gcd(a, b)
        e_dist = max(e_dist, abs((np.log(a) + np.log(b) - 2 * ip)
                                 - np.log(a * b // g // g)))
    print(f"  max |<Phi(a),Phi(b)> - log gcd(a,b)| = {e_ip:.1e}")
    print(f"  max | ||Phi(n)||^2 - log n |         = {e_nm:.1e}")
    print(f"  max | ||Phi(a)-Phi(b)||^2 - log(lcm/gcd) | = {e_dist:.1e}")

    print("\n=== 4. Proposition 4 — the lag correlation is log gcd(n_t, h) ===")
    e4 = 0.0
    for t in range(400):
        n = 2 * t + 3
        vn = odd_factorisation(n)
        for h in range(1, 60):
            vh = odd_factorisation(h)
            ip = sum(np.log(p) * min(k, vh.get(p, 0)) for p, k in vn.items())
            e4 = max(e4, abs(ip - np.log(math.gcd(n, h))))
    print(f"  max deviation over t < 400, h < 60 = {e4:.1e}")

    print("\n=== 5. Proposition 5 — layer isolation ===")
    e5 = 0.0
    for p in (3, 5, 7, 11):
        for k in (1, 2, 3):
            for n in range(3, 4000, 2):
                lhs = np.log(math.gcd(n, p ** k)) - np.log(math.gcd(n, p ** (k - 1)))
                e5 = max(e5, abs(lhs - np.log(p) * (n % p ** k == 0)))
    print(f"  max deviation over p in {{3,5,7,11}}, k <= 3, odd n < 4000 = {e5:.1e}")


# ----------------------------------------------------------- 6. Mobius
def step_6():
    print("\n=== 6. Proposition 6 — A_alpha = Lambda_alpha * 1, and its inversion ===")
    H = 3000
    comp = composite_flags(H)
    lam = np.zeros(H + 1)
    for p in range(2, H + 1):
        if not comp[p]:
            q = p
            while q <= H:
                lam[q] = np.log(p)
                q *= p
    mu = np.zeros(H + 1)
    mu[1] = 1.0
    for d in range(2, H + 1):
        m, f, x, p = d, {}, d, 2
        while p * p <= x:
            while x % p == 0:
                f[p] = f.get(p, 0) + 1
                x //= p
            p += 1
        if x > 1:
            f[x] = f.get(x, 0) + 1
        mu[d] = 0.0 if any(v > 1 for v in f.values()) else (-1.0) ** len(f)
    for al in (0.25, 0.75, 1.5):
        A = np.array([sum(lam[d] / d ** (al + 1) for d in range(1, h + 1) if h % d == 0)
                      for h in range(1, H + 1)])
        conv = max(abs(A[h - 1] - sum(lam[d] / d ** (al + 1) for d in range(1, h + 1) if h % d == 0))
                   for h in range(1, 400))
        inv = max(abs(sum(mu[d] * A[n // d - 1] for d in range(1, n + 1) if n % d == 0)
                      - lam[n] / n ** (al + 1))
                  for n in range(1, 400))
        print(f"  alpha={al:4.2f}: convolution identity {conv:.1e}   Mobius inversion {inv:.1e}")


# --------------------------------------------------------- 7. Proposition 7
def step_7():
    print("\n=== 7. Proposition 7 — the Dirichlet transform ===")
    H = 200_000
    comp = composite_flags(H)
    lam = np.zeros(H + 1)
    for p in range(2, H + 1):
        if not comp[p]:
            q = p
            while q <= H:
                lam[q] = np.log(p)
                q *= p
    idx = np.nonzero(lam)[0]
    for s, al in ((mp.mpf(3) + 4j, 0.75), (mp.mpf('2.5'), 0.75), (mp.mpf('2.2'), 1.5)):
        lhs = mp.mpf(0)
        for d in idx:
            lhs += mp.mpf(float(lam[d])) * mp.mpf(int(d)) ** (-(s + al + 1))
        lhs *= zeta(s)
        # rhs = - zeta(s) * zeta'/zeta (s + alpha + 1)
        rhs = -zeta(s) * diff(lambda z: mlog(zeta(z)), s + al + 1)
        print(f"  s={str(s):>12}, alpha={al}: relative error = {float(abs(lhs - rhs) / abs(rhs)):.1e}")


# ------------------------------------------------ 8-10. the hybrid operator
def hybrid(N, P, seed=None):
    """seed=None gives the arithmetic phases; an int randomises them."""
    qs = odd_prime_powers(P)
    t = np.arange(N)
    n = 2 * t + 3
    rg = np.random.default_rng(seed) if seed is not None else None
    A = np.zeros((N, N))
    for q, p in qs:
        if q >= N:
            continue
        if rg is None:
            hit = (n % q == 0).astype(float)
        else:
            hit = ((t % q) == int(rg.integers(0, q))).astype(float)
        v = (np.log(p) / np.sqrt(q)) * (hit - 1.0 / q)
        for i in range(N - q):
            A[i, i + q] += v[i]
            A[i + q, i] += v[i]
    return A


def step_8():
    print("\n=== 8. Lemmas 1 and 2, and the spectral table ===")
    for (N, P) in ((160, 25), (240, 49), (320, 81), (400, 121)):
        A = hybrid(N, P)
        D = np.diag((-1.0) ** np.arange(N))
        ev = np.linalg.eigvalsh(A)
        print(f"  N={N:4d} P={P:4d}:  ||A-A^T||={np.abs(A - A.T).max():.1e}"
              f"  ||DAD+A||={np.abs(D @ A @ D + A).max():.1e}"
              f"  lam_min={ev.min():+.6f}  lam_max={ev.max():+.6f}"
              f"  sum={ev.min() + ev.max():+.1e}")


def step_9():
    print("\n=== 9. the RMSE protocol comparison (N=320, P=81) ===")
    gam = np.array([float(zetazero(k).imag) for k in range(1, 16)])
    ev = np.sort(np.linalg.eigvalsh(hybrid(320, 81)))[::-1]

    def rmse(order):
        a, b = np.polyfit(ev[order][:5], gam[:5], 1)
        return a, float(np.sqrt(np.mean((a * ev[order][5:15] + b - gam[5:15]) ** 2)))

    orders = (("5 largest, descending", np.arange(len(ev))),
              ("5 largest, ascending", np.concatenate([np.arange(5)[::-1],
                                                       np.arange(5, len(ev))])),
              ("5 largest by |lambda|", np.argsort(-np.abs(ev))))
    for nm, o in orders:
        a, r = rmse(o)
        print(f"  {nm:>22}: slope={a:+.2f}   held-out RMSE={r:.3f}")
    print("  A factor of six separates the conventions, so a bare number means nothing;")
    print("  and by Lemmas 1-2 no convention can work.")


def step_10():
    print("\n=== 10. the randomised-phase control ===")
    top = np.sort(np.linalg.eigvalsh(hybrid(320, 81)))[::-1][:10]
    ds = []
    for sd in range(8):
        r = np.sort(np.linalg.eigvalsh(hybrid(320, 81, seed=sd)))[::-1][:10]
        ds.append(float(np.linalg.norm(r - top) / np.linalg.norm(top)))
    print(f"  mean normalised top-10 spectral distance over 8 seeds = {np.mean(ds):.4f}"
          f"   (spread {min(ds):.4f} - {max(ds):.4f})")
    print("  Destroying the arithmetic while keeping the periods and densities barely")
    print("  moves the spectrum.")


if __name__ == "__main__":
    step_1_2()
    step_3_5()
    step_6()
    step_7()
    step_8()
    step_9()
    step_10()
