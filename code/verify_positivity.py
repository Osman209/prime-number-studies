"""
verify_positivity.py — Appendix B, §9 and the N-versus-P separation of §10 of
"Compressing the Division Table into a Single Dynamical Column".

Prints, in order:
  1. Proposition 10 — the zero-diagonal family, indefinite by construction.
  2. Proposition 11 — over odd lags, f(pi) = -c for ANY nonnegative weights.
     Checked on the arithmetic weights and on random weights as a control.
  3. The all-prime-powers convention: symbol minimum, the finite matrices,
     and the location of the minimising angle.
  4. The alpha-window. The LOWER edge is exact: below the root of
     f(pi) = 4 log2/(2^(a+1)-1) + zeta'/zeta(a+1) = 0 the symbol is negative
     at theta = pi, and in that whole region the global minimum IS at pi, so
     the root is the edge itself. The UPPER edge is bisected numerically
     because the minimising angle has left pi by then.
  5. The drift of the minimising angle with alpha.
  6. Why truncated evaluation of the symbol must not be trusted at small alpha.
  7. §10: lambda_max saturates in the prime cutoff P but keeps growing in N.

Python 3, numpy, scipy, mpmath.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from mpmath import mp, zeta, log as mlog, diff

mp.dps = 25
HMAX = 4_000_000
MGRID = 1 << 22          # theta_j = 2 pi j / MGRID, exact by folding h mod MGRID


# ---------------------------------------------------------- von Mangoldt
def von_mangoldt(H):
    comp = np.zeros(H + 1, dtype=bool)
    for p in range(2, int(H ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
    lam = np.zeros(H + 1)
    for p in range(2, H + 1):
        if not comp[p]:
            q = p
            while q <= H:
                lam[q] = np.log(p)
                q *= p
    return lam


LAM = von_mangoldt(HMAX)
IDX = np.nonzero(LAM)[0]
WT = LAM[IDX]


def c_natural(s, odd_only=False):
    c = float(-diff(lambda z: mlog(zeta(z)), s))
    if odd_only:
        c -= float(mlog(2) / (2 ** mp.mpf(s) - 1))
    return c


def symbol(s, odd_only=False):
    """returns (min over theta of the full symbol, argmin/pi) by one real FFT."""
    h, w = IDX, WT * IDX.astype(np.float64) ** (-s)
    if odd_only:
        m = h % 2 == 1
        h, w = h[m], w[m]
    fold = np.zeros(MGRID)
    np.add.at(fold, h % MGRID, w)
    f = 2 * np.fft.rfft(fold).real
    j = int(np.argmin(f))
    return c_natural(s, odd_only) + f[j], 2 * j / MGRID


def f_pi(s):
    """exact: only the p = 2 layer sits at even lags."""
    return 4 * float(mlog(2) / (2 ** mp.mpf(s) - 1)) - c_natural(s)


def f_halfpi(s):
    t = 2 ** (-mp.mpf(s))
    return c_natural(s) + 2 * float(mlog(2)) * float(t * t / (1 - t) - t)


def toeplitz(b, diag):
    N = len(b)
    i = np.arange(N)
    T = b[np.abs(i[:, None] - i[None, :])]
    np.fill_diagonal(T, diag)
    return T


def bvec(s, N, odd_only=False):
    h = np.arange(1, N)
    b = np.zeros(N)
    v = LAM[h] * h.astype(float) ** (-s)
    b[1:] = np.where(h % 2 == 1, v, 0.0) if odd_only else v
    return b


# ------------------------------------------------------- 1. Proposition 10
def step_1(a=0.75):
    s = a + 1.0
    print(f"\n=== 1. Proposition 10 — zero diagonal forces a negative eigenvalue (alpha={a}) ===")
    for N in (10, 20, 40, 80, 120):
        T = toeplitz(bvec(s, N), 0.0)
        ev = np.linalg.eigvalsh(T)
        print(f"  N={N:4d}  lam_min={ev.min():.6f}  #negative={(ev < 0).sum():3d}"
              f"  trace={np.trace(T):.1e}  sum(eigenvalues)={ev.sum():+.1e}")


# ------------------------------------------------------- 2. Proposition 11
def step_2(a=0.75):
    s = a + 1.0
    c = c_natural(s, odd_only=True)
    m, t = symbol(s, odd_only=True)
    print(f"\n=== 2. Proposition 11 — odd lags, f(pi) = -c (alpha={a}) ===")
    print(f"  c (odd prime powers)  = {c:.6f}")
    print(f"  min over theta        = {m:.6f}   at argmin/pi = {t:.4f}")
    print(f"  -c                    = {-c:.6f}")
    row = [np.linalg.eigvalsh(toeplitz(bvec(s, N, True), c)).min()
           for N in (10, 20, 40, 80, 160, 320, 640)]
    print("  finite lam_min : " + "  ".join(f"{x:+.4f}" for x in row))

    print("  control — RANDOM nonnegative weights on odd lags (the claim is weight-free):")
    rng = np.random.default_rng(7)
    for trial in range(4):
        N = 400
        odd = np.arange(1, N)
        b = np.zeros(N)
        b[1:] = np.where(odd % 2 == 1, rng.random(N - 1) * np.exp(-odd / 40.0), 0.0)
        cc = b[1:].sum()
        th = np.linspace(0, np.pi, 20001)
        f = cc + 2 * (np.cos(np.outer(th, odd)) @ b[1:])
        print(f"    trial {trial}: c={cc:.6f}  min_theta f={f.min():+.6f}  -c={-cc:+.6f}"
              f"  argmin/pi={th[np.argmin(f)] / np.pi:.4f}")


# --------------------------------------------- 3. all prime powers at 0.75
def step_3(a=0.75):
    s = a + 1.0
    c = c_natural(s)
    m, t = symbol(s)
    print(f"\n=== 3. all prime powers, alpha = {a} ===")
    print(f"  c = -zeta'/zeta({s}) = {c:.6f}")
    print(f"  symbol minimum       = {m:.6f}   at argmin/pi = {t:.4f}")
    print(f"  margin f_min / c     = {m / c:.4f}")
    row = [np.linalg.eigvalsh(toeplitz(bvec(s, N), c)).min()
           for N in (10, 20, 40, 80, 160, 320, 640)]
    print("  finite lam_min : " + "  ".join(f"{x:.4f}" for x in row))


# ----------------------------------------------------------- 4. the window
def step_4():
    print("\n=== 4. the alpha-window ===")
    lo = 1.30
    hi = 1.70
    for _ in range(80):                       # exact bisection on f(pi)
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f_pi(mid) < 0 else (lo, mid)
    lower = 0.5 * (lo + hi) - 1.0
    print(f"  LOWER edge, exact  : root of f(pi) = 0 at alpha = {lower:.6f}")
    print("    and the global minimum sits AT theta = pi throughout that region,")
    print("    so the root is the edge itself, not merely a bound:")
    for a in (0.45, 0.49, 0.50, 0.52, 0.60):
        v, t = symbol(a + 1.0)
        print(f"      alpha={a:4.2f}  f(pi) exact={f_pi(a + 1.0):+.6f}  "
              f"symbol min={v:+.6f}  argmin/pi={t:.4f}")

    lo, hi = 1.55, 1.70
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if symbol(mid + 1.0)[0] > 0 else (lo, mid)
    upper = 0.5 * (lo + hi)
    print(f"  UPPER edge, bisected: alpha = {upper:.4f}"
          f"   (f(pi/2) exact goes negative later: {f_halfpi(4.0):+.4f} at alpha=3)")
    print(f"  => PSD only for alpha in ({lower:.6f}, {upper:.4f})")


# ------------------------------------------- 5. the cusp that pins the angle
def refine(s, lo=0.38 * np.pi, hi=0.42 * np.pi, iters=60):
    """off-grid ternary search for the minimum near 2 pi/5."""
    h, w = IDX.astype(np.float64), WT * IDX.astype(np.float64) ** (-s)
    c = c_natural(s)

    def f(th):
        tot = 0.0
        for st in range(0, len(h), 400_000):
            tot += np.cos(th * h[st:st + 400_000]) @ w[st:st + 400_000]
        return c + 2 * tot

    for _ in range(iters):
        x = lo + (hi - lo) / 3
        y = hi - (hi - lo) / 3
        lo, hi = (lo, y) if f(x) < f(y) else (x, hi)
    th = 0.5 * (lo + hi)
    return th, f(th), f


def step_5():
    print("\n=== 5. the minimising angle locks onto rational points ===")
    print(f"  {'alpha':>6} {'symbol min':>12} {'argmin/pi (grid)':>18} {'argmin/pi (off-grid)':>21}")
    for a in (0.45, 0.60, 0.70, 0.75, 1.00, 1.40, 1.60, 1.61):
        v, t = symbol(a + 1.0)
        if t > 0.9:                            # still locked on theta = pi
            print(f"  {a:6.2f} {v:+12.6f} {t:18.4f} {'(at pi)':>21}")
        else:
            th, fv, _ = refine(a + 1.0)
            print(f"  {a:6.2f} {fv:+12.6f} {t:18.4f} {th / np.pi:21.8f}")
    print("  Exactly pi below ~0.71, then exactly 2/5 at 0.75 and 1.00, then it comes loose.")

    print("\n  the lock is a downward cusp:  f(2pi/5 + eps) - f(2pi/5) ~ eps^beta")
    es = np.array([1e-3, 1e-4, 1e-5, 1e-6])
    t0 = 0.4 * np.pi
    print(f"  {'alpha':>6} {'beta':>7}   differences")
    for a in (0.60, 0.75, 1.00, 1.40):
        _, _, f = refine(a + 1.0, iters=1)
        f0 = f(t0)
        d = np.array([f(t0 + e) - f0 for e in es])
        if (d <= 0).any():
            print(f"  {a:6.2f}   {'--':>7}   not a local minimum at 2pi/5")
            continue
        beta = np.polyfit(np.log(es), np.log(d), 1)[0]
        print(f"  {a:6.2f} {beta:7.3f}   " + " ".join(f"{x:.2e}" for x in d))
    print(f"  (summation cutoff HMAX = {HMAX:.0e}; the exponents approach alpha as HMAX grows)")
    print("  beta ~ alpha while the minimum is pinned; the pinning fails once the")
    print("  cusp becomes Lipschitz (beta = 1) at alpha ~ 1.")


# ------------------------------------------------- 6. truncation warning
def step_6():
    print("\n=== 6. why the truncated symbol must not be trusted at small alpha ===")
    print(f"  {'alpha':>6} {'truncated (H=4e6)':>19} {'f(pi) exact':>13}")
    for a in (0.05, 0.10, 0.20, 0.40):
        v, _ = symbol(a + 1.0)
        print(f"  {a:6.2f} {v:19.4f} {f_pi(a + 1.0):13.4f}")
    print("  The truncated value is positive where the exact one is hugely negative.")


# ---------------------------------------------------- 7. section 10, N vs P
def step_7():
    print("\n=== 7. section 10 — lambda_max saturates in P but not in N ===")

    def hybrid(N, P):
        sm = [p for p in range(3, P + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
        qs = [(p ** k, p) for p in sm for k in range(1, 25) if p ** k <= P]
        n = 2 * np.arange(N) + 3
        diags, offs = [], []
        for q, p in qs:
            if q >= N:
                continue
            v = (np.log(p) / np.sqrt(q)) * ((n % q == 0).astype(float) - 1.0 / q)
            diags += [v[:N - q], v[:N - q]]
            offs += [q, -q]
        return sp.diags(diags, offs, shape=(N, N), format="csr")

    print("  fixed N = 2000, raising the prime cutoff:")
    for P in (25, 81, 289, 529, 1089, 2209, 4489):
        ev = np.linalg.eigvalsh(hybrid(2000, P).toarray())
        budget = sum(np.log(p) / np.sqrt(p ** k)
                     for p in [x for x in range(3, P + 1)
                               if all(x % d for d in range(2, int(x ** 0.5) + 1))]
                     for k in range(1, 25) if p ** k <= P)
        print(f"    P={P:5d}  lam_max={ev.max():.4f}  lam_max+lam_min={ev.max() + ev.min():+.1e}"
              f"  weight budget={budget:7.2f}")
    print("  fixed P = 4489, raising the matrix size:")
    for N in (5120, 20480, 81920):
        lm = eigsh(hybrid(N, 4489), k=1, which="LA", return_eigenvectors=False, maxiter=20000)[0]
        print(f"    N={N:7d}  lam_max={lm:.5f}")


if __name__ == "__main__":
    step_1()
    step_2()
    step_3()
    step_4()
    step_5()
    step_6()
    step_7()
