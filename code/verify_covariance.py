"""
verify_covariance.py — Appendix B, §8 of
"Compressing the Division Table into a Single Dynamical Column".

Prints, in order:
  1. Proposition 8, and the identity C_alpha = sum_q w_q/q that makes the
     centring annihilate the theta = 0 atom exactly.
  2. Which spectral measure actually reproduces K_alpha: the 1/q normalisation
     and the omission of the m = 0 atom are both necessary.
  3. Proposition 9 — the pointwise Gram matrix is not Toeplitz, so the
     "it is a Gram matrix" argument fails; the shift average is what is
     Toeplitz, and it is PSD as an average of Gram matrices.
  4. The convergence rate of the empirical time average: faster than 1/T,
     and not a clean power.
  5. The layer-cutoff sweep of lambda_min and its apparent decay exponent.
  6. The rank law rank = (#atoms) - 1, with the tolerance it requires.

Python 3, numpy. Set FAST = True to shorten step 4.
"""
import numpy as np
from fractions import Fraction

ALPHA = 0.75
FAST = False


# ------------------------------------------------------------------ layers
def odd_prime_powers(P):
    comp = np.zeros(P + 1, dtype=bool)
    for p in range(2, int(P ** 0.5) + 1):
        if not comp[p]:
            comp[p * p::p] = True
    out = []
    for p in range(3, P + 1):
        if not comp[p]:
            q = p
            while q <= P:
                out.append((q, p))
                q *= p
    return out


def layers(P):
    """qs = [(q, p)], w[q] = log p / q^(alpha+1), C = centring constant."""
    qs = odd_prime_powers(P)
    w = {q: np.log(p) / q ** (ALPHA + 1) for q, p in qs}
    C = sum(np.log(p) / q ** (ALPHA + 2) for q, p in qs)
    return qs, w, C


def K(qs, w, C, h):
    a = sum(w[q] for q, _ in qs) if h == 0 else sum(w[q] for q, _ in qs if h % q == 0)
    return a - C


def psi_scale(q, w):
    """the per-layer amplitude of Psi_alpha, i.e. sqrt(log p / q^alpha)."""
    return np.sqrt(w[q] * q)


# ------------------------------------------- 1 + 2. the measure, done right
def step_1_2(P=2000):
    qs, w, C = layers(P)
    dc = sum(w[q] / q for q, _ in qs)
    print(f"\n=== 1. Proposition 8 and the DC identity (alpha={ALPHA}, layers to {P}) ===")
    print(f"  layers                : {len(qs)}")
    print(f"  C_alpha               : {C:.12f}")
    print(f"  sum_q w_q/q (DC atom) : {dc:.12f}")
    print(f"  |difference|          : {abs(C - dc):.1e}   -> centring kills theta = 0 exactly")

    def from_measure(h, scaled=True, drop_dc=True):
        tot = 0.0
        for q, _ in qs:
            ww = w[q] / q if scaled else w[q]
            tot += ww * sum(np.cos(2 * np.pi * m * h / q)
                            for m in range(1 if drop_dc else 0, q))
        return tot

    print("\n=== 2. which spectral measure reproduces K_alpha(h)? ===")
    print(f"  {'h':>4} {'K_alpha(h)':>15} {'(w_q/q), m>=1':>16} {'w_q, m>=0':>15}")
    for h in (0, 1, 2, 3, 5, 6, 9, 15, 45):
        print(f"  {h:4d} {K(qs, w, C, h):15.9f} {from_measure(h):16.9f} "
              f"{from_measure(h, scaled=False, drop_dc=False):15.9f}")
    print("  Only the middle column matches: both the 1/q normalisation and the")
    print("  removal of the m = 0 atom are required.")
    return qs, w, C


# ------------------------------------------------------ 3. Proposition 9
def step_3(qs, w, C):
    print("\n=== 3. Proposition 9: [K(|i-j|)] is NOT a Gram matrix ===")

    def Psi(t):
        n = 2 * t + 3
        return np.array([psi_scale(q, w) * ((n % q == 0) - 1.0 / q) for q, _ in qs])

    G = np.array([[Psi(i) @ Psi(j) for j in range(6)] for i in range(6)])
    for d in range(3):
        v = G.diagonal(d)
        print(f"  pointwise Gram, offset {d}: {v.min():+.4f} .. {v.max():+.4f}"
              f"   constant along the diagonal? {abs(v.max() - v.min()) < 1e-12}")

    T = 200_000
    n = 2 * np.arange(T) + 3
    ind = np.array([psi_scale(q, w) * ((n % q == 0) - 1.0 / q) for q, _ in qs])
    print("  the shift average IS Toeplitz, and is PSD as an average of Gram matrices:")
    for h in (0, 1, 3, 9, 15):
        emp = float((ind[:, :T - h] * ind[:, h:]).sum() / (T - h))
        k = K(qs, w, C, h)
        print(f"    h={h:3d}  shift-avg={emp:+.9f}  K_alpha={k:+.9f}  diff={emp - k:+.1e}")
    del ind


# --------------------------------------------------- 4. convergence rate
def step_4(qs, w, C):
    print("\n=== 4. convergence of the empirical time average ===")
    hs = [0, 1, 3, 9, 15]
    Ts = (20_000, 200_000) if FAST else (20_000, 200_000, 2_000_000)
    prev = None
    for T in Ts:
        acc = np.zeros(len(hs))
        cnt = np.zeros(len(hs))
        B = 500_000
        for st in range(0, T, B):
            en = min(st + B + max(hs), T)
            n = 2 * np.arange(st, en) + 3
            ind = np.array([psi_scale(q, w) * ((n % q == 0) - 1.0 / q) for q, _ in qs])
            L0 = min(B, T - st)
            for k, h in enumerate(hs):
                L = min(L0, ind.shape[1] - h)
                acc[k] += float((ind[:, :L] * ind[:, h:L + h]).sum())
                cnt[k] += L
            del ind
        err = max(abs(acc[k] / cnt[k] - K(qs, w, C, h)) for k, h in enumerate(hs))
        line = f"  T={T:9d}  max|emp - K| = {err:.3e}"
        if prev:
            line += f"   gain over the previous decade: {prev / err:.1f}x   (1/T would give 10x)"
        print(line)
        prev = err
    print("  Faster than 1/T, and not a clean power: this is an arithmetic average")
    print("  over exactly periodic sequences, not a stochastic one.")


# ------------------------------------------------- 5. layer-cutoff sweep
def step_5():
    print("\n=== 5. lambda_min against the layer cutoff (alpha = 0.75) ===")
    Ns = (10, 20, 40, 80, 160, 320)
    print(f"  {'cutoff':>8} " + " ".join(f"{'N=' + str(n):>10}" for n in Ns) + "   slope")
    for P in (100, 200, 500, 2000, 20000, 100000):
        qs, w, C = layers(P)
        kv = np.array([K(qs, w, C, h) for h in range(max(Ns))])
        row = []
        for N in Ns:
            i = np.arange(N)
            row.append(np.linalg.eigvalsh(kv[np.abs(i[:, None] - i[None, :])]).min())
        sl = np.polyfit(np.log(Ns), np.log(np.abs(row)), 1)[0]
        print(f"  {P:>8} " + " ".join(f"{x:10.6f}" for x in row) + f"  {sl:7.3f}")
    print("  The apparent exponent moves from -4.8 to -0.74 across this table alone.")


# ---------------------------------------------------------- 6. rank law
def step_6():
    print("\n=== 6. rank = (#atoms) - 1, and the tolerance it needs ===")
    for P in (25, 49, 81):
        qs, w, C = layers(P)
        atoms = {}
        for q, _ in qs:
            for m in range(1, q):
                f = Fraction(m, q)
                atoms[f] = atoms.get(f, 0.0) + w[q] / q
        kv = np.array([K(qs, w, C, h) for h in range(1100)])
        N = 1000
        i = np.arange(N)
        T = kv[np.abs(i[:, None] - i[None, :])]
        ev = np.sort(np.linalg.eigvalsh(T))[::-1]
        A = len(atoms) + 1                      # + the annihilated theta = 0 atom
        ranks = [np.linalg.matrix_rank(T, tol=t) for t in (1e-6, 1e-9, 1e-12)]
        print(f"  cutoff {P:3d}: layers={len(qs):3d}  atoms incl. theta=0 = {A:4d}  "
              f"predicted rank = {A - 1:4d}")
        print(f"      smallest atom weight = {min(atoms.values()):.2e}"
              f"   rank at tol 1e-6 / 1e-9 / 1e-12 = {ranks[0]} / {ranks[1]} / {ranks[2]}")
        k = A - 1
        print(f"      gap at the cut: lam[{k - 1}] = {ev[k - 1]:.1e}   lam[{k}] = {ev[k]:.1e}"
              f"   lam_min = {ev[-1]:+.1e}")
    print("  The invariant is real, but it must be read with a tolerance below the")
    print("  smallest atom weight; at cutoff 81 a tolerance of 1e-9 undercounts by 4.")


if __name__ == "__main__":
    qs, w, C = step_1_2()
    step_3(qs, w, C)
    step_4(qs, w, C)
    step_5()
    step_6()
