"""
Follow-up to prime_edge_jump.py.

SUPERSEDED MECHANISM, KEPT AS A RECORD.  An earlier version of this file guessed
that as the shift t -> L the overlap region [0, L-t] shrinks to the point x = 0,
giving a rank-one limit with the single vector v_j ~ phi_j'(0) ~ j.  That guess is
WRONG and the numerics below show why: the measured alignment against j plateaus
at 1/sqrt(2) = 0.7071, not 1.

THE CORRECT ACCOUNT (A. Groskin, correspondence; see the public derivation linked
in papers/).  With delta = L - t,

    R_jk(L - delta) = (pi^2 delta^3 / 6 L^3) j k ((-1)^(j+1) + (-1)^(k+1)) + O(delta^5)

opposite parities decouple exactly, and the leading matrix is one rank-one ODD
block minus one rank-one EVEN block.  For even m the combined ratio tends to
(m-1)/(m+2).  The second endpoint vector is (-1)^(j+1) j, and j itself is a
45-degree mix of the odd-only and even-only vectors -- which is exactly why the
alignment against j saturates at 1/sqrt(2).

NUMERICAL REGIME.  All quantities here are cancellation-limited in float64 once
delta falls below about 1e-3; ratios printed below that are not reliable and the
edge limit should be taken in extended precision instead.
"""
import numpy as np, math
from sympy import primerange

def _I(w, ph, S):
    w = np.asarray(w, float); out = np.empty_like(w)
    small = np.abs(w) < 1e-13
    out[small] = S * np.cos(ph[small] if np.ndim(ph) else ph)
    ws = w[~small]; phs = ph[~small] if np.ndim(ph) else ph
    out[~small] = (np.sin(ws * S + phs) - np.sin(phs)) / ws
    return out

def overlap(t, L, m):
    if t >= L: return np.zeros((m, m))
    a = np.arange(1, m + 1) * np.pi / L
    A = a[:, None] * np.ones(m)[None, :]; B = np.ones(m)[:, None] * a[None, :]
    S = L - t
    R = (2.0 / L) * 0.5 * (_I(A - B, -B * t, S) - _I(A + B, B * t, S))
    return 0.5 * (R + R.T)

L = 4.5
print("=" * 80)
print("D  DOES THE OVERLAP BECOME RANK ONE AT THE SUPPORT EDGE?")
print("=" * 80)
print(f"   L = {L}.   eps = 1 - t/L  is the distance to the edge.")
for m in [40, 80, 160]:
    print(f"\n   basis size m = {m}")
    print(f"   {'eps':>10} {'sigma1':>12} {'sigma2/sigma1':>15} {'|cos(u1, v)|':>14}   v_j = j")
    a = np.arange(1, m + 1) * np.pi / L
    v = a / np.linalg.norm(a)
    for eps in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]:
        t = L * (1 - eps)
        R = overlap(t, L, m)
        U, s, _ = np.linalg.svd(R)
        u1 = U[:, 0]
        cos = abs(float(u1 @ v))
        print(f"   {eps:>10.0e} {s[0]:>12.3e} {s[1]/s[0]:>15.6f} {cos:>14.6f}")

print()
print("=" * 80)
print("E  THE SAME TEST ON THE ACTUAL JUMP, with the prime placed AT the edge")
print("=" * 80)
print("   choose L = log q + eps so that the entering prime sits just inside the support")
def prime_powers(cmax):
    out = []
    for p in primerange(2, int(cmax) + 2):
        k, q = 1, p
        while q <= cmax:
            out.append((q, math.log(p))); k += 1; q = p ** k
    return sorted(out)
pps = prime_powers(500)

def P(u, L, m):
    M = np.zeros((m, m))
    for q, Lam in pps:
        lq = math.log(q)
        if lq <= u and lq < L:
            M += -2.0 * Lam / math.sqrt(q) * (1.0 - lq / u) * overlap(lq, L, m)
    return M

m = 80
a = np.arange(1, m + 1) * np.pi
print(f"   m = {m}")
print(f"   {'q':>5} {'eps':>8} | {'meas/pred':>11} {'sigma2/sigma1':>14} {'|cos(u1,v)|':>13} {'eff.rank':>9}")
for q, Lam in [(13, math.log(13)), (29, math.log(29)), (49, math.log(7)), (97, math.log(97))]:
    lq = math.log(q)
    for eps in [3e-2, 3e-3, 3e-4]:
        Lx = lq / (1 - eps)
        h, gap = 1e-6, 1e-5
        J = (P(lq + gap + h, Lx, m) - P(lq + gap - h, Lx, m)) / (2 * h) \
            - (P(lq - gap + h, Lx, m) - P(lq - gap - h, Lx, m)) / (2 * h)
        pred = -2.0 * Lam / (math.sqrt(q) * lq) * overlap(lq, Lx, m)
        ratio = np.linalg.norm(J) / np.linalg.norm(pred)
        s = np.linalg.svd(J, compute_uv=False)
        U = np.linalg.svd(J)[0][:, 0]
        vv = (np.arange(1, m + 1) * np.pi / Lx); vv /= np.linalg.norm(vv)
        print(f"   {q:>5} {eps:>8.0e} | {ratio:>11.6f} {s[1]/s[0]:>14.6f} {abs(float(U@vv)):>13.6f} "
              f"{int(np.sum(s > 1e-10*s[0])):>9d}")

print()
print("=" * 80)
print("F  SUMMARY OF WHAT IS AND IS NOT REPRODUCED")
print("=" * 80)
print("""   scalar  -2 Lambda(q) / (sqrt q log q) : reproduced, but this is NOT evidence.
       With the triangular window w(q,u) = 1 - log q / u, differentiating gives
       -2 Lambda(q) q^{-1/2} (log q / u^2) R, and at u = log q that IS the scalar
       in one line.  Lambda(q) appears in the output because it was placed in the
       sum as input, prime powers included.  Reproducing it verifies this file's
       own arithmetic, not the cited theorem, whose content is the UNIVERSALITY
       of the matrix factor -- which section G measures, and which fails here.
   rank-one structure : appears only as the shift approaches the support edge,
       which is precisely the configuration L = log c.  Away from the edge the
       jump is full rank in this basis.""")

# ---------------------------------------------------------------------------
# G  IS THE MATRIX FACTOR UNIVERSAL ACROSS q?   (the content of the cited theorem)
# ---------------------------------------------------------------------------
print()
print("=" * 80)
print("G  UNIVERSALITY OF THE MATRIX FACTOR")
print("=" * 80)
Lg, mg = 4.5, 60
qs = [3, 5, 7, 9, 11, 13, 25, 27, 49]
Rn = {q: (lambda M: M / np.linalg.norm(M))(overlap(math.log(q), Lg, mg)) for q in qs}
print("   pairwise cosine between normalised matrix factors R(log q):")
print("        " + "".join(f"{q:>8}" for q in qs))
for a_ in qs:
    print(f"   {a_:>5}" + "".join(f"{float(np.sum(Rn[a_]*Rn[b])):>8.3f}" for b in qs))
off = [float(np.sum(Rn[a_]*Rn[b])) for i, a_ in enumerate(qs) for b in qs[i+1:]]
print(f"\n   off-diagonal: min {min(off):+.3f}  max {max(off):+.3f}  mean {np.mean(off):+.3f}")
print("   -> a universal (q-independent) factor would give 1.000 everywhere.")
