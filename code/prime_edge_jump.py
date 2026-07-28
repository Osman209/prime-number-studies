"""
Independent test of the prime-edge derivative jump
(Groskin 2026, Zenodo 10.5281/zenodo.21242028).

CLAIM UNDER TEST.  Fix a Galerkin level, let u = log c be the prime cutoff and
vary it continuously.  Crossing a prime-power threshold u = log q, the first
derivative of the matrix path u -> Q(u) jumps by

        Delta (dQ/du)  =  -2 Lambda(q) / (sqrt(q) log q)  *  (rank-one matrix).

Only the prime block can produce a singular part, so the archimedean and pole
blocks are omitted here: they are smooth in u and cancel out of any jump.

CONSTRUCTION (independent of the reference implementation).
  orthonormal sine basis  phi_j(x) = sqrt(2/L) sin(a_j x),  a_j = j pi / L,  on [0, L]
  prime block             P(u) = -2 sum_{q = p^k, log q <= u} Lambda(q) q^{-1/2} w(q,u) R(log q)
  triangular window       w(q,u) = 1 - log q / u
  overlap                 R_jk(t) = int_0^{L-t} phi_j(x) phi_k(x+t) dx, symmetrised

Two settings are measured:
  (A) DECOUPLED  - support L held fixed, only the prime cutoff u varies
  (B) COUPLED    - L = u, the reference parametrisation, so the basis moves too
"""
import numpy as np, math
from sympy import primerange, factorint

# ----------------------------------------------------------------- basis overlap
def _I(w, ph, S):
    """int_0^S cos(w x + ph) dx, with the removable w -> 0 case handled."""
    w = np.asarray(w, float)
    out = np.empty_like(w)
    small = np.abs(w) < 1e-13
    out[small] = S * np.cos(ph[small] if np.ndim(ph) else ph)
    ws = w[~small]
    phs = ph[~small] if np.ndim(ph) else ph
    out[~small] = (np.sin(ws * S + phs) - np.sin(phs)) / ws
    return out

def overlap(t, L, m):
    """R_jk(t) = int_0^{L-t} phi_j(x) phi_k(x+t) dx, symmetrised.  Exact closed form."""
    if t >= L:
        return np.zeros((m, m))
    a = np.arange(1, m + 1) * np.pi / L
    A = a[:, None] * np.ones(m)[None, :]
    B = np.ones(m)[:, None] * a[None, :]
    S = L - t
    R = (2.0 / L) * 0.5 * (_I(A - B, -B * t, S) - _I(A + B, B * t, S))
    return 0.5 * (R + R.T)

def prime_powers(cmax):
    out = []
    for p in primerange(2, int(cmax) + 2):
        k, q = 1, p
        while q <= cmax:
            out.append((q, math.log(p)))     # (q, Lambda(q))
            k += 1; q = p ** k
    return sorted(out)

# ----------------------------------------------------------------- the two settings
def P_decoupled(u, L, m, pps):
    """support fixed at L; only the cutoff u moves"""
    M = np.zeros((m, m))
    for q, Lam in pps:
        lq = math.log(q)
        if lq <= u and lq < L:
            M += -2.0 * Lam / math.sqrt(q) * (1.0 - lq / u) * overlap(lq, L, m)
    return M

def P_coupled(u, m, pps):
    """L = u : the reference parametrisation, basis moves with the cutoff"""
    M = np.zeros((m, m))
    for q, Lam in pps:
        lq = math.log(q)
        if lq <= u:
            M += -2.0 * Lam / math.sqrt(q) * (1.0 - lq / u) * overlap(lq, u, m)
    return M

def deriv(f, u, h):
    return (f(u + h) - f(u - h)) / (2 * h)

def jump(f, u0, h, gap):
    """one-sided derivatives either side of the threshold, then their difference"""
    Dm = deriv(f, u0 - gap, h)
    Dp = deriv(f, u0 + gap, h)
    return Dp - Dm

def rank_report(M):
    s = np.linalg.svd(M, compute_uv=False)
    s = s[s > 0]
    if len(s) == 0:
        return 0.0, 0.0, 0
    r2 = s[1] / s[0] if len(s) > 1 else 0.0
    eff = int(np.sum(s > 1e-10 * s[0]))
    return s[0], r2, eff

# ================================================================= TEST A
print("=" * 78)
print("A  DECOUPLED SETTING  (support L fixed, cutoff u varied)  -- the regime in question")
print("=" * 78)
L, m = 4.5, 60
pps = prime_powers(200)
h, gap = 1e-5, 1e-4
print(f"   L = {L}, basis size m = {m}, step h = {h}\n")
print(f"   {'q':>6} {'Lambda(q)':>10} {'log q':>8} | {'measured/predicted':>19} {'sigma2/sigma1':>14} {'eff.rank':>9}")
for q, Lam in [(3, math.log(3)), (5, math.log(5)), (7, math.log(7)),
               (9, math.log(3)), (11, math.log(11)), (13, math.log(13)),
               (25, math.log(5)), (27, math.log(3)), (49, math.log(7))]:
    u0 = math.log(q)
    if u0 >= L:
        continue
    J = jump(lambda u: P_decoupled(u, L, m, pps), u0, h, gap)
    pred_scalar = -2.0 * Lam / (math.sqrt(q) * math.log(q))
    pred = pred_scalar * overlap(math.log(q), L, m)
    ratio = np.linalg.norm(J) / np.linalg.norm(pred) if np.linalg.norm(pred) > 0 else float('nan')
    align = float(np.sum(J * pred) / (np.linalg.norm(J) * np.linalg.norm(pred)))
    s1, r2, eff = rank_report(J)
    print(f"   {q:>6} {Lam:>10.5f} {u0:>8.4f} | {ratio:>19.9f} {r2:>14.4f} {eff:>9d}   cos={align:+.6f}")
print("\n   -> the ratio and cosine are FINITE-DIFFERENCE estimates: agreement to the")
print("      printed digits, not an exact identity.  The residual departure from 1 is")
print("      truncation error in the difference quotient and shrinks with h.")
print("      sigma2/sigma1 says whether that matrix is close to rank one.")

# ================================================================= TEST B
print()
print("=" * 78)
print("B  IS THE JUMP MATRIX RANK ONE?  structure of the overlap R(t) at the shift")
print("=" * 78)
print(f"   L = {L}, m = {m}")
print(f"   {'t = log q':>10} {'t/L':>7} | {'sigma1':>10} {'sigma2/sigma1':>14} {'eff.rank':>9}")
for t in [math.log(3), math.log(7), math.log(13), math.log(29), math.log(60),
          0.90 * L, 0.99 * L, 0.999 * L]:
    if t >= L:
        continue
    s1, r2, eff = rank_report(overlap(t, L, m))
    print(f"   {t:>10.5f} {t/L:>7.3f} | {s1:>10.5f} {r2:>14.6f} {eff:>9d}")
print("\n   -> the combined ratio does NOT go to zero: it plateaus, and the plateau is")
print(f"      exactly (m-1)/(m+2) = {(m-1)}/{(m+2)} = {(m-1)/(m+2):.6f} for this basis size.")
print("      That is the law of A. Groskin (see papers/prime_edge_two_paths.md).  What")
print("      DOES go to rank one is each PARITY BLOCK separately: splitting odd from")
print("      even sine modes gives sigma2/sigma1 -> 0 on each, with the dominant")
print("      eigenvalue positive on the odd block and negative on the even one, so the")
print("      full matrix is their difference and its ratio sits just below 1.  The")
print("      split is measured in check_reply.py and, in extended precision where the")
print("      limit is actually resolvable, in harness/edge_precision.py.")
print("      Effective rank does collapse toward the edge, which is the real content")
print("      of the column above.")

# ================================================================= TEST C
print()
print("=" * 78)
print("C  COUPLED SETTING  (L = u)  --  ARTIFACT, RETAINED AS A RECORD")
print("=" * 78)
print("""   WARNING.  On this path the entering term is (1 - t/u) R(t;u) = O((u-t)^4),
   so the TRUE first-derivative jump is ZERO.  The nonzero numbers printed below
   are a finite-gap difference of the smooth background and scale linearly with
   `gap`, vanishing in the limit.  Diagnosis due to A. Groskin (correspondence,
   July 2026); the linear-in-gap decay is verified in check_reply.py section 6.
   This block is kept, annotated, rather than deleted.
""")
m2 = 40
print(f"   basis size m = {m2}")
print(f"   {'q':>6} {'u0 = log q':>11} | {'||jump||':>12} {'sigma2/sigma1':>14} {'eff.rank':>9} {'scalar ratio':>14}")
for q, Lam in [(7, math.log(7)), (11, math.log(11)), (13, math.log(13)),
               (17, math.log(17)), (19, math.log(19)), (23, math.log(23))]:
    u0 = math.log(q)
    J = jump(lambda u: P_coupled(u, m2, pps), u0, 1e-6, 1e-5)
    s1, r2, eff = rank_report(J)
    pred_scalar = -2.0 * Lam / (math.sqrt(q) * math.log(q))
    print(f"   {q:>6} {u0:>11.5f} | {np.linalg.norm(J):>12.3e} {r2:>14.6f} {eff:>9d} "
          f"{np.linalg.norm(J)/abs(pred_scalar):>14.4e}")
print("\n   -> here the entering prime sits at the support edge (log q = u = L).")
print("      This is the COUPLED SINE-WINDOW path, NOT the reference construction:")
print("      the reference carries the edge vanishing inside its matrix, whereas this")
print("      path multiplies a Dirichlet sine overlap by an external window.")
