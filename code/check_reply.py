import numpy as np, math

def _I(w, ph, S):
    w = np.asarray(w, float); out = np.empty_like(w); sm = np.abs(w) < 1e-13
    out[sm] = S * np.cos(ph[sm] if np.ndim(ph) else ph)
    ws = w[~sm]; phs = ph[~sm] if np.ndim(ph) else ph
    out[~sm] = (np.sin(ws * S + phs) - np.sin(phs)) / ws
    return out

def overlap(t, L, m):
    if t >= L: return np.zeros((m, m))
    a = np.arange(1, m + 1) * np.pi / L
    A = a[:, None] * np.ones(m)[None, :]; B = np.ones(m)[:, None] * a[None, :]; S = L - t
    R = (2.0 / L) * 0.5 * (_I(A - B, -B * t, S) - _I(A + B, B * t, S))
    return 0.5 * (R + R.T)

L = 4.5
print("=" * 78)
print("1  THE LAW  sigma2/sigma1 -> (m-1)/(m+2)")
print("=" * 78)
print("   columns are eps = 1 - t/L.  The shift itself is delta = L - t = L*eps,")
print(f"   so eps = 1e-4 here means delta = {L:.1f}e-4.")
print(f"   {'m':>5} {'predicted':>12} " + " ".join(f"eps=1e-{k}" for k in [3, 4, 5, 6]))
for m in [20, 40, 80, 160, 200]:
    pred = (m - 1) / (m + 2)
    vals = []
    for k in [3, 4, 5, 6]:
        R = overlap(L * (1 - 10.0 ** (-k)), L, m)
        s = np.linalg.svd(R, compute_uv=False)
        vals.append(f"{s[1]/s[0]:.7f}")
    print(f"   {m:>5} {pred:>12.7f} " + " ".join(f"{v:>11}" for v in vals))
print()
print("   NOTE.  The eps = 1e-6 column is float64 noise, not a worse approximation:")
print("   R is O(delta^3) formed from O(delta) pieces, so the cancellation has eaten")
print("   the answer by then.  Taken in extended precision the discrepancy keeps")
print("   falling like delta^2 -- see harness/edge_precision.py, check 2.")

print()
print("=" * 78)
print("2  DO OPPOSITE PARITIES DECOUPLE?   (the ratio at eps=1e-6 is a noise floor,")
print("   not a failure: same-parity has itself collapsed to ~1e-14 by then)")
print("=" * 78)
m = 80
idx = np.arange(1, m + 1)
odd = np.where(idx % 2 == 1)[0]; even = np.where(idx % 2 == 0)[0]
for eps in [1e-2, 1e-4, 1e-6]:
    R = overlap(L * (1 - eps), L, m)
    cross = np.abs(R[np.ix_(odd, even)]).max()
    same = max(np.abs(R[np.ix_(odd, odd)]).max(), np.abs(R[np.ix_(even, even)]).max())
    print(f"   eps={eps:.0e}:  max|cross-parity| = {cross:.3e}   max|same-parity| = {same:.3e}"
          f"   ratio = {cross/same:.3e}")

print()
print("=" * 78)
print("3  THE EXPANSION  R_jk(L-d) = (pi^2 d^3 / 6 L^3) j k ((-1)^(j+1)+(-1)^(k+1)) + O(d^5)")
print("   (here d IS delta = L - t.  The rise at d = 1e-5 is cancellation, not the")
print("   remainder: in extended precision the ratio stays at 100 -- edge_precision.py)")
print("=" * 78)
m = 40
j = np.arange(1, m + 1)
sgn = (-1.0) ** (j + 1)
Lead = np.outer(j, j) * (sgn[:, None] + sgn[None, :])
print(f"   {'delta':>10} {'||R - lead||/||R||':>20} {'scale ratio':>14}")
for d in [1e-2, 1e-3, 1e-4, 1e-5]:
    R = overlap(L - d, L, m)
    pref = math.pi ** 2 * d ** 3 / (6 * L ** 3)
    P = pref * Lead
    print(f"   {d:>10.0e} {np.linalg.norm(R - P)/np.linalg.norm(R):>20.3e} "
          f"{np.linalg.norm(R)/np.linalg.norm(P):>14.6f}")

print()
print("=" * 78)
print("4  IS THE SECOND ENDPOINT VECTOR  v_j = (-1)^(j+1) j ?")
print("=" * 78)
for m in [40, 80, 160]:
    j = np.arange(1, m + 1)
    v1 = j / np.linalg.norm(j)
    v2 = ((-1.0) ** (j + 1) * j); v2 /= np.linalg.norm(v2)
    R = overlap(L * (1 - 1e-5), L, m)
    U = np.linalg.svd(R)[0]
    print(f"   m={m:>4}:  |cos(u1, j)| = {abs(float(U[:,0]@v1)):.6f}   "
          f"|cos(u1, (-1)^(j+1) j)| = {abs(float(U[:,0]@v2)):.6f}   "
          f"|cos(u2, (-1)^(j+1) j)| = {abs(float(U[:,1]@v2)):.6f}")

print()
print("=" * 78)
print("5  SPLIT ODD AND EVEN MODES: IS EACH BLOCK RANK ONE, WITH OPPOSITE SIGN?")
print("   measured at eps = 1e-5, i.e. delta = L*1e-5")
print("=" * 78)
for m in [40, 80, 160]:
    idx = np.arange(1, m + 1)
    o = np.where(idx % 2 == 1)[0]; e = np.where(idx % 2 == 0)[0]
    R = overlap(L * (1 - 1e-5), L, m)
    Ro, Re = R[np.ix_(o, o)], R[np.ix_(e, e)]
    so = np.linalg.svd(Ro, compute_uv=False); se = np.linalg.svd(Re, compute_uv=False)
    wo = np.linalg.eigvalsh(Ro); we = np.linalg.eigvalsh(Re)
    print(f"   m={m:>4}:  odd  sigma2/sigma1 = {so[1]/so[0]:.3e}  dominant eig sign = "
          f"{'+' if abs(wo[-1])>abs(wo[0]) else '-'}")
    print(f"          even sigma2/sigma1 = {se[1]/se[0]:.3e}  dominant eig sign = "
          f"{'+' if abs(we[-1])>abs(we[0]) else '-'}")

print()
print("=" * 78)
print("6  IS THE SECTION-C 'JUMP' AN ARTIFACT THAT VANISHES WITH gap?")
print("=" * 78)
from sympy import primerange
def prime_powers(cmax):
    out = []
    for p in primerange(2, int(cmax) + 2):
        k, q = 1, p
        while q <= cmax: out.append((q, math.log(p))); k += 1; q = p ** k
    return sorted(out)
pps = prime_powers(500)

def P_coupled(u, m):
    M = np.zeros((m, m))
    for q, Lam in pps:
        lq = math.log(q)
        if lq <= u:
            M += -2.0 * Lam / math.sqrt(q) * (1.0 - lq / u) * overlap(lq, u, m)
    return M

m = 40; q = 13; u0 = math.log(q)
print(f"   coupled path, q = {q},  m = {m}")
print(f"   {'gap':>10} {'h':>10} {'||measured jump||':>20}")
for g in [1e-3, 1e-4, 1e-5, 1e-6]:
    h = g / 10
    J = (P_coupled(u0 + g + h, m) - P_coupled(u0 + g - h, m)) / (2 * h) \
        - (P_coupled(u0 - g + h, m) - P_coupled(u0 - g - h, m)) / (2 * h)
    print(f"   {g:>10.0e} {h:>10.0e} {np.linalg.norm(J):>20.6e}")
print("\n   the entering term is (1 - t/u) R(t;u) = O((u-t)^4), so the true first-derivative")
print("   jump on this path is ZERO; a nonzero reading is the smooth background over a finite gap.")
