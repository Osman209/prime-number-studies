#!/usr/bin/env python3
"""
verify_row_inheritance.py -- regenerates every number printed in section 11,
"Dynamic row inheritance", of *A Numerical Study of the Division Table*.

Odd model throughout: the candidates are the odd numbers, the columns are the
odd primes, and column p strikes the numbers whose least prime factor is p.

Checks
------
  1  A_p = p * S_p^-                            (elementwise, K = 16,666)
  2  Delta a_p(k) / (2p) = Delta m(k) / 2       (elementwise, first 1000)
  3  X_p(K) table, the ratio law p/(p'-1), and the efficiency = Mertens product
  4  X_p(K) ~ e^gamma K p log p
  5  motion table: mean step, longest step, number of distinct step sizes
  6  column 7's step histogram {14: 6249, 28: 6250, 42: 4166}
  7  collapse: a_p(k)/a_p(K) against k/K is the diagonal
  8  phi(2pP) = (p-1) phi(2P), and the deleted set is exactly p * S_p^-
  9  G_new = max(g_i + g_{i+1}) on complete cycles, 7 consecutive stages
 10  the merge-capacity multiset at struck points equals the old multiset
 11  the window table of 7.6\n 12  cycle lengths, first occurrence of the record gap, its multiplicity N, 2P/N\n 13  the window shadow of the multiset identity

Usage
-----
    python verify_row_inheritance.py           full  (windows to 10^8)
    python verify_row_inheritance.py --fast    quick (windows to 10^7)
    python verify_row_inheritance.py --deep    adds the 10^9 window column

Exits 0 only if every check passes.
"""

import sys, math
import numpy as np

FAST = "--fast" in sys.argv
DEEP = "--deep" in sys.argv
K = 16666
OPS = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
GAMMA = 0.5772156649015329
FAIL = []


def fail(tag, msg):
    FAIL.append(tag)
    print(f"    FAIL [{tag}] {msg}")


# ----------------------------------------------------------------- least prime factor
XMAX = 7_000_000
print(f"building the least-prime-factor sieve to {XMAX:,} ...", flush=True)
spf = np.zeros(XMAX + 1, dtype=np.int32)
for p in range(3, XMAX + 1, 2):
    if spf[p] == 0:
        sl = spf[p * p::2 * p]
        spf[p * p::2 * p] = np.where(sl == 0, p, sl)
spf[3::2] = np.where(spf[3::2] == 0, np.arange(3, XMAX + 1, 2), spf[3::2])
odds = np.arange(1, XMAX + 1, 2)


def contributions(p, n):
    c = odds[(spf[odds] == p) & (odds != p)]
    return c[:n]


def survivors_before(p, upto):
    s = odds[(odds > 1) & (spf[odds] >= p)]
    return s[s <= upto]


def mertens(p):
    m = 1.0
    for q in OPS:
        if q >= p:
            break
        m *= 1 - 1 / q
    return m


# ----------------------------------------------------------------- 1 and 2
print("\nCHECK 1 - A_p = p * S_p^-  (elementwise)")
for p in OPS[1:10]:
    A = contributions(p, K)
    S = survivors_before(p, int(A[-1]) // p)
    S = S[S >= p]
    if not np.array_equal(A, p * S[:len(A)]):
        fail("1", f"p={p}")
    else:
        print(f"    p={p:>3}: {len(A):,} contributions identical to p*S_p^-")

print("\nCHECK 2 - Delta a_p(k)/(2p) = Delta m(k)/2  (first 1000 steps)")
for p in (5, 7, 11, 13, 17):
    A = contributions(p, 1001)
    d = np.diff(A)
    S = survivors_before(p, int(A[-1]) // p)
    S = S[S >= p]
    g = np.diff(S[:1001]) // 2
    ok = np.array_equal(d // (2 * p), g) and bool(np.all(d % (2 * p) == 0))
    print(f"    p={p:>3}: {'ok' if ok else 'MISMATCH'}")
    if not ok:
        fail("2", f"p={p}")

# ----------------------------------------------------------------- 3 and 4
print("\nCHECK 3 - X_p(K), the ratio law, and the efficiency")
paper_X = [99999, 249995, 437507, 802043, 1042691, 1477181, 1754251, 2241419,
           2954317, 3269477, 4030891, 4588433, 4929907, 5516437, 6355813]
print(f"    {'p':>4} {'X_p':>11} {'paper':>11} {'ratio':>8} {'p/(p-1)':>9} {'eff':>7} {'prod':>7}")
Xs = {}
prevX = prevp = None
for i, p in enumerate(OPS):
    A = contributions(p, K)
    X = int(A[-1]); Xs[p] = X
    raw = int(np.count_nonzero((odds % p == 0) & (odds <= X) & (odds != p)))
    eff = K / raw
    r = X / prevX if prevX else 1.0
    rt = p / (prevp - 1) if prevp else 1.0
    print(f"    {p:>4} {X:>11,} {paper_X[i]:>11,} {r:>8.4f} {rt:>9.4f} {eff:>7.4f} {mertens(p):>7.4f}")
    if X != paper_X[i]:
        fail("3", f"p={p}: X_p = {X}, paper {paper_X[i]}")
    if prevX and abs(r - rt) / rt > 0.005:
        fail("3", f"p={p}: ratio {r:.4f} vs p/(p'-1) = {rt:.4f}")
    if abs(eff - mertens(p)) > 2e-3:
        fail("3", f"p={p}: efficiency {eff:.4f} vs Mertens {mertens(p):.4f}")
    prevX, prevp = X, p

print("\nCHECK 4 - X_p(K) / (K p log p) against e^gamma = 1.7811")
vals = [Xs[p] / (K * p * math.log(p)) for p in OPS]
print("    " + " ".join(f"{v:.3f}" for v in vals))
print(f"    min = {min(vals):.3f}, max = {max(vals):.3f}, e^gamma = {math.exp(GAMMA):.4f}   (paper: 1.808 to 1.927)")
if abs(round(min(vals),3) - 1.808) > 1e-9 or abs(round(max(vals),3) - 1.927) > 1e-9:
    fail("4", f"range [{min(vals):.3f}, {max(vals):.3f}] vs paper [1.808, 1.927]")

# ----------------------------------------------------------------- 5 and 6
print("\nCHECK 5 - the motion table")
paper = {3: (6.00, 6, 1), 5: (15.00, 20, 2), 7: (26.25, 42, 3), 11: (48.12, 110, 5),
         13: (62.56, 182, 7), 17: (88.62, 374, 10), 19: (105.24, 456, 12),
         23: (134.47, 782, 13), 29: (177.23, 986, 15), 31: (196.13, 1054, 15)}
print(f"    {'p':>4} {'mean':>9} {'2p/prod':>9} {'longest':>8} {'p*G_prev':>9} {'#sizes':>7}")
for p in OPS[:10]:
    A = contributions(p, K); d = np.diff(A)
    S = survivors_before(p, int(A[-1]) // p); S = S[S >= p]
    gprev = int(np.diff(S[:K]).max())
    mean, longest, nsizes = float(d.mean()), int(d.max()), len(set(d.tolist()))
    pm, pl, pn = paper[p]
    print(f"    {p:>4} {mean:>9.2f} {2*p/mertens(p):>9.2f} {longest:>8} {p*gprev:>9} {nsizes:>7}")
    if abs(mean - pm) > 0.01 or longest != pl or nsizes != pn:
        fail("5", f"p={p}: ({mean:.2f},{longest},{nsizes}) vs paper ({pm},{pl},{pn})")
    if longest != p * gprev:
        fail("5", f"p={p}: longest step {longest} != p*G_prev = {p*gprev}")

print("\nCHECK 6 - column 7's step histogram")
d = np.diff(contributions(7, K))
v, c = np.unique(d, return_counts=True)
got = dict(zip(v.tolist(), c.tolist()))
norm = sorted({int(x) // 14 for x in v.tolist()})
print(f"    {got}   normalised by 2p = 14: {norm}")
if got != {14: 6249, 28: 6250, 42: 4166}:
    fail("6", f"{got}")
if norm != [1, 2, 3]:
    fail("6", f"normalised step sizes {norm}, paper {{1,2,3}}")

# ----------------------------------------------------------------- 7
print("\nCHECK 7 - collapse onto the diagonal")
worst = 0.0
for p in OPS[:10]:
    A = contributions(p, K).astype(float)
    x = np.arange(1, K + 1) / K
    dev = float(np.max(np.abs(A / A[-1] - x)))
    worst = max(worst, dev)
    if p in (3, 31):
        step = float(np.max(np.diff(A))) / A[-1]
        print(f"    p={p:>3}: max deviation {dev:.5f}   (longest step)/X_p = {step:.6f}")
        if p == 31 and abs(round(step, 6) - 0.000322) > 1e-9:
            fail("7", f"(longest step)/X_p = {step:.6f}, paper 3.2e-4")
if worst > 0.002:
    fail("7", f"worst deviation {worst:.5f}")

# ----------------------------------------------------------------- 8, 9, 10
def cycle(ps):
    P = 1
    for q in ps:
        P *= q
    m = np.ones(P, bool)
    for q in ps:
        m[(q - 1) // 2::q] = False
    return P, 2 * np.nonzero(m)[0] + 1        # odds coprime to prod(ps), period 2P


print("\nCHECK 8 - phi(2pP) = (p-1) phi(2P), and the deleted set is p * S_p^-")
for k in range(1, 6):
    ps = OPS[:k]; q = OPS[k]
    P, v = cycle(ps)
    Pn, vn = cycle(ps + [q])
    if len(vn) != (q - 1) * len(v):
        fail("8", f"q={q}: {len(vn)} != {(q-1)*len(v)}")
    big = np.concatenate([v + 2 * P * j for j in range(q)])
    deleted = set(big.tolist()) - set(vn.tolist())
    if deleted != set((q * v).tolist()):
        fail("8", f"q={q}: deleted set is not q*S")
    print(f"    q={q:>3}: survivors {len(v):>5} -> {len(vn):>6} = (q-1)*{len(v)}, deleted set = q*S")

print("\nCHECK 9 - G_new = max(g_i + g_{i+1}) on complete cycles")
prevH = None
chain = []
for k in range(1, 9):
    ps = OPS[:k]
    P, v = cycle(ps)
    g = np.diff(np.concatenate([v, [v[0] + 2 * P]]))
    G = int(g.max()); H = int((g + np.roll(g, -1)).max())
    chain.append(G)
    if prevH is not None and G != prevH:
        fail("9", f"columns<= {ps[-1]}: G = {G}, predicted {prevH}")
    prevH = H
print(f"    chain of records: {' -> '.join(map(str, chain))} -> {prevH}")
if chain != [4, 6, 10, 14, 22, 26, 34, 40]:
    fail("9", f"chain {chain}")

print("\nCHECK 10 - merge-capacity multiset at struck points = old multiset")
for k in range(1, 6):
    ps = OPS[:k]; q = OPS[k]
    P, v = cycle(ps)
    gg = np.diff(np.concatenate([v, [v[0] + 2 * P]]))
    C_old = np.roll(gg, 1) + gg
    big = np.concatenate([v + 2 * P * j for j in range(q)])
    g = np.diff(np.concatenate([big, [big[0] + 2 * P * q]]))
    C_big = np.roll(g, 1) + g
    dead = big % q == 0
    same = np.array_equal(np.sort(C_big[dead]), np.sort(C_old))
    print(f"    q={q:>3}: struck {int(dead.sum()):>5}, families {len(v):>5}, multisets identical: {same}")
    if not same:
        fail("10", f"q={q}")

# ----------------------------------------------------------------- 11
print("\nCHECK 11 - the window table of 7.6")
CHECK = [10**5, 10**6, 10**7] + ([] if FAST else [10**8]) + ([10**9] if DEEP else [])
Gtrue = {19: 34, 23: 40, 29: 46, 31: 58, 37: 66, 43: 90, 53: 106}
paper_rows = {19: [34, 34, 34, 34, 34], 23: [34, 34, 36, 40, 40], 29: [34, 34, 40, 42, 46],
              31: [36, 36, 48, 48, 50], 37: [40, 46, 48, 50, 54], 43: [46, 46, 50, 54, 60],
              53: [58, 58, 58, 60, 64]}


def maxgap_curve(ps, checkpoints, seg=25_000_000):
    res = {}; cur = 0; last = None; start = 1
    for cp in checkpoints:                       # no segment crosses a checkpoint
        while start <= cp:
            end = min(start + seg - 1, cp)
            lo = start if start % 2 else start + 1
            n = (end - lo) // 2 + 1
            if n > 0:
                m = np.ones(n, bool)
                for p in ps:
                    inv2 = pow(2, p - 2, p)
                    m[(-lo * inv2) % p::p] = False
                idx = np.nonzero(m)[0]
                if len(idx):
                    vals = lo + 2 * idx
                    if last is not None:
                        cur = max(cur, int(vals[0] - last))
                    if len(vals) > 1:
                        cur = max(cur, int(np.max(np.diff(vals))))
                    last = int(vals[-1])
            start = end + 1
        res[cp] = cur
    return res


hdr = " ".join(f"{'1e%d' % int(round(math.log10(c))):>6}" for c in CHECK)
print(f"    {'rows<=':>7} {hdr} {'cycle':>7}")
for p in (19, 23, 29, 31, 37, 43, 53):
    k = OPS.index(p) + 1
    r = maxgap_curve(OPS[:k], CHECK)
    row = [r[c] for c in CHECK]
    print(f"    {p:>7} " + " ".join(f"{x:>6}" for x in row) + f" {Gtrue[p]:>7}")
    for j, c in enumerate(CHECK):
        want = paper_rows[p][j]
        if row[j] != want:
            fail("11", f"cols<={p}, X=1e{int(round(math.log10(c)))}: {row[j]} vs paper {want}")
    if row[-1] > Gtrue[p]:
        fail("11", f"cols<={p}: window value exceeds the cycle value")

# ----------------------------------------------------------------- verdict
print("\nCHECK 12 - cycle lengths quoted in the text, and the first occurrence of the record gap")
for p_, want in ((19, "9.7e6"), (53, "3.26e19")):
    P = 1
    for q in OPS[:OPS.index(p_) + 1]:
        P *= q
    print(f"    rows<= {p_:>3}: 2P = {2*P:>22,}  = {2*P:.3g}   (paper {want})")
if 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 != 9699690:
    fail("12", "2P for rows<=19 is not 9,699,690")
Pall = 1
for q in OPS:
    Pall *= q
if f"{2*Pall:.3g}" != "3.26e+19":
    fail("12", f"2P for rows<=53 is {2*Pall:.3g}, paper 3.26e19")



def first_occurrence(ps, G, limit, seg=20_000_000):
    """smallest X such that a gap of size >= G already appears in [1,X]"""
    start, last = 1, None
    while start <= limit:
        end = min(start + seg - 1, limit)
        lo = start if start % 2 else start + 1
        n = (end - lo) // 2 + 1
        if n > 0:
            m = np.ones(n, bool)
            for q in ps:
                inv2 = pow(2, q - 2, q)
                m[(-lo * inv2) % q::q] = False
            idx = np.nonzero(m)[0]
            if len(idx):
                vals = lo + 2 * idx
                if last is not None and vals[0] - last >= G:
                    return int(vals[0])
                d = np.diff(vals)
                hit = np.nonzero(d >= G)[0]
                if len(hit):
                    return int(vals[hit[0] + 1])
                last = int(vals[-1])
        start = end + 1
    return None


paper_L = {13: 9461, 17: 217153, 19: 60077, 23: 20332511}
rows = [13, 17, 19] + ([] if FAST else [23])
print(f"    {'rows<=':>7} {'G':>4} {'2P':>14} {'N':>4} {'2P/N':>14} {'first occurrence':>17} {'paper':>12}")
for p_ in rows:
    k = OPS.index(p_) + 1
    P, v = cycle(OPS[:k])
    g = np.diff(np.concatenate([v, [v[0] + 2 * P]]))
    G = int(g.max()); N = int((g == G).sum())
    L = first_occurrence(OPS[:k], G, 2 * P)
    print(f"    {p_:>7} {G:>4} {2*P:>14,} {N:>4} {2*P/N:>14,.0f} {L:>17,} {paper_L[p_]:>12,}")
    if L != paper_L[p_]:
        fail("12", f"rows<= {p_}: first occurrence {L}, paper {paper_L[p_]}")
    if N > 20:
        fail("12", f"rows<= {p_}: multiplicity {N} is not small")

print("\nCHECK 13 - the window shadow of the multiset identity (merge capacity, [1,10^5])")
W = 100_000
o = np.arange(1, W + 1, 2)
S = o.copy()
ratios = []
print(f"    {'row':>4} {'mean C, all survivors':>22} {'mean C, struck':>16} {'ratio':>8}")
for q in OPS[:9]:
    dead = (S % q == 0) & (S != q)
    g = np.diff(S)
    C_all = g[:-1] + g[1:]
    idx = np.nonzero(dead)[0]
    ii = idx[(idx > 0) & (idx < len(S) - 1)]
    C_del = (S[ii] - S[ii - 1]) + (S[ii + 1] - S[ii])
    r = float(C_del.mean() / C_all.mean())
    ratios.append(r)
    print(f"    {q:>4} {C_all.mean():>22.4f} {C_del.mean():>16.4f} {r:>8.4f}")
    S = S[~dead]
print(f"    range: {min(ratios):.4f} to {max(ratios):.4f}   (paper: 0.997 to 1.006)")
if abs(round(min(ratios), 3) - 0.997) > 1e-9 or abs(round(max(ratios), 3) - 1.006) > 1e-9:
    fail("13", f"window ratios span [{min(ratios):.4f}, {max(ratios):.4f}], paper [0.997, 1.006]")

print("\n" + "=" * 62)
if FAIL:
    print(f"FAILED checks: {sorted(set(FAIL))}")
    sys.exit(1)
print("all checks passed" + ("   (--fast: windows to 10^7)" if FAST else "")
      + ("   (--deep: 10^9 window included)" if DEEP else ""))
sys.exit(0)
