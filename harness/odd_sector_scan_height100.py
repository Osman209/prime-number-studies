"""Spurious-root check: extend the unseeded sign-change scan to a larger height,
reusing an eigenvector already cached by odd_sector_unseeded_scan.py. No matrix build,
so this runs in under a minute.

    python3 odd_sector_scan_height100.py <vector_cache/..._odd_vector.json> [max_t]

Writes the discovered root list to scan_roots.json before any comparison with
mp.zetazero, so the discovery stays unseeded. Roots are found by sign change on a
uniform sweep and refined by bisection only."""
import json, sys, time
import mpmath as mp

if len(sys.argv) < 2:
    sys.exit(__doc__)
VECTOR = sys.argv[1]
SCAN_DPS, REF_DPS = 40, 120
LO, HI, STEP = mp.mpf('0.1'), mp.mpf(sys.argv[2] if len(sys.argv) > 2 else 100), mp.mpf('0.02')

with open(VECTOR, encoding='utf-8') as fh:
    d = json.load(fh)
c, N, T = d['c'], d['N'], d['T']
mp.mp.dps = REF_DPS
L = mp.log(c)
coef = [mp.mpf(v) for v in d['vector']]
assert len(coef) == 2 * N + 1

def F(tau, dps):
    """Same reconstruction as extract_zeros, on the odd-sector coefficients."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        t = mp.mpf(tau)
        Lm = mp.mpf(L)
        e = mp.exp(-1j * t * Lm)
        tot = mp.mpc(0, 0)
        for k in range(-N, N + 1):
            ck = coef[k + N]
            if ck == 0:
                continue
            den = 2 * mp.pi * k / Lm - t
            tot += ck * (mp.mpc(Lm, 0) if abs(den) < mp.mpf('1e-30') else (e - 1) / (1j * den))
        return mp.re(mp.exp(1j * t * Lm / 2) * tot / mp.sqrt(Lm))
    finally:
        mp.mp.dps = old

print(f"c={c} N={N} T={T}; scanning ({LO}, {HI}) at step {STEP}, dps {SCAN_DPS} -> refine {REF_DPS}", flush=True)
t0 = time.time()
brackets = []
x = LO
fx = F(x, SCAN_DPS)
n = 0
while x < HI:
    y = x + STEP
    fy = F(y, SCAN_DPS)
    if fx == 0:
        brackets.append((x, x))
    elif mp.sign(fx) != mp.sign(fy):
        brackets.append((x, y))
    x, fx = y, fy
    n += 1
    if n % 1000 == 0:
        print(f"  {float(x):.1f} / {float(HI)}  ({time.time()-t0:.0f} s, {len(brackets)} brackets)", flush=True)
print(f"scan done in {time.time()-t0:.0f} s: {len(brackets)} sign changes", flush=True)

roots = []
for a, b in brackets:
    lo, hi = mp.mpf(a), mp.mpf(b)
    fl = F(lo, REF_DPS)
    for _ in range(400):
        mid = (lo + hi) / 2
        fm = F(mid, REF_DPS)
        if fm == 0:
            lo = hi = mid; break
        if mp.sign(fm) != mp.sign(fl):
            hi = mid
        else:
            lo, fl = mid, fm
        if hi - lo < mp.mpf(10) ** (-REF_DPS + 20):
            break
    roots.append((lo + hi) / 2)

mp.mp.dps = REF_DPS + 30
json.dump([mp.nstr(r, 60) for r in roots], open('scan_roots.json', 'w'))
print("root list written to scan_roots.json BEFORE any zeta comparison\n")
gam = [mp.im(mp.zetazero(k)) for k in range(1, 60)]
print(f"\n{len(roots)} roots refined. Comparing to zeta ordinates below {HI}:\n")
print(f"{'#':>3} {'root':>22} {'nearest gamma_k':>16} {'k':>4} {'|r-gamma|':>12}")
used = {}
for i, r in enumerate(roots, 1):
    k = min(range(len(gam)), key=lambda j: abs(r - gam[j]))
    err = abs(r - gam[k])
    used.setdefault(k, []).append(i)
    print(f"{i:>3} {mp.nstr(r,18):>22} {mp.nstr(gam[k],14):>16} {k+1:>4} {mp.nstr(err,4):>12}")
below = [j for j in range(len(gam)) if gam[j] < HI]
print(f"\nzeta ordinates below {HI}: {len(below)}   roots found: {len(roots)}")
missing = [j + 1 for j in below if j not in used]
dup = {k + 1: v for k, v in used.items() if len(v) > 1}
print("missing ordinates:", missing if missing else "none")
print("ordinates claimed by more than one root:", dup if dup else "none")

# --- gate: the completeness claim this script exists to support ---
problems = []
if missing:
    problems.append(f"{len(missing)} zeta ordinate(s) below {HI} were not found: {missing}")
if dup:
    problems.append(f"one ordinate claimed by several roots: {dup}")
if len(roots) != len(below):
    problems.append(f"{len(roots)} roots against {len(below)} ordinates — a spurious "
                    f"sign change, or two roots merged by the step {STEP}")
if problems:
    print("\nFAILED:")
    for p_ in problems:
        print("  " + p_)
    sys.exit(1)
print(f"\ngate passes: {len(roots)} roots, {len(below)} ordinates, none missing or doubled")
