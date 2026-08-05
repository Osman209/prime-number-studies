"""Cross-check the unseeded scan against the package's own root finder.

    python3 cross_check_extract_zeros.py <vector_cache/..._odd_vector.json> [scan_roots.json]

The reconstruction F in odd_sector_scan_height100.py is re-implemented from the
published formula; a slip there would move every root and nothing else in the pipeline
would notice. This runs `extract_zeros` from the package on the SAME eigenvector — a
different F, a different solver, and a bracket seeded at mp.zetazero — and compares the
root positions with the unseeded ones. Agreement at the working-precision floor means
the seeding changes where one looks, not what is found.
"""
import json, sys
import mpmath as mp
from connes_cvs.operator import extract_zeros

if len(sys.argv) < 2:
    sys.exit(__doc__)
VECTOR = sys.argv[1]
OURS = sys.argv[2] if len(sys.argv) > 2 else 'scan_roots.json'
DPS = 120

mp.mp.dps = DPS
with open(VECTOR, encoding='utf-8') as fh:
    d = json.load(fh)
c, N = d['c'], d['N']
L = mp.log(c)
v = mp.zeros(2 * N + 1, 1)
for i, s in enumerate(d['vector']):
    v[i, 0] = mp.mpf(s)

with open(OURS, encoding='utf-8') as fh:
    ours = [mp.mpf(x) for x in json.load(fh)]
n = len(ours)
print(f"c={c} N={N} T={d.get('T')}: comparing {n} unseeded roots against extract_zeros\n")

res = extract_zeros(v, L, n_zeros=n, dps=DPS)
mp.mp.dps = DPS + 30
g = [mp.im(mp.zetazero(k)) for k in range(1, n + 1)]

print(f"{'k':>3} {'seeded err':>13} {'unseeded err':>14} {'|difference|':>14}")
worst = mp.mpf(0)
for i in range(n):
    e_his = res[i]['error']
    e_our = abs(ours[i] - g[i])
    det = res[i]['gamma_detected']
    diff = abs(det - ours[i]) if det is not None else None
    if diff is not None and diff > worst:
        worst = diff
    print(f"{i+1:>3} {mp.nstr(e_his,4) if e_his is not None else 'None':>13} "
          f"{mp.nstr(e_our,4):>14} {mp.nstr(diff,4) if diff is not None else 'no root':>14}")
smallest_error = min(abs(ours[i] - g[i]) for i in range(n))
print(f"\nlargest disagreement between the two implementations: {mp.nstr(worst,4)}")
print(f"smallest error being reported:                        {mp.nstr(smallest_error,4)}")

# The two implementations are not bitwise identical — different near-pole cutoffs and a
# different summation order, so F itself differs at the level of the cancellation inside
# it. What matters is not that they agree exactly, but that they agree far better than
# the quantity being reported: otherwise the errors in the table would be measuring the
# implementation rather than the operator. Three orders of headroom is the condition.
if worst * 1000 > smallest_error:
    print("\nFAILED: the two implementations disagree at a level comparable to the")
    print("        errors being reported, so those errors would be implementation-")
    print("        limited. Check the reconstruction of F and the root-finder tolerances.")
    sys.exit(1)
print(f"the implementations agree {float(mp.log10(smallest_error / worst)):.0f} orders "
      f"below the smallest reported error")
