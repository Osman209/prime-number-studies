#!/usr/bin/env python3
"""
compare_families.py -- runs the fixed-support (Dirichlet sine) constructor and the
periodic Fourier (CvS) constructor on the SAME truncated Weil form and reports,
sector by sector, where the two agree.

Both constructors are reached through conventions.py section 7: builder_sine.build
and builder_periodic.build take the same (L, p_c, m) and return the same keys, so
neither family forks the harness.

WHAT IS BEING TESTED.  Section 7 predicts, from the two bases alone:

    ODD  sector -- the periodic combinations (e_k - e_{-k}) ~ sin(2 pi k x/L) vanish
                   at the centre and at both endpoints, so they ARE the Dirichlet
                   sine modes: the SAME operator, and the residual difference is
                   the archimedean truncation, which must fall with T.
    EVEN sector -- the periodic side carries the constant and cosines that do not
                   vanish at the endpoints, which the Dirichlet space cannot
                   represent at finite size: the difference must NOT fall with T.
                   (It should fall with m instead -- the odd-indexed sine modes are
                   complete in that subspace asymptotically -- but that is an
                   inference from completeness and is not what this script tests.)

Both are gated below, so the script fails if either half of that stops holding.

Usage:
    python3 compare_families.py            c = 13, N = 12, dps = 60, T = 100/200
    python3 compare_families.py --full     c = 13, N = 16, dps = 80, T = 100/200/400
Requires: pip install connes-cvs   (python-flint strongly recommended)
"""
from __future__ import annotations
import math, sys
import numpy as np

import builder_sine as BS
import builder_periodic as BP
from conventions import TAU_DEFAULT

FULL = "--full" in sys.argv
C = 13
N = 16 if FULL else 12
DPS = 80 if FULL else 60
TS = (100, 200, 400) if FULL else (100, 200)
FAILS = []


def entrywise(Wp, Ws, N, m):
    """max |G_ij - S_ij| on the ODD sector, the quantity the report quotes.

    G is the periodic matrix in the basis (e_k - e_{-k})/sqrt2; S is the sine
    matrix restricted to phi_2, phi_4, ..., phi_{2N}. The spectral gap below is a
    different quantity and the two must not be conflated."""
    _, Vo = BP.parity_projectors(N)
    G = Vo.T @ Wp @ Vo
    # sine_parity_index returns (j ODD, j EVEN) -- the first are the functions even
    # under the involution, the second the odd ones. The periodic odd sector pairs
    # with phi_2, phi_4, ..., i.e. the SECOND return value. Taking the first is the
    # wrong index set and gives a relative difference of order 1, not 1e-6.
    _, i_even_indexed = BP.sine_parity_index(m)
    S = Ws[np.ix_(i_even_indexed, i_even_indexed)]
    return float(np.abs(G - S).max()), float(np.abs(G - S).max() / np.abs(G).max())


def gap(a, b, k=6):
    """Largest of the top-k absolute eigenvalue gaps, scaled by the largest
    eigenvalue. A ratio would be dominated by the near-null end, which carries no
    information about whether the two operators agree."""
    return float(np.abs(a[-k:] - b[-k:]).max() / max(abs(a[-1]), abs(b[-1])))


def main() -> int:
    L = math.log(C)
    m = 2 * N
    print(f"c = {C}   L = log c = {L:.9f}   p_c = {C} (saturated)   N = {N}   m = {m}")

    sine = BS.build(L, C, m=m)
    print("  sine family     ", BS.inertia(sine["W"], TAU_DEFAULT))

    coup = float(np.abs(sine["W"][np.ix_(*BP.sine_parity_index(m))]).max())
    ok = coup < 1e-12
    print(f"  {'ok  ' if ok else 'FAIL'} parity block-diagonality of the sine form: "
          f"max off-block entry {coup:.2e}")
    if not ok:
        FAILS.append("parity block-diagonality")
    s = BP.sine_sectors(sine["W"], m)

    print(f"\n{'T':>6} {'build':>7}   {'ODD gap':>12}   {'EVEN gap':>12}   "
          f"{'nulls odd G/S':>14} {'even G/S':>10}")
    odd_track, even_track, entry_track = [], [], []
    for T in TS:
        import time
        t0 = time.time()
        per = BP.build(L, C, m=m, T=T, dps=DPS)
        g = BP.sectors(per["W"], per["N"])
        go, ge = gap(g["odd"], s["odd"]), gap(g["even"], s["even"])
        ew_abs, ew_rel = entrywise(per["W"], sine["W"], per["N"], m)
        entry_track.append(ew_rel)
        odd_track.append(go); even_track.append(ge)
        print(f"{T:>6} {time.time()-t0:>6.0f}s   {go:>12.3e}   {ge:>12.3e}   "
              f"[entrywise odd {ew_abs:.2e} abs, {ew_rel:.2e} rel]  "
              f"{int((np.abs(g['odd'])<TAU_DEFAULT).sum()):>7}/"
              f"{int((np.abs(s['odd'])<TAU_DEFAULT).sum()):<6}"
              f"{int((np.abs(g['even'])<TAU_DEFAULT).sum()):>5}/"
              f"{int((np.abs(s['even'])<TAU_DEFAULT).sum()):<4}")
        last = g

    print("\n  odd sector, largest six")
    print("    periodic:", " ".join(f"{v:.6f}" for v in last["odd"][-6:]))
    print("    sine    :", " ".join(f"{v:.6f}" for v in s["odd"][-6:]))
    print("  even sector, largest six")
    print("    periodic:", " ".join(f"{v:.6f}" for v in last["even"][-6:]))
    print("    sine    :", " ".join(f"{v:.6f}" for v in s["even"][-6:]))

    print()
    ok = all(odd_track[i + 1] < odd_track[i] / 3 for i in range(len(odd_track) - 1))
    print(f"  {'ok  ' if ok else 'FAIL'} the ODD-sector gap falls by at least 3x per doubling "
          f"of T: {' -> '.join(f'{v:.2e}' for v in odd_track)}")
    if not ok:
        FAILS.append("odd sector converging in T")

    ok = max(even_track) / min(even_track) < 1.5 and even_track[-1] > 10 * odd_track[-1]
    print(f"  {'ok  ' if ok else 'FAIL'} the EVEN-sector gap stays flat instead: "
          f"{' -> '.join(f'{v:.2e}' for v in even_track)}")
    if not ok:
        FAILS.append("even sector structurally different")

    print("\n  Reading: on the odd sector the two constructors agree and the residual is")
    print("  the archimedean truncation; on the even one they do not agree at this m, and")
    print("  the gap does not respond to T. Whether it responds to m is not tested here.")

    ok = all(entry_track[i + 1] < entry_track[i] / 3 for i in range(len(entry_track) - 1))
    print(f"  {'ok  ' if ok else 'FAIL'} the ENTRYWISE odd-sector difference falls with T "
          f"too: {' -> '.join(f'{v:.2e}' for v in entry_track)}")
    if not ok:
        FAILS.append("entrywise odd-sector agreement converging in T")

    if FAILS:
        print(f"\nFAILED: {len(FAILS)} check(s): " + ", ".join(FAILS))
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
