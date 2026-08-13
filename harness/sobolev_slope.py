#!/usr/bin/env python3
"""
sobolev_slope.py -- tests the Sobolev-Paley-Wiener prediction of Groskin
arXiv:2605.20224 §8.2 / §11.

THE PREDICTION.  §8.1 measures a Galerkin convergence exponent s(c) from the
pre-saturation regime N in {40, 60, 80} by regressing log10|lambda_min^even|
against log10 N; the slope is -2s.  Across six cutoffs it fits

    s(c) ~ 55 log c - 128.

§8.2 proposes a Paley-Wiener mechanism giving s = sigma_eff * T / (2 pi) with
sigma_eff ~ A log c, i.e. s(c, T) proportional to T.  §11 states the test: repeat
the measurement at T = 1600; if the slope of s against log c doubles, the
mechanism is supported; if it stays flat, it is not.

AN AMBIGUITY THE AUTHOR HAS SINCE RESOLVED.  §8.2 treats the measured slope 55 as
a T = 800 datum ("With T = 800 and the measured slope 55, this gives A = 0.432"),
while Table 8 — the c = 23 N-sweep the measurement rests on — is labelled T = 400.
Groskin has confirmed (correspondence, August 2026) that TABLE 14 IS T = 400,
dps = 150, on N in {40, 60, 80}, and that §8.2's T = 800 came from the sweep
convention; he is filing the caption as an erratum candidate.  So the T = 1600
prediction to test against is ~220, not ~110.

Stage 1 below is therefore a reproduction of his published slope on his own grid,
not a disambiguation.  Stage 2 is the test.

USAGE
    pip install connes-cvs python-flint
    python3 sobolev_slope.py --pilot     small, ~10 min, checks the machinery
    python3 sobolev_slope.py --stage1    c in {17,23}, T = 400   (~80 min, settles the T ambiguity)
    python3 sobolev_slope.py --stage2    the same at T = 1600  (~5 h, the actual test)
    python3 sobolev_slope.py --show      print what is already computed
    python3 sobolev_slope.py C N T DPS   one cell

Results append to sobolev_results.jsonl and are never recomputed, so the run can
be stopped and restarted.  Expect the T = 1600 cells to cost roughly four times
the T = 400 ones: the archimedean quadrature is the dominant term.
"""
from __future__ import annotations
import json, math, os, sys, time

try:
    import mpmath as mp
    from connes_cvs.operator import build_galerkin_matrix, HAS_FLINT
except ImportError:
    sys.exit("need: pip install connes-cvs python-flint")

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sobolev_results.jsonl")

# N values are Groskin's pre-saturation regime; c and T are what we vary
N_REGIME = (40, 60, 80)
PILOT = [(13, n, 400, 80) for n in (20, 25, 30)] + [(13, n, 1600, 80) for n in (20, 25, 30)]
# THE GRID IS THE ONE GROSKIN SPECIFIED for a cross-cutoff slope: N in {40, 60, 80},
# T = 400, dps = 150, c = 17 and c = 23.  c = 13 is deliberately EXCLUDED — it is
# already saturating on that grid, which is why its row in his Table 14 carries
# R^2 = 0.87.  Stage 1 should reproduce s(17) = 27.8 and s(23) = 46.1, a slope of
# 60.5 against log c; Stage 2 repeats it at T = 1600, where the §8.2 mechanism
# predicts ~220.  Stage 2 costs about four times as much per cell.
STAGE1 = [(c, n, 400, 150) for c in (17, 23) for n in N_REGIME]
STAGE2 = [(c, n, 1600, 250) for c in (17, 23) for n in N_REGIME]
SWEEP = STAGE1 + STAGE2


def done() -> set:
    if not os.path.exists(RESULTS):
        return set()
    out = set()
    for line in open(RESULTS):
        line = line.strip()
        if line:
            r = json.loads(line)
            out.add((r["c"], r["N"], r["T"], r["dps"]))
    return out


def even_lambda_min(c: int, N: int, T: int, dps: int):
    """Smallest-positive even-sector eigenvalue, Groskin's quantity."""
    mp.mp.dps = dps
    t0 = time.time()
    Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    build_min = (time.time() - t0) / 60
    D = 2 * N + 1
    r = 1 / mp.sqrt(2)
    Ve = mp.zeros(D, N + 1)
    Ve[N, 0] = 1
    for k in range(1, N + 1):
        Ve[N + k, k] = r
        Ve[N - k, k] = r
    E = mp.eigsy(Ve.T * Q * Ve, eigvals_only=True)
    vals = [E[i] for i in range(N + 1)]
    pos = [v for v in vals if v > 0]
    lam = min(pos) if pos else min(abs(v) for v in vals)
    return lam, build_min, sum(1 for v in vals if v < 0)


def run(cell):
    c, N, T, dps = cell
    print(f"  building c={c} N={N} T={T} dps={dps} ...", flush=True)
    lam, build_min, n_neg = even_lambda_min(c, N, T, dps)
    row = {"c": c, "N": N, "T": T, "dps": dps,
           "lambda_min_even": mp.nstr(lam, 10),
           "log10_lambda": float(mp.log10(abs(lam))),
           "n_negative": n_neg, "build_min": round(build_min, 1)}
    print(f"    log10|lambda| = {row['log10_lambda']:.3f}   "
          f"({build_min:.1f} min, {n_neg} negative)", flush=True)
    depth = -int(math.floor(row["log10_lambda"]))
    if depth + 20 > dps:
        row["WARNING"] = (f"depth {depth} digits against dps {dps}: at or below the "
                          f"precision floor, repeat at dps >= {depth + 40}")
        print("  !! " + row["WARNING"], flush=True)
    return row


def slope_from(rows):
    """s = -(slope of log10|lambda| against log10 N)/2, Groskin §8.1."""
    if len(rows) < 2:
        return None, None
    xs = [math.log10(r["N"]) for r in rows]
    ys = [r["log10_lambda"] for r in rows]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    a = sxy / sxx
    syy = sum((y - my) ** 2 for y in ys)
    r2 = (sxy ** 2) / (sxx * syy) if syy else 1.0
    return -a / 2, r2


def show():
    if not os.path.exists(RESULTS):
        print("nothing computed yet")
        return
    rows = [json.loads(l) for l in open(RESULTS) if l.strip()]
    print(f"{'c':>4} {'T':>6} {'N values':>16} {'s(c,T)':>10} {'R^2':>8}")
    by = {}
    for r in rows:
        by.setdefault((r["c"], r["T"]), []).append(r)
    slopes = {}
    for (c, T), rs in sorted(by.items()):
        rs.sort(key=lambda r: r["N"])
        s, r2 = slope_from(rs)
        slopes[(c, T)] = s
        print(f"{c:>4} {T:>6} {str([r['N'] for r in rs]):>16} "
              f"{(f'{s:.2f}' if s else '-'):>10} {(f'{r2:.5f}' if r2 else '-'):>8}"
              + ("   <-- precision warning on a cell" if any("WARNING" in r for r in rs) else ""))
    print("\nslope of s against log c, per T  (Groskin's fit gives 55; §11 predicts it "
          "should scale linearly with T):")
    for T in sorted({T for _, T in slopes}):
        pts = [(c, slopes[(c, T)]) for c, TT in slopes if TT == T and slopes[(c, TT)]]
        if len(pts) >= 2:
            (c1, s1), (c2, s2) = sorted(pts)[0], sorted(pts)[-1]
            m = (s2 - s1) / (math.log(c2) - math.log(c1))
            print(f"  T = {T:>5}:  s({c1}) = {s1:.1f}, s({c2}) = {s2:.1f}  ->  slope = {m:.1f}")
        else:
            print(f"  T = {T:>5}:  need at least two cutoffs")
    print("\nReading: if the slope roughly doubles from T=400 to 800 to 1600, the")
    print("Paley-Wiener mechanism of §8.2 is supported; if it is flat, it is not.")


def main() -> int:
    if "--show" in sys.argv:
        show(); return 0
    print(f"python-flint: {HAS_FLINT}" + ("" if HAS_FLINT else "   <-- install it"))
    if "--pilot" in sys.argv:
        todo = PILOT
    elif "--stage1" in sys.argv:
        todo = STAGE1
    elif "--stage2" in sys.argv:
        todo = STAGE2
    elif "--sweep" in sys.argv:
        todo = SWEEP
    else:
        try:
            todo = [tuple(int(x) for x in sys.argv[1:5])]
        except (ValueError, TypeError):
            sys.exit("usage: sobolev_slope.py <c> <N> <T> <dps>   or   --pilot | "
                     "--stage1 | --stage2 | --sweep | --show")
    have = done()
    bad = []
    for cell in todo:
        if cell in have:
            print(f"skipping {cell}, already computed")
            continue
        row = run(cell)
        with open(RESULTS, "a") as f:
            f.write(json.dumps(row) + "\n")
        # At c <= 67 the finite matrix is observed positive on every cell (Groskin
        # §5.1), so a negative eigenvalue here is an arithmetic failure, not a result.
        # This is the failure mode that produced log10|lambda| ~ -3.4 at T = 1600,
        # dps = 80, forty-four orders adrift, with five negatives — and a floor test
        # on digit count did NOT catch it, because the failure made lambda too
        # shallow rather than too deep.
        if row.get("c", 0) <= 67 and row.get("n_negative", 0) > 0:
            bad.append(row)
    show()
    if bad:
        print("\nFAILED: negative eigenvalues at a cutoff where the matrix should be")
        print("        positive on every cell. Rerun these at higher dps; do not fit")
        print("        a slope through them:")
        for row in bad:
            print(f"        c={row['c']} N={row['N']} T={row['T']} dps={row['dps']}"
                  f"  ({row['n_negative']} negative, log10|lambda| = {row.get('log10_lambda')})")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
