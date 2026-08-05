#!/usr/bin/env python3
"""
convergence_sweep.py -- how fast does the CvS ground state converge to the zeta
zeros, and does the ODD channel converge at the same rate as the EVEN one?

BACKGROUND.  The CvS ground state is taken in the EVEN parity sector; that is what
connes_cvs.compute_ground_state does and what the published numbers are.  The ODD
sector is the half spanned by the fixed-support (Dirichlet sine) family, and it
also reconstructs the zeros -- measured, at c = 13: even error 1.45e-55, odd
8.70e-53.  The open question this script feeds is whether those two errors follow
the same law in c, and that cannot be answered until each point is shown to be
converged in N and resolved in dps.

TWO THINGS THAT WILL FOOL YOU, both already caught once:

  * N MUST BE CONVERGED AT THAT c.  N = 60 is nearly converged at c = 13 and is
    badly unconverged at c = 19, where going to N = 80 moved the error by three
    orders AND reversed which channel was better.  Never quote a c-to-c comparison
    without an N-doubling beside it at the same c.
  * dps MUST EXCEED THE DEPTH.  The eigenvalues fall about 5.3 decades per unit of
    c.  At c = 17 they are ~1e-74 and dps = 80 is still exact (checked against
    dps = 150: every digit identical).  At c = 19, N = 80 they are ~1e-82, i.e.
    BELOW dps = 80, and that run has to be repeated at dps = 200 before anything
    is read off it.  The guard below FIRES at depth + 20 > dps, which is where a
    cell is demonstrably unsafe; the value it RECOMMENDS is depth + 40, which is
    where cells have actually been reproducible.  The two numbers are different on
    purpose: c = 13, N = 100 at dps = 80 has depth 59 and passes the guard, and it
    reproduces digit for digit at dps = 150.

USAGE
    pip install connes-cvs python-flint      # flint is ~100x; do not skip it
    python3 convergence_sweep.py 19 80 200   # one point: c, N, dps
    python3 convergence_sweep.py --sweep     # the four points that are needed
    python3 convergence_sweep.py --show      # print what is already in the file

Every finished point is appended as one JSON line to convergence_results.jsonl and
is never recomputed, so the sweep can be stopped and restarted freely.  Expect
roughly 10 min at (c=13, N=60, dps=80) and several hours at (c=23, N=80, dps=200);
cost grows with all three.
"""
from __future__ import annotations
import json, os, sys, time

try:
    import mpmath as mp
    from connes_cvs.operator import build_galerkin_matrix, extract_zeros, HAS_FLINT
except ImportError:
    sys.exit("need: pip install connes-cvs python-flint")

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "convergence_results.jsonl")
T = 400

# the points still needed, in the order they answer questions
SWEEP = [
    (19, 80, 200),   # is the channel reversal at c=19 real, or below the dps floor?
    (17, 80, 150),   # does the c=17 ratio survive an N-doubling?
    (23, 60, 150),   # a fourth c
    (23, 80, 200),   # its convergence check
]


def done() -> set:
    if not os.path.exists(RESULTS):
        return set()
    out = set()
    with open(RESULTS) as f:
        for line in f:
            line = line.strip()
            if line:
                r = json.loads(line)
                out.add((r["c"], r["N"], r["dps"]))
    return out


def run(c: int, N: int, dps: int) -> dict:
    mp.mp.dps = dps
    L = mp.log(c)
    print(f"  building c={c} N={N} dps={dps} T={T} ...", flush=True)
    t0 = time.time()
    Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    build_s = time.time() - t0
    print(f"  built in {build_s/60:.1f} min, dim {Q.rows}", flush=True)

    D = 2 * N + 1
    r = 1 / mp.sqrt(2)
    Ve = mp.zeros(D, N + 1); Ve[N, 0] = 1
    Vo = mp.zeros(D, N)
    for k in range(1, N + 1):
        Ve[N + k, k] = r; Ve[N - k, k] = r
        Vo[N + k, k - 1] = r; Vo[N - k, k - 1] = -r

    row = {"c": c, "N": N, "dps": dps, "T": T, "build_min": round(build_s / 60, 1)}
    for name, V in (("even", Ve), ("odd", Vo)):
        t0 = time.time()
        B = V.T * Q * V
        E, U = mp.eigsy(B)
        # nearest-to-zero, NOT the algebraically smallest: with a negative block the
        # two differ. At c <= 67 the matrix is observed positive on every cell so
        # they coincide, and the negative count below is recorded so that a cell
        # where they would not can never be read as if they did.
        n_neg = sum(1 for j in range(B.rows) if E[j] < 0)
        i = min(range(B.rows), key=lambda j: abs(E[j]))
        full = V * mp.matrix([U[j, i] for j in range(B.rows)])
        nrm = mp.sqrt(sum(full[j] ** 2 for j in range(D)))
        v = mp.zeros(D, 1)
        for j in range(D):
            v[j, 0] = full[j] / nrm
        z = extract_zeros(v, L, n_zeros=1, dps=dps)[0]
        row[name + "_lambda"] = mp.nstr(abs(E[i]), 8)
        row[name + "_n_negative"] = n_neg
        row[name + "_err"] = mp.nstr(z["error"], 8) if z["error"] is not None else None
        print(f"    {name:4s}: lambda = {row[name+'_lambda']}   gamma_1 error = "
              f"{row[name+'_err']}   ({time.time()-t0:.0f} s)", flush=True)

    # the warning that would otherwise be missed
    try:
        depth = -int(mp.floor(mp.log10(abs(mp.mpf(row["even_lambda"])))))
        if depth + 20 > dps:
            row["WARNING"] = (f"eigenvalue depth {depth} digits against dps {dps}: "
                              f"this point is at or below the precision floor, repeat "
                              f"it at dps >= {depth + 40}")
            print("  !! " + row["WARNING"], flush=True)
    except Exception as exc:                      # the guard must never fail silently
        row["WARNING"] = f"precision-floor check could not run: {exc!r}"
        print("  !! " + row["WARNING"], flush=True)
    return row


def show():
    if not os.path.exists(RESULTS):
        print("nothing computed yet")
        return
    rows = [json.loads(l) for l in open(RESULTS) if l.strip()]
    rows.sort(key=lambda r: (r["c"], r["N"], r["dps"]))
    print(f"{'c':>4} {'N':>4} {'dps':>5} {'even lambda':>14} {'even err':>14} "
          f"{'odd lambda':>14} {'odd err':>14} {'odd/even':>10}")
    for r in rows:
        try:
            ratio = f"{float(mp.mpf(r['odd_err']) / mp.mpf(r['even_err'])):.3g}"
        except Exception:
            ratio = "-"
        print(f"{r['c']:>4} {r['N']:>4} {r['dps']:>5} {r['even_lambda']:>14} "
              f"{r['even_err'] or '-':>14} {r['odd_lambda']:>14} "
              f"{r['odd_err'] or '-':>14} {ratio:>10}"
              + ("   <-- " + r["WARNING"][:40] if "WARNING" in r else ""))
    print("\nA ratio is meaningless unless the same c appears twice with different N")
    print("and the two agree. Check that before reading anything off this table.")


def main() -> int:
    if "--show" in sys.argv:
        show(); return 0
    print(f"python-flint: {HAS_FLINT}" + ("" if HAS_FLINT else "   <-- install it, this will crawl"))
    todo = SWEEP if "--sweep" in sys.argv else [tuple(int(x) for x in sys.argv[1:4])]
    have = done()
    warned = []
    for c, N, dps in todo:
        if (c, N, dps) in have:
            print(f"skipping c={c} N={N} dps={dps}, already in {os.path.basename(RESULTS)}")
            continue
        row = run(c, N, dps)
        with open(RESULTS, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"  saved.\n", flush=True)
        if "WARNING" in row:
            warned.append((c, N, dps))
    show()
    if warned:
        print("\nFAILED: these cells are at or below the precision floor and must be")
        print("        rerun at higher dps before any number from them is quoted:")
        for c, N, dps in warned:
            print(f"        c={c} N={N} dps={dps}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
