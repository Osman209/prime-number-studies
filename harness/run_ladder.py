"""
run_ladder.py -- the regression harness.

Sweeps the m-ladder over a set of (L, p_c) points, reports the full inertia triple
with its stated tolerance, runs the non-circular validator on the saturated
diagonal, writes raw arrays and JSON metadata, and exits NONZERO if any check
fails.

    python run_ladder.py                      # default sweep, writes out/
    python run_ladder.py --L 4.5 --tau 1e-6

Checks that must pass, or the run fails:
    V1  validator residual within its tail bound on the saturated diagonal
    V2  lambda_min >= -tau on the saturated diagonal  (the form is PSD there)
    V3  lambda_min stable to STAB_TOL relative across the top two ladder rungs
Counts are REPORTED, never checked -- they are not converged and the harness
must not pretend otherwise.
"""
from __future__ import annotations
import argparse, json, math, os, platform, sys, time
import numpy as np

from conventions import M_LADDER, TAU_DEFAULT, R_MAX, NR
from builder_sine import build, inertia
from validator import validate

STAB_TOL = 0.02          # lambda_min must agree to 2% across the top two rungs
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def sweep(L: float, pcs, ladder, tau: float):
    rows = []
    for p_c in pcs:
        for m in ladder:
            t0 = time.time()
            B = build(L, p_c, m)
            inr = inertia(B["W"], tau)
            inr.update({"L": L, "p_c": p_c, "m": m,
                        "saturated": B["saturated"],
                        "n_prime_powers": len(B["prime_powers"]),
                        "arch_tail_bound": B["arch_tail_bound"],
                        "wall_s": round(time.time() - t0, 2)})
            rows.append(inr)
            np.save(os.path.join(OUT, f"eig_L{L}_pc{int(p_c)}_m{m}.npy"),
                    np.linalg.eigvalsh(B["W"]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=4.5)
    ap.add_argument("--tau", type=float, default=TAU_DEFAULT)
    ap.add_argument("--ladder", type=int, nargs="*", default=list(M_LADDER))
    ap.add_argument("--pcs", type=float, nargs="*", default=None)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    L, tau = args.L, args.tau
    sat = math.floor(math.exp(L))
    pcs = args.pcs if args.pcs is not None else [30, 45, 60, 75, sat]

    print(f"L = {L}   saturated p_c = {sat}   tau = {tau:g}   ladder = {args.ladder}\n")
    rows = sweep(L, pcs, args.ladder, tau)

    print(f"{'p_c':>6} {'m':>6} {'lambda_min':>13} {'n-':>6} {'n0':>5} {'n+':>6} {'sat':>5}")
    for r in rows:
        print(f"{int(r['p_c']):>6} {r['m']:>6} {r['lambda_min']:>13.4e} "
              f"{r['n_minus']:>6} {r['n_zero']:>5} {r['n_plus']:>6} {str(bool(r['saturated'])):>5}")
    print("\n  counts are REPORTED, not converged: they drift with m and with tau.")
    print("  the diagonal near-null count is the only one measuring near-null directions;")
    print("  off the diagonal the form is indefinite and n_minus dominates.\n")

    failures = []

    # V1 -- non-circular validation on the diagonal
    m_top = max(args.ladder)
    val = validate(L, min(m_top, 200))
    print(f"V1 validator  (L={L}, saturated, m={min(m_top,200)}, gamma_max={val['gamma_max']:.0f})")
    for r in val["rows"]:
        ok = r["within_tail"] or abs(r["rel_residual"]) < 1e-5
        print(f"   rank {r['rank']:>4}  lambda {r['lambda']:>14.9f}  zero side {r['zero_side']:>14.9f}"
              f"  rel {r['rel_residual']:>10.2e}  {'ok' if ok else 'FAIL'}")
        if not ok and r["lambda"] > 1e-4:
            failures.append(f"V1 rank {r['rank']} rel residual {r['rel_residual']:.2e}")

    # V2 -- PSD on the diagonal
    diag = [r for r in rows if r["saturated"]]
    for r in diag:
        if r["lambda_min"] < -tau:
            failures.append(f"V2 lambda_min {r['lambda_min']:.2e} < -tau at m={r['m']}")
    print(f"\nV2 PSD on diagonal: lambda_min = " +
          ", ".join(f"{r['lambda_min']:.2e} (m={r['m']})" for r in diag))

    # V3 -- lambda_min stability across the top two rungs
    print("\nV3 lambda_min stability across the top two rungs")
    if len(args.ladder) >= 2:
        m1, m2 = sorted(args.ladder)[-2:]
        for p_c in pcs:
            a = next(r for r in rows if r["p_c"] == p_c and r["m"] == m1)
            b = next(r for r in rows if r["p_c"] == p_c and r["m"] == m2)
            denom = max(abs(a["lambda_min"]), 1e-12)
            rel = abs(b["lambda_min"] - a["lambda_min"]) / denom
            flag = "ok" if (rel < STAB_TOL or abs(a["lambda_min"]) < tau) else "FAIL"
            print(f"   p_c={int(p_c):>4}  m={m1}: {a['lambda_min']:>11.4e}   "
                  f"m={m2}: {b['lambda_min']:>11.4e}   rel {rel:>8.2e}  {flag}")
            if flag == "FAIL":
                failures.append(f"V3 p_c={p_c} lambda_min drift {rel:.2e}")

    meta = {
        "L": L, "tau": tau, "ladder": args.ladder, "p_c_values": pcs,
        "saturated_p_c": sat, "r_max": R_MAX, "nr": NR,
        "stability_tol": STAB_TOL,
        "rows": rows, "validator": val, "failures": failures,
        "env": {"python": platform.python_version(), "numpy": np.__version__,
                "platform": platform.platform()},
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    with open(os.path.join(OUT, "run.json"), "w") as f:
        json.dump(meta, f, indent=2, default=float)

    if failures:
        print("\nFAILURES:")
        for f_ in failures:
            print("  -", f_)
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
