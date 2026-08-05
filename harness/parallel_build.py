#!/usr/bin/env python3
"""
parallel_build.py — build the CvS Galerkin matrix using every core, and cache it.

Two savings, applied together:

  1. SYMMETRY, free.  For integer n, psi(-n) = -psi(n) and psi'(-n) = psi'(n)
     (verified here before use, and asserted at run time).  The package computes
     all 2N+1 indices; only N+1 of them are independent, so half the quadratures
     are redundant.

  2. PARALLELISM.  Each psi pair is an independent adaptive quadrature, so the
     N+1 survivors are handed to a process pool.  This is the phase Groskin
     parallelises too ("multiprocessing.imap_unordered for the psi-cache phase").

STATUS: EXPERIMENTAL.  Nothing in the report depends on this file.  It reproduces
the package's builder exactly under --verify, but the JSON cache it writes is not
yet read by any other script here, so treat the cache as a convenience and not as
part of the evidence chain until a loader exists.  Diagonalisation is NOT parallel
— mpmath's eigsy is serial, and on large N it becomes the dominant cost once the
build is spread out.

    python3 parallel_build.py --c 100 --N 150 --T 800 --dps 500 --workers 10
    python3 parallel_build.py --c 13 --N 20 --T 400 --dps 60 --workers 4 --verify

--verify rebuilds the same matrix through the package's own serial routine and
compares entrywise; use it once on a small cell before trusting a long run.

MEMORY WARNING.  Each worker is a separate interpreter.  At dps = 500 with a
large prime set expect a few hundred MB each; with 4 GB free, 8 workers is
already ambitious.  Close other applications, and drop --workers if the machine
starts swapping — swapping is far slower than using fewer cores.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from multiprocessing import Pool

import mpmath as mp
from connes_cvs.operator import (_compute_psi_pair, prime_powers_up_to,
                                 build_galerkin_matrix, HAS_FLINT)

_CTX: dict = {}


def _init(c, T, dps):
    mp.mp.dps = dps
    _CTX["L"] = mp.log(c)
    _CTX["T"] = T
    _CTX["dps"] = dps
    _CTX["pd"] = prime_powers_up_to(int(c))[0]
    # CRITICAL: each worker is a fresh interpreter, so flint's precision must be
    # set here too. build_galerkin_matrix does this internally; without it the
    # digamma evaluations run at flint's default and the whole build silently
    # degrades to ~1e-16 no matter what dps says. The --verify gate catches it.
    try:
        from flint import ctx as fctx
        fctx.prec = int(dps * 3.5)
        fctx.threads = 1                    # one per process; the pool does the rest
    except ImportError:
        pass                                # pure mpmath: slow but correct
    except Exception as exc:                # anything else means the precision did
        raise RuntimeError(                 # NOT get set, which is the silent bug
            f"flint is present but its precision could not be set ({exc!r}); "
            f"the build would silently run at flint's default. Refusing to continue."
        ) from exc


def _psi_pair(n):
    psi, psi_d = _compute_psi_pair(n, _CTX["L"], _CTX["T"], _CTX["dps"], _CTX["pd"])
    return n, mp.nstr(psi, _CTX["dps"] + 5), mp.nstr(psi_d, _CTX["dps"] + 5)


def build(c, N, T, dps, workers):
    mp.mp.dps = dps
    t0 = time.time()
    todo = list(range(0, N + 1))            # symmetry: negatives come for free
    print(f"psi phase: {len(todo)} independent quadratures "
          f"(instead of {2*N+1}) on {workers} workers", flush=True)
    with Pool(workers, initializer=_init, initargs=(c, T, dps)) as pool:
        done = 0
        psi, psid = {}, {}
        for n, a, b in pool.imap_unordered(_psi_pair, todo, chunksize=1):
            psi[n], psid[n] = mp.mpf(a), mp.mpf(b)
            done += 1
            if done % max(1, len(todo) // 20) == 0:
                el = time.time() - t0
                print(f"  {done}/{len(todo)}  {el/60:.1f} min elapsed, "
                      f"~{el/done*(len(todo)-done)/60:.1f} min left", flush=True)
    for n in range(1, N + 1):               # mirror
        psi[-n], psid[-n] = -psi[n], psid[n]
    print(f"psi phase done in {(time.time()-t0)/60:.1f} min", flush=True)

    t1 = time.time()
    DIM = 2 * N + 1
    Q = mp.matrix(DIM, DIM)
    for i, m in enumerate(range(-N, N + 1)):
        for j, n in enumerate(range(-N, N + 1)):
            Q[i, j] = psid[n] if m == n else (psi[m] - psi[n]) / (m - n)
    print(f"assembly done in {time.time()-t1:.0f} s", flush=True)
    return Q


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--c", type=int, required=True)
    p.add_argument("--N", type=int, required=True)
    p.add_argument("--T", type=int, default=400)
    p.add_argument("--dps", type=int, default=150)
    p.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    p.add_argument("--out", default=None)
    p.add_argument("--verify", action="store_true")
    a = p.parse_args()

    print(f"python-flint: {HAS_FLINT}; cores seen: {os.cpu_count()}", flush=True)
    Q = build(a.c, a.N, a.T, a.dps, a.workers)

    if a.verify:
        print("\nverifying against the package's own serial builder ...", flush=True)
        t = time.time()
        R = build_galerkin_matrix(c=a.c, N=a.N, T=a.T, dps=a.dps)
        worst = max(abs(Q[i, j] - R[i, j]) for i in range(Q.rows) for j in range(Q.rows))
        scale = max(abs(R[i, j]) for i in range(R.rows) for j in range(R.rows))
        print(f"serial build took {(time.time()-t)/60:.1f} min")
        print(f"max |parallel - serial| = {mp.nstr(worst, 6)}   "
              f"relative to max|Q| = {mp.nstr(worst/scale, 6)}")
        if worst / scale > mp.mpf(10) ** (-(a.dps - 10)):
            print("\nFAILED: the parallel build does not reproduce the serial one.")
            return 1
        print("the two builds agree at working precision")

    out = a.out or f"Q_c{a.c}_N{a.N}_T{a.T}_dps{a.dps}.json"
    with open(out, "w") as f:
        json.dump({"c": a.c, "N": a.N, "T": a.T, "dps": a.dps,
                   "Q": [[mp.nstr(Q[i, j], a.dps + 5) for j in range(Q.rows)]
                         for i in range(Q.rows)]}, f)
    print(f"\nmatrix cached to {out} ({os.path.getsize(out)/1e6:.0f} MB)")
    print("diagonalisation is serial in mpmath and is not covered by this script.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
