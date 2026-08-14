#!/usr/bin/env python3
"""
negative_branch.py — build the CvS Galerkin matrix if needed, then probe the
NEGATIVE-sign block of its even sector.  One command, two cached stages.

arXiv:2605.20224 §6.6 reports that at c = 100 the even block carries a handful of
eigenvalues with a negative sign and that the count grows with N.

READ THE CURRENT VERSION OF THAT PAPER BEFORE READING ANYTHING INTO THIS SCRIPT'S
OUTPUT.  The June 2026 revision attributes those blocks to the FINITE ARCHIMEDEAN
CUTOFF and reports that they disappear as T grows and in cutoff-free evaluation, so
the smallest-positive branch is the genuine lowest eigenvalue of the cutoff-free
problem.  An earlier draft of this file quoted the pre-revision "this still requires
certification" as though the question were open; it is not, and that mistake cost a
full rewrite of the accompanying note.

What remains worth measuring is narrower: how the FINITE-T branches behave — which
is what this script does.  It says nothing about the cutoff-free operator.

THIS SCRIPT DOES NOT SETTLE THAT.  It applies one calibrated test: a vector whose
reconstruction lands on the zeta ordinates does so by a very large margin (measured
at c = 13, N = 100: the ground state reaches every one of the first ten ordinates to
between 1e-52 and 1e-35, while a random vector in the same basis gets no closer than
0.144 to any of them), so each branch falls on one side of that gap or the other.
Either outcome is a measurement, not a certification, and a positive one would not
separate a structural feature from an artifact that happens to preserve the
reconstruction.

TWO CACHES, AND WHY THE STAGES STAY SEPARATE INSIDE ONE FILE.  The build costs ~55
min at c = 100, N = 150; the scan costs minutes.  Every question so far has been
answered by rescanning the SAME matrix at a different --step, --ref-dps or
--n-positive.  So:

    Q_c<c>_N<N>_T<T>_dps<dps>.json           the matrix        — built once
    Q_..._<sector>_eig.json                  the eigenvectors  — diagonalised once
    Q_..._<sector>_branch_scan.json          the roots         — rewritten each scan

The matrix is shared by both parity sectors, so `--sector odd` costs a diagonalisation
and a scan and no rebuild at all.  Both later caches carry the sector in their names.

A stage is skipped whenever its cache is present.  --rebuild and --rediagonalise
force one stage; nothing forces a rebuild by accident.

    python3 negative_branch.py --c 100 --N 150 --T 400 --dps 500 --workers 6
    python3 negative_branch.py --c 100 --N 150 --T 400 --dps 500 --step 0.005
    python3 negative_branch.py --c 13 --N 20 --T 400 --dps 60 --workers 4 --verify

The psi phase writes each quadrature to a JSONL as it lands, so an interrupted build
resumes where it stopped — which matters because memory, not cores, is the binding
constraint on a long high-precision run.
"""
from __future__ import annotations
import argparse, hashlib, json, os, random, sys, time
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

    # RESUME. Each psi pair is appended to a JSONL cache the moment it lands, so a
    # crash, a swap-death or a Ctrl-C costs only the quadratures still in flight.
    # Re-running the same command picks up exactly where it stopped. This matters
    # because memory, not cores, is the binding constraint on a long high-dps run.
    cache = f"psi_c{c}_N{N}_T{T}_dps{dps}.jsonl"
    psi, psid = {}, {}
    if os.path.exists(cache):
        for line in open(cache, encoding="utf-8"):
            if line.strip():
                r = json.loads(line)
                psi[r["n"]], psid[r["n"]] = mp.mpf(r["psi"]), mp.mpf(r["psid"])
        print(f"resuming: {len(psi)} psi pairs already in {cache}", flush=True)

    todo = [n for n in range(0, N + 1) if n not in psi]   # symmetry: negatives free
    print(f"psi phase: {len(todo)} quadratures left of {N+1} independent "
          f"(instead of {2*N+1}) on {workers} workers", flush=True)
    if todo:
        with Pool(workers, initializer=_init, initargs=(c, T, dps)) as pool, \
             open(cache, "a", encoding="utf-8") as fh:
            done = 0
            for n, a, b in pool.imap_unordered(_psi_pair, todo, chunksize=1):
                psi[n], psid[n] = mp.mpf(a), mp.mpf(b)
                fh.write(json.dumps({"n": n, "psi": a, "psid": b}) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
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


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_matrix(path):
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    mp.mp.dps = d["dps"]
    N, DIM = d["N"], 2 * d["N"] + 1
    Q = mp.matrix(DIM, DIM)
    for i, row in enumerate(d["Q"]):
        for j, v in enumerate(row):
            Q[i, j] = mp.mpf(v)
    return d, Q


def odd_projector(N):
    """V_odd: columns (e_k - e_{-k})/sqrt2, k = 1..N.  Dimension N, not N+1.

    The real projected combination differs from phi_{2k} by the overall phase i,
    which leaves the Gram matrix and its zeros unchanged — so F built from these
    real coefficients is the same function the odd-sector work already scans.
    """
    D, r = 2 * N + 1, 1 / mp.sqrt(2)
    V = mp.zeros(D, N)
    for k in range(1, N + 1):
        V[N + k, k - 1] = r
        V[N - k, k - 1] = -r
    return V


def even_projector(N):
    """V_even: columns e_0 and (e_k + e_{-k})/sqrt2, k = 1..N."""
    D, r = 2 * N + 1, 1 / mp.sqrt(2)
    V = mp.zeros(D, N + 1)
    V[N, 0] = 1
    for k in range(1, N + 1):
        V[N + k, k] = r
        V[N - k, k] = r
    return V


def lift(V, col, D):
    full = V * col
    nrm = mp.sqrt(sum(full[j] ** 2 for j in range(D)))
    v = mp.zeros(D, 1)
    for j in range(D):
        v[j, 0] = full[j] / nrm
    return v


def make_F(coef, N, L):
    """The reconstruction of extract_zeros, on a full coefficient vector."""
    def F(tau, dps):
        old = mp.mp.dps
        mp.mp.dps = dps
        try:
            t, Lm = mp.mpf(tau), mp.mpf(L)
            e = mp.exp(-1j * t * Lm)
            tot = mp.mpc(0, 0)
            for k in range(-N, N + 1):
                ck = coef[k + N]
                if ck == 0:
                    continue
                den = 2 * mp.pi * k / Lm - t
                tot += ck * (mp.mpc(Lm, 0) if abs(den) < mp.mpf(10) ** (-(dps - 10))
                             else (e - 1) / (1j * den))
            return mp.re(mp.exp(1j * t * Lm / 2) * tot / mp.sqrt(Lm))
        finally:
            mp.mp.dps = old
    return F


def bisection_plan(step, ref_dps):
    """Return the requested root tolerance and a non-binding iteration cap.

    The cap is derived from the initial bracket width, with ten safety
    halvings, so the tolerance — rather than a hard-coded loop count — is the
    stopping condition.
    """
    step = mp.mpf(step)
    root_tol = mp.mpf(10) ** (-(ref_dps - 20))
    ratio = step / root_tol
    required = 0 if ratio <= 1 else int(mp.ceil(mp.log(ratio, 2)))
    return root_tol, max(1, required + 10)


def unseeded_scan(F, lo, hi, step, scan_dps, ref_dps,
                  root_tol, max_bisections):
    """Sign changes on a uniform sweep, refined by bisection ONLY.

    No zeta ordinate is used anywhere in here; the caller compares afterwards."""
    brackets, x, fx = [], mp.mpf(lo), None
    fx = F(x, scan_dps)
    hi = mp.mpf(hi)
    while x < hi:
        y = min(x + step, hi)          # never sweep past the requested interval
        fy = F(y, scan_dps)
        if fx == 0:
            brackets.append((x, x))
        elif mp.sign(fx) != mp.sign(fy):
            brackets.append((x, y))
        x, fx = y, fy
    roots, rejected = [], []
    for a, b in brackets:
        lo_, hi_ = mp.mpf(a), mp.mpf(b)
        # RE-CONFIRM THE BRACKET AT THE WORKING PRECISION. The sign change was found at
        # scan_dps, which is deliberately low; if it was an artefact of that precision
        # the two ends may agree in sign here, and bisecting a non-bracket returns a
        # confident number that is not a root. Only the left end used to be re-evaluated.
        fl, fh = F(lo_, ref_dps), F(hi_, ref_dps)
        if fl != 0 and fh != 0 and mp.sign(fl) == mp.sign(fh):
            rejected.append(mp.nstr((lo_ + hi_) / 2, 8))
            continue
        for _ in range(max_bisections):
            mid = (lo_ + hi_) / 2
            fm = F(mid, ref_dps)
            if fm == 0:
                lo_ = hi_ = mid
                break
            if mp.sign(fm) != mp.sign(fl):
                hi_ = mid
            else:
                lo_, fl = mid, fm
            if hi_ - lo_ < root_tol:
                break
        roots.append((lo_ + hi_) / 2)
    if rejected:
        print(f"    {len(rejected)} bracket(s) rejected — the sign change did not "
              f"survive re-evaluation at dps {ref_dps}: {', '.join(rejected[:4])}"
              + (" ..." if len(rejected) > 4 else ""), flush=True)
    return roots




def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default=None,
                    help="path to a cached matrix; default is the name implied by "
                         "--c/--N/--T/--dps, which is also where a new build is written")
    # NOT required when --matrix is given: c and N are recorded IN the cache, and
    # asking for them again invites the worst kind of mistake — passing --N 200
    # alongside an N = 150 matrix and reading the result as the wrong cell.
    ap.add_argument("--c", type=int, default=None)
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--T", type=int, default=None,
                    help="read from --matrix when it exists; otherwise 400")
    ap.add_argument("--dps", type=int, default=None,
                    help="read from --matrix when it exists; otherwise 500")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--verify", action="store_true",
                    help="check a fresh build against the package's serial builder")
    ap.add_argument("--rebuild", action="store_true",
                    help="rebuild the matrix even if its cache exists (~55 min at c=100, N=150)")
    ap.add_argument("--rediagonalise", action="store_true",
                    help="discard the eigenvector cache and diagonalise again")
    ap.add_argument("--max-t", type=float, default=60.0)
    ap.add_argument("--min-t", type=float, default=0.1)
    ap.add_argument("--step", type=float, default=0.02)
    ap.add_argument("--scan-dps", type=int, default=60)
    ap.add_argument("--ref-dps", type=int, default=150)
    ap.add_argument("--random-seed", type=int, default=20260804)
    ap.add_argument("--sector", choices=["even", "odd"], default="even",
                    help="which parity sector to diagonalise and scan. The negative "
                         "block and everything published about it is the EVEN sector; "
                         "--sector odd asks whether the ordering is a property of the "
                         "whole near-null space or of the even sector only. Both use "
                         "the SAME cached matrix — no rebuild.")
    ap.add_argument("--n-positive", type=int, default=3,
                    help="how many of the smallest POSITIVE eigenvectors to scan. The "
                         "monotonicity claim is that the spurious-root count depends on "
                         "|lambda| alone with the sign playing no role; with only one "
                         "positive per cell that claim rests on a single interleaved "
                         "point. Three puts positives throughout the ordering.")
    a = ap.parse_args()
    # ALL FOUR PARAMETERS COME FROM AN EXISTING MATRIX, not just c and N. Reading only
    # c and N meant `--matrix Q_..._T800_... --rebuild` would rebuild at the DEFAULT
    # T = 400 and write the result under a filename that still says T800 — a wrong
    # matrix under a right name, which no later check could catch.
    _exists = a.matrix is not None and os.path.exists(a.matrix)
    if _exists:
        _d = json.load(open(a.matrix, encoding="utf-8"))
        for _k in ("c", "N", "T", "dps"):
            _v = getattr(a, _k)
            if _v is not None and _v != _d[_k]:
                ap.error(f"--{_k} {_v} contradicts {a.matrix}, which holds "
                         f"{_k} = {_d[_k]}. Remove the flag or point at another file.")
            setattr(a, _k, _d[_k])
    else:
        # no matrix to read from, so the geometry must be given and the rest defaulted
        if a.c is None or a.N is None:
            ap.error(f"{a.matrix} does not exist, so --c and --N are required"
                     if a.matrix else "give either --matrix or both --c and --N")
        if a.T is None:
            a.T = 400
        if a.dps is None:
            a.dps = 500

    c, N, T, dps = a.c, a.N, a.T, a.dps
    # RAISE dps BEFORE CREATING L. mpmath fixes an mpf's precision AT CREATION, so
    # mp.log(c) evaluated at the default 15 digits stays a 15-digit number however
    # high dps goes afterwards. F depends on L through exp(i tau L/2) and 2 pi k / L,
    # so that caps every root at ~1e-16 — the same failure as passing math.log(c),
    # which is already on record from the odd-sector work. It also made the two code
    # paths disagree: the fresh build kept the low-precision L while a rerun from the
    # cache recomputed it after load_matrix had raised dps.
    mp.mp.dps = dps
    matrix = a.matrix or f"Q_c{c}_N{N}_T{T}_dps{dps}.json"
    DIM, L = 2 * N + 1, mp.log(mp.mpf(c))

    # ---- STAGE 1: the matrix. Built only if its cache is absent. -------------
    if os.path.exists(matrix) and not a.rebuild:
        log(f"reusing {matrix}")
        d, Q = load_matrix(matrix)
        c, N, T, dps = d["c"], d["N"], d["T"], d["dps"]
        mp.mp.dps = dps
        DIM, L = 2 * N + 1, mp.log(mp.mpf(c))
    else:
        log(f"building c={c} N={N} T={T} dps={dps} on {a.workers} workers")
        print(f"python-flint: {HAS_FLINT}; cores seen: {os.cpu_count()}", flush=True)
        Q = build(c, N, T, dps, a.workers)
        if a.verify:
            print("\nverifying against the package's own serial builder ...", flush=True)
            R = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
            worst = max(abs(Q[i, j] - R[i, j]) for i in range(Q.rows) for j in range(Q.rows))
            scale = max(abs(R[i, j]) for i in range(R.rows) for j in range(R.rows))
            print(f"max |parallel - serial| = {mp.nstr(worst, 6)}   "
                  f"relative to max|Q| = {mp.nstr(worst/scale, 6)}")
            if worst / scale > mp.mpf(10) ** (-(dps - 10)):
                print("\nFAILED: the parallel build does not reproduce the serial one.")
                return 1
            print("the two builds agree at working precision")
        mp.mp.dps = dps
        with open(matrix, "w") as fh:
            json.dump({"c": c, "N": N, "T": T, "dps": dps,
                       "Q": [[mp.nstr(Q[i, j], dps + 5) for j in range(Q.rows)]
                             for i in range(Q.rows)]}, fh)
        log(f"matrix cached to {matrix} ({os.path.getsize(matrix)/1e6:.0f} MB)")

    # ---- STAGE 2: the eigenvectors. Diagonalised only if absent. -------------
    # THE SECTOR GOES IN BOTH CACHE NAMES. Without it an odd run overwrites the even
    # eigendecomposition and the even scan results, and the two are not interchangeable.
    # A REBUILT MATRIX INVALIDATES EVERYTHING DOWNSTREAM. Without this, --rebuild
    # produces a new matrix and then scans it with the PREVIOUS matrix's eigenvectors,
    # which is worse than either result alone.
    if a.rebuild:
        for _sec in ("even", "odd"):
            for _suf in ("_eig.json", "_branch_scan.json"):
                _stale = matrix.replace(".json", f"_{_sec}{_suf}")
                if os.path.exists(_stale):
                    os.remove(_stale)
                    log(f"--rebuild: removed the now-stale {_stale}")
    eig_cache = matrix.replace(".json", f"_{a.sector}_eig.json")
    scan_out = matrix.replace(".json", f"_{a.sector}_branch_scan.json")
    if a.rediagonalise and os.path.exists(eig_cache):
        os.remove(eig_cache)
    V = even_projector(N) if a.sector == "even" else odd_projector(N)
    # THE CACHE MUST NAME THE MATRIX IT CAME FROM. --rebuild clears it from inside this
    # script, but a matrix replaced from outside — a different build under the same
    # filename — would leave a stale eigendecomposition that nothing here could detect.
    with open(matrix, "rb") as _fh:
        matrix_sha = hashlib.sha256(_fh.read()).hexdigest()
    if os.path.exists(eig_cache):
        ec = json.load(open(eig_cache, encoding="utf-8"))
        _was = ec.get("matrix_sha256")
        if _was is not None and _was != matrix_sha:
            print(f"\nFAILED: {eig_cache} was diagonalised from a DIFFERENT matrix "
                  f"(sha {_was[:12]}) than the one now at {matrix} (sha "
                  f"{matrix_sha[:12]}).\n        Rerun with --rediagonalise.")
            return 1
        if _was is None:
            # A cache with no sha was written by an earlier version of this script and
            # therefore carries none of the checks added since — no ||Qv - lambda v||,
            # no parity figure, no matrix identity. It is reused rather than refused,
            # because refusing would invalidate finished work, but silence would let a
            # run inherit a predecessor's blind spots without saying so.
            log(f"NOTE: {eig_cache} predates the matrix-identity and residual checks. "
                f"It is being reused as-is; pass --rediagonalise (about a minute here) "
                f"to regenerate it under the current checks.")
        log(f"reusing {eig_cache}")
        E = [mp.mpf(x) for x in ec["eigenvalues"]]
        vecs = {int(k): [mp.mpf(x) for x in v] for k, v in ec["vectors"].items()}
    else:
        log(f"projecting onto the {a.sector} sector "
            f"({V.cols}-dimensional) and diagonalising (this is the slow part)")
        t0 = time.time()
        B = V.T * Q * V
        Ev, U = mp.eigsy(B)
        E = [Ev[j] for j in range(B.rows)]
        log(f"diagonalised in {(time.time()-t0)/60:.1f} min")
        order = sorted(range(len(E)), key=lambda j: E[j])       # ALGEBRAIC order
        keep = [j for j in order if E[j] < 0]                    # the negative block
        pos = sorted((j for j in order if E[j] > 0), key=lambda j: E[j])
        keep += pos[:max(1, a.n_positive)]        # the smallest few POSITIVE ones
        vecs = {j: [x for x in lift(V, mp.matrix([U[i, j] for i in range(B.rows)]), DIM)]
                for j in keep}
        # THE CHECKS THE WRITE-UP CLAIMS — computed here rather than asserted in prose.
        # For each kept vector: the eigenpair residual, the norm, and the parity of the
        # lifted vector (even: v_{-n} = v_n; odd: v_{-n} = -v_n and v_0 = 0).
        sgn = 1 if a.sector == "even" else -1
        checks = {}
        for j in keep:
            # AGAINST Q ON THE LIFTED VECTOR, not against B on the projected one.
            # ||Bu - lambda u|| only certifies the diagonalisation of B; ||Qv - lambda v||
            # certifies that AND that the sector is genuinely invariant under Q, which
            # is the assumption the whole parity split rests on.
            vcol = mp.matrix(vecs[j])
            res = Q * vcol - E[j] * vcol
            rn = mp.sqrt(sum(res[i] ** 2 for i in range(DIM)))
            v = vecs[j]
            # NOT an independent test: v was built by multiplying by V, so this is
            # necessarily near zero. It is a serialisation check on the cache, and it
            # is labelled as one rather than counted as evidence.
            par = max(abs(v[N + k] - sgn * v[N - k]) for k in range(1, N + 1))
            if a.sector == "odd":
                par = max(par, abs(v[N]))
            checks[str(j)] = {"residual_vs_Q": mp.nstr(rn, 6),
                              "norm": mp.nstr(mp.sqrt(sum(x ** 2 for x in v)), 12),
                              "parity_defect_serialisation_only": mp.nstr(par, 6)}
        # THE TWO CHECKS THE RESIDUAL DOES NOT COVER.
        # ||Qv - lambda v|| verifies the retained eigenpairs in the full space, but it
        # does not measure the coupling BETWEEN the sectors. Full parity invariance is
        # ||P_even Q P_odd||; and the symmetry of Q is asserted everywhere and had never
        # been measured. Both are cheap and both belong in the record rather than in a
        # sentence hedging what the residual does not show.
        Ve, Vo = even_projector(N), odd_projector(N)
        cross = Ve.T * Q * Vo
        cross_norm = mp.sqrt(sum(cross[i, j] ** 2
                                 for i in range(cross.rows) for j in range(cross.cols)))
        qnorm = mp.sqrt(sum(Q[i, j] ** 2 for i in range(DIM) for j in range(DIM)))
        asym = max(abs(Q[i, j] - Q[j, i]) for i in range(DIM) for j in range(i + 1, DIM))
        log(f"parity coupling ||P_even Q P_odd||_F = {mp.nstr(cross_norm, 4)} "
            f"({mp.nstr(cross_norm / qnorm, 4)} relative to ||Q||_F); "
            f"max |Q_ij - Q_ji| = {mp.nstr(asym, 4)}")
        checks["_global"] = {"cross_sector_norm": mp.nstr(cross_norm, 8),
                             "cross_sector_relative": mp.nstr(cross_norm / qnorm, 8),
                             "max_asymmetry": mp.nstr(asym, 8)}
        worst_res = max(mp.mpf(x["residual_vs_Q"]) for x in checks.values()
                        if "residual_vs_Q" in x)
        worst_par = max(mp.mpf(x["parity_defect_serialisation_only"])
                        for x in checks.values()
                        if "parity_defect_serialisation_only" in x)
        log(f"||Qv - lambda v|| <= {mp.nstr(worst_res, 4)} over {len(keep)} vectors "
            f"({a.sector} sector); parity defect <= {mp.nstr(worst_par, 4)} "
            f"(serialisation check, not independent)")
        json.dump({"matrix_sha256": matrix_sha,
                   "eigenvalues": [mp.nstr(x, dps + 5) for x in E],
                   "checks": checks,
                   "vectors": {str(k): [mp.nstr(x, dps + 5) for x in v]
                               for k, v in vecs.items()}},
                  open(eig_cache, "w"))
        log(f"cached eigendecomposition to {eig_cache}")

    neg = [j for j in range(len(E)) if E[j] < 0]
    pos_min = min((j for j in range(len(E)) if E[j] > 0), key=lambda j: E[j])
    # THE CACHE MUST SATISFY THE REQUEST OR THE RUN IS NOT THE RUN THAT WAS ASKED FOR.
    # A warning is not enough: the printed table would silently cover fewer branches
    # than the command specified, and the operator would have no reason to notice.
    # COMPARE THE KEY SET, NOT THE COUNT. Equal counts can still be different vectors
    # — a cache holding one extra negative and one fewer positive would pass a count
    # test and scan a set the command never asked for.
    _pos_sorted = sorted((j for j in range(len(E)) if E[j] > 0), key=lambda j: E[j])
    want_keys = set(neg) | set(_pos_sorted[:max(1, a.n_positive)])
    have_keys = set(int(k) for k in vecs)
    if have_keys != want_keys:
        print(f"\nFAILED: {eig_cache} holds {len(have_keys)} vectors but this command "
              f"needs {len(want_keys)} ({len(neg)} negative + "
              f"{max(1, a.n_positive)} positive), and the sets differ by "
              f"{len(want_keys ^ have_keys)} branch(es).\n"
              f"        Rerun with --rediagonalise, or with --n-positive "
              f"{len(have_keys) - len(neg)} to match what is cached.")
        return 1
    print(f"\n{a.sector} sector at c={c}, N={N}, T={T}, dps={dps}")
    print(f"  eigenvalues            : {len(E)}")
    print(f"  with a NEGATIVE sign   : {len(neg)}")
    for j in sorted(neg, key=lambda j: E[j]):
        print(f"      lambda = {mp.nstr(E[j], 8)}")
    print(f"  smallest POSITIVE      : {mp.nstr(E[pos_min], 8)}")

    # THE PRECISION-FLOOR CHECK, from an external criterion rather than a guess.
    # Andrews (Zenodo 10.5281/zenodo.21272569) documents an artifact in the same
    # family of computations: once |log10 lambda| approaches the working digits, the
    # smallest eigenvector STOPS LOOKING EVEN and the natural and forced-even
    # eigenvalues appear to split by ~24% -- and the whole effect collapses when the
    # precision is raised. He finds a ratio dps/|log10 lambda| near 0.79 broken and
    # near 1.6 safe. Any negative eigenvalue in that regime would be indistinguishable
    # from that artifact, so the ratio is printed here rather than assumed.
    _deep = min((abs(E[j]) for j in range(len(E)) if E[j] != 0), default=None)
    if _deep is not None and _deep > 0:
        _ratio = mp.mpf(dps) / abs(mp.log(_deep, 10))
        _verdict = ("above the ~1.6 Andrews finds safe" if _ratio >= mp.mpf("1.6")
                    else "*** BELOW 1.6 -- rerun at higher dps before reading the "
                         "sign of anything here ***")
        print(f"  precision-floor ratio  : {mp.nstr(_ratio, 4)} "
              f"(dps {dps} against the deepest |lambda|) — {_verdict}")

    # the calibration ruler, measured in this same basis
    random.seed(a.random_seed)
    # V.cols, not N+1: the odd projector has N columns and the even one N+1, and a
    # hardcoded N+1 makes the random control the one thing in the run that cannot be
    # built in the odd sector.
    rnd = mp.matrix([mp.mpf(random.gauss(0, 1)) for _ in range(V.cols)])
    vecs = dict(vecs)
    vecs["random"] = [x for x in lift(V, rnd, DIM)]

    mp.mp.dps = a.ref_dps + 30
    gam = [mp.im(mp.zetazero(k)) for k in range(1, 60)]
    below = [g for g in gam if g < a.max_t]

    # THE RESOLUTION FLOOR.  The iteration cap is derived from the bracket width
    # and carries ten safety halvings, so it cannot bind before ROOT_TOL under the
    # ordinary bisection path.  The floor is therefore the requested break
    # tolerance, not step/2^a hard-coded iteration count.
    ROOT_TOL, MAX_BISECTIONS = bisection_plan(str(a.step), a.ref_dps)
    FLOOR = ROOT_TOL
    FLOOR_TOL = mp.mpf("1e-3")
    print(f"\nresolution floor: {mp.nstr(FLOOR, 4)} — set by the break tolerance; "
          f"the derived cap is {MAX_BISECTIONS} halvings for a {a.step}-wide bracket "
          f"(10 safety halvings). Any approach at or below this is saturated, not measured.")

    print(f"\nunseeded scan of ({a.min_t}, {a.max_t}) at step {a.step} — "
          f"roots are frozen before any zeta comparison")
    print(f"{'vector':<24} {'roots':>5} {'err min':>13} {'median':>13} {'max':>13}   reading")
    out, fingerprints, tight = {}, {}, []
    # SORTED BY |lambda|, which is the ordering the monotonicity claim is about. Printing
    # them in eigenvalue order instead makes a monotone sequence look scattered.
    for key in sorted(vecs, key=lambda k: mp.inf if k == "random" else abs(E[int(k)])):
        label = ("random control" if key == "random" else
                 f"negative branch lambda={mp.nstr(E[int(key)],4)}" if int(key) in neg else
                 f"positive branch lambda={mp.nstr(E[int(key)],4)}")
        F = make_F(vecs[key], N, L)
        t0 = time.time()
        roots = unseeded_scan(F, a.min_t, a.max_t, mp.mpf(str(a.step)),
                              a.scan_dps, a.ref_dps, ROOT_TOL, MAX_BISECTIONS)
        best = min((min(abs(r - g) for g in gam) for r in roots), default=None)
        out[str(key)] = {"n_roots": len(roots),
                         "lambda": (mp.nstr(E[int(key)], 10) if key != "random" else None),
                         # FULL precision, not 40 digits. At 40 the stored strings for
                         # different branches come out identical and every per-branch
                         # accuracy comparison silently collapses to the truncation error.
                         "roots": [mp.nstr(r, a.ref_dps) for r in roots],
                         "closest_approach": mp.nstr(best, 8) if best is not None else None}
        # COMPLETENESS is the test, not the best of many shots. Match each root to
        # its nearest ordinate and count how many ordinates are hit and how many
        # roots are left over. A vector with 41 roots and 13 targets is not doing
        # what a vector with 13 roots and 13 targets is doing, however close its
        # single best root happens to fall.
        # ONE ROOT PER ORDINATE. Assigning every root to its nearest ordinate lets two
        # roots near the SAME zero both count as matched while the set records the zero
        # once, so a branch with 14 roots could print "13/13, 0 extras". Each ordinate
        # takes its closest root and no other; every remaining root is an extra.
        claims = {}
        for i, r in enumerate(roots):
            k = min(range(len(gam)), key=lambda j: abs(r - gam[j]))
            if gam[k] < a.max_t and abs(r - gam[k]) < FLOOR_TOL:
                if k not in claims or abs(r - gam[k]) < abs(roots[claims[k]] - gam[k]):
                    claims[k] = i
        matched = set(claims)
        extra = len(roots) - len(claims)
        # PER-ORDINATE, sorted. `best` alone is the minimum over thirteen and always
        # flatters the branch; the max is the honest figure and the two differ here by
        # up to four orders of magnitude.
        # THE SCAN'S OWN RESOLUTION. A uniform sweep cannot see two roots closer
        # together than `step`: they merge into one bracket, or the two sign changes
        # cancel and both vanish. So the extras count has a CEILING set by the sweep,
        # and a branch whose comb is dense enough will be UNDERCOUNTED without any
        # sign of it. Report the tightest adjacent gap so the ceiling is visible.
        srt = sorted(roots)
        gap = min((srt[i + 1] - srt[i] for i in range(len(srt) - 1)), default=None)
        per = sorted(min(abs(r - g) for r in roots) for g in below) if roots else []
        med = per[len(per) // 2] if per else None
        worst = per[-1] if per else None
        at_floor = best is not None and best < FLOOR
        out[str(key)].update({"ordinates_matched": len(matched),
                              "ordinates_available": len(below),
                              "unmatched_roots": extra,
                              "err_min": mp.nstr(best, 8) if best is not None else None,
                              "err_median": mp.nstr(med, 8) if med is not None else None,
                              "err_max": mp.nstr(worst, 8) if worst is not None else None,
                              "per_ordinate": [mp.nstr(e, 8) for e in per],
                              "min_root_gap": mp.nstr(gap, 6) if gap is not None else None,
                              "at_resolution_floor": bool(at_floor)})
        if best is None:
            reading = "no roots"
        elif best > mp.mpf("1e-3"):
            reading = "random-like"
        elif len(matched) == len(below) and extra == 0:
            reading = "matches every ordinate, nothing left over"
        else:
            reading = (f"{len(matched)}/{len(below)} ordinates, {extra} roots unmatched "
                       f"— NOT a clean reconstruction")
        if gap is not None and gap < 5 * mp.mpf(str(a.step)):
            reading += (f"  [ROOTS {mp.nstr(gap,3)} APART, only "
                        f"{float(gap)/a.step:.1f}x the step — the count may be an "
                        f"UNDERCOUNT; rerun with a smaller --step]")
            tight.append(str(key))
        if at_floor:
            reading += "  [SATURATED at the resolution floor, not a measurement]"
        # If two branches produce IDENTICAL root lists, their accuracies are not being
        # measured — the solver has saturated and both are reporting the same point.
        # lambda differing by forty decades with identical roots is the tell.
        fp = tuple(mp.nstr(r, 40) for r in roots)
        for other, ofp in fingerprints.items():
            if ofp == fp:
                reading += f"  [ROOT LIST IDENTICAL TO {other} — saturated, not resolved]"
                break
        fingerprints[str(key)] = fp
        print(f"{label:<24} {len(roots):>5} "
              f"{mp.nstr(best,5) if best is not None else '-':>13} "
              f"{mp.nstr(med,5) if med is not None else '-':>13} "
              f"{mp.nstr(worst,5) if worst is not None else '-':>13}   {reading}")
    # PROVENANCE. Without it two scan files are indistinguishable and every step or
    # precision comparison rests on remembering which file came from which command.
    json.dump({"provenance": {"sector": a.sector, "c": c, "N": N, "T": T, "dps": dps,
                              "step": a.step, "scan_dps": a.scan_dps,
                              "ref_dps": a.ref_dps, "min_t": a.min_t, "max_t": a.max_t,
                              "match_tol": mp.nstr(FLOOR_TOL, 3),
                              "n_positive": a.n_positive,
                              "random_seed": a.random_seed,
                              "matrix": os.path.basename(matrix),
                              "eig_cache": os.path.basename(eig_cache)},
               "branches": out},
              open(scan_out, "w"))
    if tight:
        print(f"\n*** {len(tight)} branch(es) have roots within five steps of each other. "
              f"The sweep cannot resolve a pair closer than {a.step}, so those extras "
              f"counts are lower bounds. Rerun with --step {a.step/4:g} and check they "
              f"do not move before reading anything into the ordering. ***")
    print(f"\nzeta ordinates below {a.max_t}: {len(below)}")
    print("root lists written to " + scan_out)
    print("\nThis is a measurement of where each branch's roots fall. It is not a")
    print("certification of the negative block in either direction.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
