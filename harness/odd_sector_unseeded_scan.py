#!/usr/bin/env python3
"""Long-running, resumable unseeded zero scan for the odd CvS sector.

The discovery stage never calls ``mp.zetazero``.  It builds the CvS matrix,
extracts the lowest odd-sector eigenvector, scans a user-selected interval for
sign changes of its centred Fourier transform, and refines every bracket by
plain bisection.  Only after that root list is frozen does the optional
comparison stage evaluate known zeta ordinates.

Install:
    python -m pip install connes-cvs python-flint

Recommended overnight run:
    python odd_sector_unseeded_scan.py \
        --Ns 20 30 40 60 80 100 --c 13 --T 400 --dps 120 \
        --min-t 0.1 --max-t 50 --step 0.02 \
        --output odd_scan_c13 2>&1 | tee odd_scan_c13.log

Resume after interruption by running the same command again.  Completed cases
and cached eigenvectors are reused unless ``--force`` is supplied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
import time
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Callable

import mpmath as mp

try:
    from connes_cvs import build_galerkin_matrix
except ImportError as exc:
    raise SystemExit(
        "Missing dependency. Run: python -m pip install connes-cvs python-flint"
    ) from exc


@dataclass
class RootRecord:
    index: int
    root: str
    residual: str
    bracket_left: str
    bracket_right: str


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(temporary, path)


def odd_ground_state(Q: mp.matrix) -> tuple[mp.mpf, mp.matrix]:
    """Return the lowest eigenpair of the odd parity block."""
    dim = Q.rows
    if Q.cols != dim or dim % 2 != 1:
        raise ValueError("Q must be a square (2N+1)-dimensional matrix")
    N = (dim - 1) // 2
    projector = mp.matrix(dim, N)
    inv_sqrt_two = 1 / mp.sqrt(2)
    for k in range(1, N + 1):
        projector[N + k, k - 1] = inv_sqrt_two
        projector[N - k, k - 1] = -inv_sqrt_two

    odd_block = projector.T * Q * projector
    eigenvalues, eigenvectors = mp.eigsy(odd_block)
    eigenvalue = eigenvalues[0]
    vector = projector * eigenvectors[:, 0]
    vector /= mp.norm(vector)

    parity_error = max(
        abs(vector[N + k] + vector[N - k]) for k in range(1, N + 1)
    )
    if abs(vector[N]) > mp.eps or parity_error > mp.sqrt(mp.eps):
        raise ArithmeticError(f"odd-parity lift failed: {parity_error=}")
    return eigenvalue, vector


def save_vector_cache(
    path: Path,
    *,
    c: int,
    N: int,
    T: int,
    dps: int,
    eigenvalue: mp.mpf,
    vector: mp.matrix,
) -> None:
    atomic_json(
        path,
        {
            "schema": 1,
            "c": c,
            "N": N,
            "T": T,
            "dps": dps,
            "eigenvalue": mp.nstr(eigenvalue, dps),
            "vector": [mp.nstr(vector[i], dps) for i in range(vector.rows)],
        },
    )


def load_vector_cache(
    path: Path, *, c: int, N: int, T: int, dps: int
) -> tuple[mp.mpf, mp.matrix] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    expected = {"schema": 1, "c": c, "N": N, "T": T, "dps": dps}
    if any(data.get(key) != value for key, value in expected.items()):
        return None
    values = data.get("vector", [])
    if len(values) != 2 * N + 1:
        return None
    vector = mp.matrix([mp.mpf(value) for value in values])
    vector /= mp.norm(vector)
    return mp.mpf(data["eigenvalue"]), vector


def centred_transform(vector: mp.matrix, L: mp.mpf) -> Callable[[mp.mpf], mp.mpf]:
    """Construct the real centred Fourier transform used for root discovery."""
    N = (vector.rows - 1) // 2
    pi = mp.pi
    coefficients = [vector[k + N] for k in range(-N, N + 1)]

    def transform(tau: mp.mpf) -> mp.mpf:
        tau = mp.mpf(tau)
        exp_minus = mp.exp(-1j * tau * L)
        total = mp.mpc(0)
        for offset, coefficient in enumerate(coefficients):
            if not coefficient:
                continue
            k = offset - N
            denominator = 2 * pi * k / L - tau
            threshold = mp.power(10, -(mp.mp.dps - 10))
            if abs(denominator) < threshold:
                term = mp.mpc(L)
            else:
                term = (exp_minus - 1) / (1j * denominator)
            total += coefficient * term
        centred = mp.exp(1j * tau * L / 2) * total / mp.sqrt(L)
        return mp.re(centred)

    return transform


def refine_bracket(
    function: Callable[[mp.mpf], mp.mpf],
    left: mp.mpf,
    right: mp.mpf,
    f_left: mp.mpf,
    f_right: mp.mpf,
    digits: int,
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Bracket-preserving bisection, stopped by interval width."""
    if f_left == 0:
        return left, left, left
    if f_right == 0:
        return right, right, right
    if f_left * f_right > 0:
        raise ValueError("refine_bracket requires a sign-changing bracket")

    # log2(10) bisections per decimal digit, plus a small safety margin.
    iterations = math.ceil((digits + 8) * math.log2(10))
    a, b, fa = left, right, f_left
    for _ in range(iterations):
        midpoint = (a + b) / 2
        f_midpoint = function(midpoint)
        if fa * f_midpoint <= 0:
            b = midpoint
        else:
            a, fa = midpoint, f_midpoint
    return (a + b) / 2, a, b


def scan_sign_changing_roots(
    function: Callable[[mp.mpf], mp.mpf],
    minimum: mp.mpf,
    maximum: mp.mpf,
    step: mp.mpf,
    digits: int,
) -> list[RootRecord]:
    """Discover roots only from sign-changing adjacent grid values."""
    if not minimum < maximum or step <= 0:
        raise ValueError("Require minimum < maximum and step > 0")

    roots: list[RootRecord] = []
    left = minimum
    f_left = function(left)
    while left < maximum:
        right = min(left + step, maximum)
        f_right = function(right)
        candidate = None
        bracket_left = left
        bracket_right = right
        if f_left == 0:
            candidate = left
            bracket_right = left
        elif f_right == 0:
            candidate = right
            bracket_left = right
        elif f_left * f_right < 0:
            candidate, bracket_left, bracket_right = refine_bracket(
                function, left, right, f_left, f_right, digits
            )

        if candidate is not None:
            duplicate_tolerance = mp.power(10, -(digits // 2))
            if not roots or abs(candidate - mp.mpf(roots[-1].root)) > duplicate_tolerance:
                roots.append(
                    RootRecord(
                        index=len(roots) + 1,
                        root=mp.nstr(candidate, digits),
                        residual=mp.nstr(abs(function(candidate)), digits),
                        bracket_left=mp.nstr(bracket_left, digits),
                        bracket_right=mp.nstr(bracket_right, digits),
                    )
                )
        left, f_left = right, f_right
    return roots


def compare_after_discovery(roots: list[RootRecord], count: int) -> list[dict]:
    """Compare a frozen root list to zeta ordinates; never used for discovery."""
    if not roots or count <= 0:
        return []
    root_values = [mp.mpf(record.root) for record in roots]
    comparison = []
    for index in range(1, count + 1):
        gamma = mp.im(mp.zetazero(index))
        nearest = min(root_values, key=lambda root: abs(root - gamma))
        comparison.append(
            {
                "zeta_index": index,
                "gamma": mp.nstr(gamma, mp.mp.dps),
                "nearest_discovered_root": mp.nstr(nearest, mp.mp.dps),
                "absolute_error": mp.nstr(abs(nearest - gamma), mp.mp.dps),
            }
        )
    return comparison


def random_odd_vector(N: int, seed: int) -> mp.matrix:
    rng = random.Random(seed)
    vector = mp.matrix(2 * N + 1, 1)
    for k in range(1, N + 1):
        value = mp.mpf(str(rng.gauss(0, 1)))
        vector[N + k] = value
        vector[N - k] = -value
    vector /= mp.norm(vector)
    return vector


def write_summary_csv(path: Path, completed: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "N",
                "odd_lambda",
                "discovered_roots",
                "random_roots",
                "first_zeta_error",
                "elapsed_seconds",
            ]
        )
        for key in sorted(completed, key=int):
            result = completed[key]
            comparison = result.get("comparison", [])
            writer.writerow(
                [
                    key,
                    result["odd_lambda"],
                    len(result["roots"]),
                    len(result.get("random_control", {}).get("roots", [])),
                    comparison[0]["absolute_error"] if comparison else "",
                    result["elapsed_seconds"],
                ]
            )
    os.replace(temporary, path)


def package_version() -> str:
    try:
        return version("connes-cvs")
    except PackageNotFoundError:
        return "unknown"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resumable unseeded odd-sector CvS zero scan"
    )
    parser.add_argument("--Ns", nargs="+", type=int, default=[20, 30, 40, 60])
    parser.add_argument("--c", type=int, default=13)
    parser.add_argument("--T", type=int, default=400)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--min-t", type=str, default="0.1")
    parser.add_argument("--max-t", type=str, default="50")
    parser.add_argument("--step", type=str, default="0.02")
    parser.add_argument("--compare-zeros", type=int, default=10)
    parser.add_argument("--random-seed", type=int, default=20260804)
    parser.add_argument("--no-random-control", action="store_true")
    parser.add_argument("--no-compare", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("odd_scan_results"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.c < 2 or args.T < 1 or args.dps < 30 or any(N < 2 for N in args.Ns):
        raise SystemExit("Require c>=2, T>=1, dps>=30, and every N>=2")

    args.output.mkdir(parents=True, exist_ok=True)
    checkpoint_path = args.output / "checkpoint.json"
    summary_path = args.output / "summary.csv"
    roots_directory = args.output / "roots"
    cache_directory = args.output / "vector_cache"

    configuration = {
        "schema": 1,
        "c": args.c,
        "T": args.T,
        "dps": args.dps,
        "min_t": args.min_t,
        "max_t": args.max_t,
        "step": args.step,
        "compare_zeros": 0 if args.no_compare else args.compare_zeros,
        "random_control": not args.no_random_control,
        "random_seed": args.random_seed,
        "connes_cvs_version": package_version(),
        "python": sys.version,
    }

    completed: dict[str, dict] = {}
    if checkpoint_path.exists() and not args.force:
        with checkpoint_path.open("r", encoding="utf-8") as handle:
            previous = json.load(handle)
        if previous.get("configuration") == configuration:
            completed = previous.get("completed", {})
        else:
            log("Existing checkpoint has different settings; it will not be reused.")

    for N in args.Ns:
        key = str(N)
        if key in completed and not args.force:
            log(f"N={N}: already completed; skipping")
            continue

        started = time.time()
        mp.mp.dps = args.dps
        L = mp.log(args.c)
        cache_path = cache_directory / (
            f"c{args.c}_N{N}_T{args.T}_dps{args.dps}_odd_vector.json"
        )
        cached = None if args.force else load_vector_cache(
            cache_path, c=args.c, N=N, T=args.T, dps=args.dps
        )
        if cached is None:
            log(f"N={N}: building Q(c={args.c}, T={args.T}, dps={args.dps})")
            Q = build_galerkin_matrix(c=args.c, N=N, T=args.T, dps=args.dps)
            log(f"N={N}: diagonalising odd block")
            odd_lambda, odd_vector = odd_ground_state(Q)
            save_vector_cache(
                cache_path,
                c=args.c,
                N=N,
                T=args.T,
                dps=args.dps,
                eigenvalue=odd_lambda,
                vector=odd_vector,
            )
            del Q
            log(f"N={N}: eigenvector cached")
        else:
            odd_lambda, odd_vector = cached
            log(f"N={N}: loaded cached odd eigenvector")

        minimum = mp.mpf(args.min_t)
        maximum = mp.mpf(args.max_t)
        step = mp.mpf(args.step)
        transform = centred_transform(odd_vector, L)
        log(f"N={N}: unseeded scan on [{minimum}, {maximum}] with step {step}")
        roots = scan_sign_changing_roots(
            transform, minimum, maximum, step, args.dps
        )
        frozen_roots = [asdict(record) for record in roots]
        atomic_json(roots_directory / f"N{N}_discovered_roots.json", frozen_roots)
        log(f"N={N}: froze {len(roots)} discovered sign-changing roots")

        # This is deliberately after discovery and after the root list is saved.
        comparison = []
        if not args.no_compare:
            comparison = compare_after_discovery(roots, args.compare_zeros)
            log(f"N={N}: post-discovery zeta comparison complete")

        random_control = {}
        if not args.no_random_control:
            control_vector = random_odd_vector(N, args.random_seed + N)
            control_transform = centred_transform(control_vector, L)
            control_roots = scan_sign_changing_roots(
                control_transform, minimum, maximum, step, args.dps
            )
            random_control = {
                "seed": args.random_seed + N,
                "roots": [asdict(record) for record in control_roots],
                "comparison": (
                    compare_after_discovery(control_roots, args.compare_zeros)
                    if not args.no_compare
                    else []
                ),
            }
            log(f"N={N}: random odd control found {len(control_roots)} roots")

        result = {
            "N": N,
            "odd_lambda": mp.nstr(odd_lambda, args.dps),
            "roots": frozen_roots,
            "comparison": comparison,
            "random_control": random_control,
            "elapsed_seconds": round(time.time() - started, 3),
        }
        completed[key] = result
        atomic_json(
            checkpoint_path,
            {"configuration": configuration, "completed": completed},
        )
        write_summary_csv(summary_path, completed)
        log(
            f"N={N}: checkpoint saved; elapsed={result['elapsed_seconds']} s; "
            f"lambda_odd={mp.nstr(odd_lambda, 12)}"
        )

    log(f"All requested cases complete. Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
