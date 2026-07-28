# Prime Number Studies

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21638887.svg)](https://doi.org/10.5281/zenodo.21638887)

Notes, measurements and verified code from an independent, self-directed study of the
multiplicative structure of the integers — starting from an elementary object, the
**division table**, and following it as far as it goes.

---

## What this is, and what it is not

**This repository contains independent re-derivations of known mathematics.** Every
headline result here already exists in the literature, and each paper says so explicitly
and names the source. Specifically:

| what the papers reach | where it already lives |
|---|---|
| the map `Φ(j,t) = 2jt + j + t − 1` characterising odd composites | the **sieve of Sundaram** (1934) |
| the wheel pair-count `W(H)` and its normalization | the **Hardy–Littlewood** singular series (1923), truncated |
| peaks of a smoothed prime signal landing on the zeta ordinates | the **explicit formula** — the object is a truncated `ξ'/ξ` |
| reconstruction of `τ, φ, μ, Λ` from valuation vectors | the fundamental theorem of arithmetic |
| lag functions of `gcd(n,h)` and their transforms | **Ramanujan–Fourier** analysis of `r`-even functions |
| the positive `log gcd` Gram kernel | the classical **gcd-matrix** literature |
| `λ_min → 0` for a purely atomic spectral measure | the **Szegő limit theorem** (Grenander–Szegő, 1958) |
| the critical balance of the Weil form at `β* = 1` | **Bombieri**, Rend. Lincei 11 (2000) |

**There is no claim here about the Riemann Hypothesis.** Where a construction touches it,
the paper states precisely why the numerics certify nothing about the critical line.

What the repository *does* offer:

1. **A route.** An elementary division table, followed honestly, arrives at the explicit
   formula and at the truncated Weil quadratic form. Seeing that path laid out end to end
   is the main educational value here.
2. **Negative results and controls**, stated plainly rather than buried — including
   several that overturn earlier readings of the same data.
3. **Runnable code for every number printed.** No table appears in a paper without a
   script in `code/` that regenerates it.

---

## Contents

```
papers/     four self-contained papers, plus one short structural note (Markdown)
figures/    the figures they reference
code/       verification scripts for every table and number, plus the figure generators
```

### `papers/division_table_prime_towers.md`
From the odd division table to prime-power towers. Builds the remainder field
`R(H,j) = 2(H−j+1) mod (2j+1)`, shows every column is one arithmetic progression, derives
the factor map `Φ`, identifies it as Sundaram's sieve, proves the preimage multiplicity is
`τ(n) − 2`, and reduces the whole table to intersections of prime-power layers. Includes a
from-scratch explanation of the division table for readers meeting it for the first time.
Deliberately stops before the zeta function.

*Also records that the compression over trial division is a factor ≈ 3, asymptotically only
`(log N)/2` — i.e. the prime number theorem, and nothing better.*

### `papers/prime_gap_wheel_study.md`
Counting residue pairs in the sieve cycle `W = 30030`, and the frequency of prime gaps.
The pair count is an exact theorem with no free parameters. The empirical half is a worked
example of a fitted parameter that looks real and is not: a stretched exponent `β ≈ 1.07`
fits better *and* transfers better out of sample, yet moves about seven times more with the
arbitrary fitting window than with six decades of height, and even its drift direction
depends on that window. What it absorbs is identified and measured.

### `papers/zeta_dynamical_model.md`
A dynamical model `D_X(t) = A(t) − P_X(t)` built from divisibility via Möbius inversion.
Its peaks match the zeta ordinates (29/29 in `10 < t < 100`, reproduced here with the full
protocol). The paper then shows why that match is structurally forced — `D_X` is a
truncated `ξ'/ξ`, so its peaks are its poles — and that two features which looked like
independent structure resolve into truncation error and the Rayleigh resolution criterion.

### `papers/division_table_single_column.md`
The odd division table rewritten as **one time-ordered series**, with nothing lost. Each
odd prime is a clock starting at `p²` and stepping by `2p`, so a prime is the *absence* of
a strike; each integer becomes a vector of its prime-power layers whose inner product is
exactly `log gcd`; and the table's columns reappear as lag diagonals of the single column.
Time averages give a Dirichlet convolution of the von Mangoldt function, Möbius inversion
isolates `Λ`, and the Dirichlet transform gives the prime side of `−ζ'/ζ`. The chain is
exact at every step — and stops, provably, at its half-plane of convergence, which is where
the construction stops being a re-derivation and starts being a citation of `ζ`.

Three spectral objects built from it are then tested, and all three are closed for
**structural** reasons rather than numerical ones. The covariance kernel is positive
because it is a stationary autocovariance (*not* because it is a Gram matrix — the
pointwise Gram matrix is not even Toeplitz), and its decay to zero is forced by a purely
atomic spectral measure. Isolating `Λ` gives a Toeplitz family whose positivity is a
**convention artifact**: over odd prime powers the symbol satisfies `f(π) = −c < 0` for
*any* nonnegative weights, while admitting the prime 2 restores positivity but only inside
a bounded window of the parameter. And the hybrid prime-power operator is bipartite, so its
spectrum is exactly symmetric about zero and its top does not move when primes are added,
while the zeta ordinates are positive and unbounded.

*One open item is flagged as possibly having content beyond bookkeeping: the minimising
angle of the symbol locks onto rational points — exactly `π`, then exactly `2π/5` to eight
decimals — held there by a downward cusp whose exponent tracks the parameter, and comes
loose only when that cusp becomes Lipschitz.*

### `papers/prime_edge_two_paths.md`
A short note, not a paper of its own. It measures the prime-edge derivative jump of the
truncated Weil form in a Dirichlet sine construction and compares it with the matrix-valued
von Mangoldt measure of Groskin (2026), where the jump is rank one. The two differ in where
the edge vanishing lives: inside the matrix in the reference path, outside it as a scalar
window here — so they are distinct families rather than the same object in different
coordinates. The note carries a closed form for the measured plateau, `(m-1)/(m+2)`, due to
A. Groskin in correspondence and verified here to seven digits, and a corrections section
recording what earlier versions of the note got wrong.

---

## Reproducibility

Every table and figure is regenerated by a script in `code/`.

```bash
pip install numpy sympy mpmath scipy matplotlib

python code/verify_paper2.py          # division table → prime towers
python code/verify_paper3.py          # wheel pair count and gap frequencies
python code/robust_paper3.py          # the window sweeps that decide the β question
python code/w_fast2.py                # D_X vs ξ'/ξ, peak shape
python code/w_stab2.py                # peak extraction, prominence and tolerance sweeps
python code/w_fast3.py                # angular criterion and the resolution law
python code/w_alg.py                  # divisibility matrix and operator algebra
python code/verify_core.py            # single column: Propositions 1–7, the operator tables, controls
python code/verify_covariance.py      # single column: the spectral measure, decay, and rank saturation
python code/verify_positivity.py      # single column: both conventions, the α-window, N versus P
python code/prime_edge_jump.py        # prime-edge jump; section C retained with its artifact annotated
python code/prime_edge_rank.py        # rank, edge analysis, and factor universality
python code/check_reply.py            # closed form for the plateau, and the split blocks

cd code
python fig.py       # figure for the division-table paper
python robust_paper3.py && python fig3.py   # fig3 needs the gap data robust_paper3 caches
python fig4.py      # figure for the dynamical-model paper
cd ..
```

Tested with Python 3.12, NumPy 2.4, SymPy 1.14, mpmath 1.3.
The gap-frequency scripts segment-sieve windows of width `4×10⁷` up to `10¹³` and take a
few minutes. `verify_covariance.py` and `verify_positivity.py` also take a few minutes —
each has a long-average or large-matrix step, and `verify_covariance.py` has a `FAST` flag
that shortens the slowest one. Everything else runs in seconds.

---

## Three method rules these papers exist to illustrate

All three were learned by getting them wrong first.

**Test any threshold law on the bare kernel before calling it arithmetic.** A resolution
law `log X_crit ≍ C/g` looked like a property of the zeta zeros. Two bare Fejér bumps at
separation `g` give `L_crit · g = 5.21233`, constant to `0.0000%` over a fifteenfold range
of `g`. It was the Rayleigh criterion the whole time.

**Match the number of transported quantities before comparing models out of sample.** A
two-parameter model beating a zero-parameter model across a six-decade extrapolation says
almost nothing until the zero-parameter model is given one calibrated constant too. Doing
that closed most of the gap.

**Print the truncation next to every number that depends on it.** The same covariance
kernel, at the same parameter, has an apparent decay exponent of `−4.8`, `−2.2`, `−1.2` or
`−0.74` depending only on where the layers are cut off. Two exponents in this literature
that look like different phenomena are the same object at two truncations — and at any
*fixed* cutoff there is no exponent at all, because the spectral measure has finitely many
atoms and the rank simply saturates.

---

## Status

Four papers and one note are here.

The fourth paper closes the division-table route as a direct approach to RH, and says so in
its own conclusion. What it leaves genuinely open is listed in its Appendix A — most
interestingly the rational-angle locking of the symbol's minimum, which is the one measured
structure in it that does not obviously reduce to bookkeeping.

Corrections, counterexamples and pointers to prior art are all welcome; open an issue.
Being told that something here is already known is a useful outcome, not an unwelcome one.

---

## License

- Text, papers and figures: **CC BY 4.0** — see `LICENSE-CONTENT`.
- Code: **MIT** — see `LICENSE`.

## Citation

```
Osman, M. (2026). Prime Number Studies: independent re-derivations, measurements
and negative results around the division table. Zenodo.
https://doi.org/10.5281/zenodo.21638887
```

Repository: https://github.com/Osman209/prime-number-studies
ORCID: https://orcid.org/0009-0004-5912-999X
