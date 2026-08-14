# Prime Number Studies

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21638887.svg)](https://doi.org/10.5281/zenodo.21638887)

**Read the papers online — <https://osman209.github.io/prime-number-studies/>**

Every paper and note has a page there with an abstract and a typeset PDF.

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
| the shifted Möbius map, its derivative object, and the equality between them | the **Li–Sekatskii coefficients** — Sekatskii (2013, 2014), on **Bombieri–Lagarias** (1999) and **Li** (1997) |
| explicit indices at which a Li-type coefficient certifies an off-line zero | **Brown** (2005), **Bucur–Ernvall-Hytönen–Odžak–Smajlović** (2016), **Palojärvi** (2019) |
| lag functions of `gcd(n,h)` and their transforms | **Ramanujan–Fourier** analysis of `r`-even functions |
| the positive `log gcd` Gram kernel | the classical **gcd-matrix** literature |
| `λ_min → 0` for a purely atomic spectral measure | the **Szegő limit theorem** (Grenander–Szegő, 1958) |
| the critical balance of the Weil form at `β* = 1` | **Bombieri**, Rend. Lincei 11 (2000) |
| the arithmetic weight `μ(q)/φ(q)` on a major arc, and the cusp at rational angles | **Hardy–Littlewood** (1923) major arcs; **Ramanujan** sums (1918); the **Lerch–Wood** polylogarithm expansion at `z = 1` |
| the determinant of the division table after one column is falsified | **Redheffer** (1977); Barrett–Forcade–Pollington (1988); Vaughan (1993, 1996) |
| the density of the numbers a truncated sieve has not yet struck, as a function of `log x / log z` | **Buchstab's** function `ω` (1937), the dual of **Dickman's** `ρ` |
| repeat the cycle, delete one member of each family, merge the two neighbouring gaps | the **wheel sieve** — Pritchard (1982) — and cycles of gaps in Eratosthenes' sieve |

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
papers/     six self-contained papers, plus four short structural notes (Markdown)
figures/    the papers' figures, and two more drawn for the division-table paper
code/       verification scripts for every table and number, plus the figure generators
harness/    the fixed-support builder, the zero-side validator, the CvS comparison
            harness reporting inertia across an m-ladder with JSON metadata
docs/       the GitHub Pages site: one page and one PDF per paper
```

### `papers/a_numerical_study_of_the_division_table.md`
**Start here if you are new to the repository.** The most elementary paper of the set, and
deliberately expository: it starts from the quotients `n/d` for `n ≤ 100`, written as mixed
numbers, and reports at every step which classical object the sheet has become. Three
features before any theory — the whole cells are `D(x)`, the column sums are the harmonic
numbers, and the determinant of the `0–1` shadow, after one column is falsified, is the
**Mertens function**, so `det Rₙ = O(n^{1/2+ε})` is the Riemann hypothesis wearing a matrix.
Then a walk across the table traces an **arch**, the arch is **Fermat's method** (1643), its
failure is repaired by **Kraitchik** and the **quadratic sieve**, and rotating the picture
forty-five degrees lands on the **Dirichlet divisor** and **Gauss circle** problems,
Voronoï's series, and the `x^{1/4}` that both share. Two closing sections truncate the sheet
instead of walking it, which turns it into the sieve of Eratosthenes and reaches **Buchstab's
function** and the **wheel sieve**; the second of them measures that which row owns a struck
point is independent of the additive geometry around it — the additive/multiplicative divide
seen from the sieve side rather than argued.

Its ledger section lists twenty-one items and decomposes all twenty-one. **No result in it
is new**, and it says so in the abstract; what it offers is the route and an honest map of
where each step already lives.

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

*One open item was flagged as possibly having content beyond bookkeeping: the minimising
angle of the symbol locks onto rational points — exactly `π`, then exactly `2π/5` to eight
decimals — held there by a downward cusp whose exponent tracks the parameter, and comes
loose only when that cusp becomes Lipschitz. It has since been closed, and dissolved, in
`papers/rational_angle_lock.md`; the crossing point stated in this paper is corrected there.*

### `papers/odd_parity_sector.md`
A probe of the **odd parity sector** of the Connes–van Suijlekom Galerkin construction
(A. Groskin, `connes-cvs`; arXiv:2605.20224), which §10 of that paper lists as unprobed:
*we compute the even sector only; the odd-sector zero, if any, is not probed*.

The odd sector of the periodic construction is shown to coincide **exactly, at finite
size**, with a specific finite sine space — an identification that holds there and not on
the even sector, where only asymptotic completeness is available. On that sector, at
`c = 13, 17, 19, 23`, the eigenvector belonging to the smallest eigenvalue yields a
function whose **unseeded** sign-change scan returns the first ten zeta ordinates in
order, none missing and no additional sign-changing root detected. At `c = 13` the scan
extends to height 100 and returns exactly the twenty-nine ordinates below it. The scan
builds its function from the eigenvector coefficients alone, refines by bisection, writes
the root list to disk, and only then calls `mp.zetazero` to compare — no ordinate enters
the discovery. A random odd vector in the same basis, put through the identical scan,
comes no closer to any ordinate than order `1e-1`.

Two further measurements. The first-zero error tracks the smallest eigenvalue through
`C = |γ̂₁ − γ₁| / λ_min`, which grows with `c` in both parity sectors while the ratio
between them moves only a few percent. And the Galerkin exponent is nearly insensitive to
the archimedean cutoff: `46.14, 46.03, 45.93` at `T = 400, 800, 1600` on the published
grid, where the source paper's §8.2 mechanism predicted proportionality and would have
put the last value near `184.6`.

A later section extends the same scan from the ground state to **every near-null branch of
both parity blocks** at `c = 100`. Across twelve cell–sector combinations the number of
roots a branch produces beyond the thirteen ordinates is non-decreasing in `|λ|`, with the
sign of the eigenvalue interleaved throughout and no separate sign effect detected — in one
cell the count steps by one through four alternating sign changes. The negative
eigenvalues of those blocks, which are a finite-cutoff artifact, are counted against the
archimedean cutoff and do not reach zero at the cutoffs tested; a parity-count pattern that
held in five successive cells is reported together with the sixth cell that breaks it.

The source paper's §8.2 mechanism has since been withdrawn and the provenance of its
Table 14 corrected; both are recorded in that package's `ERRATA.md`, along with a later
correction to the cutoff at which the negative block is said to vanish. The `extract_zeros`
precision trap noted here is fixed in `connes-cvs` 0.3.0.

*Everything in it is numerical, at four cutoffs, on the first ten ordinates below height
50 — twenty-nine below 100 at `c = 13` only. It says nothing about whether these roots
converge to the zeta zeros as `c → ∞`, which is the open question the construction exists
to raise, and nothing about RH.*

### `papers/prime_edge_two_paths.md`
A short note, not a paper of its own. It measures the prime-edge derivative jump of the
truncated Weil form in a Dirichlet sine construction and compares it with the matrix-valued
von Mangoldt measure of Groskin (2026), where the jump is rank one. The two differ in where
the edge vanishing lives: inside the matrix in the reference path, outside it as a scalar
window here — so they are distinct families rather than the same object in different
coordinates. The note carries a closed form for the measured plateau, `(m-1)/(m+2)`, due to
A. Groskin in correspondence and verified here to seven digits, and a corrections section
recording what earlier versions of the note got wrong.

### `papers/li_lens_law.md`
A short note on the **Li–Sekatskii coefficients** — the one-parameter family of Li
coefficients obtained by shifting the Möbius map, whose unit circle is the critical line
for every admissible shift. The whole construction is Sekatskii's and the note says so in
a prior-art table; what it adds is that the shift behaves as a lens with an **exact optimal
width**: for a hypothetical zero at `ρ = 1/2 + β + iT` the detection rate is maximal at
`d = |ρ − 1/2|`, where it equals `artanh(β/|ρ − 1/2|)`, so matching the lens to the zero
improves the detection index by a factor `T`. A previously measured empirical `30.9×`
improvement for a planted zero at height `30` turns out to be that `T`.

The note also records a truncation artifact that manufactures a false detection — its rate
is `(1 + 1/|a|)`, which is Sekatskii's own exponential term — and asks the same optimisation
question of the truncated Weil form, where the answer is the opposite: the support length
has no interior optimum, being bounded by a Paley–Wiener price rather than a matching
condition. It also states how far its own prior-art search went, and which sources it did
not consult.

*Its own conclusion, and the reason it is a note rather than a paper: a more sensitive
detector is not a shorter proof.*

### `papers/rational_angle_lock.md`
A short note that closes the open item flagged in `papers/division_table_single_column.md`
above, by dissolving it. The cusp law is
`f(2πa/q + ε) − f(2πa/q) ~ 2 (μ(q)/φ(q)) Γ(−α) cos(πα/2) |ε|^α`: the **Hardy–Littlewood
singular series** multiplying the **Lerch–Wood expansion of the polylogarithm** at `z = 1`.
That gives the exponent (`β = α`, no arithmetic in it), the sign (down exactly where
`μ(q) = −1`), the strength (`1/φ(q)`) and the transition at `α = 1` (the cusp becomes
Lipschitz). Rescaling every measured edge by `φ(q)/μ(q)` collapses them onto one prime-free
universal cusp — mean `0.9944`, sd `0.0096`. The phenomenon is the major-arc approximation
of the circle method drawn as a graph.

It also **corrects a number in the single-column paper**: the handover from `π` to `2π/5` is at
`α* = 0.74005083`, not "about `α = 0.71`", and which rational wins is decided by the value
at the angle, not by the cusp strength — which is why `q = 3`, twice as strong as `q = 5`,
never wins.

*What survives as its own is small and labelled as such: for `0 < α < 1` a sub-Lipschitz
cusp cannot be dislodged by a smooth background, so the minimiser is confined to rational
angles with `μ(q) = −1` and can only ever jump.*

### `papers/phase_and_masking.md`
A companion to the note above, turning the same lens around to ask when it is **blind**. At
the matched width the phase of the image is exactly `−π/2` independent of `β` and `T`, so a
planted off-line quartet contributes on a strict `n mod 4` cycle: detection at `n ≡ 0`,
active masking at `n ≡ 2`, and no exponential term at all for odd `n`. Following the rate
`artanh(2βd/(d² + β² + T²))` to the far end recovers **Sekatskii's Theorem 3** as the other
zero of one curve — a shift taken too far does not merely fail to help, it actively delays —
and the detection index is measured to grow **linearly** in the shift, `N(d)·r(d) ≈ log M`,
so hiding a violation past order `m` costs a shift of order `mR`.

Two natural refinements are tested and **fail**: phase-tuned lenses do not beat `d = R` on
the first negative index, and the far detection pockets at `d ≈ nT/πk`, though real and
exactly where predicted, are swamped by the background by seven orders of magnitude.

*Its usable consequence: a numerical search that finds the first `m` coefficients positive
has established nothing unless it states its shift.*

### `harness/`
`conventions.py` is the single source of truth for the support, basis, block definitions,
tolerance and reporting requirements. `builder_sine.py` assembles the truncated Weil form
with archimedean, pole and prime blocks exposed separately. `validator.py` checks that
assembly non-circularly against ordinates from `mpmath.zetazero` — the zeros enter no part
of the construction — reporting residuals against explicit tail bounds rather than
asserting agreement. `run_ladder.py` sweeps the `m`-ladder, prints the full inertia triple
with its stated tolerance, writes raw arrays and JSON metadata, and exits nonzero on
regression failure.

`compare_families.py` sets the periodic and sine constructions against each other sector
by sector, with two deliberately mispaired controls so that a small agreement figure is
evidence rather than arithmetic. `odd_sector_unseeded_scan.py`,
`odd_sector_scan_height100.py`, `cross_check_extract_zeros.py`,
`reconstruction_controls.py`, `odd_vector_reality.py`, `odd_sector_zeros.py`,
`convergence_sweep.py` and `sobolev_slope.py` are the odd-parity-sector campaign; each
carries a failure gate on the claim it supports and exits nonzero when that claim does
not hold. `negative_branch.py` scans every near-null branch of either parity block —
`--sector even|odd` reuses one matrix, and the matrix, the eigendecomposition and the scan
are cached separately so a control is a rescan rather than a rebuild — and
`root_precision_probe.py` recomputes a root at several working precisions and compares the
roots to each other rather than to zeta, which is what makes a displacement figure a
measurement rather than a guess.

Built to a specification supplied by A. Groskin so that a second constructor can be run
against the same harness.

---

## Reproducibility

Every table and number is regenerated by a script in `code/` — or, for the odd-parity-sector
paper, in `harness/` — and the figures of three of the papers by the generators listed below.

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
python code/verify_li_lens.py         # the lens law; eleven checks, the proof verified line by line
python code/verify_phase_masking.py   # the phase cycle and the masking law; 13 checks, --fast
python code/detection_index.py        # the detection index at a few heights
python code/detection_index_sweep.py  # the 24-run gain table of the lens note's section 3.3
python code/verify_angle_lock.py      # the rational-angle cusp law; --fast for a 5 s run
python code/verify_division_table.py  # the division-table paper, sections 1-9; 105 checks
python code/verify_redheffer.py       # its section 1 again, in exact rational arithmetic
python code/verify_knowledge_clock.py # its section 10, Buchstab's omega and the clock u
python code/verify_row_inheritance.py # its section 11, A_p = p S_p^- and the window limit

cd code
python fig.py       # figure for the division-table paper
python robust_paper3.py && python fig3.py   # fig3 needs the gap data robust_paper3 caches
python fig4.py      # figure for the dynamical-model paper
cd ..

cd harness
python validator.py 4.5 100            # non-circular check against the zeta zeros
python run_ladder.py                   # inertia ladder + JSON; nonzero exit on failure

# the odd-parity-sector paper; these need `pip install connes-cvs` and, above N ≈ 40,
# `python-flint`. Each exits nonzero when the claim it supports does not hold.
python compare_families.py             # the sector comparison and its two mispaired controls
python odd_sector_unseeded_scan.py     # the unseeded scan, frozen root list, random control
python odd_sector_scan_height100.py    # the extension to height 100 from a cached eigenvector
python cross_check_extract_zeros.py    # the two-implementation cross-check on all 29 roots
python reconstruction_controls.py      # the four control vectors
python convergence_sweep.py            # the lag table between the sectors
python odd_vector_reality.py           # the real-rootedness ensemble, experiments A and B
python sobolev_slope.py --stage1       # the T-scaling grid
python negative_branch.py --c 100 --N 100 --T 400 --dps 500 --sector even   # section 6
python negative_branch.py --c 100 --N 100 --T 400 --dps 500 --sector odd    # same matrix
python edge_precision.py --m 20 40     # the edge limit in extended precision
cd ..
```

Four of the scripts — `verify_li_lens.py`, `detection_index.py`, `detection_index_sweep.py`
and `verify_phase_masking.py` — read a list of critical-line ordinates from the environment
variable `ZETA_ZEROS`: one ordinate per line, increasing. The tables in the two Li-lens notes
use the first 100,000 (`γ` from `14.1347` to `74920.8275`). Such a list can be taken from
Odlyzko's tables (`www.dtc.umn.edu/~odlyzko/zeta_tables`) or the LMFDB, or generated with
`mpmath.zetazero`, which is far slower. Without the variable each script says so and skips
only the tables that need it. A list **shorter** than a table needs is never silently used in
its place: the table is skipped and reported at the end.

Tested with Python 3.12, NumPy 2.4, SymPy 1.14, mpmath 1.3.
The four division-table scripts take about ten, one, four and five minutes respectively on a
full run; every one of them accepts `--fast`, and `verify_row_inheritance.py` also accepts
`--deep` for the `10⁹` window row of section 11.6.
The gap-frequency scripts segment-sieve windows of width `4×10⁷` up to `10¹³` and take a
few minutes. `verify_covariance.py` and `verify_positivity.py` also take a few minutes —
each has a long-average or large-matrix step, and `verify_covariance.py` has a `FAST` flag
that shortens the slowest one. `harness/run_ladder.py` sweeps five prime cutoffs
across the full `m`-ladder and validates at every rung, so the default run takes
about ten minutes; pass `--pcs 90 --ladder 100 200` for a quick check. Everything
else runs in seconds.

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

## On the use of AI assistance

The verification scripts in `code/` and `harness/`, and much of the prose in `papers/`,
were written with the assistance of **Claude (Anthropic)**, used as a working collaborator
throughout: drafting and rewriting code, running the computations, drafting and editing
text, searching the literature, and — most usefully — auditing the papers against their own
scripts.

The research direction, the questions asked, the decisions about what to publish and what
to withdraw, and the final responsibility for every claim are the author's.

A great many of the corrections recorded in these papers were found by that auditing: a
claim stated more precisely than the computation supported, a script that had not caught up
with a correction to the text, a formula quoted but never tested. Several were errors the
assistant had itself introduced and then found on a later pass. Where a result is reported
here, it is because a script regenerates it and the script has been read; that discipline,
rather than any assurance about the tool, is what the reader is asked to rely on.

---

## Status

Six papers and four notes are here.

The newest paper, `a_numerical_study_of_the_division_table.md`, is the odd one out and the
natural entry point: it is expository, it claims nothing new anywhere, and its ledger
decomposes all twenty-one of its items into named classical results. It is also the only
paper here that could reasonably be submitted to a journal, because an expository paper is
accepted for its route rather than for a result.

`papers/division_table_single_column.md` closes the division-table route as a direct
approach to RH, and says so in
its own conclusion. What it leaves genuinely open is listed in its Appendix A. The most
interesting of those items — the rational-angle locking of the symbol's minimum, the one
measured structure that did not obviously reduce to bookkeeping — has since been closed in
`papers/rational_angle_lock.md`, and it did reduce to bookkeeping: the singular series
times a polylogarithm cusp, which is the circle method's major arcs seen as a graph. That
note also corrects the crossing point that paper reported.

The four notes are smaller and were written the same way: an elementary route is followed
until it lands on something already published, the prior art is named, and what survives is
a measurement or two. The first ends on two constructions that look like one object and are
not. The second ends on the finding this repository keeps arriving at — a sharper instrument
is not a shorter proof. The third goes one step further along, on a structure that looked
like it might be an exception and was not. The fourth turns the same instrument around and
asks when it is blind, which is the one question of the four whose answer is usable. That is by now the
pattern the repository is really a record of.

The sixth paper, `papers/odd_parity_sector.md`, is the only one here that is not about the
division table. It probes a sector of someone else's construction that that paper listed as
unprobed, and it was written in correspondence with its author, who has reproduced its
measurements independently. Nothing in it is new mathematics either; what it adds is a
sector that had not been looked at and a set of numbers taken carefully enough to be worth
checking.

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
