# Two Paths Through the Prime Edge

### Why the derivative jump is rank one in one truncation of the Weil form and not in another

July 2026

---

## Summary

Groskin (2026), *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path* (Zenodo [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028)), shows that on the finite matrix path `u ↦ Q(u)`, with `u = log c` the prime cutoff, crossing a prime-power threshold produces a first-derivative jump

```
Δ (dQ/du)  =  − 2 Λ(q) / ( √q · log q )  ×  (universal rank-one all-ones matrix).
```

This note records a measurement of the corresponding jump in a **different construction** — a Dirichlet sine space with the cutoff applied as a separate scalar window — where the matrix factor is instead `q`-dependent and full rank. It then gives the reason, which is not a discrepancy but a genuine difference of path.

**The two constructions differ in where the edge vanishing lives.** In the Connes–van Suijlekom parametrisation the periodic Fourier modes `U_n(x) = L^{−1/2} e^{2πinx/L}` have symmetrised overlap exactly `A_N(1 − t/L)` for a source at `t`. At `L = u = t` this gives `A_N(0) = 0` and `A_N′(0) = 2·𝟙𝟙ᵀ` — the vanishing at the edge is carried **inside the matrix**, and the derivative jump is rank one for that reason. There is no additional scalar window outside it.

The construction measured here instead takes a Dirichlet sine overlap and multiplies it by an external window `1 − t/u`. The sine space is not an invertible change of the finite periodic Fourier space — the endpoint maps and the frequencies differ — so this is a different family, and its first-derivative factor is the overlap `R(t)`, which depends on where the prime power sits.

**The account of the numerical plateau, and of the second endpoint vector, is due to A. Groskin** (correspondence, July 2026); it is reproduced and verified here with his permission and attribution. `L = u = log c` is structural for the path in his paper. A decoupled `(L, u)` family is mathematically legitimate but is a different object.

---

## Construction measured here

Only the prime block can contribute a singular part; the archimedean and pole blocks are analytic at a prime-power threshold and cancel out of any jump.

```
basis        φ_j(x) = √(2/L) · sin(a_j x),   a_j = jπ/L,   on [0, L]
overlap      R_jk(t) = ∫₀^{L−t} φ_j(x) φ_k(x+t) dx,  symmetrised
window       w(q, u) = 1 − log q / u          ← applied outside the matrix
prime block  P(u) = − 2 Σ_{q = p^k, log q ≤ u} Λ(q) q^{−1/2} w(q, u) R(log q)
```

Holding `L` fixed and varying `u`, differentiation gives

```
Δ (dP/du)  =  − 2 Λ(q) / ( √q · log q )  ×  R(log q) .
```

**The scalar here is not evidence of anything.** It drops out of the window in one line, and `Λ(q)` appears in the output because it was placed in the sum as input. What distinguishes the two paths is the matrix factor.

---

## Measurement 1 — the matrix factor is not universal here

Pairwise cosines between the normalised factors `R(log q) / ‖R(log q)‖`, at `L = 4.5`, `m = 60`, over `q = 3, 5, 7, 9, 11, 13, 25, 27, 49`:

```
off-diagonal range   [−0.083, +0.100]
mean                 +0.006
```

A universal factor would give `1.000` throughout. Proportionality to a common matrix is basis-independent — if `A_q = λ_q J` in one basis then `Sᵀ A_q S = λ_q (Sᵀ J S)` in any other — and proportional matrices have pairwise `|cos| = 1`. So the two families are not the same object in different coordinates.

Consistently, the jump is full rank at every shift tested: `σ₂/σ₁ ≈ 1.0`, effective rank `= m`.

## Measurement 2 — the edge behaviour, and its closed form

Pushing the shift to the support edge, with `ε = 1 − t/L`, the leading singular value collapses while the ratio plateaus:

| ε | 10⁻¹ | 10⁻² | 10⁻³ | 10⁻⁴ |
|---|---|---|---|---|
| `σ₁` (m = 80) | 5.00×10⁻¹ | 1.96×10⁻¹ | 2.90×10⁻⁴ | 2.91×10⁻⁷ |
| `σ₂/σ₁` | 1.000 | 0.973 | 0.9635 | 0.9634 |

**The plateau has an exact value.** Writing `δ = L − t`,

```
R_jk(L − δ)  =  (π² δ³ / 6L³) · j k · ( (−1)^{j+1} + (−1)^{k+1} )  +  O(δ⁵),
```

opposite parities decouple exactly, and the leading matrix is one rank-one odd block minus one rank-one even block. For even `m` this gives

```
σ₂ / σ₁  ⟶  (m − 1) / (m + 2) .
```

*(Expansion and law: A. Groskin, correspondence.)* Verified here at `δ = 10⁻⁴`:

| m | 20 | 40 | 80 | 160 | 200 |
|---|---|---|---|---|---|
| `(m−1)/(m+2)` | 0.8636364 | 0.9285714 | 0.9634146 | 0.9814815 | 0.9851485 |
| measured | 0.8636366 | 0.9285719 | 0.9634156 | 0.9814833 | 0.9851509 |

Agreement to seven digits. At `δ = 10⁻⁶` the measurement degrades — floating-point cancellation on an `O(δ³)` quantity, not a departure from the law.

The expansion itself checks out: `‖R − leading‖ / ‖R‖ = 5.0×10⁻³, 5.0×10⁻⁵, 7.8×10⁻⁶` at `δ = 10⁻², 10⁻³, 10⁻⁴`, falling like `δ²` as an `O(δ⁵)` remainder requires, with scale ratio `1.000001` confirming the `π²/6L³` prefactor. Parity decoupling is exact to machine precision: `max|cross-parity| = 4×10⁻¹⁶` against `max|same-parity| = 1×10⁻²`.

## Measurement 3 — the split blocks are each rank one

Separating odd from even sine modes, at `δ` corresponding to `ε = 10⁻⁵`:

| m | 40 | 80 | 160 |
|---|---|---|---|
| odd block `σ₂/σ₁` | 8.4×10⁻⁶ | 2.6×10⁻⁵ | 4.0×10⁻⁶ |
| even block `σ₂/σ₁` | 8.1×10⁻⁶ | 2.7×10⁻⁵ | 4.0×10⁻⁶ |
| dominant eigenvalue sign | odd `+`, even `−` | odd `+`, even `−` | odd `+`, even `−` |

Each block is rank one, with opposite signs. The full matrix is their difference, which is why the combined ratio sits just below `1` rather than near `0`.

**The second endpoint vector is `(−1)^{j+1} j`,** not `j`. Testing the leading singular vector against `j` alone gives `0.720083, 0.713664, 0.710403` at `m = 40, 80, 160` — converging to `1/√2 = 0.70711`, because `j` is a 45° mix of the odd-only and even-only vectors. That plateau was a mixing angle, not a failed match.

---

## Corrections to earlier versions of this note

Recorded rather than silently edited.

**1. The scalar was presented as an independent verification. It is not.** With the external window `w(q,u) = 1 − log q/u`, differentiating gives `−2Λ(q)q^{−1/2}(log q/u²)R`, which at `u = log q` is the stated scalar in one line. Reproducing it confirms that this code matches this arithmetic, nothing more. The earlier claim that it "tests the von Mangoldt weight itself, not merely a normalising coefficient" was wrong and is withdrawn.

**2. The Galerkin nodes were misdescribed.** An earlier version attributed "squared-integer nodes `0, 1, 4, …, N²`" to the reference Galerkin construction. Those nodes belong to the source-to-jet dictionary in that paper; the Galerkin modes are the **integer modes `{−N, …, N}`**.

**3. The coupled-path figure in `prime_edge_jump.py` §C was an artifact.** On that path the entering term is `(1 − t/u) R(t;u) = O((u−t)⁴)`, so the true first-derivative jump is **zero**; the nonzero quantity printed there is a finite-gap difference of the smooth background. It scales linearly with the gap — `1.881×10⁰, 1.881×10⁻¹, 1.881×10⁻², 1.881×10⁻³` at `gap = 10⁻³ … 10⁻⁶` — and vanishes in the limit. *(Diagnosis: A. Groskin.)* The section is retained with this annotation rather than deleted.

**4. The plateau was recorded as an open question.** It has a closed form; see Measurement 2.

---

## What remains open

The decoupled `(L, p_c)` family is well-defined and distinct. Whether it is *arithmetically* interesting is a separate matter, and worth stating plainly: the coupling `L = log c` is what makes the reference path the Weil quadratic form, so a decoupled family may well turn out to be well-defined and arithmetically inert.

The concrete quantity available for that question is the near-null count `d = #{λ < 10⁻⁶}`. On the saturated diagonal `p_c = ⌊e^L⌋` it reproduces `d = 2c − 6`, i.e. the Connes–Consani `2λ²` count (arXiv:2106.01715, ultimately Slepian 1961) — a rediscovery, not a result. Off the diagonal, at `L = 4.5` fixed:

| `p_c` | 30 | 45 | 60 | 75 | 90 |
|---|---|---|---|---|---|
| `d` | 94 | 103 | 122 | 156 | 173 |
| `2·p_c` | 60 | 90 | 120 | 150 | 180 |

with a non-monotonic per-unit slope. Whether the off-diagonal region is positive semidefinite at all is not yet established here, and that is the first thing to settle before the surface itself is worth fitting.

Nothing in this note bears on the Riemann Hypothesis. No claim of priority is made and no error in the cited work is alleged.

---

## Acknowledgement

The closed form `(m−1)/(m+2)`, the `δ³` expansion, the identification of the second endpoint vector, the account of where the two paths separate, and the diagnosis of the `§C` artifact are all due to **A. Groskin**, in correspondence of July 2026. They are reproduced here with attribution; the verification of each against the numerics is the contribution of this note.

## Reproducibility

```bash
python code/prime_edge_jump.py    # the jump; §C retained with its artifact annotated
python code/prime_edge_rank.py    # rank, edge analysis, factor universality
python code/check_reply.py        # verification of the closed form and the split blocks
```

Python 3.12, NumPy 2.4, SymPy 1.14.

## References

1. A. Groskin (2026). *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path.* Zenodo, [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028).
2. A. Groskin (2026). *High-Precision Approximation of Riemann Zeros via the Truncated Weil Form.* arXiv:2605.20224.
3. A. Connes and W. van Suijlekom (2026). *Quadratic forms, real zeros and echoes of the spectral action.* arXiv:2511.23257; Comm. Math. Phys. **406**. *(Proposition 4.1.)*
4. A. Connes, C. Consani and H. Moscovici (2025). *Zeta spectral triples.* arXiv:2511.22755.
5. A. Connes and C. Consani (2021). *Spectral triples and ζ-cycles.* arXiv:2106.01715; Enseign. Math. **69**.
6. D. Slepian and H. O. Pollak (1961). *Prolate spheroidal wave functions, Fourier analysis and uncertainty — I.* Bell System Tech. J. **40**, 43–63.
