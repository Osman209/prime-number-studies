# Independent Check of the Prime-Edge Derivative Jump

### A verification note on the matrix-valued von Mangoldt measure of Groskin (2026)

July 2026

---

## What is being checked

Groskin (2026), *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path* (Zenodo [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028)), states the following. Fix a Galerkin level, let `u = log c` be the prime cutoff and vary it continuously. Crossing a prime-power threshold `u = log q`, the first derivative of the matrix path `u ↦ Q(u)` jumps by

```
Δ (dQ/du)  =  − 2 Λ(q) / ( √q · log q )  ×  (rank-one all-ones matrix).
```

This note reports an independent numerical check of that statement in a **different construction**: an orthonormal sine basis rather than the divided-difference form on squared-integer nodes used in the reference implementation.

The result splits. **The scalar is reproduced exactly. The rank-one structure is not.** Since matrix rank is invariant under change of basis, the second half is read here as a difference of convention between the two constructions, and is stated as an open question rather than a disagreement.

---

## Construction

Only the prime block can contribute a singular part — the archimedean and pole blocks are smooth in `u` and cancel out of any jump — so only the prime block is built.

```
basis        φ_j(x) = √(2/L) · sin(a_j x),   a_j = jπ/L,   on [0, L]
overlap      R_jk(t) = ∫₀^{L−t} φ_j(x) φ_k(x+t) dx,  symmetrised
window       w(q, u) = 1 − log q / u
prime block  P(u) = − 2 Σ_{q = p^k, log q ≤ u} Λ(q) q^{−1/2} w(q, u) R(log q)
```

Differentiating and isolating the discontinuity at `u = log q` gives the prediction

```
Δ (dP/du)  =  − 2 Λ(q) / ( √q · log q )  ×  R(log q) ,
```

which is what is measured below, by one-sided numerical derivatives either side of the threshold.

**Setting.** The support `L` is held fixed and only the cutoff `u` varies. This decouples the two parameters, which the reference parametrisation `L = log c` ties together.

---

## Result 1 — the scalar is exact

`L = 4.5`, basis size `m = 60`, step `h = 10⁻⁵`. "Ratio" is `‖measured jump‖ / ‖predicted jump‖`; "cos" is the cosine between them as matrices.

| q | Λ(q) | log q | ratio | cos |
|---|---|---|---|---|
| 3 | 1.09861 | 1.0986 | 0.999828 | +1.000000 |
| 5 | 1.60944 | 1.6094 | 0.999869 | +1.000000 |
| 7 | 1.94591 | 1.9459 | 0.999891 | +1.000000 |
| **9** | **1.09861** | 2.1972 | 0.999960 | +1.000000 |
| 11 | 2.39790 | 2.3979 | 0.999904 | +1.000000 |
| 13 | 2.56495 | 2.5649 | 0.999902 | +1.000000 |
| **25** | **1.60944** | 3.2189 | 0.999982 | +1.000000 |
| **27** | **1.09861** | 3.2958 | 0.999894 | +0.999999 |
| **49** | **1.94591** | 3.8918 | 0.999839 | +0.999999 |

The residual departure from `1` is the finite-difference truncation error and shrinks with `h`.

**The prime powers are the sharp part of this test.** For `q = 9, 25, 27, 49` the correct weight is `Λ(q) = log p`, not `log q` — `Λ(9) = log 3`, not `log 9`. The measured jumps carry `log p`. The identity therefore tests the von Mangoldt weight itself, not merely a normalising coefficient, and it is confirmed in a construction that shares no code with the reference implementation.

---

## Result 2 — the rank-one structure does not appear in this basis

In the sine basis the jump matrix is **full rank at every shift tested**: `σ₂/σ₁ ≈ 1.0` and effective rank `= m`.

The natural place to look for the rank-one form is the support edge, since with `L = log c` an entering prime sits exactly there. Writing `ε = 1 − t/L`, at `m = 80`:

| ε | 10⁻¹ | 10⁻² | 10⁻³ | 10⁻⁴ |
|---|---|---|---|---|
| `σ₁` | 5.00×10⁻¹ | 1.96×10⁻¹ | 2.90×10⁻⁴ | 2.91×10⁻⁷ |
| `σ₂/σ₁` | 1.000 | 0.973 | 0.9635 | 0.9634 |

The leading singular value collapses by six orders of magnitude, but the ratio `σ₂/σ₁` **plateaus near `0.96` and never approaches zero**. The same plateau appears at `m = 40, 80, 160`, so it is not a basis-size artifact.

A candidate mechanism was also tested and does not hold. As `t → L` the overlap region `[0, L−t]` shrinks to the single point `x = 0`, which would suggest a rank-one limit with vector `v_j ∝ φ_j′(0) ∝ j`. The measured alignment `|cos(u₁, v)|` plateaus at `0.71`, not `1`.

*Numerical caveat.* At `ε = 3×10⁻⁴` the measured/predicted ratio jumps to `6–17`. This is the finite-difference step interacting with a rapidly varying overlap, not a feature; the `ε = 3×10⁻²` and `3×10⁻³` rows are the reliable ones.

---

## Reading

Matrix rank is invariant under change of basis, so a genuinely rank-one operator would present as rank one in any basis. The discrepancy is therefore located in the constructions rather than in the theorem.

The most likely source is that the two objects are not the same matrix. The reference implementation builds entries of divided-difference type on squared-integer nodes, with no archimedean cutoff, whereas the construction here is a Gram-type overlap on a sine basis with the support tied to the shift range. The all-ones form is natural in the first parametrisation and would not be expected to survive verbatim into the second — but which convention carries the rank-one structure, and whether the decoupled setting used here is a legitimate object at all, is not something this note can settle.

**Summary.**

| statement | status here |
|---|---|
| jump scalar `−2Λ(q)/(√q log q)` | reproduced exactly, independent basis, incl. prime powers |
| jump `=` scalar `×` basis overlap at the shift | reproduced exactly |
| jump is rank one | not reproduced; full rank at all shifts tested |

Nothing here bears on the Riemann Hypothesis, and no claim of priority or of error is intended.

---

## Reproducibility

```bash
python code/prime_edge_jump.py    # Result 1, and the coupled-setting comparison
python code/prime_edge_rank.py    # Result 2, the edge and rank analysis
```

Python 3.12, NumPy 2.4, SymPy 1.14.

## References

1. A. Groskin (2026). *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path.* Zenodo, [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028).
2. A. Groskin (2026). *High-Precision Approximation of Riemann Zeros via the Truncated Weil Form.* arXiv:2605.20224; Zenodo [10.5281/zenodo.19546514](https://doi.org/10.5281/zenodo.19546514).
3. A. Connes and W. van Suijlekom (2026). *Quadratic forms, real zeros and echoes of the spectral action.* arXiv:2511.23257; Comm. Math. Phys. **406**.
4. A. Connes, C. Consani and H. Moscovici (2025). *Zeta spectral triples.* arXiv:2511.22755.
