# A Structural Difference in the Prime-Edge Jump

### On the matrix factor accompanying the von Mangoldt weight in two constructions of the truncated Weil form

July 2026

---

## Summary

Groskin (2026), *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path* (Zenodo [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028)), states that on the finite matrix path `u ↦ Q(u)`, with `u = log c` the prime cutoff, crossing a prime-power threshold produces a first-derivative jump

```
Δ (dQ/du)  =  − 2 Λ(q) / ( √q · log q )  ×  (universal rank-one all-ones matrix).
```

This note reports a measurement of the same jump in a **different construction** — an orthonormal sine basis with a Gram-type overlap, rather than divided-difference entries on squared-integer nodes. The result is a clean structural difference, and it is stated here as a question about conventions, not as a disagreement with the theorem.

**In the construction used here the matrix factor is not universal.** It depends on `q`, and across nine prime powers the normalised matrix factors are close to mutually orthogonal, with pairwise cosines in `[−0.083, +0.100]` and mean `+0.006`. A universal factor would give `1.000` throughout.

**A note on what is not being claimed.** The scalar `−2Λ(q)/(√q log q)` also appears here and matches to four decimal places, but this is *not* offered as independent evidence for anything. With the triangular window `w(q,u) = 1 − log q/u`, the scalar follows from one line of differentiation, and `Λ(q)` appears in the output because it was placed in the sum as input. Reproducing it confirms that the code here matches the arithmetic here. It is the *matrix factor*, not the scalar, that carries the content of the cited theorem — and that is where the two constructions part company.

---

## Construction

Only the prime block can contribute a singular part; the archimedean and pole blocks are smooth in `u` and cancel out of any jump. So only the prime block is built.

```
basis        φ_j(x) = √(2/L) · sin(a_j x),   a_j = jπ/L,   on [0, L]
overlap      R_jk(t) = ∫₀^{L−t} φ_j(x) φ_k(x+t) dx,  symmetrised
window       w(q, u) = 1 − log q / u
prime block  P(u) = − 2 Σ_{q = p^k, log q ≤ u} Λ(q) q^{−1/2} w(q, u) R(log q)
```

Differentiating and isolating the discontinuity at `u = log q` gives

```
Δ (dP/du)  =  − 2 Λ(q) / ( √q · log q )  ×  R(log q) ,
```

so in this construction **the matrix factor is `R(log q)`, the basis overlap at that particular shift.** Since different prime powers sit at different shifts, the factor is manifestly `q`-dependent. The question is how far apart the factors actually are, and whether anything in the setup could make them collapse onto a common matrix.

The support `L` is held fixed and only the cutoff `u` varies, decoupling two parameters that the reference parametrisation `L = log c` ties together.

---

## Measurement 1 — the matrix factors are near-orthogonal across `q`

Pairwise cosines between the normalised matrix factors `R(log q) / ‖R(log q)‖`, at `L = 4.5`, `m = 60`:

| | 3 | 5 | 7 | 9 | 11 | 13 | 25 | 27 | 49 |
|---|---|---|---|---|---|---|---|---|---|
| **3** | 1.000 | 0.027 | −0.020 | 0.017 | −0.013 | −0.013 | 0.007 | −0.007 | −0.006 |
| **5** | 0.027 | 1.000 | 0.064 | −0.017 | 0.025 | 0.013 | −0.012 | 0.011 | 0.008 |
| **7** | −0.020 | 0.064 | 1.000 | −0.083 | 0.005 | 0.025 | −0.001 | 0.002 | −0.000 |
| **9** | 0.017 | −0.017 | −0.083 | 1.000 | 0.100 | 0.020 | −0.019 | 0.017 | 0.011 |
| **11** | −0.013 | 0.025 | 0.005 | 0.100 | 1.000 | 0.087 | 0.004 | −0.002 | −0.003 |
| **13** | −0.013 | 0.013 | 0.025 | 0.020 | 0.087 | 1.000 | 0.027 | −0.023 | −0.013 |
| **25** | 0.007 | −0.012 | −0.001 | −0.019 | 0.004 | 0.027 | 1.000 | −0.008 | 0.009 |
| **27** | −0.007 | 0.011 | 0.002 | 0.017 | −0.002 | −0.023 | −0.008 | 1.000 | −0.014 |
| **49** | −0.006 | 0.008 | −0.000 | 0.011 | −0.003 | −0.013 | 0.009 | −0.014 | 1.000 |

Off-diagonal range `[−0.083, +0.100]`, mean `+0.006`. The factors are essentially orthogonal.

## Measurement 2 — the factors are also not rank one

At every shift tested the jump matrix is full rank: `σ₂/σ₁ ≈ 1.0`, effective rank `= m`. Pushing the shift toward the support edge — where an entering prime sits when `L = log c` — does not produce a rank-one limit either. With `ε = 1 − t/L` at `m = 80`:

| ε | 10⁻¹ | 10⁻² | 10⁻³ | 10⁻⁴ |
|---|---|---|---|---|
| `σ₁` | 5.00×10⁻¹ | 1.96×10⁻¹ | 2.90×10⁻⁴ | 2.91×10⁻⁷ |
| `σ₂/σ₁` | 1.000 | 0.973 | 0.9635 | 0.9634 |

The leading singular value collapses by six orders of magnitude while `σ₂/σ₁` plateaus near `0.96`. The same plateau appears at `m = 40, 80, 160`, so it is not an artifact of basis size.

A candidate mechanism was tested and fails. As `t → L` the overlap region `[0, L−t]` shrinks to the point `x = 0`, which would suggest a rank-one limit with vector `v_j ∝ φ_j′(0) ∝ j`; the measured alignment `|cos(u₁, v)|` plateaus at `0.71`, not `1`.

*Numerical caveat.* At `ε = 3×10⁻⁴` the measured/predicted ratio jumps to `6–17`; this is the finite-difference step interacting with a rapidly varying overlap, not a feature. The `ε = 3×10⁻²` and `3×10⁻³` rows are the reliable ones.

---

## What this does and does not establish

Rank is invariant under change of basis, and so is the mutual orthogonality of a family of matrices under a fixed congruence. So the two constructions cannot be describing the same operator family in different coordinates; something in the setup differs.

The likely candidates, none of which this note can adjudicate:

- **Node placement.** The reference builds divided-difference entries on squared-integer nodes `0, 1, 4, …, N²`; this note uses a Gram-type overlap on a sine basis. Universality of the matrix factor may be a property of the former.
- **The archimedean cutoff.** The reference construction carries none; the form here is built with the support and shift range tied together.
- **The decoupling itself.** Holding `L` fixed while `u` varies may not be an admissible object in the reference framework at all, in which case the comparison is not like for like.

**Open question.** Which feature of the construction makes the matrix factor universal? Equivalently: is `L = log c` structural, or a parametrisation choice?

Nothing here bears on the Riemann Hypothesis. No claim of priority is made and no error is alleged; the cited theorem is a statement about a specific construction, and this note only records that a differently built object behaves differently, together with the numbers.

---

## Reproducibility

```bash
python code/prime_edge_jump.py    # the jump, and the coupled-setting comparison
python code/prime_edge_rank.py    # rank and edge analysis
```

Python 3.12, NumPy 2.4, SymPy 1.14.

## References

1. A. Groskin (2026). *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path.* Zenodo, [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028).
2. A. Groskin (2026). *High-Precision Approximation of Riemann Zeros via the Truncated Weil Form.* arXiv:2605.20224.
3. A. Connes and W. van Suijlekom (2026). *Quadratic forms, real zeros and echoes of the spectral action.* arXiv:2511.23257; Comm. Math. Phys. **406**.
4. A. Connes, C. Consani and H. Moscovici (2025). *Zeta spectral triples.* arXiv:2511.22755.
