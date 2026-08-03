# Two Paths Through the Prime Edge

### Why the derivative jump is rank one in one truncation of the Weil form and not in another

July 2026

---

## Summary

Groskin (2026), *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path* (Zenodo [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028)), shows that on the finite matrix path $u \mapsto Q(u)$, with $u = \log c$ the prime cutoff, crossing a prime-power threshold produces a first-derivative jump equal to $-2\Lambda(q)/(\sqrt q \log q)$ times a universal rank-one matrix.

This note measures the corresponding jump in a **different construction** — a Dirichlet sine space with the cutoff applied as a separate external window — where the matrix factor is instead $q$-dependent. It then gives the reason, which is a genuine difference of path rather than a discrepancy.

**The two constructions differ in where the edge vanishing lives.** In the Connes–van Suijlekom parametrisation the periodic Fourier modes $U_n(x) = L^{-1/2} e^{2\pi inx/L}$ have symmetrised overlap exactly $A_N(1 - t/L)$ for a source at $t$. At $L = u = t$ this gives $A_N(0) = 0$ and $A_N'(0) = 2\cdot \mathbf{1}\mathbf{1}^{\mathsf{T}}$, so the vanishing at the edge is carried **inside the matrix** and the jump is rank one for that reason. There is no additional scalar window outside it.

The construction here takes a Dirichlet sine overlap and multiplies it by an external window $1 - t/u$. The sine space is not an invertible change of the finite periodic Fourier space — the endpoint maps and the frequency sets differ — so this is a different family, and its first-derivative factor is the overlap $R(t)$, which depends on where the prime power sits.

$L = u = \log c$ is structural for the reference path. A decoupled $(L, u)$ family is mathematically legitimate but is a different object.

> **Attribution.** The edge expansion, the $(m-1)/(m+2)$ law, the alternating endpoint vector, the path separation, and the section-C diagnosis are due to Akiva Groskin; see the [public derivation](https://github.com/akivag613/connes-cvs-/issues/2#issuecomment-5102029641). They are independently checked numerically here.

---

## Construction measured here

Only the prime block can contribute a singular part; the archimedean and pole blocks are analytic at a prime-power threshold and cancel out of any jump.

```
basis        φ_j(x) = √(2/L) · sin(a_j x),   a_j = jπ/L,   on [0, L]
overlap      R_jk(t) = ∫₀^{L−t} φ_j(x) φ_k(x+t) dx,  symmetrised
window       w(q, u) = 1 − log q / u          ← applied outside the matrix
prime block  P(u) = − 2 Σ_{q = p^k, log q ≤ u} Λ(q) q^{−1/2} w(q, u) R(log q)
```

Holding $L$ fixed and varying $u$, differentiation gives $\Delta(dP/du) = -2\Lambda(q)/(\sqrt q\cdot \log q) \times R(\log q)$.

**The scalar here is not evidence of anything.** It drops out of the window in one line, and $\Lambda(q)$ appears in the output because it was placed in the sum as input. What distinguishes the two paths is the matrix factor.

---

## Measurement 1 — the matrix factor is not universal here

Pairwise cosines between the normalised factors $R(\log q)/\Vert R(\log q)\Vert$, at $L = 4.5$, $m = 60$, over $q = 3, 5, 7, 9, 11, 13, 25, 27, 49$:

```
off-diagonal range   [−0.083, +0.100]
mean                 +0.006
```

A universal factor would give $1.000$ throughout. Proportionality to a common matrix is basis-independent — if $A_q = \lambda_q J$ in one basis then $S^{\mathsf{T}}A_qS = \lambda_q(S^{\mathsf{T}}JS)$ in any other — and proportional matrices have pairwise $|\cos | = 1$. So the two families are not the same object in different coordinates.

## Measurement 2 — the edge behaviour, and its closed form

Write $\varepsilon = 1 - t/L$ for the fractional distance to the edge, so $\delta = L - t = L\varepsilon$. At $L = 4.5$ these differ by a factor $4.5$; the tables below are indexed by $\varepsilon$.

| ε | 10⁻¹ | 10⁻² | 10⁻³ | 10⁻⁴ |
|---|---|---|---|---|
| $\sigma_{1}$ (m = 80) | 5.00×10⁻¹ | 1.96×10⁻¹ | 2.90×10⁻⁴ | 2.91×10⁻⁷ |
| $\sigma_{2}/\sigma_{1}$ | 1.000 | 0.973 | 0.9635 | 0.9634 |

**The plateau has an exact value.** With $\delta = L - t$,

```
R_jk(L − δ)  =  (π² δ³ / 6L³) · j k · ( (−1)^{j+1} + (−1)^{k+1} )  +  O(δ⁵),
```

opposite parities decouple exactly, and the leading matrix is one rank-one odd block minus one rank-one even block. For even $m$,

```
σ₂ / σ₁  ⟶  (m − 1) / (m + 2) .
```

Checked at $\varepsilon = 10^{-4}$, i.e. $\delta = 4.5\times 10^{-4}$:

| m | 20 | 40 | 80 | 160 | 200 |
|---|---|---|---|---|---|
| $(m-1)/(m+2)$ | 0.8636364 | 0.9285714 | 0.9634146 | 0.9814815 | 0.9851485 |
| measured | 0.8636366 | 0.9285719 | 0.9634156 | 0.9814833 | 0.9851509 |
| \|difference\| | 2.3×10⁻⁷ | 4.5×10⁻⁷ | 9.2×10⁻⁷ | 1.9×10⁻⁶ | 2.3×10⁻⁶ |

**Agreement to within $2.4\times 10^{-6}$.** That is a statement about the absolute discrepancy at one value of $\varepsilon$; it is not a digit count.

**In extended precision the law converges, and it converges like $\delta^{2}$.** The discrepancy from $(m-1)/(m+2)$:

| $\varepsilon$ | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ |
|---|---|---|---|---|
| $m = 20$ | 2.1×10⁻⁵ | 2.1×10⁻⁷ | 2.1×10⁻⁹ | 2.1×10⁻¹¹ |
| $m = 40$ | 4.45×10⁻⁵ | 4.45×10⁻⁷ | 4.45×10⁻⁹ | 4.45×10⁻¹¹ |

Two decades of $\varepsilon$ per two decades of discrepancy, monotone throughout — which is what a $\delta \to 0$ limit is supposed to do, and what float64 could not show.

In float64 the expansion checks out only over the range where float64 survives. Residuals $\Vert R - leading\Vert/\Vert R\Vert$ at $\delta = 10^{-2}, 10^{-3}, 10^{-4}$ are $5.0\times 10^{-3}, 5.0\times 10^{-5}, 7.8\times 10^{-6}$, with successive ratios $100.3$ and $6.4$; only the first step shows the $\delta^{2}$ scaling an $O(\delta^{5})$ remainder predicts, and at $\delta = 10^{-5}$ the residual rises again to $6.5\times 10^{-3}$.

**In extended precision the remainder is confirmed to be $O(\delta^{5})$ throughout** (`harness/edge_precision.py`, $m = 20$, precision set per $\delta$):

| $\delta$ | 10⁻² | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ |
|---|---|---|---|---|---|
| $\Vert R - leading\Vert/\Vert R\Vert$ | 1.290×10⁻³ | 1.289×10⁻⁵ | 1.289×10⁻⁷ | 1.289×10⁻⁹ | 1.289×10⁻¹¹ |
| ratio | — | 100.07 | 100.001 | 100.0 | 100.0 |

The turnaround was the arithmetic, not the mathematics. Parity decoupling is exact to machine precision in either regime: $\max |cross-parity| = 4\times 10^{-16}$ against $\max |same-parity| = 1\times 10^{-2}$.

## Measurement 3 — the split blocks approach rank one

*(This check was suggested by A. Groskin; the numbers below are the result of running it.)*

Separating odd from even sine modes, the odd-block ratio $\sigma_{2}/\sigma_{1}$ falls as the edge is approached and then **rises again** once cancellation dominates:

| ε | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ |
|---|---|---|---|---|
| odd block $\sigma_{2}/\sigma_{1}$ (m = 80) | 3.6×10⁻⁷ | 9.8×10⁻⁹ | 2.6×10⁻⁵ | 1.4×10⁻² |

Each block therefore **approaches** rank one, and the non-monotone tail is the float64 limit rather than the mathematics.

**Taken in extended precision the ratio falls monotonically, and its rate is $\varepsilon^{4}$** ($m = 20$):

| $\varepsilon$ | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ | 10⁻⁸ |
|---|---|---|---|---|---|
| odd block | 1.3546×10⁻⁹ | 1.35434×10⁻¹³ | 1.354341×10⁻¹⁷ | 1.354341×10⁻²¹ | 1.354341×10⁻²⁹ |
| even block | 1.6516×10⁻⁹ | 1.65129×10⁻¹³ | 1.651282×10⁻¹⁷ | 1.651282×10⁻²¹ | 1.651282×10⁻²⁹ |

so $\sigma_{2}/\sigma_{1} = C(m) \varepsilon^{4} + O(\varepsilon^{6})$ with $C$ converged to ten figures — $1.354341014$ on the odd block and $1.651282092$ on the even. A finite ratio at finite $\varepsilon$ still does not *prove* exact rank one, but the limit is now exhibited rather than asserted. The dominant eigenvalue is positive on the odd block and negative on the even block, so the full matrix is their difference — which is why the combined ratio sits just below $1$ rather than near $0$.

**The second endpoint vector is $(-1)^{j+1} j$,** not $j$. Testing the leading singular vector against $j$ alone gives $0.720083, 0.713664, 0.710403$ at $m = 40, 80, 160$ — converging to $1/\sqrt 2 = 0.70711$, because $j$ is a 45° mix of the odd-only and even-only vectors. That plateau was a mixing angle, not a failed match.

---

## The fixed-support family: what is measured, and how it is reported

The decoupled $(L, p_c)$ family is well-defined and distinct from the reference path. Whether it is *arithmetically* interesting is a separate matter, and worth stating plainly: the coupling $L = \log c$ is what makes the reference path the Weil quadratic form, so a decoupled family may well turn out to be well-defined and arithmetically inert.

**A bare count `#{λ < τ}` is not reported here, and earlier versions of this note were wrong to report one.** Such a count is a near-null count only where the form is positive semidefinite. Off the saturated diagonal this family is strongly indefinite, so the count silently sums negative and near-null directions. Every figure below carries the full inertia triple with its stated tolerance.

At $L = 4.5$, $\tau = 10^{-6}$, saturated $p_c = \lfloor e^{4.5}\rfloor = 90$:

| $p_c$ | m | $\lambda_{\min}$ | $n_{-}$ | $n_{0}$ | $n_{+}$ |
|---|---|---|---|---|---|
| 30 | 100 | −4.6621 | 25 | 17 | 58 |
| 30 | 200 | −4.7059 | 45 | 16 | 139 |
| 60 | 100 | −2.0688 | 10 | 44 | 46 |
| 60 | 200 | −2.1182 | 19 | 65 | 116 |
| 60 | 300 | −2.1341 | 27 | 76 | 197 |
| 60 | 400 | −2.1425 | 35 | 79 | 286 |
| 90 | 100 | −1.1×10⁻¹³ | 0 | 59 | 41 |
| 90 | 200 | −1.2×10⁻¹² | 0 | 96 | 104 |

On the saturated diagonal $n_{-} = 0$ and the form is positive semidefinite to precision. Off it, $\lambda_{\min}$ is $O(1)$ negative and $n_{-}$ dominates.

**$\lambda_{\min}$ converges; the counts do not.** These are two different situations and the note keeps them apart, because only one of them is a usable number.

At $p_c = 60$ the successive $\lambda_{\min}$ drifts are $2.39$%, $0.75$%, $0.39$% across $m = 100 \to 200 \to 300 \to 400$, roughly halving at each rung and extrapolating to about $-2.151$. The harness gate compares the top two rungs of the ladder, so on the default ladder $(100, 200, 300, 400, 500)$ it passes; it is only the $100 \to 200$ step that exceeds the $2$% tolerance, and a run restricted to those two rungs is correctly rejected.

The counts show no such behaviour. Over the same rungs $n_{-}$ runs $10, 19, 27, 35$ and $n_{0}$ runs $44, 65, 76, 79$, both still climbing with no sign of settling, while $n_{+}$ merely absorbs the remainder of $m$. They drift with $\tau$ as well as with $m$.

**This is why the harness gates on $\lambda_{\min}$ and reports the triple without gating on it.** $\lambda_{\min}$ is a converged quantity of the form and may be quoted; the inertia counts at a fixed absolute $\tau$ are basis-size artefacts as much as they are properties of the form, and nothing about the shape of the off-diagonal surface should be inferred from them.

**Prior art does not transfer here.** Bombieri's negative-eigenvalue count theorem is stated for his zero-indexed matrix $H(\Gamma ; t)$, not for this sine-window family. Inertia remains the natural invariant, but that theorem is not invoked.

**No relation to the diagonal law is asserted.** On the saturated diagonal the near-null count is the object counted by Connes–Consani as $1 + \nu(\lambda^{2}) \sim 2\lambda^{2}$ (arXiv:2106.01715, ultimately Slepian–Pollak 1961). Earlier versions asserted $d = 2c - 6$ as an equality; at $L = 4.5$ the measured value is $173$ against $2\cdot 90 - 6 = 174$, and the offset is $6$ at some $L$ and $7$ at others. **That equality is withdrawn** pending a precise statement of the counting convention, which is not settled here.

---

## Corrections to earlier versions

Recorded rather than silently edited. Items 5–10 are due to A. Groskin.

1. **The scalar was presented as an independent verification. It is not.** It follows from the external window in one line, and $\Lambda(q)$ appears in the output because it was input. The claim that it "tests the von Mangoldt weight itself" was wrong and is withdrawn.
2. **The Galerkin nodes were misdescribed** as squared integers $0, 1, 4, \ldots, N^{2}$. Those belong to the source-to-jet dictionary; the Galerkin modes are the integer modes ${-N, \ldots, N}$.
3. **The coupled-path figure in `prime_edge_jump.py` §C was an artifact.** The entering term is $(1 - t/u)R(t;u) = O((u-t)^{4})$, so the true first-derivative jump is zero; the printed value is a finite-gap difference of the smooth background, scaling linearly with the gap ($1.881\times 10^{0} \ldots 1.881\times 10^{-3}$ at $gap = 10^{-3} \ldots 10^{-6}$). The section is retained with a warning rather than deleted.
4. **A stability argument was stated as invariance of orthogonality under congruence.** The Frobenius inner product is not congruence-invariant; the correct invariant is proportionality to a common matrix, which the argument now uses.
5. **$\varepsilon$ was labelled $\delta$.** The law table is taken at $\varepsilon = 1 - t/L = 10^{-4}$, hence $\delta = L - t = 4.5\times 10^{-4}$.
6. **"Seven digits" was wrong.** The correct statement is agreement within $2.4\times 10^{-6}$ at that $\varepsilon$.
7. **The $\delta^{2}$ claim was over-extended *as float64 evidence*.** In double precision only the first two residuals scale that way and the third is cancellation-limited, so the original claim was not supported by the computation behind it. Extended precision (§ Measurement 2) since confirms the scaling itself, out to $\delta = 10^{-6}$ — the correction was to the evidence, not to the mathematics.
8. **"Is rank one" and "full rank = $m$" were both overstated.** The parity blocks *approach* rank one; effective ranks near the edge fall well below $m$. Extended precision now exhibits that approach with a rate, $\sigma_{2}/\sigma_{1} = C(m) \varepsilon^{4}$, but a limit at finite $\varepsilon$ still does not prove exact rank one.
9. **Bombieri's count theorem was invoked where it does not apply** — it is stated for $H(\Gamma ; t)$, not this family.
10. **$d = 2c - 6$ was asserted as an equality** despite an offset that is $6$ or $7$ depending on $L$; withdrawn above.


---

## Reproducibility

```bash
python code/prime_edge_jump.py       # the jump; §C retained with its artifact annotated
python code/prime_edge_rank.py       # edge analysis; superseded mechanism kept as a record
python code/check_reply.py           # the closed form and the split blocks

cd harness
python validator.py 4.5 100          # non-circular check against the zeta zeros
python run_ladder.py --L 4.5         # inertia ladder, JSON metadata, nonzero exit on failure
python edge_precision.py --m 20 40   # the edge limit in extended precision
```

`harness/conventions.py` is the single source of truth for the support, basis, block definitions, tolerance and reporting requirements. Python 3.12, NumPy 2.4, SymPy 1.14, mpmath 1.3.

Edge-limit quantities are cancellation-limited in float64. They are recomputed in extended precision by `harness/edge_precision.py`, which sets the working precision from the target $\delta$ and evaluates the overlap through a reduction with no trigonometric cancellation left in its structure; the float64 figures are retained above for contrast, with the extended-precision result beside each.

## AI assistance

The verification script and much of the prose in this note were written with the assistance of Claude (Anthropic); the research direction, the decisions about what to publish, and responsibility for every claim are the author's. See the repository README for a fuller statement.

---
## References

1. A. Groskin (2026). *A matrix-valued von Mangoldt measure in the finite Connes–van Suijlekom path.* Zenodo, [10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028).
2. A. Groskin (2026). *High-Precision Approximation of Riemann Zeros via the Truncated Weil Form.* arXiv:2605.20224.
3. A. Groskin (2026). *A finite Guinand–Weil dictionary and archimedean tail order for the truncated Weil quadratic form.* arXiv:2607.02828.
4. A. Connes and W. D. van Suijlekom (2025). *Quadratic forms, real zeros and echoes of the spectral action.* Comm. Math. Phys. **406**, 312. [doi:10.1007/s00220-025-05493-1](https://doi.org/10.1007/s00220-025-05493-1)
5. A. Connes, C. Consani and H. Moscovici (2025). *Zeta spectral triples.* arXiv:2511.22755.
6. A. Connes and C. Consani (2021). *Spectral triples and ζ-cycles.* arXiv:2106.01715; Enseign. Math. **69**.
7. D. Slepian and H. O. Pollak (1961). *Prolate spheroidal wave functions, Fourier analysis and uncertainty — I.* Bell System Tech. J. **40**, 43–63.
8. E. Bombieri (2000). *Remarks on Weil's quadratic functional in the theory of prime numbers.* Rend. Lincei **11**, 183–233.
