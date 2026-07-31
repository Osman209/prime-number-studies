# The Rational-Angle Lock

### Why the minimum of a von Mangoldt symbol sits on rational angles, and what that costs

---

## Summary

The fourth paper in this repository, *The Odd Division Table as a Single Column*,
ends with three open items. The second of them is the only place in that document
where a measured structure looked as though it might carry content beyond
bookkeeping: the minimising angle of the Toeplitz symbol

```
f_α(θ) = 2 Σ_{q = p^k} Λ(q) q^{−(1+α)} cos(q θ)
```

**locks onto rational angles**. It sits exactly at `θ = π`, then jumps — without
sliding — to exactly `2π/5`, and only comes loose above `α ≈ 1`. The paper measured
a downward cusp `f(2π/5 + ε) − f(2π/5) ~ ε^β` with `β ≈ α`, and left four questions
open: why `β = α`, why the transition is at `α = 1`, which rational angle takes over
and when, and what the upper edge of the `α`-window is.

This note answers the first three and dissolves the phenomenon. The edge law is

```
f_α(2πa/q + ε) − f_α(2πa/q)  ~  2 (μ(q)/φ(q)) Γ(−α) cos(πα/2) |ε|^α ,   (a,q) = 1
```

— a product of two entirely classical ingredients. The arithmetic factor
`μ(q)/φ(q)` is the residue at `s = 1` of `Σ Λ(n) e(an/q) n^{−s}`, which is the
**Hardy–Littlewood singular series**; the shape `Γ(−α)|ε|^α` is the **Lerch–Wood
expansion of the polylogarithm** at `z = 1`, and contains no arithmetic at all.
Their product, evaluated near rationals, is the **major-arc approximation of the
circle method drawn as a graph**. Rescaling every measured edge by `φ(q)/μ(q)`
collapses all of them onto one prime-free universal cusp: mean ratio `0.9944`,
standard deviation `0.0096`, over seven values of `q` and four values of `ε`.

So: `β = α` because the cusp is a polylogarithm; the sign is `−μ(q)` (downward at
`μ(q) = −1`, upward at `μ(q) = +1`, absent at `μ(q) = 0` — verified, the amplitude
falls by four orders of magnitude at `q = 4, 9, 12`); the strength is `1/φ(q)`; and
the transition at `α = 1` is where the cusp becomes Lipschitz and can no longer
hold a minimum against a linear background.

Which rational wins is *not* decided by cusp strength. Every downward cusp sits at
`μ(q) = −1`, so the global minimum is at one of those angles, and the competition is
between the **values** `f_α(2πa/q)`. That is why `q = 3`, with strength `1/2`, never
wins although it is twice as strong as `q = 5`. Enumerating all `97` candidates with
`μ(q) = −1, q ≤ 40`: `π` wins below the crossing, `2π/5` above it, and nothing else
ever comes within `3%`.

**This corrects a number in the fourth paper.** That paper puts the handover at
"about `α = 0.71`". Bisecting `f_α(π) − f_α(2π/5)` gives

```
α* = 0.74005083     (cutoff 2×10⁷; 0.740038 at cutoff 4×10⁶)
```

At `α = 0.71` the gap is still `−4.56×10⁻²`, i.e. `π` is winning comfortably.

What survives as this note's own is small and elementary, and is stated as such in
§6: for `0 < α < 1` the cusp exponent beats any smooth linear background, so the
minimum **cannot** slide and must instead jump when two locked wells cross. That is
a correct observation about a classical object, not a new object.

---

## Status legend

| tag | meaning |
|---|---|
| **[K]** | Known. Re-derived here for self-containedness; source named in §1. No novelty claimed. |
| **[P]** | Proved here. |
| **[M]** | Measured, with parameters stated. |
| **[C]** | Caveat. |

---

## 1. Prior art

| component used here | where it already lives |
|---|---|
| `c_q(n) = Σ_{(a,q)=1} e(an/q)`, and `c_q(1) = μ(q)` | **Ramanujan** (1918); Hardy & Wright, *An Introduction to the Theory of Numbers*, §16.6 |
| `Σ_{n≤x} Λ(n) e(an/q) ~ (μ(q)/φ(q)) x` for fixed `q`, `(a,q)=1` | the prime number theorem in arithmetic progressions; **Siegel–Walfisz**. Davenport, *Multiplicative Number Theory*, §22 |
| `μ(q)/φ(q)` as *the* arithmetic weight on a major arc | **Hardy–Littlewood** (1923); Vaughan, *The Hardy–Littlewood Method*, §2 |
| `Li_s(e^{iε}) − Li_s(1) ~ Γ(1−s)(−iε)^{s−1}` | **Lerch**; **Wood** (1992), *The computation of polylogarithms*; Erdélyi et al., *Higher Transcendental Functions* I §1.11 |
| the Dirichlet-character decomposition `e(an/q) = φ(q)^{−1} Σ_χ χ̄(a) τ(χ) χ(n)` | standard; Davenport §9 |
| non-continuability of `Σ Λ(n) e(nθ) n^{−s}` past `Re s = 1`, and the density of its singularities on rationals | classical, and the reason the circle method splits into major and minor arcs |

**Nothing in §3 is new.** The edge law is the major-arc expansion, and this note's
contribution to it is a verification and an identification, not a theorem.

### How far the prior-art search went, and what it did not cover [C]

Searched: the circle method literature for major-arc expansions of `ψ(x,θ)`;
polylogarithm asymptotics at the unit circle; Ramanujan-sum expansions of `Λ`.
**Not** searched, and where an exact prior statement of the edge law is most likely
to sit: the literature on **boundary behaviour of lacunary and arithmetic Dirichlet
series on their line of convergence**, and on `Σ Λ(n)n^{−s}e(nθ)` as a function of
`θ` for fixed `s`. A reader who knows that literature should assume the law is in it
and tell us where.

---

## 2. The object, and where the question comes from

`f_α` is the Toeplitz symbol of §9 of the fourth paper: `f_α(θ) = c + 2 Σ_{h≥1} b(h) cos(hθ)`
with `b(h) = Λ(h) h^{−(α+1)}`, `c = Σ_q Λ(q) q^{−(α+1)}`. Its positivity for all `N`
is equivalent to `f_α ≥ 0`, so the location and value of `min_θ f_α` is what decides
the `α`-window in that paper. Two conventions occur there: all prime powers, and odd
prime powers only (§4 of that paper excludes the prime 2). Everything below is
verified in **both**.

---

## 3. The edge law [P for the residue, K for the transfer, M throughout]

### 3.1 Statement

For `(a,q) = 1`, `0 < α < 1`, and `ε → 0`:

```
f_α(2πa/q + ε) − f_α(2πa/q)  =  2 (μ(q)/φ(q)) Γ(−α) cos(πα/2) |ε|^α  +  o(|ε|^α).
```

Three consequences, all of them the whole content of the phenomenon:

- **sign** — `Γ(−α)cos(πα/2) < 0` for `0 < α < 1`, so `μ(q) = −1` gives a *downward*
  cusp (a local minimum), `μ(q) = +1` an upward one, and `μ(q) = 0` no leading edge at all;
- **strength** — proportional to `1/φ(q)`;
- **exponent** — exactly `α`, with no arithmetic in it.

### 3.2 Where `μ(q)/φ(q)` comes from [P]

Write `D_{a,q}(s) = Σ_n Λ(n) e(an/q) n^{−s}`. For `(n,q) = 1`,
`e(an/q) = φ(q)^{−1} Σ_χ χ(a) τ(χ̄) χ(n)`, so

```
D_{a,q}(s) = φ(q)^{−1} Σ_χ χ(a) τ(χ̄) (−L′/L)(s, χ) + E_q(s),
```

with `E_q` collecting the finitely many `p | q` and holomorphic for `Re s > 0`. Only
`χ_0` has a pole at `s = 1`, simple, of residue `1`, inherited from `ζ`; and
`τ(χ_0) = Σ_{(b,q)=1} e(b/q) = c_q(1) = μ(q)`. Hence

```
Res_{s=1} D_{a,q}(s) = μ(q)/φ(q).
```

*Verified [M]:* `c_q(1) = μ(q)` to machine precision for `q = 2 … 40`; and
`Ψ_{a,q}(x)/x` lands on `μ(q)/φ(q)` at `x = 10⁵, 10⁶, 3×10⁶` for
`q = 2, 3, 5, 6, 7, 10, 11, 15, 30`, and on `0` for `q = 4, 9, 12`.

### 3.3 The transfer to the boundary, and the one step that is not a one-liner [C]

The residue must be carried from `s = 1` to the edge `|ε|^α`. The route is Mellin:
`G(z) = (2πi)^{−1} ∫ Γ(w) D(1+α+w) z^{−w} dw`, shifted left past `w = −α`, picking up
`Γ(−α)(μ/φ) z^α`.

This is not free, and the obstruction should be stated rather than hidden. One
cannot shift to a fixed line `Re w = −α−δ`, because `−L′/L` has poles at every zero
of `L(s,χ)` and zeros approach `Re s = 1`. The contour must hug the classical
zero-free region `σ > 1 − c/log(q(|t|+2))`, which costs `exp(−c√(log(1/|ε|)))`
rather than `|ε|^δ`. For **fixed** `q` this is Siegel–Walfisz, and the statement is
unconditional with error `O(|ε|^α exp(−c√(log(1/|ε|))))`. A zero `ρ = β + iγ`
contributes `z^{1+α−ρ}`, of size `|z|^{1+α−β} = o(|z|^α)` since `β < 1`; under GRH
the second layer is exactly `|z|^{α+1/2}`.

**Partial summation does not close this** [C]. With `|R(x)| ≪ x E(x)` it leaves
`|ε| ∫_{1/|ε|}^∞ E(x) x^{−α} dx`, which **diverges** for `α < 1` unless the
oscillation of `e^{ixε}` is used. The `Γ`-factor decay on the Mellin contour is
exactly what supplies it. So Mellin is not a stylistic preference here; it is the
only route that closes.

*The transfer lemma isolated and tested prime-free* [M], on `a_n ≡ 1`, at `α = 0.5`,
summed to `2×10⁷`: measured `/|ε|^α` gives `−2.492, −2.481, −2.462, −2.425` against
`Γ(−α)cos(πα/2) = −2.5066` at `ε = 10⁻³ … 3×10⁻⁵`. The approach is from above and
slow, because the remainder is only `O(1)` relative to a vanishing leading term. It
confirms that the shape and its constant are prime-free, and that all arithmetic
enters solely through `μ(q)/φ(q)`.

### 3.4 Measurements [M]

Symmetric second difference `½[f(θ₀+ε) + f(θ₀−ε)] − f(θ₀)`, which kills the smooth
linear background. `Λ` sieved to `2×10⁷` (`1 271 339` prime powers). Fitted exponent
`β` from `ε = 10⁻³` to `3×10⁻⁵`; amplitude `C_q` at `ε = 3×10⁻⁵`, reported as the
normalisation-free ratio `C_q/C_2`, whose prediction is `−μ(q)/φ(q)` because
`μ(2)/φ(2) = −1`.

**All prime powers, `α = 0.75`:**

| `a/q` | `μ(q)/φ(q)` | `β` | `C_q/C_2` | predicted |
|---|---|---|---|---|
| `1/2` | `−1.0000` | `0.751` | `1.0000` | `1.0000` |
| `1/3` | `−0.5000` | `0.751` | `0.5000` | `0.5000` |
| `1/5` | `−0.2500` | `0.750` | `0.2499` | `0.2500` |
| `1/6` | `+0.5000` | `0.750` | `−0.4999` | `−0.5000` |
| `1/7` | `−0.1667` | `0.748` | `0.1668` | `0.1667` |
| `1/10` | `+0.2500` | `0.750` | `−0.2500` | `−0.2500` |
| `1/11` | `−0.1000` | `0.737` | `0.1001` | `0.1000` |
| `1/15` | `+0.1250` | `0.760` | `−0.1242` | `−0.1250` |
| `1/30` | `−0.1250` | `0.749` | `0.1249` | `0.1250` |
| `1/4`, `1/9`, `1/12` | `0` | — | `≤ 10⁻⁴` | `0` |

**Odd prime powers only, `α = 0.75`:** the same table to the same digits
(`β = 0.737 … 0.760`, ratios agreeing to `10⁻³`). The reason is one line: dropping
the `2`-tower subtracts `log 2 · Σ_k e(a2^k/q) 2^{−ks}`, which is holomorphic at
`s = 1`. The residue, and hence the entire leading edge, is untouched; only the
smooth background moves.

**The paper's measured exponents were converging, and to `α`.** That paper reported
`β = 0.624, 0.769, 1.004, 1.063` at `α = 0.60, 0.75, 1.00, 1.40` at cutoff
`1.6×10⁷` (those four, and the drift values quoted in §4, are regenerated by
`code/verify_positivity.py`, not by this note's script), and noted that they moved *towards* `α` as the cutoff grew. They were, and
the drift is checked directly here: at `θ = π`, `α = 0.75`, the fitted exponent is
`0.7548` at cutoff `4×10⁶` and `0.7511` at `2×10⁷`. There is no separate exponent to
find; there is one exponent, equal to `α`, seen through a truncation.

**One universal cusp** [M]. Rescaling each measured edge by `φ(q)/μ(q)` and dividing
by the prime-free form `2Γ(−α)cos(πα/2)|ε|^α`, over `q = 2, 3, 5, 6, 7, 11, 30` and
`ε = 10⁻³ … 3×10⁻⁵`: mean `0.9944`, sd `0.0096`, range `0.949 … 1.008`. **The shape
carries no prime information.** All of it is in one arithmetic number per angle.

---

## 4. Why the minimum locks, and why it stops at `α = 1` [P, elementary]

Near a downward cusp the symbol behaves as `f(θ₀) − C|ε|^α + b₁ε + O(ε²)` with
`C > 0`. For `α < 1` the cusp term dominates the linear background for all small `ε`,
so `θ₀` is a strict local minimum **for every** `b₁`: the minimum cannot be moved by
a perturbation of the background, and therefore cannot slide. The global minimiser
must be one of the countably many downward-cusp points, and it changes only by
**jumping** when two of them cross in value.

At `α = 1` the cusp becomes Lipschitz and the two terms are of the same order. Above
it the minimum detaches, at distance

```
|ε*| ~ (|b₁| / (Cα))^{1/(α−1)},
```

which is the slow drift the fourth paper measured (`0.39950π` at `α = 1.40`,
`0.39892π` at `α = 1.60`) and could not account for.

---

## 5. Which rational angle wins, and where the handover is [M]

Since every downward cusp is at `μ(q) = −1`, the global minimum for `α < 1` is at one
of those angles, and is decided by the **value** there, not by the cusp strength.
Enumerating all `97` candidates with `μ(q) = −1`, `q ≤ 40`:

| `α` | winner | `f_min` | runners-up |
|---|---|---|---|
| `0.50` | `π` | `−1.493197` | `2π/5: −0.8783`, `2π/7: −0.7540` |
| `0.60` | `π` | `−1.007147` | `2π/5: −0.7319`, `2π/7: −0.6338` |
| `0.70` | `π` | `−0.690138` | `2π/5: −0.6280`, `12π/31: −0.5475` |
| `0.735` | `π` | `−0.605498` | `2π/5: −0.5982`, `12π/31: −0.5299` |
| `0.745` | `2π/5` | `−0.590249` | `π: −0.5833`, `12π/31: −0.5250` |
| `0.75` | `2π/5` | `−0.586331` | `π: −0.5725`, `12π/31: −0.5226` |
| `0.90` | `2π/5` | `−0.488253` | `12π/31: −0.4562`, `14π/37: −0.4430` |
| `1.00` | `2π/5` | `−0.438174` | `12π/31: −0.4177`, `14π/37: −0.4078` |

Note `q = 3`: cusp strength `1/2`, twice that of `q = 5`, and it never wins — the
handover is a competition of values, not of strengths.

**The crossing.** Bisecting `f_α(π) − f_α(2π/5)`:

```
α* = 0.74005083   at cutoff 2×10⁷
α* = 0.740038     at cutoff 4×10⁶      (shift 1.3×10⁻⁵)
```

### 5.1 Correction to the fourth paper [C]

Appendix A item 2 of *The Odd Division Table as a Single Column* states that the
minimum is at `π` "up to about `α = 0.71`". **That is wrong**, and the correct value
is `α* = 0.74005`. At `α = 0.71` the gap `f(π) − f(2π/5)` is still `−4.56×10⁻²`.
The two conventions of that paper are not the cause: the `0.71` does not correspond
to the odd-convention crossing either. Over `α = 0.30 … 0.99` the odd-convention
minimum stays at `π` throughout, and comfortably — the deepest rival reaches at most
`0.65` of the winner's depth — so there is no handover in that convention to confuse
it with.

*Reproduced on the way, to confirm the convention is the paper's own:* at `α = 0.75`,
`f(2π/5) = −0.586331` against the paper's `−0.58633` (all prime powers), and
`f(π) = −1.159010` against its `−1.158989` (odd convention).

---

## 6. What is offered as new, and what is not

**Not new:** everything in §3. The law is the major-arc expansion of the circle
method; `μ(q)/φ(q)` is the singular series; the cusp is the polylogarithm at `z = 1`.
The identification is this note's contribution to §3, and the identification is a
dissolution — the phenomenon that looked like it might carry content beyond
bookkeeping is bookkeeping, of the most classical kind available.

**Offered as new, and it is small:** §4 — that for `0 < α < 1` a sub-Lipschitz cusp
cannot be dislodged by a smooth background, so the minimiser of such a symbol is
confined to the rational angles with `μ(q) = −1` and can only ever *jump*; together
with §5, the resulting handover structure and its crossing point. This is an
elementary consequence, correctly observed, about an object that was already fully
understood. It is offered at that weight.

**What it does not settle** [C]: the upper edge of the `α`-window, `α = 1.6096`,
remains numerical — but now for a stated reason rather than as an unexplained gap.
Above `α = 1` the cusp is Lipschitz and the minimiser leaves every rational angle,
so there is no closed form of the kind that made the *lower* edge exact. Appendix A
item 2 should therefore be read as answered in its first three parts and
**structurally closed** in its fourth: not "we have not found the closed form" but
"the minimiser is not at a special point, so there is none to find".

---

## 7. Reproducibility

```bash
python code/verify_angle_lock.py          # cutoff 2×10⁷, about 30 s
python code/verify_angle_lock.py --fast   # cutoff 4×10⁶, about 5 s
```

Eleven checks, in the order in which the sections above use them: the Ramanujan sum,
the residue `Ψ_{a,q}(x)/x → μ(q)/φ(q)`, the prime-free transfer lemma, the edge law's
exponent and amplitude in both conventions, the universal cusp, the angle
enumeration, the crossing, the fourth paper's two numbers, and the cutoff stability
of `α*` and of a fitted exponent. The script regenerates every table and number
above, exits nonzero on failure, and prints its summation cutoff next to each result
— none of these numbers is cutoff-free, and the third of the repository's method
rules requires the cutoff to be printed beside anything that depends on it.

---

*The verification script and much of the prose in this note were written with the
assistance of Claude (Anthropic); the research direction, the decisions about what to
publish, and responsibility for every claim are the author's. See the repository README
for a fuller statement.*

---

## References

1. S. Ramanujan, *On certain trigonometrical sums and their applications in the theory of numbers*, Trans. Camb. Phil. Soc. **22** (1918) 259–276.
2. G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford, §16.6.
3. H. Davenport, *Multiplicative Number Theory*, 3rd ed., Springer GTM 74, §§9, 20, 22.
4. G. H. Hardy and J. E. Littlewood, *Some problems of 'Partitio Numerorum' III*, Acta Math. **44** (1923) 1–70.
5. R. C. Vaughan, *The Hardy–Littlewood Method*, 2nd ed., Cambridge, §2.
6. D. C. Wood, *The computation of polylogarithms*, Technical Report 15-92, University of Kent, 1992.
7. A. Erdélyi et al., *Higher Transcendental Functions*, vol. I, McGraw–Hill 1953, §1.11.
8. H. Iwaniec and E. Kowalski, *Analytic Number Theory*, AMS Colloquium Publications 53, §§4, 5.
9. M. Osman, *The Odd Division Table as a Single Column*, in this repository — Appendix A, item 2, which this note answers.

**No claim is made about the Riemann Hypothesis.** The object studied here is a
Toeplitz symbol whose positivity was already shown, in the paper this note serves, to
be a convention artifact carrying a comfortable margin; nothing about its minimum
bears on the critical line.
