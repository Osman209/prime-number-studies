# Matching the Lens to the Zero

### An exact optimal shift for the Li–Sekatskii coefficients, measured against the Weil window

July 2026

---

## Summary

The generalized Li criterion of **Sekatskii (2013, 2014)** — whose coefficients are now called the **Li–Sekatskii coefficients** in the literature — replaces Li's map by

```
q_a(ρ) = (ρ − a) / (ρ − (1 − a)),        a real, a ≠ 1/2,
```

whose unit circle is exactly the critical line, and shows that the non-negativity of the associated sums — equivalently, of certain derivatives of `log ξ` evaluated at `s = 1 − a` — is again equivalent to the Riemann Hypothesis, for every admissible `a`. The construction below is his; §1 names the source for each part.

This note adds a measurement and the small piece of analysis around it: **the shift is a lens with an optimal width, and the optimum has a closed form.**

For a hypothetical zero off the line at `ρ = 1/2 + β + iT`, write `d = 1/2 − a` for the lens width. The exponential rate at which that zero announces itself in the shifted sums is

```
log |q_a(ρ)|  =  ½ · log[ ((β+d)² + T²) / ((β−d)² + T²) ]
              =  2βd / (d² + T²)  +  O(β³),
```

This has an **exact** maximum. Writing `R = |ρ − 1/2| = √(β² + T²)`, the rate is stationary at `d = R` and only there, and its value is `artanh(β/R)`. Li's classical choice `a = 0` is the lens `d = 1/2`; for `β ≪ T` the two rates are `β/T` and `β/T²`.

> **The lens is optimal when its foci sit at the same distance from `1/2` as the zero does. Matching it improves the exponential detection rate — and hence the index at which the zero first shows — by a factor `T`.**

§3.4 shows why the detection index should scale like `1/log R`, so that maximising the rate minimises the index: the on-line terms are `1 − cos(nθ)` with `θ = 2arctan(γ/d) − π`, and the far zeros contribute a positive `n²d²` term. Those are the two ingredients a Palojärvi-type interval needs; the uniform estimates that would make it a certified interval are not supplied here.

§6 asks the same question of the truncated Weil form and gets the opposite answer: its support length has **no** interior optimum, the ceiling `e^{βL/2}` being monotone by Paley–Wiener. The Li parameter is a position and admits a matching condition; the Weil parameter is a bandwidth and is limited by a price rather than an optimum.

The gain is real and it is measurable. It is also, on its own, worth very little: §5 explains why a better detector is not a shorter proof, and §4 records a truncation artifact — already identified analytically by Sekatskii — that will manufacture a false detection if the prime side is cut before the cancellation is done.

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
| positivity of `λ_n` ⟺ RH | **Li**, *J. Number Theory* **65** (1997) 325–333 |
| the general multiset theorem behind it | **Bombieri & Lagarias**, *J. Number Theory* **77** (1999) 274–287 |
| the shifted map `q_a` and the derivative object `L_n(a)` | **Sekatskii**, arXiv:1304.7895 = *Ukrainian Math. J.* **66** (2014) 415–431 — both appear verbatim in the abstract |
| the equality between the zero sum and the derivative | same paper |
| the prime side as a generalized Laguerre transform of `Λ` | **Sekatskii**, arXiv:1404.7276, Lemma 1 and eq. (10); Laguerre polynomials in this setting go back to **Coffey** (2005, 2007, 2010) |
| the arithmetic interpretation | **Sekatskii**, arXiv:1305.1421 |
| asymptotics of the generalized sums as `n → ∞` | **Sekatskii**, arXiv:1403.4484; and **Voros**, *Math. Phys. Anal. Geom.* **9** (2006) 53–63 |
| the spurious exponential term `(1 + 1/b)^n` and its cancellation | **Sekatskii**, arXiv:1404.7276, eq. (7) and §4 |
| Li coefficients beyond `ζ` | **Lagarias** (automorphic `L`), **Smajlović** (Selberg class) |
| the *other* one-parameter family, `ρ ↦ ρ/(ρ−τ)` | **Freitas** (2006); **Droll** (2012) — the `τ`-Li coefficients |
| explicit detection indices: which `n` certifies a zero outside a region | **Brown** (2005); **Bucur, Ernvall-Hytönen, Odžak & Smajlović** (2016); **Palojärvi** (2019), whose Theorems 3.1 and 3.3 give explicit intervals `[n₀,n₁]` with tables |
| numerical behaviour of `τ`-Li coefficients for functions that violate RH | **Bucur, Ernvall-Hytönen, Odžak & Smajlović** (2016) |
| the ongoing programme | **Ernvall-Hytönen, Odžak, Smajlović & Zubača** (2024), *Variants of the Li-Type Criteria*; **Voros** (2020), discretized Keiper/Li |
| the shifted family carried to the Selberg class, with arithmetic and asymptotic representations | **Mazhouda & Sodaïgui** (2022), *The Li–Sekatskii coefficients for the Selberg class* — the name this paper adopts |
| the shifted criterion restated and extended by its author | **Sekatskii** (2021), in *Schrödinger Operators, Spectral Analysis and Number Theory*, Springer Proc. Math. Stat. **348** |

**Nothing in §2 is new.** It is written out because §3 needs it in one notation.

### 1.1 Two different one-parameter families, and which question each answers

There are two parameterised families of Li coefficients in circulation and they are **not** the same object.

| | `τ`-Li: `ρ ↦ ρ/(ρ−τ)` | shifted Li: `ρ ↦ (ρ−a)/(ρ−(1−a))` |
|---|---|---|
| foci of the map | `0` and `τ` | `a` and `1−a` |
| the line it tests | `Re s = τ/2` — **moves with `τ`** | `Re s = 1/2` — **fixed** |
| what the parameter does | selects which half-plane is being certified | changes how the same line is viewed |
| what it is used for | zero-free regions, explicit `[n₀,n₁]` certificates | — |
| source | Freitas 2006; Droll 2012 | Sekatskii 2013 |

**This paper is about the second family only.** The question here — which shift best reveals a zero at a given height, with the tested line held at `1/2` — is a detector question, not a zero-free-region question, and it is not the question the `τ`-literature answers.

**[C] The genre, however, is established, and the paper should not be read as opening it.** Explicit detection indices for Li-type coefficients are a developed subject: Brown (2005) related finitely many non-negative coefficients to zero-free regions; Palojärvi (2019) proved two-sided theorems giving explicit intervals `[n₀,n₁]` such that a negative coefficient in that range certifies a zero outside a region, with numerical tables; and Bucur, Ernvall-Hytönen, Odžak and Smajlović (2016) studied the numerical behaviour of `τ`-Li coefficients for functions that violate RH — which is, numerically, the same kind of experiment as §3.3 here. Those results are rigorous with explicit constants; what follows is a closed-form optimisation and a measurement, and is much smaller in scope.

**In particular, the heuristic used in §3.3 — that the crossing index scales inversely with the exponential rate — is the informal version of what Palojärvi's Theorem 3.3 makes precise**, her `n₀` and `n₁` depending on `log R`. The estimate here is not a substitute for that and does not carry its constants.

**[C] How far the prior-art search went.** Four passes: the classical sources; the `τ`-Li and zero-free-region literature of 2005–2026; Sekatskii's own continuation, which reached 2021 and gave the family its name; and a targeted sweep for parameter optimisation in any Li-type family, which returned nothing. Not consulted: MathSciNet or zbMATH, Droll's thesis, and the full texts of Mazhouda–Sodaïgui (2022) and the Ernvall-Hytönen–Odžak–Smajlović–Zubača chapter (2024), all three of which are the most likely places for a statement of the optimum to be sitting. A reader who finds it there is asked to say so.

**[O] That also names the one edit that would make this paper more than a measurement.** Palojärvi's argument produces explicit intervals `[n₀,n₁]` for the `τ` family from the zero-counting bounds and a bound on the far-zero contribution. The same argument structure should carry to the Li–Sekatskii family, since the two differ only in where the foci sit; the far-zero estimate would be taken at `q_a` rather than `ρ/(ρ−τ)`, and the rate `artanh(β/R)` of Proposition 2 would enter where `log R` enters hers. Carrying that through would turn "the matched lens reaches the evidence at an index smaller by `T`" from a scaling statement into a certified interval. It is not attempted here.

### What is offered as new

One statement and one measurement, both in §3:

1. **[P]** within Sekatskii's shifted family, the closed form for the detection rate `log|q_a(ρ)|` of an off-line zero, its exact maximum at `d = |ρ − 1/2|` with value `artanh(β/|ρ − 1/2|)`, and the resulting factor-`T` gain over the classical lens. The optimisation of the shift against a target height does not appear in the sources of §1, but the surrounding genre does — see §1.1;
2. **[M]** the numerical confirmation of that law, including the observation that a previously measured empirical improvement of `30.9×` for a planted zero at height `30` is the law's `T = 30`;
3. **[P]** the two ingredients a Palojärvi-type detection interval needs, worked out for this family (§3.4): the on-line terms are `1 − cos(nθ)` with `θ = 2arctan(γ/d) − π`, and the far zeros contribute a positive `n²d²` term. The uniform estimates that would make it a certified interval are not supplied.

---

## 2. The construction [K]

### 2.1 The map

For real `a ≠ 1/2` the Möbius map `q_a(ρ) = (ρ−a)/(ρ−(1−a))` satisfies

```
|q_a(ρ)| = 1   ⟺   Re ρ = 1/2,
```

because `|ρ−a| = |ρ−(1−a)|` says `ρ` is equidistant from the real points `a` and `1−a`, whose perpendicular bisector is the line `Re s = 1/2`. The functional equation gives `q_a(1−ρ) = 1/q_a(ρ)`, so the image multiset is closed under `z ↦ 1/z` and **the zero sum must be read with symmetric summation**, never term by term.

### 2.2 The two objects

```
S_n(a) = Σ_ρ [ 1 − q_a(ρ)^n ]
L_n(a) = 1/(n−1)! · dⁿ/dsⁿ [ (s−a)^(n−1) log ξ(s) ] |_(s = 1−a)
```

At `a = 0` these are Li's `λ_n`.

**Proposition 1 [K].** `S_n(a) = (1 − 2a) L_n(a)`.

*Verified.* At `n = 1` this is the Hadamard identity: `1 − q_a(ρ) = (2a−1)/(ρ−1+a)`, so `S_1(a) = (2a−1) Σ_ρ 1/(ρ−1+a) = (1−2a) · ξ′/ξ(1−a) = (1−2a) L_1(a)`. Numerically at `n = 1, 2` and `a = 0, −0.5, −1.5`, the zero sum converges monotonically to the derivative side as the zero count grows through 100, 200, 400; at `a = 0, n = 1` the target is `0.02309569`, against the classical `λ₁ = 1 + γ/2 − ½log(4π) = 0.0230957089661`.

**[C] A numerical trap in evaluating `L_n`.** `log ξ` is regular at `s = 1`, but its natural expression is not: `log(s−1)` and `log ζ(s)` are each singular there and cancel. Computing `(s−1)ζ(s)` as one quantity fixes the value but not the derivative, because an automatic differentiator's default step lands at `s − 1 ∼ 10⁻⁴⁰`, where `ζ ∼ 10⁴⁰` and the product loses every digit. With `mpmath.diff` at default settings this returns `−0.554` for `λ₁`. An explicit step of `10⁻²` recovers `0.0230957089658`, ten digits.

### 2.3 The prime side

For `Re s > 1` the derivative object is a generalized Laguerre transform of the von Mangoldt function,

```
S_n^prime(a) = −(1 − 2a) Σ_(m≥2) Λ(m) m^(a−1) L¹_(n−1)( (1−2a) log m ),
```

which is Sekatskii's Lemma 1.

**[M] Checked here rather than quoted.** Comparing the `ζ`-part of `L_n(a)` against the sum directly, for `a < 0` so that `1−a > 1` and the Dirichlet series converges, the two sides agree to the sign shown, and the residue is truncation rather than error:

| relative residue `\|LHS + RHS(M)\|/\|LHS\|` | `M = 500` | `2000` | `8000` | `32000` |
|---|---|---|---|---|
| `a = −0.5`, `n = 1` | 5.99×10⁻² | 3.04×10⁻² | 1.56×10⁻² | 8.14×10⁻³ |
| `a = −1.5`, `n = 1` | 3.32×10⁻⁴ | 1.54×10⁻⁴ | 1.32×10⁻⁴ | 1.29×10⁻⁴ |

At `a = −0.5` the residue falls steadily with the prime cutoff. At `a = −1.5` the sum converges fast and the residue settles on `1.3×10⁻⁴`, which is the finite-difference floor of the derivative side, not a disagreement.

**The discriminating test is the behaviour in `M`, not the size at one `M`.** A variant carrying an extra `1/log m` in the coefficient — an easy misreading of the source — gives residues that *stall* instead of converging:

| variant residue | `M = 500` | `2000` | `8000` |
|---|---|---|---|
| `a = −0.5`, `n = 1` | 3.700×10⁻¹ | 3.657×10⁻¹ | 3.639×10⁻¹ |
| `a = −1.5`, `n = 1` | 1.729×10⁻² | 1.732×10⁻² | 1.732×10⁻² |

At `a = −1.5` the wrong coefficient is only `1.7%` off and would pass a single-cutoff check; it is the flat column, against a correct residue falling to `1.3×10⁻⁴`, that excludes it.

The tail criterion `limsup_n |S_n(a)|^(1/n) ≤ 1` excludes off-line zeros, and obtaining it from this side requires cancellation inside the Laguerre sum of essentially RH strength.

---

## 3. The lens law

### 3.1 The rate

**Proposition 2 [P].** Let `ρ = 1/2 + β + iT` with `β > 0`, let `d = 1/2 − a` be the lens width, and put `R = |ρ − 1/2| = √(β² + T²)`. Then

```
log |q_a(ρ)| = ½ · log[ ((β+d)² + T²) / ((β−d)² + T²) ],
```

this has a unique stationary point in `d > 0`, at

```
d = R,
```

and it is a maximum, with value

```
max_d log |q_a(ρ)| = ½ · log[ (R+β)/(R−β) ] = artanh(β/R).
```

*Proof.* Substituting `a = 1/2 − d` gives `ρ − a = (β+d) + iT` and `ρ − (1−a) = (β−d) + iT`, which is the closed form. Write `N = (β+d)²+T²` and `D = (β−d)²+T²`. Then `(N/D)′ = 2[(β+d)D + (β−d)N]/D²`, and with `u = β+d`, `v = β−d` the bracket is `u(v²+T²) + v(u²+T²) = (u+v)(uv + T²) = 2β(β² − d² + T²)`. Since `β > 0` this vanishes exactly at `d² = β² + T²`, changing sign from `+` to `−`. At `d = R` one has `N = 2R(R+β)` and `D = 2R(R−β)`, giving the stated value. ∎

**[M] Every line of that proof is machine-checked.** `verify_li_lens.py` check 3a verifies the substitution, the form of `|q_a|²`, the derivative's numerator, the factorisation through `u = β+d`, `v = β−d`, its collapse to `2β(β²−d²+T²)`, the values of `N` and `D` at `d = R`, the uniqueness of the root and the sign change — eleven symbolic identities, each with residual exactly zero, plus the sign change numerically. The final `artanh` step is checked as three separate algebraic facts rather than one transcendental identity, because `artanh(z) = ½log((1+z)/(1−z))` holds only for `|z| < 1`: the script verifies that `N/D` at `d = R` is `(R+β)/(R−β)`, that `(x−1)/(x+1) = β/R` for that `x`, and that `R² − β² = T² > 0` so the hypothesis `β/R < 1` is met. The domain condition was missing from the first draft of the proof and the check is what surfaced it.

**The geometric reading.** The map `q_a` has its two foci at `a` and `1−a`, placed symmetrically about `1/2` at distance `d`. The proposition says the detector is sharpest when **the foci stand at the same distance from `1/2` as the zero does**.

**Corollary 3 [P].** The gain over the classical lens is exactly

```
G(β, T)  =  artanh(β/R)  /  ½ log[ ((β+½)² + T²) / ((β−½)² + T²) ],       R = √(β²+T²),
```

and `G → T` as `β → 0` with `T` fixed: the numerator is `β/T + O(β³/T³)` and the denominator is `β/(¼+T²) + O(β³) = β/T² + O(T⁻⁴)`.

Since the index at which an off-line zero first drives the sum negative scales inversely with the exponential rate, **the detection threshold improves by a factor of `T`.**

### 3.2 Measured

**[M]** The maximiser, located by direct search over `d ∈ [0.05, 5T]` at 4×10⁵ points:

| `T`, `β` | 10, 0.1 | 30, 0.1 | 100, 0.1 | 30, 0.01 | 500, 0.01 |
|---|---|---|---|---|---|
| `argmax_d` | 10.00049 | 30.00011 | 100.00013 | 29.99993 | 500.00000 |
| `√(β²+T²)` | 10.00050 | 30.00017 | 100.00005 | 30.00000 | 500.00000 |
| max rate | 9.999833e−3 | 3.333327e−3 | 9.999998e−4 | 3.333333e−4 | 2.000000e−5 |
| `artanh(β/R)` | 9.999833e−3 | 3.333327e−3 | 9.999998e−4 | 3.333333e−4 | 2.000000e−5 |

**[M]** The gain, as the ratio of rates:

| `T` | 10 | 30 | 100 | 500 |
|---|---|---|---|---|
| rate at `d = T` | 1.000×10⁻² | 3.333×10⁻³ | 1.000×10⁻³ | 2.000×10⁻⁴ |
| rate at `d = 1/2` | 9.974×10⁻⁴ | 1.111×10⁻⁴ | 1.000×10⁻⁵ | 4.000×10⁻⁷ |
| gain | **10.026** | **30.009** | **100.003** | **500.001** |

with `β = 0.1`; the `β = 0.01` row gives the same gains to three decimals, as the `O(β³)` term requires.

### 3.3 A previously measured improvement, explained

**[M]** An earlier experiment planted a symmetric quartet at `ρ = 1/2 ± 0.1 ± 30i` among verified critical-line zeros and located the first negative index: `n ≈ 41646` at `a = 0`, and `n ≈ 1348` with a lens of width `d ≈ 30`, an improvement of `30.89×`. Proposition 2 predicts the ratio of rates,

```
(3.333327×10⁻³) / (1.110790×10⁻⁴) = 30.01,
```

against the measured `30.89`.

**[C] The `3%` residue is not explained here, only attributed.** The crossing index is not exactly inversely proportional to the rate, because the on-line zeros contribute a background growing like `(n/2) log n` (Sekatskii, arXiv:1404.7276, Theorem 6) rather than a constant. That accounts for a residue of this sign and rough size, but no model of the crossing was fitted and the `3%` was not predicted in advance.

**The empirical `30.9` was the height of the planted zero.**


### 3.4 Why the threshold scales like `1/log R`, and what is still missing

§3.3 used the rule that the crossing index scales inversely with the rate. That rule is not assumed here; the two ingredients Palojärvi's argument needs are available for this family too.

**Proposition 4 [P].** For a zero **on** the critical line, `ρ = 1/2 + iγ`, and lens width `d`,

```
q_a(ρ) = e^{iθ},    θ = 2 arctan(γ/d) − π,    so   Re(1 − q_a(ρ)ⁿ) = 1 − cos(nθ) ≥ 0,
```

and `θ = −2 arctan(d/γ) = −2d/γ + O(d³/γ³)`. Hence for `γ ≫ nd` the term is `2n²d²/γ² + O(·)`.

*Checked:* `|q_a| − 1 = 1.1×10⁻¹⁶` on the line, and `|Re(1 − q_a^n) − (1 − cos nθ)| ≤ 1.7×10⁻¹⁴` over `d = 0.5, 30`, `γ = 50, 500, 5000`, `n = 3, 100`.

**[M] The far zeros therefore contribute a positive term growing like `n²d²`.** Summing `2n²d²/γ²` against the Riemann–von Mangoldt density gives

```
Σ_{γ > T} Re(1 − q_aⁿ)  ≈  (2n²d²/π) · (log(T/2π) + 1) / T,
```

measured against the exact sum over 2×10⁶ ordinates:

| `d`, `n` | `T = 5×10³` | `2×10⁴` | `10⁵` |
|---|---|---|---|
| `0.5`, `20` | 0.993 | 0.974 | 0.892 |
| `0.5`, `50` | 0.993 | 0.974 | 0.892 |
| `30`, `50` | — | — | 0.892 |

as ratios of exact to model. The shortfall at `T = 10⁵` is the finite zero list — the model integrates to infinity while the list stops at `1.13×10⁶`.

**The structure of the argument then closes.** The far term is the analogue of Palojärvi's `K_{F,3}(T,τ) n(n−1)`; for the near zeros, an off-line zero with `|q_a(ρ)| ≥ R` gives, by the same lemma of Montgomery she uses,

```
Σ_{|γ| ≤ T} Re(1 − q_aⁿ)  <  N(T) − (1/20) Rⁿ
```

for some `n` in an explicit range. Combining, `S_n(a)` is driven negative once

```
Rⁿ  ≳  20 [ N(T) + (2n²d²/π)(log(T/2π)+1)/T ],       i.e.   n  ≳  log(background) / log R.
```

**So the detection index is inversely proportional to `log R` up to a logarithmic correction, and Proposition 2 — which maximises `log R` over the lens — minimises it.** That is the derivation §3.3 was missing.

**[C] What this is not.** It is the argument's *structure*, not a theorem. Palojärvi's intervals come with explicit constants obtained from uniform estimates on every term, and none of that uniformity is established here: the `O(d³/γ³)` remainder is not bounded uniformly in `n`, the range of `n` for which Montgomery's lemma applies is not made explicit, and the background's own dependence on `d` — which enters through `n²d²` — is not carried through. That last omission is very likely where the `3%` residue of §3.3 lives, since the two lenses there differ in `d` by a factor of `60`. Supplying the uniform estimates would turn the factor-`T` gain from a scaling law into a certified interval, and is the concrete next step named in §1.1.


---

## 4. The truncation artifact [K]

**[C]** Truncating the prime-Laguerre sum at `m ≤ M` produces apparent exponential growth in `n` even when no off-line zero is present. Measured `n`-th root growth: `3` at `a = −0.5`, `5/3` at `a = −1.5`.

These are not empirical constants. With `d = 1/2 − a`,

```
(d + 1/2)/(d − 1/2) = (a−1)/a = (b+1)/b,      b = −a,
```

giving `3`, `5/3`, `9/7` at `a = −0.5, −1.5, −3.5` — exactly the term `(1 + 1/b)ⁿ` in Sekatskii's eq. (7), which he traces to the formal zero of `ξ` at `s = 1` and shows is cancelled by the pole of `ζ` once the prime sum is treated to first order.

**So the artifact has a known analytic source and a known cancellation.** Raising the cutoff postpones it and does not remove it; only doing the regrouping before truncating does. Any numerical search for off-line zeros through this criterion must subtract this term first, or it will find one.

---

## 5. What the lens does not buy

The lens law is a statement about **detector sensitivity**, not about provability. Three things follow, and none of them is a route to a proof.

1. **The criterion is unchanged.** For every admissible `a`, non-negativity of `S_n(a)` for all `n` is equivalent to RH — that is Sekatskii's theorem. Choosing `a` well changes which `n` first reveals a violation; it does not change what has to be proved.
2. **The hard direction is untouched.** Establishing `limsup |S_n(a)|^{1/n} ≤ 1` from the prime side requires cancellation in the Laguerre transform of von Mangoldt of essentially RH strength, at every `a`. Moving the lens changes the arithmetic weights `m^{a−1}` and the polynomial argument `(1−2a) log m`; it does not make the cancellation easier to exhibit.
3. **The gain is bounded by what one is looking for.** A factor of `T` helps only if one already knows roughly where to look. Scanning `T` costs back what matching to `T` gains.

**Where it is genuinely useful:** as a targeted test. Given a specific height at which one suspects a violation, the matched lens reaches the same evidence at an index smaller by `T`, which at `T ∼ 10³` is the difference between a feasible computation and an infeasible one. That is a numerical instrument, and it is worth having as one.

> **A more sensitive detector is not a shorter proof.**


---

## 6. The same question for the Weil form, and why the answer differs

The shifted Li criterion is not the only finite window on the explicit formula. The truncated Weil quadratic form — Connes–van Suijlekom, and the Galerkin truncations studied numerically by Groskin — is another, and Sekatskii notes the connection outright: his arithmetic interpretation is obtained from Weil's explicit formula, and he cites Bombieri's *Remarks on Weil's quadratic functional*. So the two constructions are two finite windows on one identity, and it is fair to ask the same question of both.

**In the Weil setting the free parameter is the support length `L`.** A test function supported in `[−L/2, L/2]` has a transform that extends to an entire function, and an off-line zero at `ρ = 1/2 + β + iT` is probed at the complex point `T − iβ`.

**[K]** By Paley–Wiener, for any `f` supported in `[−L/2, L/2]`,

```
|f̂(T − iβ)|  ≤  e^(βL/2) · ‖f‖₁,
```

and the factor is sharp. So the greatest amplification the support length permits is `βL/2` in the log.

**[P] That ceiling has no interior maximum.** `βL/2` is strictly increasing in `L`. There is no optimal support length; longer is always better.

**[M]** The Dirichlet sine modes sit comfortably inside the ceiling and do not saturate it: at `T = 30`, the ratio `|φ̂_j(T−iβ)| / (‖φ_j‖₁ e^{βL/2})` for the mode nearest resonance is

| `L` | 4.5 | 9 | 18 | 36 |
|---|---|---|---|---|
| `β = 0.1` | 0.632 | 0.517 | 0.363 | 0.210 |
| `β = 0.5` | 0.312 | 0.173 | 0.087 | 0.044 |

against ceilings `e^{βL/2}` of `1.25, 1.57, 2.46, 6.05` and `3.08, 9.49, 90.0, 8103`.

### The contrast

| | shifted Li | truncated Weil |
|---|---|---|
| free parameter | lens width `d = 1/2 − a` | support length `L` |
| detector rate | `artanh(β/R)`, `R = \|ρ − 1/2\|` | ceiling `βL/2` |
| behaviour in the parameter | **interior maximum at `d = R`** | **monotone, no maximum** |
| cost of moving it | none — `a` is free | on the CvS path `L = log c`, so a longer support forces a larger prime sum |

**[C] The two rows are not the same kind of object and should not be read as one.** The Li entry is the *actual* rate for a given zero, determined by the parameter alone. The Weil entry is a *ceiling* over all test functions of that support, not the rate of any particular one. The comparison is between what each parameter can buy, not between two measured rates.

**The structural reading.** The Li parameter is a *position*: `q_a` has foci at `a` and `1−a`, and there is a matching condition because a position can be right or wrong relative to the target. The Weil parameter is a *bandwidth*: more is monotonically better, and it is bounded not by an optimum but by a price. On the coupled path that price is the prime sum, which is exactly the coupling `L = log c` that Groskin identifies as structural for the reference construction.

**[O]** Whether the Weil form has an optimum in a parameter that *is* a position — the prime cutoff held independent of the support, say — is not settled here. Off the saturated diagonal that family is strongly indefinite, so "detection" would first have to be given a meaning that does not reduce to prime deletion.

---

## 7. Reproducibility

```bash
python code/verify_li_lens.py        # Propositions 1-3, the prime formula, and the tables of §3, §4, §6
```

`§2.2` needs an explicit differentiation step (see the caveat there); the script sets it and fails loudly if the recovered `λ₁` does not match the closed form. Tested with Python 3.12, NumPy 2.4, SymPy 1.14, mpmath 1.3.

---

## AI assistance

The verification script and much of the prose in this note were written with the assistance of Claude (Anthropic); the research direction, the decisions about what to publish, and responsibility for every claim are the author's. See the repository README for a fuller statement.

---

## References

1. X.-J. Li, *The positivity of a sequence of numbers and the Riemann hypothesis*, J. Number Theory **65** (1997) 325–333.
2. E. Bombieri and J. C. Lagarias, *Complements to Li's criterion for the Riemann hypothesis*, J. Number Theory **77** (1999) 274–287.
3. S. K. Sekatskii, *Generalized Bombieri–Lagarias' theorem and generalized Li's criterion with its arithmetic interpretation*, Ukrainian Math. J. **66** (2014) 415–431; arXiv:1304.7895.
4. S. K. Sekatskii, *An arithmetic interpretation of generalized Li's criterion*, arXiv:1305.1421.
5. S. K. Sekatskii, *Asymptotic of the generalized Li's sums which non-negativity is equivalent to the Riemann Hypothesis*, arXiv:1403.4484.
6. S. K. Sekatskii, *First applications of generalized Li's criterion to study the Riemann zeta-function zeroes location*, arXiv:1404.7276.
7. A. Voros, *Sharpenings of Li's criterion for the Riemann Hypothesis*, Math. Phys. Anal. Geom. **9** (2006) 53–63.
8. M. W. Coffey, *Toward verification of the Riemann Hypothesis: application of the Li criterion*, Math. Phys. Anal. Geom. **8** (2005) 211–255.
9. J. C. Lagarias, *Li coefficients for automorphic L-functions*, Ann. Inst. Fourier **57** (2007) 1689–1740.
10. K. Maślanka, *Li's criterion for the Riemann hypothesis — numerical approach*, Opuscula Math. **24** (2004).
11. L. Smajlović, *On Li's criterion for the Riemann hypothesis for the Selberg class*, J. Number Theory **130** (2010) 828–851.
12. P. Freitas, *A Li-type criterion for zero-free half planes of Riemann's zeta function*, J. London Math. Soc. **73** (2006) 399–414.
13. A. D. Droll, *Variations of Li's criterion for an extension of the Selberg class*, PhD thesis, Queen's University, Ontario, 2012.
14. F. C. Brown, *Li's criterion and zero-free regions of L-functions*, J. Number Theory **111** (2005) 1–32.
15. A. Bucur, A.-M. Ernvall-Hytönen, A. Odžak and L. Smajlović, *On a Li-type criterion for zero-free regions of certain Dirichlet series with real coefficients*, LMS J. Comput. Math. **19** (2016) 259–280.
16. N. Palojärvi, *Explicit zero-free regions and a τ-Li-type criterion*, arXiv:1807.01506.
17. A.-M. Ernvall-Hytönen, A. Odžak, L. Smajlović and M. Zubača, *Variants of the Li-Type Criteria for the Generalized Riemann Hypothesis*, in *Women in Numbers Europe IV*, Springer (2024).
18. A. Voros, *Discretized Keiper/Li approach to the Riemann hypothesis*, Exp. Math. **29** (2020) 452–469.
19. K. Mazhouda and B. Sodaïgui, *The Li–Sekatskii coefficients for the Selberg class*, Int. J. Math. **33** (2022), no. 12, Paper No. 2250075.
20. S. K. Sekatskii, *On the Generalized Li's Criterion Equivalent to the Riemann Hypothesis and Its First Applications*, in *Schrödinger Operators, Spectral Analysis and Number Theory*, Springer Proc. Math. Stat. **348** (2021).

**No claim is made about the Riemann Hypothesis.** The numerical results are exploratory measurements of a detector, not evidence about the location of any zero.
