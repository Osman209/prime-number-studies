# When the Lens Goes Blind

### The phase cycle of the Li–Sekatskii detector, and the price of a shift taken too far

July 2026

---

## Summary

A companion note, [`li_lens_law.md`](li_lens_law.md), found that the shift in the
**Li–Sekatskii coefficients** behaves as a lens with an exact optimal width: for a
hypothetical zero at $\rho = 1/2 + \beta + iT$, writing $d = 1/2 - a$ and $R = |\rho - 1/2|$,
the exponential rate $r(d) = \mathrm{log}|q_a(\rho )|$ is maximal at $d = R$.

That is a statement about the **modulus** of a complex number. This note takes the
**phase**, and then follows the rate to the far end where **Sekatskii's Theorem 3**
lives.

Three things come out of the phase, and two of them are negative.

**The rate has a closed form and a symmetry.**

```
tanh r(d) = 2βd / (d² + β² + T²),
```

which is invariant under $d \mapsto R^{2}/d$. So every sub-maximal rate is attained at *two*
lens widths whose product is $R^{2}$, and $d = R$ is that map's only fixed point — the
optimum is a centre of symmetry, not merely a stationary point. In the coordinate
$x = \mathrm{log}(d/R)$ the rate is $\mathrm{artanh}((\beta /R) \mathrm{sech} x)$: a **sech dome on the logarithmic
axis**, so the natural distance from the optimum is the ratio $d/R$, not $d - R$.

**At the optimum the phase is exactly $- \pi /2$,** independent of $\beta$ and $T$. The
quartet's contribution therefore runs on a strict cycle of period four — strong
detection at $n \equiv 0 (\mathrm{mod} 4)$, strong **masking** at $n \equiv 2$, and no exponential term
at all for odd $n$. **In 24 measured runs the first negative index of the full
coefficient was $\equiv 0 (\mathrm{mod} 4)$ every time.**

**And the far end is Sekatskii's.** As $d \to \infty$ the image $q_a(\rho ) \to - 1$, every $C_n$
turns positive, and the detection index grows **linearly** in the lens width. That
is his Theorem 3 seen from the rate curve: $r(d)$ vanishes at both ends and peaks at
$R$, so a shift taken far past the optimum does not merely fail to help — it hides.
Temporarily: for any fixed shift the rate stays positive and the violation returns.

**The two negative results.** Lenses tuned to align the phase at a chosen order do
**not** improve the first detection index; and the far detection pockets, though
real, are swamped by the background by seven orders of magnitude. Both look
promising and neither works, which is why they are here.

---

## Status legend

| tag | meaning |
|---|---|
| **[K]** | Known. Source named in §1. No novelty claimed. |
| **[P]** | Proved here, and machine-checked step by step. |
| **[M]** | Measured, with parameters stated. |
| **[C]** | Caveat. |
| **[O]** | Open. |

---

## 1. Prior art, and what this note is not

| component | where it lives |
|---|---|
| the shifted map, the derivative object, and the equality between them | the **Li–Sekatskii coefficients** — Sekatskii, arXiv:1304.7895 = *Ukrainian Math. J.* **66** (2014) 415–431 |
| positivity of $\lambda _n$ ⟺ RH, and the multiset theorem behind it | **Li** (1997); **Bombieri & Lagarias** (1999) |
| **the masking statement, proved** | **Sekatskii**, *On the Generalized Li's Criterion…*, Springer Proc. Math. Stat. **348** (2021) 241–254, **Theorem 3**: for any $m$ there is $c$ such that $\lambda _{n,b,1/2} = \Sigma _\rho (1 - ((\rho +b)/(\rho - b- 1))^{n}) \gt 0$ for all $n \le m$ and all $b \ge c$ |
| the exact optimum $d = R$ and the factor-$T$ gain | the companion note, [`li_lens_law.md`](li_lens_law.md) |
| explicit detection indices for Li-type coefficients | **Brown** (2005); **Bucur, Ernvall-Hytönen, Odžak & Smajlović** (2016); **Palojärvi** (2019) |

**[C] The relation to Theorem 3 must be stated carefully, because it is easy to
overstate.** Theorem 3 is proved for $\xi$ itself, and Sekatskii confirms
(correspondence, August 2026) that it is obtained through the analytic formula for
$\xi ′/\xi$ — so the proof uses what $\xi$ is, not only the shape of the sum. What he
describes as *apparently* claimable is the corresponding statement for general
multisets containing off-line elements, which is a separate matter.
**Everything measured below is the multiset case** — a quartet planted at
$1/2 \pm \beta \pm iT$ among genuine critical-line ordinates — so it bears on that
extension, not on the theorem as proved. With $b = d - 1/2$ his $\lambda _{n,b,1/2}$ is
exactly the sum computed here.

### What is offered as new

1. **[P]** the closed form $\mathrm{tanh} r(d) = 2\beta d/(d^{2} + \beta ^{2} + T^{2})$, its exact invariance
   under $d \mapsto R^{2}/d$, and the `sech` profile on the logarithmic axis;
2. **[P]** the phase identity $\varphi (R) = - \pi /2$ and the resulting $n \mathrm{mod} 4$ cycle;
3. **[M]** that the cycle governs the **full** coefficient, not only the isolated
   quartet — 24 of 24 first-negative indices $\equiv 0 (\mathrm{mod} 4)$;
4. **[M]** the scaling $N_M(d)\cdot r(d) = \mathrm{log} M + O(1)$ on 100,000 ordinates, the U
   shape with its minimum at $d = R$, and linear growth of the index in the lens
   width;
5. **[M]** two negative results, §6.

---

## 2. The rate, closed

Write $\rho = 1/2 + \beta + iT$ with $\beta \gt 0$, $d = 1/2 - a$, and $R = \sqrt (\beta ^{2} + T^{2}) = |\rho - 1/2|$.

**Proposition 1 [P].** $\mathrm{tanh} r(d) = 2\beta d / (d^{2} + \beta ^{2} + T^{2})$.

*Proof.* $e^{2r} = N/D$ with $N = (\beta +d)^{2}+T^{2}$ and $D = (\beta - d)^{2}+T^{2}$, so
$\mathrm{tanh} r = (N- D)/(N+D)$. Now $N - D = 4\beta d$ and $N + D = 2(\beta ^{2} + d^{2} + T^{2})$. ∎

**[C] The $\beta \ll T$ form is different and is not the identity.** Dropping $\beta ^{2}$ gives
$2\beta d/(d^{2}+T^{2})$, which agrees only to the order shown: at $d = 30$, $\beta = 0.1$,
$T = 30$, the exact value is $0.0033333148$ against $0.0033333333$.

**Proposition 2 [P].** The rate is invariant under $d \mapsto R^{2}/d$, whose only fixed
point in $d \gt 0$ is $d = R$.

*Proof.* Substituting $R^{2}/d$ into $2\beta d/(d^{2}+R^{2})$ returns it unchanged; $d = R^{2}/d$
gives $d = R$. ∎

So the sensitivity landscape is symmetric about the optimum **multiplicatively**.
A lens $k$ times too narrow and one $k$ times too wide have the *same* rate.

**Corollary 3 [P].** With $x = \mathrm{log}(d/R)$,

```
r(x) = artanh( (β/R) · sech x ),
```

even in $x$, peaked at $x = 0$. **The peak is broad on the ratio scale**: for
$\beta \ll T$ one keeps half the maximal rate over $0.268 R \lesssim d \lesssim 3.73 R$. Matching the
lens matters; matching it precisely does not.

---

## 3. The phase

### 3.1 The value at the optimum

Rationalising, $q_a(\rho ) = [(R^{2} - d^{2}) - 2idT] / [(\beta - d)^{2} + T^{2}]$.

**Proposition 4 [P].** At $d = R$ the real part of the numerator vanishes
identically, so $q_a(\rho ) = - i|q_a(\rho )|$ and

```
φ(R) = −π/2,
```

**independent of $\beta$ and $T$.**

**[M]** Checked at $(\beta ,T) = (0.1, 30), (0.01, 500), (0.4, 10), (0.001, 3)$: the
computed phase equals $- \pi /2$ to $10^{-14}$ in every case.

### 3.2 The quartet, and the cycle of four

The symmetric quartet $\rho , \bar{\rho}, 1- \rho , 1- \bar{\rho}$ maps to $e^{\pm r \pm i\varphi }$, and since
$q_a(1- \rho ) = 1/q_a(\rho )$,

**Proposition 5 [P].** $C_n(d) = 4 - 4 \mathrm{cosh}(n r(d)) \mathrm{cos}(n \varphi (d))$.

*Checked* against the direct sum over the four zeros to $5\times 10^{-14}$.

At $d = R$ this becomes $4 - 4 \mathrm{cosh}(n r) \mathrm{cos}(n\pi /2)$, and the cycle is exact:

| $n \mathrm{mod} 4$ | $\mathrm{cos}(n\varphi )$ | $C_n(R)$ | |
|---|---|---|---|
| $0$ | $+1$ | $4 - 4\mathrm{cosh}(nr) \lt 0$ | **detection** |
| $1$ | $0$ | $= 4$ exactly | no exponential term |
| $2$ | $- 1$ | $4 + 4\mathrm{cosh}(nr) \gt 0$ | **masking** |
| $3$ | $0$ | $= 4$ exactly | no exponential term |

**[M]** At $\beta = 0.1$, $T = 30$: $C_{100}(R) = - 0.224287$, $C_{102}(R) = +8.233435$, and
$C_{101}(R) = C_{103}(R) = 4.000000$.

**So the lens that maximises the envelope fails at three orders in four, in two
different ways.** At $n \equiv 2 (\mathrm{mod} 4)$ — one order in four — the off-line quartet
contributes a large *positive* amount and actively conceals itself, which is worse
than useless. At odd $n$ — half the orders — its exponential growth simply does not
enter the real part, so the quartet is invisible rather than misleading. Only
$n \equiv 0 (\mathrm{mod} 4)$ detects.

### 3.3 The cycle survives into the full coefficient

The above is arithmetic on four planted points. The question is whether it governs
$S_n = \Sigma _\rho (1 - q_a(\rho )^{n})$ with real ordinates present.

**[M] It does.** Over 24 runs — six heights $T = 10, 20, 30, 50, 80, 120$, two
displacements $\beta = 0.1, 0.05$, and two ordinate counts $2000, 8000$ — the first
index at which the full coefficient goes negative was $\equiv 0 (\mathrm{mod} 4)$ in **24 of 24**.
At the classical lens $d = 1/2$, where $\varphi \ne - \pi /2$, the same indices fall across all
four residues.

Under the null hypothesis of no phase structure the probability is $4^{-24} \approx 3\times 10^{-15}$.

---

## 4. The far end, and Theorem 3

### 4.1 Why a large lens hides

**[P]** With $z = \beta + iT$, $q_a(\rho ) = - (1 + z/d)/(1 - z/d) = - \mathrm{exp}(2z/d + O(d^{-3}))$, so

```
r(d) = 2β/d + O(d⁻³),        φ(d) = −π + 2T/d + O(d⁻³),
```

and $q_a(\rho ) \to - 1$. Hence $C_n \to 4(1 - (- 1)^{n})$: $0$ for even $n$, $8$ for odd. The
next order decides the sign:

```
n even:  C_n = 8n²(T² − β²)/d² + O(d⁻⁴)          → positive when T > |β|
n odd:   C_n = 8 + 8n²(β² − T²)/d² + O(d⁻⁴)      → 8
```

**[M]** Ratios of exact to model at $n = 2, 4, 10, 50$: at $d = 2\times 10^{5}$ they are
$1.0000$ throughout; at $d = 2\times 10^{4}$ they are $1.0000$, $1.0000$, $0.9999$,
$0.9981$ — the model is an expansion in $nd^{-1}$, so the departure grows with $n$ and
shrinks with $d$, which is what the two columns show. At $d = 10^{7}$: $C_{4} = 1.15\times 10^{-9}$, $C_{5} = 8.000000$.

**So for $T \gt |\beta |$ every $C_n$ is positive once $d$ is large enough**, and since
on-line zeros contribute $2[1 - \mathrm{cos} n\theta ] \ge 0$, the whole coefficient is non-negative
for the first $m$ orders. That is the mechanism of Theorem 3.

### 4.2 The order of the quantifiers

**[C]** Fix $d$ and let $n \to \infty$: $|q_a(\rho )| \gt 1$ whenever $\beta \gt 0$, so the violation
appears. Fix $n$ and let $d \to \infty$: the contribution is non-negative. The limits do
not commute, and the statement is

```
∀m ∃d_m ,      not      ∃d ∀m .
```

**Masking is temporary.** No finite shift hides a zero forever, which is what the
criterion's equivalence to RH requires.

### 4.3 The cost, measured

**[M]** With 2000 ordinates and the quartet at $T = 30$:

| $d/R$ | 1 | 2 | 10 | 50 | 200 |
|---|---|---|---|---|---|
| $N(d)$ | 2276 | 2863 | 11536 | 56947 | **227462** |
| $N(d)/(d/R)$ | 2276 | 1432 | 1154 | 1139 | **1137** |

**The index grows linearly in the lens width**, the ratio settling by $d = 50R$. So
pushing detection past order $m$ costs a shift proportional to $m$, which is the
scale $O(mR)$.

### 4.4 The U, and the scaling law

**[M]** Over $d/R = 0.005$ to $20$ at $M = 2000$:

| $d/R$ | 0.005 | 0.1 | 0.6 | **1** | 2 | 5 | 20 |
|---|---|---|---|---|---|---|---|
| $N(d)$ | 229947 | 11568 | 2581 | **2276** | 2863 | 5936 | 22857 |

a clean U with its floor at $d = R$. And at the matched lens, over ten values of
$M$ from 250 to 100,000:

| $M$ | 250 | 1000 | 4000 | 16000 | 64000 | 100000 |
|---|---|---|---|---|---|---|
| $N_M(R)$ | 1660 | 2072 | 2500 | 2880 | 3296 | 3520 |
| $N\cdot r - \mathrm{log} M$ | +0.012 | −0.001 | +0.039 | −0.080 | −0.080 | **+0.220** |

```
N_M(d) · r(d)  =  log M + O(1).
```

**[C] The $O(1)$ is not $o(1)$.** The residual stays bounded — the full range over
the ten values is $[- 0.080, +0.220]$ — but it does not settle toward zero, and the
last point is three times any earlier one. The law is an envelope with a bounded,
non-vanishing correction, and should not be written with $o(1)$.

---

## 5. The symmetry, in the index

Proposition 2 is an exact statement about the rate. The detection index also
depends on the background, which has no reason to respect $d \mapsto R^{2}/d$.

**[M] It nearly does.** Comparing $N(0.1R)$ with $N(10R)$:

| $M$ | 2000 | 4000 | 8000 | 16000 | 32000 | 48000 | 64000 | 100000 |
|---|---|---|---|---|---|---|---|---|
| ratio | 1.0028 | 0.9936 | 0.9965 | 0.9871 | 1.0029 | 1.0212 | 1.0291 | 1.0217 |

A lens $k$ times too narrow and one $k$ times too wide cost within a few percent of
the same number of coefficients, across a 50-fold range of $M$.

**[C] But the departure is not noise at large $M$.** The sign is positive at
$M = 2000$, negative at $4000$, $8000$ and $16000$, and then positive at **every**
value from $M = 32000$ on — four consecutive values, reaching the $2- 3$% level
from $M = 48000$: the narrow lens costs *more*. A plausible reading is that the two
lenses interact differently with the truncation of the ordinate list, since at
$d = 0.1R$ the topmost ordinates are still phase-coherent while at $d = 10R$ they
have long since begun to oscillate.

**[O] That reading is not supported by the data as it stands.** The natural
candidate variable, $\Lambda = 2N(d)d/\gamma _M$, falls monotonically from $27.6$ to $1.4$
across the range at $d = 0.1R$ and never crosses $1$, while the sign changes twice.
A monotone quantity cannot produce a non-monotone sign, so whatever fixes the sign
is not $\Lambda$ alone.

---

## 6. Two things that do not work

Both look like they should, which is the reason for recording them.

### 6.1 Phase-tuned lenses do not improve detection

Since $r′(R) = 0$ while $\varphi ′(R) = - 1/T$, a small displacement changes the phase to
first order and the rate only to second. It is therefore tempting to displace the
lens by $O(T/n)$ so that $\mathrm{cos}(n\varphi ) = 1$ at a chosen order, losing only $O(n^{-2})$ of
the rate. For a *fixed* $n$ that is correct.

**[M] It does not help the first index.** Scanning $d$ over $[R- 3, R+3]$ with 2000
ordinates at $T = 30$, the eight displaced lenses at
$d - R = - 3, - 1.5, - 0.75, - 0.25, +0.25, +0.75, +1.5, +3$ give, in that order,
$2298, 2295, 2289, 2300, 2276, 2284, 2289, 2293$, against $N = 2276$ at $d = R$
itself. **None is below.** One — at
$d - R = +0.25$ — ties, so the minimum is attained twice and never beaten; every
other lens costs between 8 and 24 coefficients.

**Why.** Minimising $C_n$ at a fixed $n$ is not the problem. The detection index is
a **first passage over all $n$**, and displacing the lens buys alignment at one
order while costing envelope at every order. The refinement also needs $n$ in
advance, which is what a search does not have.

### 6.2 The far detection pockets are swamped

**[P]** At large $d$, $\mathrm{cos}(n\varphi ) = (- 1)^{n} \mathrm{cos}(2nT/d)$, so for even $n$ the phase
realigns whenever $2nT/d = 2\pi k$. This gives a family of narrow **detection pockets**
at

```
d ≈ nT/(πk),        k = 1, 2, 3, …
```

of relative full width $2\beta /T$. **[M]** The measured centre approaches the prediction
as $n$ grows — relative departure $2.1\times 10^{-3}$, $3.4\times 10^{-4}$,
$3.0\times 10^{-5}$ at $n = 40, 100, 400$ — because the minimum of $C_n$ sits where the
$\cosh$ and the $\cos$ balance, not exactly where the phase realigns; and the width is $0.00667$, $0.00667$, $0.00666$ at
$(n,k) = (100,1), (400,1), (400,2)$ against $2\beta /T = 0.00667$.

**[M] Not one of them produces a negative coefficient.**

| $n$, $k$ | 40, 1 | 100, 1 | 400, 1 | 400, 2 |
|---|---|---|---|---|
| $C_n$ in the pocket | −8.70×10⁻⁴ | −8.76×10⁻⁴ | −8.77×10⁻⁴ | −3.51×10⁻³ |
| background $B_n$ | 3744.3029 | 3983.2184 | 3989.5495 | 3988.6781 |
| $S_n = B_n + C_n$ | 3744.3021 | 3983.2175 | 3989.5487 | 3988.6746 |

**Why they never deepen.** Along a pocket line the quartet's contribution is
$C_n \approx - 8n^{2}\beta ^{2}/d^{2} = - 8(\pi k\beta /T)^{2}$, which **does not depend on $n$ at all**. However far
out one goes and however large the order, the pocket has the same depth — about
$8.8\times 10^{-4}$ here — against a background of order the number of ordinates. The ratio
is $2\times 10^{-7}$.

A reader of the far-field analysis would reasonably expect a second family of usable
lenses. There is none.

---

## 7. What this changes, and what it does not

**It changes how a numerical search should be reported.** A search that computes
$m$ coefficients, finds them all positive, and concludes that no violation is
present at that height has established nothing unless it also states its shift.
Theorem 3 says a large enough shift makes the first $m$ positive regardless, §4.3
measures the cost as linear in the shift with threshold $O(mR)$, and §3.2 says that
even at the optimal lens only one order in four carries the signal. **Three separate
reasons why a finite check can come back clean and mean nothing.**

**It does not establish that the optimal shift is worth choosing.** Sekatskii's own
position (correspondence, August 2026) is that the maximum of the modulus, and hence
the choice of $a$, **might** be rather than is important for numerical work, and that
qualification is adopted here. The masking direction is the firmer of the two: it
rests on his Theorem 3, which is proved.

**It does not touch the Riemann Hypothesis**, and it does not make the criterion
easier to prove. Every statement here is about a detector's behaviour on a planted
configuration, in a finite truncation.

> **Knowing when an instrument is blind is not the same as improving it.**

---

## 8. Reproducibility

```bash
ZETA_ZEROS=zeros.txt python code/verify_phase_masking.py          # 13 checks, ~85 s
ZETA_ZEROS=zeros.txt python code/verify_phase_masking.py --fast   # ~3 s
```

Every number above comes from that script, which exits nonzero if any check that
could be run fails.

`ZETA_ZEROS` names a file of critical-line ordinates, one per line, increasing.
The measurements above use the first **100,000**, $\gamma$ from $14.1347$ to
$74920.8275$; a list of that size can be taken from Odlyzko's tables
(`www.dtc.umn.edu/~odlyzko/zeta_tables`) or the LMFDB, or generated with
$mpmath.zetazero$, which is far slower. Without the variable the symbolic and
quartet checks still run and the sweeps are skipped.

**A shorter list does not silently become a longer one.** `numpy` slicing is
quiet — $gam[:100000]$ on a 2000-element list returns 2000 — so every slice in
the script is guarded, and a table needing more ordinates than the list holds is
skipped and listed under $REDUCED RUN$ at the end rather than computed on fewer
and printed under the count it asked for. A short list is not a failure and does
not set the exit code; it is reported so that no printed number carries a label
it did not earn. In particular §5 needs $M = 32000$ before it can say anything
about the sign, and §4.4's reading of the residual needs the full ten-point run.

Tested with Python 3.12, NumPy 2.4, SymPy 1.14.

---

## AI assistance

The verification script and much of the prose in this note were written with the
assistance of Claude (Anthropic); the research direction, the decisions about what
to publish, and responsibility for every claim are the author's. See the repository
README for a fuller statement.

---

## Acknowledgements

I am grateful to **S. K. Sekatskii** for correspondence on the relation between the
measurements here and his Theorem 3, and in particular for confirming that the
theorem is obtained through the analytic formula for $\xi ′/\xi$ and so is a statement
about $\xi$, not about general multisets. Any error that remains is mine.

## References

1. X.-J. Li, *The positivity of a sequence of numbers and the Riemann hypothesis*,
   J. Number Theory **65** (1997) 325–333.
2. E. Bombieri and J. C. Lagarias, *Complements to Li's criterion for the Riemann
   hypothesis*, J. Number Theory **77** (1999) 274–287.
3. S. K. Sekatskii, *Generalized Bombieri–Lagarias' theorem and generalized Li's
   criterion with its arithmetic interpretation*, Ukrainian Math. J. **66** (2014)
   415–431; arXiv:1304.7895.
4. S. K. Sekatskii, *On the Generalized Li's Criterion Equivalent to the Riemann
   Hypothesis and Its First Applications*, in *Schrödinger Operators, Spectral
   Analysis and Number Theory: In Memory of Erik Balslev*, Springer Proc. Math.
   Stat. **348** (2021) 241–254; doi:10.1007/978-3-030-68490-7.
5. S. K. Sekatskii, *First applications of generalized Li's criterion to study the
   Riemann zeta-function zeroes location*, arXiv:1404.7276.
6. F. C. Brown, *Li's criterion and zero-free regions of L-functions*, J. Number
   Theory **111** (2005) 1–32.
7. A. Bucur, A.-M. Ernvall-Hytönen, A. Odžak and L. Smajlović, *On a Li-type
   criterion for zero-free regions of certain Dirichlet series with real
   coefficients*, LMS J. Comput. Math. **19** (2016) 259–280.
8. N. Palojärvi, *Explicit zero-free regions and a τ-Li-type criterion*,
   arXiv:1807.01506.
9. K. Mazhouda and B. Sodaïgui, *The Li–Sekatskii coefficients for the Selberg
   class*, Int. J. Math. **33** (2022), no. 12, Paper No. 2250075.

**No claim is made about the Riemann Hypothesis.** Everything here is a measurement
of a detector's behaviour on a planted configuration in a finite truncation, not
evidence about the location of any zero.
