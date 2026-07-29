# Compressing the Division Table into a Single Dynamical Column

### Prime clocks, gcd kernels, lag correlations, the von Mangoldt function, and a spectral test

Date: July 2026.

**Research status:** a correct dynamical reconstruction of the prime side of `−ζ'/ζ`. Not a proof, not a new criterion for RH, and — see §§8–10 — not a demonstration that positivity fails either.

---

## 0. What the division table is

Everything below is built on one very simple object, so it is worth stating plainly before any machinery arrives.

Write the odd numbers in a row, left to right: 3, 5, 7, 9, 11, 13, … Under each of them, keep one line per odd prime, and put a mark whenever that prime divides that number. The result is an infinite grid — rows indexed by primes, columns indexed by odd integers — and it records the entire multiplicative structure of the integers in the crudest possible way: *who divides whom*. That grid is the **division table**.

| n → | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 25 | 27 | 35 | 49 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **p = 3** | ● | · | · | ● | · | · | ● | · | · | ● | · | ● | · | · |
| **p = 5** | · | ● | · | · | · | · | ● | · | · | · | ● | · | ● | · |
| **p = 7** | · | · | ● | · | · | · | · | · | · | ● | · | · | ● | ● |

*Figure 0. The division table, first rows. A column with no mark anywhere is a prime.*

Three features of this grid drive the whole document.

1. **Each row is perfectly periodic.** The marks in row `p` repeat with period `p` — no exceptions, no drift. A row is a clock. Nothing about a single row is hard.
2. **Each row can be started late.** Row `p` only needs to be drawn from `p²` onwards; every earlier multiple of `p` already carries a mark from a smaller prime. This is the sieve of Eratosthenes, read column-wise.
3. **The primes are the empty columns.** Primality is not a property you compute; it is the *absence* of any mark. All the difficulty of prime distribution lives in how the periodic rows interleave, never inside any one of them.

The table is therefore an object of two natures at once: row by row it is a bank of trivial periodic signals; column by column it is the sequence of primes. The construction in this document consists of pushing on that tension. We collapse the two-dimensional table into a **single time series** — one value per odd integer, read in order — in such a way that nothing is lost: the rows reappear as lag correlations of that single column (§5). Averaging those lags produces a Dirichlet convolution of the von Mangoldt function (§6), Möbius inversion isolates `Λ` itself, and its Dirichlet transform is the prime side of `−ζ'/ζ` (§7). The remaining sections ask whether the operators produced along the way say anything about the zeros of `ζ`. They do not, and §§8–10 establish precisely why.

**In one sentence.** The division table is the record of which primes divide which odd numbers; this document rewrites that record as a single time series, verifies that the rewriting is exact, and then shows that the rewriting alone carries no information about the Riemann zeros.

---

## Status legend

Every numbered claim carries one of:

| tag | meaning |
|---|---|
| **[P]** | Proved here, and verified numerically. |
| **[V]** | Verified numerically to the stated precision; the identity is classical. |
| **[K]** | Correct, but classical — no novelty claimed. |
| **[I]** | Interpretation / language. Carries no mathematical content on its own. |
| **[C]** | Caveat. A point where the natural or commonly circulated statement is wrong, overstated, or parameter-dependent. |
| **[N]** | New here — a lemma or measurement obtained in this work. |

All numerics below were produced by the scripts in Appendix B and are reproducible. **Every table in this document states the parameters it was computed with**; several of the quantities involved are sensitive to a truncation that is easy to leave implicit, and §8 shows what happens when it is.

---

## Abstract

We compress the odd division table into a single time-ordered column. Each odd prime `p` generates a periodic clock whose composite trajectory begins at `p²` and advances by `2p`. A Hilbert-state representation records the prime-power layers dividing each integer and yields the exact Gram kernel `log gcd`. Lag diagonals recover divisibility indicators; their time averages give a Dirichlet convolution of the von Mangoldt function. Möbius inversion recovers `Λ`, and the Dirichlet transform produces a shifted logarithmic derivative of `ζ`. The chain is exact at every step.

We then test three spectral objects. **(i)** The centred covariance kernel is positive semidefinite, though not — as is easily assumed — because it is a Gram matrix; it is a *time average* of Gram matrices, and the distinction matters, since the pointwise Gram matrix is not even Toeplitz. Its smallest eigenvalue decays to zero, and we identify the exact structural reason: the spectral measure is purely atomic, so Szegő's limit is zero. At any *fixed* layer cutoff the measure has finitely many atoms, the Toeplitz rank saturates, and no decay exponent exists at all. **(ii)** Isolating `Λ` gives a Toeplitz family whose positivity is **a convention artifact**: over odd prime powers — the convention the rest of the construction requires — the symbol satisfies `f(π) = −c < 0` identically, for every parameter and for any nonnegative weights; if the prime 2 is admitted, positivity holds instead, but only for `α` in a bounded window, with a maximum margin of 33% near `α = 0.75`. Neither sign reads any arithmetic. **(iii)** The hybrid prime-power operator fails to predict held-out zeros. The failure is structural, not numerical: the operator is bipartite, so its spectrum is exactly symmetric about zero, and its top is insensitive to the prime cutoff while the zeta ordinates are positive and unbounded.

**Keywords:** division table, sieve, prime clocks, von Mangoldt function, Möbius inversion, gcd kernels, Ramanujan sums, divisibility operators, Riemann zeta function.

---

## Reader's introduction

This section is for a reader meeting the construction for the first time. It contains no new mathematics and can be skipped.

### Why compress the table into a single column

Two dimensions are inconvenient: the table has one row per prime, so it never ends. We therefore read it as a time series — one entry per odd number, in order — and ask whether anything is lost. Nothing is: the rows come back as lag correlations of that single column (§5). What was a two-dimensional bookkeeping problem becomes one signal with structured memory, and the tools of signals — correlations, spectra, positive-definite kernels — become available.

### What is being asked, and what is not

- **Asked, and answered yes:** can divisibility be re-encoded as one dynamical signal, exactly, with the classical prime-counting weight `Λ` recoverable from it? (§§1–7.)
- **Asked, and answered no:** does that re-encoding contain information about the zeros of `ζ`? (§§8–10.)
- **Not asked at all:** a proof of the Riemann Hypothesis. The construction is closed as a route to RH in §12, for reasons that are structural and stated precisely — not for lack of computation.

### The whole idea in one page

```
division table → prime clocks → single time column → Hilbert states → log gcd
→ lag correlations → time average → Möbius inversion → Λ → −ζ'/ζ ⟂ zeros
```

| step | what it does | where |
|---|---|---|
| Prime clocks | each prime becomes a periodic trajectory starting at `p²`, step `2p` | §1 |
| Single column | an event queue emits one entry per odd number; empty entry = prime | §2 |
| Hilbert states | each integer becomes a vector of its prime-power layers | §4 |
| `log gcd` kernel | inner products of those vectors are exactly `log gcd` | §4 |
| Lag correlations | the table's columns reappear as diagonals of the single column | §5 |
| Möbius inversion | time averages give `Λ ∗ 1`; inversion isolates `Λ` | §6 |
| Dirichlet transform | the prime side of `−ζ'/ζ`, valid for `Re s > 1` | §7 |

And then it stops. Three separate walls, each structural:

1. **Convergence (§7).** The series lives in `Re s > 1`. The zeros live in the analytic continuation, which `ζ` supplies and the clocks do not.
2. **Content-free positivity (§8–§9).** The covariance kernel's positivity is automatic, and its decay to zero is forced by the atomic spectral measure. The `Λ`-Toeplitz family's positivity, in either direction, is decided by which primes one admits — not by arithmetic.
3. **Wrong shape of spectrum (§10).** The hybrid operator is bipartite, so its eigenvalues are symmetric about zero, and its top does not move when the prime cutoff is raised. The zeta ordinates are positive and unbounded.

### Prime clocks in practice

A composite is marked by every prime factor at or below its square root (Proposition 1), so several clocks can strike the same number. Two examples make the rule concrete — and the second shows why the `√n` bound matters.

| n | factors | √n | clocks that strike n | why |
|---|---|---|---|---|
| 45 | 3²·5 | 6.7 | `p = 3` and `p = 5` | 45 = 9 + 6·6 on clock 3; 45 = 25 + 10·2 on clock 5 |
| 77 | 7·11 | 8.8 | `p = 7` only | 77 = 49 + 14·2; clock 11 has not started, since 11² = 121 > 77 |

*Table 0. Each clock begins at `p²`, so only primes `p ≤ √n` can strike `n`. 45 is the first odd number struck twice.*

**The idea, plainly.** A composite is never missed: its smallest prime factor is always `≤ √n`, so at least one clock is already running and already pointing at it. A prime is struck by nothing, because a clock that could strike it would have to have started at a square below it.

---

## 1. Starting point: the division table

Index the odd integers by `n_t = 2t + 3`. For every odd prime `p`, the new composites generated by `p` begin at `p²` and recur with increment `2p`:

```
M_{p,r} = p² + 2pr = p(p + 2r),   r ≥ 0.
```

**Proposition 1 [P].** Let `n` be an odd composite and let `p` be **any** prime factor of `n` with `p ≤ √n`. Then `n = p² + 2pr` for some integer `r ≥ 0`.

*Proof.* Write `n = pq`. Since `n` is odd, `q` is odd; since `p ≤ √n`, we have `q ≥ p`. Both are odd, so `q − p = 2r` with `r ≥ 0`, giving `n = p(p + 2r) = p² + 2pr`. ∎

**[C] Scope of the proposition.** The familiar version of this statement covers only the *least* prime factor. That is a special case: the statement holds for every prime factor `p ≤ √n`, and the least prime factor is simply the one that always qualifies. This matters for §2, where a composite legitimately appears on several clocks at once. The bound must be `≤`, not `<`: `n = p²` is the boundary case, with `r = 0`.

*Verified:* 82,016 odd composites below 2×10⁵, checking **every** qualifying prime factor — zero violations.

| prime | start | step | first values |
|---|---|---|---|
| 3 | 9 | 6 | 9, 15, 21, 27, … |
| 5 | 25 | 10 | 25, 35, 45, 55, … |
| 7 | 49 | 14 | 49, 63, 77, 91, … |

**[K]** The `p²` starting point is the standard optimisation in the sieve of Eratosthenes. **[I]** The "birth and motion" reading is explanatory language, not a new sieve theorem.

---

## 2. Event queue and the single column

Maintain a next-event value `c_p` for each discovered prime. Let `Z(n) = { p : c_p = n }`. Process odd `n` in increasing order: if `Z(n) = ∅`, declare `n` prime and insert the event `c_n = n²`; otherwise advance each matching clock by `2p`.

**Proposition 2 [P].** For odd `n ≥ 3`, `Z(n) = ∅` if and only if `n` is prime.

*Proof.* If `n` is composite it has a prime factor `p ≤ √n`, and by Proposition 1 lies on `p`'s trajectory, which was inserted when `p` was declared prime (at `p < n`) and has not passed `n`. Conversely all trajectory points `p(p+2r)` are composite. ∎

*Verified:* 99,997 odd `n` up to 2×10⁵, zero mismatches.

**[I]** The profinite / adding-machine language is an interpretation of the residue bookkeeping.

**Worked example — the queue from 9 to 51.** Each entry below is one step of the single column: the set `Z(n)` of clocks striking `n`.

| n | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 | 27 | 29 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `Z(n)` | {3} | ∅ | ∅ | {3} | ∅ | ∅ | {3} | ∅ | {5} | {3} | ∅ |

| n | 31 | 33 | 35 | 37 | 39 | 41 | 43 | 45 | 47 | 49 | 51 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `Z(n)` | ∅ | {3} | {5} | ∅ | {3} | ∅ | ∅ | **{3,5}** | ∅ | **{7}** | {3} |

*Figure 2. The single column. The primes are exactly the empty cells; 45 is the first double strike; 49 is where clock 7 is born.*

**The idea, plainly.** Take `n = 45`. Two clocks read 45 at that moment: clock 3 (which had read 39) and clock 5 (which had read 35). Both are advanced — to 51 and 55 — and the queue moves on. Nothing is searched, nothing is factored: the clocks arrive on their own, and primality is the absence of an arrival.

---

## 3. Periodic clocks and Fourier coordinates

Let `I_p(t) = 1[ p | (2t+3) ]`. It is periodic in `t` with period `p`.

**The phase, explicitly.** Since `2` is invertible mod `p` for odd `p`,

```
τ_p ≡ −3 · 2⁻¹ (mod p),        I_p(t) = 1[ t ≡ τ_p (mod p) ],
I_p(t) = (1/p) Σ_{m=0}^{p−1} exp( 2πi m (t − τ_p) / p ).                    [V]
```

*Verified:* max deviation `1.4 × 10⁻¹⁴` over all odd `p < 60` and three full periods.

For `q = p^k` the clock density is `q⁻¹`.

**[I] A tempting overstatement.** One might say that the coordinate `−log(q⁻¹) = k log p` "explains the multiplicative frequencies appearing in Euler products." That would be an analogy, not a derivation: nothing in §3 produces an Euler product. The honest statement is that the additive frequencies `m/p^k` and the multiplicative coordinate `k log p` are two different labellings of the same layer, and the passage between them is supplied later, in §7, by the Dirichlet transform — not here.

**[K]** Ramanujan–Fourier analysis of `r`-even functions is the established language for this section.

---

## 4. Hilbert representation and the log-gcd kernel

Define, over **odd** primes `p` and integers `k ≥ 1`,

```
Φ(n)_{p,k} = √(log p) · 1[ p^k | n ].
```

**Proposition 3 [P].** For odd `a, b`:

```
⟨Φ(a), Φ(b)⟩ = log gcd(a,b),   ‖Φ(n)‖² = log n,   ‖Φ(a) − Φ(b)‖² = log( lcm(a,b) / gcd(a,b) ).
```

*Proof.* `Σ_k 1[p^k|a] 1[p^k|b] = min(v_p(a), v_p(b)) = v_p(gcd(a,b))`, so `⟨Φ(a),Φ(b)⟩ = Σ_p log p · v_p(gcd) = log gcd(a,b)`. Setting `b = a` gives the norm. The third identity is `log a + log b − 2 log g = log(ab/g²) = log(lcm/gcd)`. ∎

*Verified:* 4,000 random odd pairs below 2×10⁵; maximum errors `8.9e−16`, `1.8e−15` and `5.3e−15` for the three identities in order.

**[C] The odd convention is binding, and it is used again in §9.** `Φ` is indexed by *odd* primes only. For odd arguments the `p = 2` coordinate vanishes identically, so Proposition 3 would be unaffected by admitting it — but two later quantities are not: the constant `C_α` of §8, and, decisively, the sign of the positivity test in §9. **Every statement below is in the odd convention unless the contrary is stated explicitly**, and §9 shows what the alternative does.

**Worked example — `gcd(45, 75)`.** Each integer becomes a checklist of prime-power layers; the inner product counts the layers they share, weighted by `log p`.

| layer `p^k` | weight | 45 = 3²·5 | 75 = 3·5² | shared |
|---|---|---|---|---|
| 3 | log 3 | ✓ | ✓ | ✓ |
| 9 | log 3 | ✓ | — | — |
| 5 | log 5 | ✓ | ✓ | ✓ |
| 25 | log 5 | — | ✓ | — |

`⟨Φ(45), Φ(75)⟩ = log 3 + log 5 = log 15 = log gcd(45,75)`. Norms: `‖Φ(45)‖² = 2 log 3 + log 5 = log 45`.

**The idea, plainly.** Every integer carries a prime-power fingerprint, and the inner product of two fingerprints is `log gcd`. Divisibility is therefore a geometry: sharing factors means pointing in the same direction, and the distance between two integers is `log(lcm/gcd)`.

**[K]** Positive gcd matrices and kernels are classical.

---

## 5. Lag diagonals recover divisibility

Since `n_{t+h} = n_t + 2h` and `n_t` is odd:

**Proposition 4 [P].** `⟨Φ(n_t), Φ(n_{t+h})⟩ = log gcd(2t+3, h)`.

*Proof.* `gcd(n_t, n_t + 2h) = gcd(n_t, 2h) = gcd(n_t, h)`, the last step because `n_t` is odd. Apply Proposition 3. ∎

**Proposition 5 [P] (layer isolation).** For every prime power,

```
log gcd(n, p^k) − log gcd(n, p^{k−1}) = log p · 1[ p^k | n ].
```

*Proof.* `log gcd(n,p^j) = log p · min(v_p(n), j)`, and `min(v,k) − min(v,k−1) = 1[v ≥ k]`. ∎

*Verified:* `4.4e−16` and `6.7e−16` respectively, over `t < 400`, `h < 60`, `p ∈ {3,5,7,11}`, `k ≤ 3`.

So the original columns of the division table are stored as lag diagonals of the single time column: the information is **reorganised, not compressed**. **[I]** No column is lost and none is cheaper to read; this is a change of coordinates.

**Worked example — one lag.** Take `n_t = 21` and step forward by `h = 3` places: `n_{t+3} = 27`. Then `gcd(21, 27) = 3`, and Proposition 4 predicts the same value from `h` alone: `log gcd(21, 3) = log 3`.

| lag `h` | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| `n_{t+h}` | 23 | 25 | 27 | 29 | 31 | 33 | 35 |
| `⟨Φ(21), Φ(n_{t+h})⟩` | 0 | 0 | log 3 | 0 | 0 | log 3 | log 7 |

The correlation at lag `h` depends only on `gcd(n_t, h)` — so column 3 of the table is the lag-3 diagonal, column 7 the lag-7 diagonal, and so on.

**[K]** Functions of `gcd(n,h)` are `r`-even; their DFT is expressed by Ramanujan sums.

---

## 6. Time average and Möbius inversion

For `α > 0` set `A_α(h) = Σ_{p^k | h} log p / p^{k(α+1)}`, over odd `p`.

**Proposition 6 [P].** `A_α = Λ_α ∗ 1` where `Λ_α(d) = Λ(d)/d^{α+1}` supported on odd prime powers; equivalently `A_α(h) = Σ_{d|h} Λ_α(d)`, and by Möbius inversion

```
Λ(n)/n^{α+1} = Σ_{d|n} μ(d) A_α(n/d)          (n odd).
```

*Proof.* `Λ(d)` is supported on prime powers `d = p^k` with value `log p`, so the two sums are the same sum. Inversion is Möbius over the divisor lattice. ∎

*Verified:* `h ≤ 3000`, `α ∈ {0.25, 0.75, 1.5}`; convolution identity exact, inversion to `≤ 1.7e−16`.

**[I]** "A composite is inherited memory, a prime power introduces a new layer" is a restatement of `μ ∗ 1 = δ`. It is accurate but carries no content beyond Möbius inversion.

**Worked example — recovering one layer.** At `α = 0.75`, with `h = 9` (divisors 1, 3, 9):

```
A_α(3) = log 3 / 3^1.75 = 1.0986 / 6.8385 = 0.160664
A_α(9) = log 3 / 3^1.75 + log 3 / 9^1.75 = 0.160664 + 1.0986/46.77 = 0.184152
Σ_{d|9} μ(d) A_α(9/d) = A_α(9) − A_α(3) = 0.023488 = log 3 / 9^1.75   ✓
```

**Why inversion matters.** The time average `A_α` mixes every layer dividing `h` — a composite inherits the layers of its divisors. Möbius inversion subtracts exactly that inheritance, leaving only what is *new* at `h`. What is new is nonzero only at prime powers, and that is precisely `Λ`.

---

## 7. Dirichlet transform and zeta

**Proposition 7 [V].** For `Re s > 1` and `α > 0`, summing over all prime powers,

```
Σ_{h≥1} A_α(h) h^{−s} = ζ(s) · Σ_d Λ(d) d^{−(s+α+1)} = − ζ(s) · ζ'(s+α+1) / ζ(s+α+1).
```

*Verified* by direct summation to `h = 2×10⁵`: relative error `3.9e−17` at `s = 3+4i, α = 0.75`; `4.6e−17` at `s = 2.5, α = 0.75`; `1.3e−17` at `s = 2.2, α = 1.5`.

**[C] Convergence caveat.** The series converges absolutely only for `Re s > 1`. The poles at `s = ρ − α − 1` sit at `Re s = 1/2 − α − 1 < 0`, far outside that half-plane. They exist **only in the analytic continuation**, which is supplied by `ζ` and is *not* produced by the clock construction. It is tempting to gesture at "external analytic continuation" and still list the poles as an output of the chain. They are not: the chain outputs a Dirichlet series in its half-plane of convergence, and every statement about zeros is imported from outside.

**This is the exact point at which the construction stops being a reconstruction and starts being a citation of `ζ`.**

**Why the chain stops here.** Three facts, in order.

1. The series converges absolutely only for `Re s > 1` — that is the whole domain the clocks can reach.
2. Everything interesting about the zeros happens at `Re s = 1/2`, i.e. in the analytic continuation of the same function.
3. The continuation is a property of `ζ`, imported from outside. No amount of extra clock bookkeeping produces it, because a Dirichlet series does not continue itself.

**[K]** The identity itself is the standard `−ζ'/ζ` Dirichlet series, reached here through lag averages.

---

## 8. The centred covariance kernel

Centre each layer: `Ψ_α(t)_{p,k} = √(log p / p^{kα}) · (1[p^k | (2t+3)] − 1/p^k)`.

**Proposition 8 [P].** With `C_α = Σ_{p,k} log p / p^{k(α+2)}` (odd `p`), the lag covariance is

```
K_α(h) = A_α(h) − C_α      (h ≠ 0),      K_α(0) = Σ_{p,k} log p / p^{k(α+1)} − C_α.
```

*Proof.* For a single odd layer `q = p^k`, the events `q | n_t` and `q | n_{t+h}` require `q | 2h`, i.e. `q | h`. Hence the layer covariance is `1/q − 1/q²` if `q | h` and `−1/q²` otherwise. Weighting by `log p / q^α` and summing gives the statement. ∎

*Verified:* the empirical time average matches `K_α(h)` at `h = 0,1,2,3,6,9,15` (`α = 0.75`, layers to 2000) with maximum residual `4.0e−5`, `3.0e−6`, `1.2e−7` at `T = 2×10⁴, 2×10⁵, 2×10⁶` — falling by factors of 13 and 26 per decade in `T`, i.e. faster than `1/T` and not a clean power, as befits an arithmetic rather than a stochastic average.

### Positivity, and why the obvious proof of it is wrong

**Proposition 9 [P].** The Toeplitz matrices `[K_α(|i−j|)]` are positive semidefinite for every size.

**[C] The one-line proof "they are Gram matrices of the vectors `Ψ_α(t)`" is false.** The pointwise Gram matrix `G_{ij} = ⟨Ψ_α(i), Ψ_α(j)⟩` is **not even Toeplitz**: at `α = 0.75` with layers to 2000, its main diagonal ranges over `0.2555 … 0.4200` and its first off-diagonal over `−0.1621 … +0.0300` within the first six indices alone. Divisibility of `n_t` depends on `t`, not only on the lag.

*Correct proof.* For each shift `s`, the matrix `G^{(s)}_{ij} = ⟨Ψ_α(s+i), Ψ_α(s+j)⟩` **is** a Gram matrix, hence PSD. Averaging over `s` gives `lim_S (1/S) Σ_{s<S} G^{(s)}_{ij} = K_α(i−j)`, and an average of PSD matrices is PSD. Equivalently: `K_α` is the autocovariance of a stationary sequence, so it is a positive-definite function by Bochner's theorem. ∎

*Verified:* the shift average converges to `K_α(h)` as tabulated above.

Positivity is therefore **automatic** — it holds for the covariance of any stationary sequence whatsoever — and cannot be sensitive to RH. This deserves to be stated as strongly as possible.

### The decay of `λ_min`, and the truncation that controls it

**[C] Any table of `λ_min` here is meaningless without its parameters.** Both `α` and the **layer cutoff** used for `A_α` and `C_α` must be stated. At `α = 0.75`, over odd prime powers up to the stated cutoff:

| cutoff | N=10 | N=20 | N=40 | N=80 | N=160 | N=320 | log–log slope |
|---|---|---|---|---|---|---|---|
| 100 | 0.099776 | 0.047711 | 0.016656 | 0.002423 | 0.000026 | 0.000000 | −4.780 |
| 200 | 0.117479 | 0.064086 | 0.030432 | 0.011165 | 0.001903 | 0.000025 | −2.222 |
| 500 | 0.128831 | 0.075076 | 0.040709 | 0.020021 | 0.008002 | 0.001880 | −1.177 |
| 2000 | 0.136795 | 0.082951 | 0.048406 | 0.027363 | 0.014641 | 0.007130 | −0.847 |
| 10⁵ | 0.141011 | 0.087158 | 0.052596 | 0.031518 | 0.018726 | 0.011075 | −0.736 |

The decay is monotone and its *shape* is robust, but **no individual entry and no exponent is meaningful without the cutoff**: the apparent slope moves from `−4.8` to `−0.74` across this table alone. Quoted exponents differing by a factor of three in the literature of this construction are the same kernel at different truncations, not different phenomena.

**[N] At a fixed cutoff there is no exponent at all — the rank saturates.** With finitely many layers the spectral measure has finitely many atoms, so the Toeplitz family has bounded rank and `λ_min` is machine zero beyond it:

| layer cutoff | layers | atoms (incl. `θ = 0`) | rank at N=1000 | `λ_min` at N=1000 |
|---|---|---|---|---|
| 25 | 10 | 117 | 116 | −1.9e−14 |
| 49 | 18 | 399 | 398 | −1.7e−14 |
| 81 | 26 | 909 | 908 | −1.3e−15 |

In each case `rank = (#atoms) − 1`, and the `−1` is not a coincidence: by the identity
`C_α = Σ_q w_q / q` above, the centring annihilates exactly the atom at `θ = 0`, leaving
`(#atoms) − 1` atoms with strictly positive weight. **The rank must be read with a
tolerance below the smallest atom weight.** At cutoff 81 the smallest atom weight is
`6.2e−6` and the rank reads 904 at tolerance `1e−9` but 908 at `1e−12`; the true cut is
unambiguous, with `λ₉₀₇ = 1.0e−10` against `λ₉₀₈ = 5.0e−15`. A decay exponent therefore
exists only in a **joint** limit `cutoff, N → ∞`, and any slope measured at one cutoff is
an artifact of where on the collapse curve one happened to sample.

**[N] The limit itself has an exact structural cause, and it is not arithmetic.** By the Szegő limit theorem, `λ_min` of a Toeplitz family converges to the essential infimum of the **absolutely continuous part** of the spectral measure. Here the measure is

```
dσ(θ) = Σ_q (w_q / q) Σ_{m=1}^{q−1} δ( θ − 2π m / q ),        w_q = log p / q^{α+1},
```

with **no atom at `θ = 0`**: the constant `C_α` subtracted in Proposition 8 is exactly the weight that atom would carry, since

```
C_α = Σ_q log p / q^{α+2} = Σ_q w_q / q,
```

so centring removes the DC atom identically rather than approximately. The measure is **purely atomic**, with no absolutely continuous part at all. Hence `ess-inf(a.c.) = 0` identically, and `λ_min → 0` is forced for *any* weights `w_q ≥ 0`, arithmetic or not. The finite-rank saturation above is precisely the finite-cutoff shadow of this statement.

**Consequence.** The decay in §8 is not a fact about primes. It would be reproduced by any set of periodic clocks whatsoever. It therefore carries no information about the zeros, and no positivity route can be built on it.

**[K]** Centred divisibility indicators also appear in the Helfgott–Radziwiłł divisibility-graph operator.

### The idea, plainly

- **Positivity is free.** It is the positivity of a stationary autocovariance. A property that holds by construction cannot be equivalent to a hard conjecture.
- **The decay to zero is free too.** The clocks are periodic, so the spectral measure sits on rational angles only — it is purely atomic. Szegő's theorem then forces `λ_min → 0`.
- **Neither fact is about primes.** Replace the primes by any set of periods and any nonnegative weights and the same two conclusions hold verbatim.

---

## 9. The positivity test

The natural next test isolates `Λ` itself: set `b(h) = Λ(h)/h^{α+1}` and study the Toeplitz family `[b(|i−j|)]`. The result is the sharpest negative in this document, because the answer flips sign depending on a bookkeeping choice that carries no mathematics.

### (a) With zero diagonal, indefiniteness is forced

**Proposition 10 [P].** Any real symmetric matrix with zero trace and at least one nonzero entry has a strictly negative eigenvalue.

*Proof.* The eigenvalues sum to the trace, which is `0`; they are not all `0` because the matrix is nonzero. ∎

Setting `b(0) = 0` sets the trace to `0` by construction, so negative eigenvalues are guaranteed before any arithmetic is inserted, and the same table would appear with `Λ` replaced by any nonzero sequence.

*Verified (all prime powers, `α = 0.75`):* `λ_min = −0.421561, −0.461791, −0.507879, −0.533619, −0.546082` with `5, 12, 23, 47, 70` negative eigenvalues at `N = 10, 20, 40, 80, 120` — the trace exactly `0` by construction in every case, and the computed eigenvalues summing to at most `2.2e−15` in magnitude.

**The smallest possible illustration.** `M = [[0, b], [b, 0]]` has eigenvalues `±b` and trace `0`, so a negative one is unavoidable — for **any** `b ≠ 0`, arithmetic or not.

### (b) With the natural diagonal, the answer depends on which primes are admitted

The natural diagonal is the value the lag-0 entry has in §8's construction, `c = A_α(0) = Σ_q Λ(q) q^{−(α+1)}`. The Toeplitz family with symbol `f(θ) = c + 2 Σ_{h≥1} b(h) cos(hθ)` is PSD for all `N` iff `f ≥ 0`.

**Proposition 11 [P] (the odd convention never gives positivity).** Let `(w_h)` be any nonnegative weights supported on **odd** `h ≥ 3` with `c = Σ_h w_h < ∞`, and let `f(θ) = c + 2 Σ_h w_h cos(hθ)`. Then

```
f(π) = c − 2c = −c ≤ 0,
```

with equality only if all weights vanish. Hence the Toeplitz family is not positive semidefinite.

*Proof.* Every `h` in the support is odd, so `cos(hπ) = −1` for every term. ∎

Since §4 fixes the odd convention, this is the statement that applies to the construction of this document. Measured at `α = 0.75` over odd prime powers:

```
c = 0.579509,     min_θ 2 Σ b(h) cos(hθ) = −1.158989,   attained at θ = π,
symbol minimum = c − 1.158989 = −0.579480,   i.e. −c to six digits.
```

| N | 10 | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|
| `λ_min` | +0.2104 | −0.0264 | −0.2184 | −0.3532 | −0.4409 | −0.4958 | −0.5293 |

negative from `N = 20` onward, converging to `−c`.

**[N] If the prime 2 is admitted instead, positivity holds — but only in a window of `α`.** With all prime powers, at `α = 0.75`:

```
min_θ 2 Σ_{h≥1} Λ(h) h^{−1.75} cos(hθ) = −0.58633   (attained at θ = 2π/5, exactly to 8 decimals)
natural diagonal: c = −ζ'/ζ(1.75) = 0.872770
⇒ symbol minimum 0.28644   (32.8% of the diagonal)
```

| N | 10 | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|
| `λ_min` | 0.4512 | 0.4110 | 0.3649 | 0.3392 | 0.3194 | 0.3067 | 0.2988 |

positive at every size, converging to the symbol minimum.

**But the margin is not generic.** Sweeping `α`, the family is PSD **only for**

```
α ∈ (0.495075…, 1.6096),
```

and the margin `f_min / c` rises to a maximum of `0.328` near `α = 0.75` and falls to `0` at both edges. Both failures have closed forms, so they need no numerics:

- **Below the window**, `f(π) = 4 log 2 / (2^{α+1} − 1) + ζ'/ζ(α+1)` goes negative. Its root is `α = 0.495075…`, and this is the lower edge **exactly**, not approximately: throughout that region the symbol's minimum sits at `θ = π` itself (measured `argmin/π = 1.0000` at `α = 0.45, 0.49, 0.50, 0.52, 0.60`), so `min_θ f = f(π)` there and the two vanish together. Values of `f(π)`: `−0.2988` at `α = 0.40`, `−7.016` at `α = 0.10`, `−16.842` at `α = 0.05`.
- **Above the window**, `f(π/2) = c + 2 log 2 ( 2^{−2(α+1)}/(1 − 2^{−(α+1)}) − 2^{−(α+1)} )` goes negative: `−0.0172` at `α = 3`, `−0.0087` at `α = 5`.

**[C] Numerical warning.** Direct truncated evaluation of the symbol is unusable below `α ≈ 0.45`: at cutoff `4×10⁶` it reports `+1.86` at `α = 0.05`, where the exact closed form gives `−16.84`. The two closed forms above must always be used as a cross-check. Efficient evaluation of the symbol: fold the prime-power weights modulo `M = 2²¹` and take one real FFT — `f(2πj/M)` is then exact on the grid.

### Corrected conclusion for §9

**The sign of the answer is a convention artifact.** Over odd prime powers — the convention this construction requires — positivity fails identically, by Proposition 11, for *any* nonnegative weights; admit the prime 2 and positivity holds instead, but only inside a bounded window of `α`. Neither outcome reads any arithmetic: the negative one is forced by parity, the positive one by the size of a single 2-adic layer.

**[N] And even at its best, the margin closes the route.** A quantity equivalent to RH must be *critically* balanced — the Weil quadratic form is positive with `β* = 1` as its unique positive scaling point, so it has no margin at all in the scaling of its prime block. An object whose best case carries a 33% margin, at a parameter value chosen to maximise it, cannot be equivalent to a zero-margin condition: it is asking whether a comfortably positive symbol is positive, not whether two large terms cancel exactly.

---

## 10. The hybrid operator — the negative result is structural

The finite operator is

```
(Hf)(t) = Σ_{q = p^k ≤ P} ( log p / √q ) · ( 1[q | n_t] − 1/q ) · [ f(t+q) + f(t−q) ].
```

**Lemma 1 [N] (self-adjointness).** `H = Hᵀ`.

*Proof.* The `(t, t+q)` entry is `w_q(1[q|n_t] − 1/q)` and the `(t+q, t)` entry is `w_q(1[q|n_{t+q}] − 1/q)`. Since `n_{t+q} = n_t + 2q` and `q | 2q`, we have `q | n_{t+q} ⇔ q | n_t`, so the two coincide. ∎

*Verified:* `‖H − Hᵀ‖_max = 0` exactly at `(N,P) = (160,25), (240,49), (320,81)`. Self-adjointness is **not** automatic here: the shift `2q` paired with the modulus `q` is precisely what makes it true.

**Lemma 2 [N] (bipartite ⇒ symmetric spectrum).** Let `D = diag((−1)^t)`. Then `D H D = −H`, so the spectrum of `H` is symmetric about `0`.

*Proof.* Every shift `q = p^k` with `p` odd is odd, so `H_{t,t'} ≠ 0` forces `t − t'` odd, hence `(−1)^t (−1)^{t'} = −1`. ∎

*Verified:* `‖DHD + H‖_max = 0` exactly, and `λ_max + λ_min = O(10⁻¹⁴)`:

| N, P | 160, 25 | 240, 49 | 320, 81 | 400, 121 |
|---|---|---|---|---|
| `λ_min` | −2.206051 | −2.271354 | −2.321267 | −2.360540 |
| `λ_max` | +2.206051 | +2.271354 | +2.321267 | +2.360540 |

**[C] A consistency check that catches errors.** By Lemma 2 the two rows must be exact negatives. Any table in which the `λ_max` column is not `−λ_min` entry by entry is wrong — the usual cause is a row offset between the two columns.

### [N] The zero-matching test is dead before it is run

**1. Symmetry.** The spectrum is symmetric about `0`; the ordinates `γ_n` are all positive. Any affine map `λ ↦ aλ + b` sending part of a symmetric spectrum onto `{γ_n}` must fail on the mirror half.

**2. The top of the spectrum does not read the primes.** This must be stated carefully, because a table that raises `P` and `N` together attributes to the prime cutoff what belongs to the matrix size. Separating them:

*Fixed `N = 2000`, raising the prime cutoff:*

| P | 25 | 49 | 81 | 121 | 169 | 289 | 529 | 1089 | 2209 | 4489 |
|---|---|---|---|---|---|---|---|---|---|---|
| `λ_max` | 2.3314 | 2.4339 | 2.4684 | 2.4852 | 2.4891 | 2.4956 | 2.4973 | 2.4977 | 2.4977 | 2.4977 |
| `Σ_q log p/√q` | 6.23 | 10.27 | 14.01 | 17.91 | 21.96 | 30.12 | 41.87 | 61.71 | 89.56 | 129.67 |

`λ_max` **saturates completely**: it is constant to four decimals from `P = 1089` on, while the raw weight budget grows by a further factor of two and by a factor of 21 overall.

*Fixed `P = 4489`, raising the matrix size:*

| N | 5120 | 20480 | 81920 |
|---|---|---|---|
| `λ_max` | 2.57867 | 2.62714 | 2.69060 |

Here it does grow, slowly, with no sign of saturation. **So the correct statement is not that `λ_max` is bounded — it is not, on this evidence — but that it is completely insensitive to the prime cutoff.** Adding primes adds nothing to the top of the spectrum; only enlarging the matrix does. Meanwhile the number of ordinates below `T` grows like `(T/2π) log(T/2π)`. There is nowhere for an unbounded set of ordinates to go, and the growth that does occur is not arithmetic in origin.

**[C] Any RMSE table here is protocol-dependent, so the protocol must be stated.** At `N = 320, P = 81`, fitting an affine map on the five largest eigenvalues against `γ₁…γ₅` and evaluating held-out RMSE on `γ₆…γ₁₅`:

| pairing convention | slope | held-out RMSE |
|---|---|---|
| 5 largest eigenvalues, descending ↦ `γ₁…γ₅` | −71.45 | **7.890** |
| 5 largest eigenvalues, ascending ↦ `γ₁…γ₅` | +67.48 | **47.233** |
| 5 largest by `\|λ\|` ↦ `γ₁…γ₅` | −0.43 | **27.768** |

A factor of six separates the conventions, so a bare number is uninterpretable. Nor does the smallest of them mean anything: by Lemmas 1–2 **no** convention can work.

**[V] Randomised-phase control.** Replacing each layer's true phase `τ_q` by a uniform random residue, preserving all periods and densities, gives a mean normalised top-10 spectral distance of `0.0309` over 8 seeds. Destroying the arithmetic while keeping the periodic structure barely moves the spectrum: the spectrum is set by the periodic density structure, not by the arithmetic phases.

**[K]** The construction is structurally close to Helfgott–Radziwiłł with prime-power shifts and von-Mangoldt weights added.

**Why the spectrum is symmetric.** Every shift `q = p^k` is odd, so the operator only ever connects an even index to an odd one. The index set splits in two with no edge inside either half — a bipartite graph. Flipping the sign of one half conjugates `H` into `−H`, so `λ` and `−λ` always come in pairs.

**The idea, plainly.** The operator's eigenvalues are symmetric about zero and its top does not move when primes are added. The zeta ordinates are all positive and grow without bound. Fitting one onto the other is not a question of tuning: half the spectrum has the wrong sign, and the part that has the right sign is not listening to the primes.

---

## 11. Literature map

| component | status |
|---|---|
| Periodic prime clock | **[K]** residue-class / sieve structure |
| Start at `p²` | **[K]** classical sieve optimisation |
| Profinite clock bank | **[I]** interpretation |
| Frequencies `m/p^k` | **[K]** finite Fourier / Ramanujan sums |
| Positive gcd kernel | **[K]** classical gcd-matrix literature |
| Lag functions `log gcd(n,h)` | **[K]** `r`-even functions |
| Möbius inversion to `Λ` | **[K]** standard |
| Dirichlet transform to `−ζ'/ζ` | **[K]** standard |
| Stationary-average (not Gram) proof of Proposition 9 | **[N]** this work |
| Finite-rank saturation of the covariance Toeplitz family | **[N]** this work |
| Atomic-measure explanation of `λ_min ↓ 0` | **[N]** this work |
| Parity obstruction to positivity in the odd convention (Prop. 11) | **[N]** this work |
| The `α`-window for positivity when `p = 2` is admitted | **[N]** this work |
| Self-adjointness and bipartite symmetry of `H` | **[N]** this work |
| Insensitivity of `λ_max(H)` to the prime cutoff | **[N]** this work |
| Single-column narrative | **[I]** uncommon synthesis, no new theorem |
| New restriction on zeros | **not obtained** |

---

## 12. Conclusion

The chain

```
division table → prime clocks → Hilbert states → log gcd → lag correlations
→ Möbius inversion → Λ → −ζ'/ζ
```

is exact at every step, and every step has been derived and verified here. It is a genuine reorganisation of divisibility into a single time series, and it recovers the prime side of `ζ` correctly.

It does not yield a proof or a new criterion for RH, and the reasons are structural:

1. **§7 is where the construction ends.** The Dirichlet series converges only for `Re s > 1`. Every statement about zeros lives in the analytic continuation, which is imported from `ζ`, not produced by the clocks.
2. **§8's positivity is automatic** — it is the positivity of a stationary autocovariance — **and its decay to zero is forced** by the purely atomic spectral measure via Szegő. At any finite layer cutoff the family even has bounded rank. Both facts would hold for any periodic clock bank whatsoever.
3. **§9's answer is a convention artifact.** In the odd convention positivity fails by parity alone (Proposition 11), for any nonnegative weights; admitting the prime 2 restores it, but only for `α` in a bounded window and with a margin that is maximised, not typical, at 33%. A quantity equivalent to RH is critically balanced; this one is decided by bookkeeping.
4. **§10 fails structurally, not numerically.** `H` is bipartite, so its spectrum is exactly symmetric about zero, and its top is insensitive to the prime cutoff while the ordinates are positive and unbounded. The randomised-phase control confirms the spectrum is a density statistic.

**Decision.** Close this route as a direct approach to RH. Retain it as a correct and unusually clean dynamical model of divisibility and of the prime side of `ζ`, and retain Propositions 9 (as proved here), 10 and 11, Lemmas 1 and 2, the finite-rank saturation, and the atomic-measure argument, which are the substantive contributions.

---

## 13. What the model teaches

The negative results are the useful part, because each one is a general lesson about a class of attempts rather than a fact about this particular construction.

| observation | the lesson |
|---|---|
| The chain to `−ζ'/ζ` is exact | Recovering the prime side of `ζ` is easy and proves nothing on its own; it is bookkeeping, not analysis. |
| A stationary autocovariance is positive | Check whether a positivity property is automatic before reading meaning into it — and check *why* it is automatic, since the obvious reason may be the wrong one. |
| `λ_min ↓ 0` from an atomic measure | Decay produced by periodicity alone carries no arithmetic content. Test it against a random-period control. |
| The rank saturates at a fixed cutoff | An exponent measured at one truncation may not exist in the limit at all. Vary the truncation before quoting a slope. |
| Zero diagonal ⇒ negative eigenvalue | Normalisation choices can manufacture a dramatic result. Always ask what the trace was forced to be. |
| The sign of §9 flips with one prime | If a conclusion changes when a bookkeeping convention changes, the conclusion was never about the mathematics. |
| A 33% margin | RH-equivalent statements are critically balanced. A comfortable margin means the statement is weaker, not that it is stronger. |
| Bipartite ⇒ symmetric spectrum | Compare the shape of a spectrum (sign, growth, symmetry) with the target before fitting anything numerically. |
| `λ_max` saturates in `P` but not in `N` | Sweep one parameter at a time. A table that moves two at once will attribute the effect to the wrong one. |
| Random phases reproduce the spectrum | A control that destroys the arithmetic but keeps the structure tells you which of the two your result depended on. |

### Things to try on the model

1. **Change the weights.** Replace `log p / p^{kα}` by any nonnegative weights and re-measure `λ_min`. It still decays to zero — the prediction of §8 that the decay is not arithmetic.
2. **Break the periodicity.** Use non-integer or irrational periods. The spectral measure gains an absolutely continuous part and `λ_min` stops going to zero — the sharpest way to see what the atomic measure was doing.
3. **Scramble the phases.** Keep every period and density, randomise `τ_q`. The top of the spectrum barely moves (§10) — the model is reading density, not arithmetic.
4. **Add an even shift or a diagonal term** to `H` and watch the spectral symmetry break. Whether anything useful survives is exactly the question §10 leaves behind.

**Closing note for the reader.** The construction does one thing well: it turns divisibility into a single signal, exactly, with nothing hidden. It then fails to say anything about the zeros — three times, for three different structural reasons, each identifiable without a single additional computation. Knowing where a route ends, and why the ending is structural rather than numerical, is a usable result: it tells the next attempt which walls it must climb rather than walk into.

---

## Appendix A. Open items that are actually open

These are the questions this document did **not** settle, stated so that a later reader does not mistake them for closed.

1. **The joint limit in §8.** At a fixed layer cutoff the Toeplitz rank saturates and no decay exponent exists; as the cutoff grows the measured slope drifts (`−4.78` at cutoff 100 to `−0.736` at 10⁵) and has not stabilised. What the correct law is in a joint limit `cutoff, N → ∞` — and at what relative rate the two must be taken — is unmeasured. Szegő already forces the limit to be `0`, so this is a rate question, not an existence question.
2. **The upper edge of the `α`-window in §9, and the cusp that pins the minimum.** The lower edge is exact (the root of `4 log 2/(2^{α+1} − 1) + ζ'/ζ(α+1) = 0`, at `α = 0.495075…`), because the symbol's minimum sits at `θ = π` throughout that region. The upper edge, `α = 1.6096`, is only numerical, and the reason is a structure this document measures but does not explain. The minimising angle **locks onto rational points**: it is exactly `π` up to about `α = 0.71`, then exactly `2π/5` — to 8 decimals at `α = 0.75` and 7 at `α = 1.00` — and only then does it come loose, drifting to `0.39950π` at `α = 1.40` and `0.39892π` at `α = 1.60`. The locking is caused by a downward cusp: measuring `f(2π/5 + ε) − f(2π/5) ~ ε^β` gives

   | `α` | 0.60 | 0.75 | 1.00 | 1.40 |
   |---|---|---|---|---|
   | `β` | 0.624 | 0.769 | 1.004 | 1.063 |

   (summation cutoff `1.6×10⁷`; the exponents move *towards* `α` as the cutoff grows — at `4×10⁶` the same measurement reads `0.654, 0.793, 1.016, 1.063`) — so `β ≈ α` while the minimum is pinned, and the pinning fails exactly when the cusp becomes Lipschitz (`β = 1`, at `α ≈ 1`) and can no longer hold a minimum against the surrounding slope. Why the exponent should equal `α`, why the transition is at `α = 1`, and which rational angle takes over and when, are all open. This is the one place in the document where a measured structure looks like it might have content beyond bookkeeping.
3. **The growth law of `λ_max(H)` in `N`.** Measured `2.57867 → 2.69060` over `N = 5120 … 81920` at fixed `P = 4489`, with the increment per decade *increasing*, so it is super-logarithmic; no fitted law is stable over this range. Since §10's argument rests on the insensitivity to `P` and on Lemma 2, this does not affect any conclusion — but the law is not identified.

## Appendix B. Reproducibility

Three scripts accompany this document; each prints the tables above.

- `verify_core.py` — Propositions 1–7, the §10 spectral table with Lemmas 1 and 2, the RMSE protocol comparison, and the randomised-phase control.
- `verify_covariance.py` — §8: the identity `C_α = Σ_q w_q/q` and the spectral measure it forces, the non-Toeplitz pointwise Gram matrix and the shift-average identity that replaces it, the convergence rate of the time average, the layer-cutoff sweep, and the finite-rank saturation with its tolerance column.
- `verify_positivity.py` — §9: both conventions, Proposition 11 with a random-weight control, the symbol minimum by FFT folding, the exact lower edge and bisected upper edge of the `α`-window, the off-grid search that locates the minimising angle, the cusp exponents, and the `N` versus `P` separation in §10.

Environment: Python 3, NumPy 2.4, SciPy, mpmath 1.3. Zeta ordinates from `mpmath.zetazero` at 20 digits.

## AI assistance

The verification script accompanying this paper, and much of its prose, were
written with the assistance of Claude (Anthropic). The research direction, the
decisions, and the responsibility for every claim are the author's. See the
repository README for a fuller statement.

---

## References

1. H. A. Helfgott and M. Radziwiłł. *Expansion, Divisibility and Parity.* arXiv:2103.06853.
2. H. A. Helfgott. *Expansion, Divisibility and Parity: An Explanation.* arXiv:2201.00799.
3. H. G. Gadiyar and R. Padma. *Linking the Circle and the Sieve: Ramanujan–Fourier Series.* arXiv:math/0601574.
4. P. Haukkanen. *Discrete Ramanujan–Fourier Transform of Even Functions.* arXiv:1210.0295.
5. L. Tóth. *The Discrete Fourier Transform of r-even Functions.* arXiv:1009.5281.
6. D. Guillot and X. Wu. *Total Nonnegativity of GCD Matrices and Kernels.* arXiv:1901.01947.
7. U. Grenander and G. Szegő. *Toeplitz Forms and Their Applications.* University of California Press, 1958. *(Szegő limit theorem, used in §8.)*
8. E. Bombieri. *Remarks on Weil's quadratic functional in the theory of prime numbers.* Rend. Lincei 11 (2000), 183–233. *(Critical balance of the Weil form, used in §9.)*
