# A probe of the odd parity sector of the CvS Galerkin matrix

**Status.** Numerical report. No theorem is claimed and no claim is made about the
Riemann hypothesis. Every newly reported numerical output is regenerated from
`harness/`, except the **direct-quadrature completeness table in §1.3** — which no script
produces and whose recipe is written out beside it — and the **derived $K(m)$ table in
§5**, which is one division away from `edge_precision.py`'s output, also written out.
Four principal measurements and one follow-up correction are reported — the sector correspondence, the unseeded zero scan,
the random-ensemble test of real-rootedness, and the $T$-scaling of the Galerkin
exponent — each with its limits stated beside it, and four open questions are
collected separately at the end.

**Novelty, stated plainly.** Nothing here is proved. Every section reports a
measurement, and the one thing that would matter theoretically — whether these roots
converge to the zeta zeros as $c \to \infty$ — is exactly the open question this
framework was built around, and is untouched by anything below. Whether the odd-sector
behaviour is already known is not something this report can settle either; it is asked
in §6 rather than assumed. What is offered is a probe of a sector the published work
lists as unprobed, with its controls and its limits attached.

**What is and is not claimed, in one paragraph.** The odd sector of the periodic
family coincides exactly with a specific finite sine space. Its near-null vector yields
a function whose **unseeded** sign-change scan of $(0.1, 100)$ returns twenty-nine
sign-changing roots, each matching one of the twenty-nine zeta ordinates below 100, in
order, with no additional sign-changing root detected — errors from $10^{-52}$ at
$\gamma_1$ to $10^{-15}$ at $\gamma_{29}$ —
although the published CvS theorems assume the even sector. Random-ensemble experiments
indicate that real-rootedness may extend to a global odd *kernel* vector, but not to the
lowest vector of the odd sector as such. Nothing here is claimed to be new, proved, or
explained.

**Provenance of the numbers.** Three different pieces of code produce what is below,
and it matters which is which.

- My own fixed-support implementation is double precision and cannot speak to the deep
  eigenvalue regime; that was stated when this thread opened and is unchanged. It
  supplies exactly two things: the entrywise and bulk-spectral sector comparisons of §1.2 and the
  completeness table of §1.3.
- **Every matrix and every eigenvector in §2 comes from A. Groskin's `connes-cvs`
  package in mpmath**, with a parity projection applied to the matrix it returns. The
  reconstruction $F$ used for the unseeded scan is re-implemented from the published
  formula rather than called from the package — which is why §2 carries a cross-check
  against the package's own `extract_zeros` on the same eigenvector.
- §3 uses neither: it is a synthetic ensemble drawn from CvS eq. (11) in numpy and
  mpmath, with no Weil form and no prime sum anywhere in it.

**Attribution.** The parity split of the sine modes, the exact parity decoupling,
the edge expansion and the $(m-1)/(m+2)$ law are due to A. Groskin (issue #2 of
`akivag613/connes-cvs-`). The answer that $L = u = \log c$ is structural for the
CvS path is his. CvS Theorem 5.6, Theorem 6.1, eq. (11), eq. (21) and Appendix B.3
are Connes and van Suijlekom's, arXiv:2511.23257.

---

## 1. The correspondence between the two constructions

### 1.1 The parameter map

`build_galerkin_matrix` sets $L = \log c$ and sums prime powers $n \le c$;
`extract_zeros` integrates against $\exp(2\pi ikx/L)$ on $[0,L]$, i.e. the basis
$U_n(x) = L^{-1/2}\exp(2\pi inx/L)$ of CvS Proposition 4.1. Against the sine basis

$$\varphi_j(x)  =  \sqrt{\tfrac{2}{L}} \sin\Big(\frac{j\pi (x + L/2)}{L}\Big),
\qquad j = 1,\dots,m,$$

on $x \in [-L/2, L/2]$, this fixes

$$L = \log c, \qquad p_c = c, \qquad m = 2N,$$

which is the saturated diagonal. Everything below sits on it.

### 1.2 The odd sector: an exact identification on a subspace

The statement that the sine space is not an invertible basis change of the finite
periodic Fourier space is correct as a statement about the whole spaces, and is not
disputed here. What holds is narrower, and the variable matters: the periodic modes
$U_n$ live on $[0,L]$, so both sides are written in $x' = x + L/2 \in [0,L]$, in which
$\varphi_j(x') = \sqrt{2/L} \sin(j\pi x'/L)$. Then

$$\frac{U_k - U_{-k}}{i\sqrt{2}}  =  \sqrt{\tfrac{2}{L}} \sin\Big(\frac{2\pi k x'}{L}\Big)  =  \varphi_{2k}(x'),$$

with no sign. The real projected combination used numerically,
$(e_k - e_{-k})/\sqrt2$, differs from $\varphi_{2k}$ by the overall phase $i$; that
leaves the Gram matrix and its zeros unchanged, and is noted only because the sine
coefficients are imaginary in the complex-exponential basis. (Written instead in the
centred variable $x$ one has
$\varphi_{2k}(x) = (-1)^k\sqrt{2/L}\sin(2\pi kx/L)$, so mixing the two variables
introduces a spurious $\operatorname{diag}((-1)^k)$; the algebra fixes the pairing used
below.) **The two finite odd basis spaces coincide exactly. The observed spectral discrepancy
decreases rapidly with $T$, consistently with archimedean truncation, though I have
not established that this is its only source.**

Measured two ways, and the distinction matters. **Entrywise**, with $G$ the periodic
matrix restricted to its odd sector and $S$ the sine matrix restricted to
$\varphi_2, \varphi_4, \dots, \varphi_{2N}$, `compare_families.py` reports
$\max_{i,j}|G_{ij} - S_{ij}| / \max|G|$ at $c = 13$, $N = 12$, $\mathrm{dps} = 60$:

$$2.04\times10^{-4}  \longrightarrow  2.87\times10^{-5} \qquad (T = 100,\ 200),$$

falling with the archimedean truncation. **Two deliberately mispaired controls on the
same quantity** show that the small figure is an identification and not a coincidence:
pairing the periodic odd sector with the *odd*-indexed sine modes gives $0.856$, and
conjugating the correct pairing by $D = \mathrm{diag}((-1)^k)$ — the sign pattern that
appears if the two sides are read in different variables — gives $0.765$. Both are of
order one; the correct pairing is four orders below them and falling.

**Spectrally**, the same script reports the bulk metric

$$\Delta_{\mathrm{odd}}
=\frac{\max_{j\ \text{in the largest six}}|\lambda_j(G)-\lambda_j(S)|}
{\max(|\lambda_{\max}(G)|,|\lambda_{\max}(S)|)},$$

which at $c = 13$, $N = 16$, $\mathrm{dps} = 80$ falls
$1.71\times10^{-4} \to 2.42\times10^{-5} \to 3.22\times10^{-6}$ at $T = 100, 200, 400$.
That is a comparison of the top six eigenvalues — neither an entrywise comparison nor a
statement about the near-null eigenvalue. **I do not know whether either residual is
entirely explained by the archimedean truncation.**

### 1.3 The even sector: completeness without a finite identification

The odd-indexed modes satisfy $\varphi_j(L - x') = (-1)^{j+1}\varphi_j(x')$, so for
$j$ odd they carry the same parity as $\lbrace 1, \cos(2\pi k x'/L)\rbrace $ and are complete in
the same subspace of $L^2[0,L]$ — but only asymptotically. Expanding the constant
and $\cos(2\pi kx'/L)$ in them gives coefficients decaying like $1/j$ and an $L^2$
residual like $1/m$:

| modes up to $j$ | 10 | 50 | 200 | 400 |
|---|---|---|---|---|
| residual, constant | $4.04\times10^{-2}$ | $8.10\times10^{-3}$ | $2.03\times10^{-3}$ | $1.01\times10^{-3}$ |
| residual, $\cos(2\pi x'/L)$ | $8.30\times10^{-2}$ | $1.62\times10^{-2}$ | $4.05\times10^{-3}$ | $2.03\times10^{-3}$ |

No script in `harness/` produces this table; it is a direct quadrature, quoted here so
that it is reproducible without one. With $c_j = \int_0^L f \varphi_j$ over odd $j$ on a
uniform grid, the entry is $1 - \sum_{j \le M} c_j^2 / \lVert f\rVert^2$; the values above
are at $c = 13$ on $4\times10^5$ points, and the leading coefficients of the constant go
$1.4419,\ 0.4806,\ 0.2884,\ 0.2060,\ 0.1602$ at $j = 1,3,5,7,9$, i.e. $1.442/j$.

So the finite even blocks of the two constructions disagree. In the same top-six
spectral metric, that disagreement was measured against $T$ and is flat there —
$5.0\times10^{-2}$, $4.7\times10^{-2}$, $4.9\times10^{-2}$ at $T = 100, 200, 400$ —
while the completeness above is a statement
about the basis and not about the matrices. **That the block disagreement converges in
$m$ is therefore an inference from completeness, not something measured here**; the
$m$-dependence of the even-block discrepancy was not swept. A naive $1/m$ extrapolation to a $10^{-55}$
function-space residual would require $m$ of order $10^{55}$; **that is an
extrapolation of an $L^2$ residual and not a conditioning bound on the matrix or its
eigenvalues**, and no claim is made that the two quantities are related by that
rate. The practical consequence is only that this route was not attempted.

---

## 2. An unseeded scan of the odd-sector function

Groskin's §10 lists among its limitations: *"We compute the even sector only. The
odd-sector zero, if any, is not probed."* This section probes it.

**Protocol, stated before the numbers.** The test function $F(\tau)$ is built from the
odd-sector eigenvector coefficients alone. The interval $(0.1, 50)$ is then swept at
step $0.02$ for **sign changes**; each is refined by **bisection only**; the resulting
root list is **written to disk**; and only after that is `mp.zetazero` called, to
compare. No zeta ordinate enters the discovery at any point. (This differs from
`extract_zeros`, which brackets `findroot` on $(\gamma_k - 0.005,\ \gamma_k + 0.005)$
and so is seeded at the answer; the two protocols agree where they overlap, see below.)

Four properties of `odd_sector_unseeded_scan.py` are checkable in the source rather than
taken on trust: `scan_sign_changing_roots` admits a candidate only where two adjacent
grid values differ in sign; `refine_bracket` is bracket-preserving bisection with
$\lceil(\text{dps}+8)\log_2 10\rceil$ halvings and no Newton or secant step; the root
list is written by `atomic_json` **before** `compare_after_discovery` is reached, and the
only call to `mp.zetazero` in the file is inside that later function; and the `residual`
recorded with each root is $|F|$ evaluated at the located point. The near-pole cutoff in
`centred_transform` is $10^{-(\text{dps}-10)}$ rather than the package's fixed
$10^{-130}$ — the cross-check below covers that difference.

**Result at $c = 13$, $T = 400$, $\mathrm{dps} = 120$.** The scan returns **exactly ten
roots**, and they are the first ten zeta ordinates **in order** — no extra
sign-changing root in the interval, and none of the ten missing. Errors at $N = 100$:

| $k$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| $\lvert r_k - \gamma_k\rvert$ | 8.70e−53 | 6.67e−50 | 4.28e−48 | 2.50e−45 | 4.14e−44 | 4.70e−42 | 3.74e−40 | 3.93e−39 | 1.10e−36 | 9.87e−36 |

and $\log_{10}$ of the first-zero error falls with $N$ as
$-32.9,\ -41.3,\ -46.9,\ -51.7,\ -52.1$ at $N = 20, 30, 40, 60, 100$.

**Extended to height 100.** Rerunning the sweep over $(0.1, 100)$ at the same step, on
the cached $N = 100$ eigenvector from the $T = 600$ ladder rather than the $T = 400$ one
above — which is why $\gamma_1$ reads $1.065\times10^{-52}$ here against
$8.70\times10^{-53}$ in the table, the same $T$-shift described below — returns **twenty-nine sign
changes, and there are exactly twenty-nine zeta ordinates below 100.** None missing, none claimed by two roots, and no
spurious sign change anywhere in the interval. The errors continue smoothly:
$1.065\times10^{-52}$ at $\gamma_1$, $1.199\times10^{-35}$ at $\gamma_{10}$,
$4.161\times10^{-23}$ at $\gamma_{20}$, $2.711\times10^{-15}$ at
$\gamma_{29} = 98.83$.

**The height law.** Fitting $\log_{10}$ of the error against the ordinate over all
twenty-nine gives

$$\log_{10}\lvert r_k - \gamma_k\rvert  \approx  0.453 \gamma_k - 58.0,$$

i.e. **about $0.453$ decades lost per unit of height**, so at this cutoff and basis size
the error reaches $O(1)$ near $\gamma \approx 128$. The fit is a straight line through
points spanning 37 decades, and its residuals reach a few tenths of a decade at either
end, so the crossing point is an extrapolation and not a measurement.

(The height-law paragraph below quotes $0.482$ and $\gamma \approx 123$ for the same
cutoff. That fit uses the ten zeros below 50 on the $T = 400$ vector; this one uses all
twenty-nine below 100 on the $T = 600$ vector. The slopes differ by $6.4$% and the
$O(1)$ crossings by $4.1$% — the disagreement between a ten-point and a twenty-nine-point
fit of one line, not two different laws. The longer baseline is the better estimate.)
That is the honest statement of the range: the agreement is not unlimited, it degrades
linearly in the height, and it runs out a little past 100.

**Consistency with the seeded route.** At $N = 100$ the seeded `extract_zeros` gives
8.70e−53, 6.67e−50, 4.28e−48 for the first three — the same values. The two protocols
agree; the unseeded one simply does not need the answer to find them.

**The random control.** A random odd vector in the same basis does produce
sign-changing roots in the interval — 9 to 14 of them across the $c = 13$ runs, and 11
at $c = 17$. In the $c = 13$, $N = 100$, $T = 600$ run its **closest approach to any zeta
ordinate is $0.144$**, against $1.199\times10^{-35}$ for the worst of the ten in that
same run — the two figures must come from one cell, and these do. Thirty-four orders of
magnitude separate them, so the agreement is not a generic feature of the tested random
odd controls.

**Cross-check against the package's own root finder.** The reconstruction $F$ above is
re-implemented here from the published formula, so an error in that re-implementation
would move every root and nothing else in this pipeline would notice. As a control, the
same twenty-nine roots were recomputed with `extract_zeros` from the package itself,
acting on the same cached eigenvector — an independent evaluation of the same
mathematical reconstruction $F$, using a different solver
(`findroot`/Anderson rather than bisection), and a bracket seeded at `mp.zetazero`.
**Every one of the twenty-nine error figures agrees to all printed digits, and the root
positions differ by between $10^{-62}$ and $4.4\times10^{-59}$** — six orders below the
smallest error in the table, which rules out the root-location implementation as the
source of the observed errors. Seeding therefore changes where one looks, not what is found.

**Step control.** The sweep was repeated at $c = 13$, $N = 60$, $T = 400$,
$\mathrm{dps} = 180$ with the step halved to $0.01$: the same ten roots, and
$\lambda_{\min}^{\text{odd}} = 5.07520120029\times10^{-55}$ identical to every printed
digit. So the list is not an artefact of the sampling interval at that resolution.

**Three further cutoffs.** The same unseeded protocol was run at $c = 17$ ($T = 400$,
$\mathrm{dps} = 180$, $N = 60, 80, 100$), at $c = 19$ ($\mathrm{dps} = 220$,
$N = 80, 100, 120$) and at $c = 23$ ($\mathrm{dps} = 300$, $N = 100, 120, 140$).
**All nine runs return ten roots, each closest to its own $\gamma_k$, sorted, none
missing.** The random control returned 11, 11, 14 at $c = 17$; 11, 15, 12 at $c = 19$;
15, 13, 14 at $c = 23$ — again unstable for random vectors and stable for the
eigenvector. At $c = 17$, $N = 60$ the errors are

$$1.240\times10^{-67},\ 1.227\times10^{-64},\ 9.547\times10^{-63},\ \dots,\
1.504\times10^{-50},\ 1.633\times10^{-49}.$$

and the eigenvalue at each rung reproduces the earlier sweep to every printed digit —
$2.41746920332\times10^{-70}$, $1.98674205805\times10^{-75}$,
$2.13985905178\times10^{-76}$ at $N = 60, 80, 100$, against $2.41747\times10^{-70}$,
$1.9867421\times10^{-75}$, $2.1398591\times10^{-76}$ from the earlier run at
$\mathrm{dps} = 80$.

**What that agreement does and does not test.** $\lambda_{\min}^{\text{odd}}$ comes from
building the matrix and diagonalising the odd block; **no root finder touches it**, so
the agreement says nothing about seeded against unseeded search. What it does establish
at $c = 17$ is that the build, the parity projection and the eigensolver reproduce
between the $\mathrm{dps} = 80$ and higher-precision runs — useful, and a different
question. The
protocol itself is tested by the random control, the step control, the completeness
count, and the cross-check against `extract_zeros` on a fixed eigenvector; those four are
what carry it, not this. At $c = 19$ and $c = 23$ the eigenvalues likewise reproduce the
earlier precision-resolved sweep; the contaminated $c = 19$, $N = 80$,
$\mathrm{dps} = 80$ cell is excluded as documented in §2.2. The reproduced values are
$2.57039874544\times10^{-83}$ and $4.92683492454\times10^{-87}$ at $N = 80, 120$, and
$6.42778802803\times10^{-103}$ and $7.43680210541\times10^{-107}$ at $N = 100, 120$.

**The error tracks the eigenvalue, with a constant that grows in $c$.** Writing

$$C=\frac{|\widehat{\gamma}_1-\gamma_1|}{\lambda_{\min}},$$

**each read at the deepest common rung used for the paired comparison**:

| $c$ | $N$ (both sectors) | $C_{\text{odd}}$ | $C_{\text{even}}$ | ratio |
|---|---|---|---|---|
| 13 | 100 | 432.425 | 7005.193 | **16.200** |
| 17 | 100 | 489.384 | 8052.033 | **16.453** |
| 19 | 120 | 510.275 | 8441.326 | **16.543** |
| 23 | 120 | 545.112 | 9097.646 | **16.689** |

Every row is one cell: the same $c$, the same $N$, the same $T = 400$ for both
sectors. An earlier draft mixed the deepest available rung per sector, which makes a
ratio of two numbers that were not measured together.

**The $N$-dependence must be reported with $C$, and the fourth cutoff shows why.** At
$c = 23$, $C_{\text{odd}}$ runs $553.7,\ 545.1,\ 542.9$ at
$N = 100, 120, 140$ — so a table taken uniformly at $N = 100$ would put $c = 23$ two
percent high while the other three are already settled there. At $c = 19$ it runs
$520.0,\ 511.2,\ 510.3$ at $N = 80, 100, 120$; at $c = 17$, $512.8,\ 491.8,\ 489.4$ in
the odd sector and $8514,\ 8103,\ 8052,\ 8049$ in the even one at $N = 60, 80, 100, 120$,
so the $N = 60$ figures overstate both by four to five percent. The even-sector value at
$c = 17$ sits beside the published $8047$ at the same cutoff, measured there at $T = 800$
against $T = 400$ here.

Fitting the four:

$$C_{\text{odd}} \approx 197.9 \log c - 73.7\ (R^2 = 0.9982), \qquad
C_{\text{even}} \approx 3673.6 \log c - 2392.3\ (R^2 = 0.9987),$$

the same functional form as the published even-sector fit
$C \approx 6730\log c - 11268$ (fifteen cutoffs, $R^2 = 0.956$), with different
coefficients. Four points fit two parameters, so neither $R^2$ is evidence of much; the
form is quoted because it matches, not because four points establish it.

**The ratio varies substantially less than either constant over the four tested
cutoffs.** Across $c = 13$ to $23$ the odd constant moves $26.1$% and the even one
$29.9$% — they do not even grow at the same rate — while their ratio moves $3.0$%,
monotonically: $16.200,\ 16.453,\ 16.543,\ 16.689$. Four points is four points.

$C_{\text{odd}}$ is also stable in $N$ from $N = 60$ onward at $c = 13$ (806.4, 534.0,
463.1, 433.1, 432.2, 432.2 at $N = 20, 30, 40, 60, 80, 100$, all at $T = 600$), and it
**was stable under the one change of $T$ that was tested**: at $c = 13$, $N = 60$ it is
433.4 at $T = 400$ and 433.1 at $T = 600$, although $\lambda$ itself moves 24% between
those two. That is a single test at a single cell, not a demonstration of
$T$-independence. What it does show is that at that cell the absolute error is not
$T$-converged because $\lambda$ is not, while the ratio is.

**Height law across the cutoffs.** Fitting $\log_{10}$ of the error against the ordinate
at $N = 100$: **$0.482$, $0.504$, $0.526$ decades lost per unit of $\gamma$ at
$c = 13, 17, 23$**, with the error reaching $O(1)$ near $\gamma \approx 123,\ 160,\ 204$.
The slope creeps up by about $9$% across that range while the reach grows by $66$%:
**in these three fits, most of the observed gain is in the reachable height rather than
in an improvement of the decay rate.** At $N = 60$ the slopes are almost the
same ($0.483$ and $0.513$ at $c = 13, 17$), so the observed trend does not appear to be
driven primarily by $N$ on the rows tested.

**Multi-zero behaviour across four cutoffs.** With $c = 13, 17, 19, 23$ all at
$N = 100$, the per-zero slope is a fit with two degrees of freedom. Regressing
$\log_{10}|r_k - \gamma_k|$ on $\log c$:

Since the error falls as $c$ grows, the fitted slope is negative; the positive quantity
below is the **decay exponent**
$\alpha_k = - d\log_{10}|r_k - \gamma_k| / d\log c$.

| $k$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| $\alpha_k$ | 83.14 | 82.86 | 82.65 | 82.30 | 82.12 | 81.74 | 81.44 | 81.20 | 80.69 | 80.49 |
| ratio to $k=1$ | 1.000 | 0.997 | 0.994 | 0.990 | 0.988 | 0.983 | 0.980 | 0.977 | 0.971 | 0.968 |

with mean $R^2 = 0.9989$. **The spread across the ten is $3.2$%**, against the published
even-sector figure of all ten within $3.8$% over fifteen cutoffs. Read on two cutoffs the same statistic gave $3.6$%, on three $3.4$%: **the spread
decreases modestly as points are added, an encouraging but still short sequence.**

Two checks on that number. The $c = 19$ errors are not converged at $N = 100$ — going to
$N = 120$ moves them by a factor $3.1$ — but that is a nearly uniform shift in $\log$
across all ten $k$, so it moves the absolute slopes (to $82.16 \dots 79.46$) and leaves
the **spread at $3.3$%**. And read at $N = 60$, where none of the cutoffs has settled,
the same statistic gives $6.9$% — nearly double. The spread is robust to the first and
not to the second, which is worth stating rather than quoting one figure.

Four cutoffs are still four, not fifteen, and the fit has two degrees of freedom.

**What is not established.** Four cutoffs only, $c = 13, 17, 19, 23$; heights below 100
at $c = 13$ and below 50 at the other three. The scan detects **sign changes**, so a root of
even order would be invisible to it, and two roots closer together than the step would
merge — halving the step to $0.01$ leaves the list unchanged, which bounds that risk at
this resolution but does not exclude it in general. The interval begins at $0.1$, so
anything near the origin is excluded by construction. Nothing here is a theorem, and
none of it bears on the Riemann hypothesis.

And $\lambda_{\min}^{\text{odd}}$ itself is **not converged in $T$**: at
$\mathrm{dps} = 180$ it is identical to the $\mathrm{dps} = 120$ value to every printed
digit, but moving $T = 400 \to 600$ shifts it by $+11.7\%,\ +7.3\%,\ +7.3\%,\ +24.1\%,\
+23.0\%,\ +22.5\%$ at $N = 20, 30, 40, 60, 80, 100$. So at the $10^{-55}$ level the
eigenvalue sits on a $T$-limited plateau rather than an $N$-limited one. The root
positions move with it: at $c = 13$, $N = 60$ every one of the ten errors grows by the
same factor $1.24$ that $\lambda$ does, which is why the ratio above is the stable
quantity and the absolute error is not.

### 2.1 Controls

Four vectors from each sector, all from one cell — $c = 13$, $N = 60$, $T = 400$,
$\mathrm{dps} = 80$ — put through the **seeded** `extract_zeros`, which is the cheaper
route when the question is only whether a vector reconstructs at all. The near-null
figures are therefore the $N = 60$ ones and not the deeper $N = 100$ values quoted
above.

| vector | $|\widehat{\gamma}_1-\gamma_1|$, even sector | $|\widehat{\gamma}_1-\gamma_1|$, odd sector |
|---|---|---|
| near-null | $4.22\times10^{-55}$ | $2.20\times10^{-52}$ |
| largest eigenvalue | $1.27$ | $0.78$ |
| random (Gaussian coefficients) | did not converge in the seeded search | $0.34$ |
| second smallest eigenvalue | $1.32\times10^{-49}$ | $7.35\times10^{-47}$ |

The even random vector did not converge to a root in the seeded search, while the
odd random vector returned one $0.34$ away; a failure of `findroot` does not establish
that no root exists in the window. **Close reconstruction occurred at the low end of the
spectrum, not only at the very bottom**: the second-lowest vector also lands within
$1.32\times10^{-49}$ and $7.35\times10^{-47}$, five to six orders worse than the ground
vector but nothing like the controls. The largest-eigenvalue and random vectors did not
closely reconstruct in these trials. That is an observation over four vectors per sector, not a law: no
functional relation between the eigenvalue and the root error is claimed here, and
establishing one would need a sweep over the whole low end of the spectrum with
$\lambda$ tabulated beside the error.

### 2.2 The lag between the sectors

Every cell run, at $T = 400$. The lag is $\log_{10}$ of the ratio of the two
first-zero errors.

| $c$ | $N$ | dps | $\lambda_{\min}^{\text{even}}$ | $|\widehat{\gamma}_1-\gamma_1|$ even | $|\widehat{\gamma}_1-\gamma_1|$ odd | lag |
|---|---|---|---|---|---|---|
| 13 | 40 | 80 | $6.79298\times10^{-54}$ | $5.17506\times10^{-50}$ | $1.19574\times10^{-47}$ | 2.364 |
| 13 | 60 | 80 | $5.99929\times10^{-59}$ | $4.21519\times10^{-55}$ | $2.19967\times10^{-52}$ | 2.718 |
| **13** | **100** | **80** | $2.077\times10^{-59}$ | $1.4550\times10^{-55}$ | $8.6968\times10^{-53}$ | **2.776** |
| 17 | 60 | 80 | $3.52427\times10^{-74}$ | $3.00057\times10^{-70}$ | $1.23960\times10^{-67}$ | 2.616 |
| 17 | 60 | 150 | $3.52427\times10^{-74}$ | $3.00057\times10^{-70}$ | $1.23960\times10^{-67}$ | 2.616 |
| 17 | 80 | 150 | $1.71325\times10^{-79}$ | $1.38825\times10^{-75}$ | $9.77031\times10^{-73}$ | 2.847 |
| 17 | 100 | 250 | $1.40108\times10^{-80}$ | $1.12815\times10^{-76}$ | $1.04721\times10^{-73}$ | 2.968 |
| **17** | **120** | **250** | $1.01444\times10^{-80}$ | $8.16509\times10^{-77}$ | $7.80056\times10^{-74}$ | **2.980** |
| 19 | 60 | 80 | $1.14409\times10^{-79}$ | $1.04719\times10^{-75}$ | $4.54572\times10^{-73}$ | 2.638 |
| 19 | 80 | 80 | $6.89173\times10^{-82}$ | $1.00758\times10^{-78}$ | $1.23528\times10^{-79}$ | −0.912 (discarded) |
| 19 | 80 | 200 | $2.02675\times10^{-87}$ | $1.75041\times10^{-83}$ | $1.33669\times10^{-80}$ | 2.883 |
| 19 | 120 | 250 | $2.42209\times10^{-91}$ | $2.04456\times10^{-87}$ | $2.51404\times10^{-84}$ | 3.090 |
| 23 | 60 | 150 | $2.76142\times10^{-88}$ | $3.04435\times10^{-84}$ | $1.18794\times10^{-81}$ | 2.591 |
| 23 | 80 | 200 | $1.14292\times10^{-99}$ | $1.11320\times10^{-95}$ | $8.04765\times10^{-93}$ | 2.859 |
| 23 | 100 | 250 | $4.18176\times10^{-107}$ | $3.87555\times10^{-103}$ | $3.55935\times10^{-100}$ | 2.963 |
| 23 | 120 | 300 | $3.13658\times10^{-111}$ | $2.85355\times10^{-107}$ | $4.05389\times10^{-104}$ | 3.152 |
| 23 | 160 | 350 | $5.83185\times10^{-113}$ | $5.27672\times10^{-109}$ | $9.60204\times10^{-106}$ | 3.260 |

Bold rows are the deepest $N$ reached at each of $c = 13$ and $c = 17$; they are not
equally settled, and the convergence paragraph below this table says how far each is.

**Discarded row.** At $c = 19$, $N = 80$, $\mathrm{dps} = 80$ the eigenvalue depth is
82 digits and the two sectors returned $6.89\times10^{-82}$ and $4.40\times10^{-82}$,
reading as the odd channel winning. At $\mathrm{dps} = 200$ the same cell gives
$2.03\times10^{-87}$ and $2.57\times10^{-83}$. This is the backward-error
contamination analysed in Groskin §7, at a different $N$.

**Precision control.** The two $c = 17$, $N = 60$ rows are the same cell at
$\mathrm{dps} = 80$ and $150$; they agree to every printed digit.

**Convergence.** The lag climbs with $N$ at fixed $c$ — 2.364, 2.718, 2.776 at
$c = 13$; 2.616 through 2.980 at $c = 17$ — so no lag figure is meaningful until its
$N$-step has settled. The last $N$-step at each cutoff moves the error by

$$2.9\ (c = 13,\ N = 60 \to 100), \quad 1.38\ (c = 17,\ 100 \to 120), \quad
8561\ (c = 19,\ 80 \to 120), \quad 54\ (c = 23,\ 120 \to 160),$$

and none of them is small enough to call the underlying error converged. The lag is the
steadier quantity, and every cutoff does have a last step to read it on:

$$0.058\ (c = 13),\quad 0.012\ (c = 17),\quad 0.207\ (c = 19),\quad 0.108\ (c = 23)$$

decades over each cutoff's deepest available $N$-step ($2.718 \to 2.776$,
$2.968 \to 2.980$, $2.883 \to 3.090$, $3.152 \to 3.260$). So $c = 17$ is much the
steadiest and $c = 19$ and $c = 23$ the least, which matches the error steps above.
**None of the four is treated here as a converged bound.** The lag increases monotonically with $N$ in all four
series measured here, which suggests the unsettled rows understate it, but no
monotonicity in $N$ is established and none of them is quoted as a bound.

Across the three intervals the even channel gains $5.31$, $5.30$, $5.40$ decades per
unit of $c$; the lag grows by $0.051$, $0.055$, $0.043$ per unit, about one percent
of that. **Only the first of those three intervals rests on two of the better-settled
rows**; the second and third use the $c = 19$ and $c = 23$ figures, which are
provisional in the sense just described, so the near-constancy across the three should
be read as suggestive rather than measured.

---

## 3. Is real-rootedness forced by the matrix class?

Evenness is a hypothesis in every statement this rests on: CvS Theorem 5.6 requires
the one-dimensional kernel to be even with respect to $\gamma$; Theorem 6.1 the
lowest eigenfunction; CCM Theorem 1.1 the smallest eigenvector; and Groskin's finite
Guinand–Weil dictionary is stated for real even Galerkin coefficient vectors. The
vector of §2 satisfies none of them: the even sector is always lower, so it is the
lowest eigenvector **of the odd sector**, and it is not in the kernel.

To see whether that matters, random matrices were drawn from the **algebraic matrix
class of CvS eq. (11)** — an ensemble respecting that algebraic structure, not
matrices arising from the Weil form itself —

$$q_{i,i} = a_i, \qquad q_{i,j} = \frac{b_i - b_j}{i-j}\ (j \ne i), \qquad
a_{-i} = a_i, \quad b_{-i} = -b_i,$$

and $Q - \lambda_{\min}I$ formed — positive semidefinite of the same form, with the
ground eigenvector in its kernel. The roots of eq. (21),

$$P(s)  =  \sum_{k \in \lbrace -N,\dots,N\rbrace } \xi_k
\prod_{\substack{j \in \lbrace -N,\dots,N\rbrace  \cr  j \ne k}} (j - s),$$

were then tested. **The root test was validated against CvS Appendix B.3 first**: at
$N = 2$ with $u = 1$ the roots are real exactly when $v \notin (-2,-\tfrac12)$,
agreeing with their criterion at all 24 values of $v$ tried.

**Experiment A — the global kernel vector is odd.** This is Appendix B.3's
configuration, carried to higher $N$. Roots real in **693 of 693** trials at
$N = 3, 4, 5, 6, 8$. The float64 failures are root-finding noise on the degree-$2N$
polynomial; recomputed in mpmath, **no nonzero imaginary part was detected at the working
precision**.

**Experiment B — the global ground state is even, and the vector tested is the
lowest eigenvector of the odd sector.** This is the configuration of §2. Roots real in
**686 of 1307 trials, about 52%**, with genuine failures reaching
$|\mathrm{Im}| \approx 1.7\times10^{2}$.

Both experiments re-test every float64 failure in mpmath before counting it: a
degree-$2N$ polynomial is badly conditioned in double precision, and an apparent
complex root there is not evidence of one. Without that recheck the raw float64 counts
read 683/693 and 633/1307; the figures above are the corrected ones.

Both experiments, the mpmath recheck and the Appendix B.3 validation are in
`harness/odd_vector_reality.py`; it refuses to report unless A comes out all-real and B
lands between 30% and 70%, so a silent regression in either would fail rather than
print.

**Reading.** Real-rootedness in configuration B is **not forced in the sampled
ensemble**: close to half the draws break it, against none in configuration A. No
explicit certified counterexample is exhibited here, so this is a statement about the
ensemble and not a theorem about every matrix of the form. The CvS matrices possess some
additional structure not captured by this random ensemble; **whether that structure
is arithmetic is precisely what is not known here.** The ensemble was not designed
to isolate any particular constraint, and a non-arithmetic constraint absent from it
would produce the same contrast.

---

## 4. The $T$-scaling of the Galerkin exponent

Groskin §8.2 proposes $s = \sigma_{\text{eff}}T/2\pi$ and §11 asks for the
$N$-convergence measurement to be repeated at larger $T$. Regressing
$\log_{10}|\lambda_{\min}^{\text{even}}|$ on $\log_{10} N$ as in §8.1, at a fixed $N$
grid, $T = 400$ against $T = 1600$:

| $c$ | $N$ grid | $T = 400$ | $T = 800$ | $T = 1600$ |
|---|---|---|---|---|
| 13 | {20, 25, 30} | 24.23 | — | 24.14 |
| 23 | {40, 60, 80} | **46.14** | **46.03** | **45.93** |

with $R^2 \ge 0.9990$ on every fit. At $c = 23$ the three values were all measured here at
$\mathrm{dps} = 250$ on the published $N$ grid, so they form a sequence rather than a
pair: **each doubling of $T$ moves $s$ by $0.23$% and then $0.21$%**, and the
quadrupling from $400$ to $1600$ by $0.46$%. At $c = 13$ the quadrupling moves it by
$0.37$%. The mechanism $s = \sigma_{\text{eff}}T/2\pi$ predicts that each doubling
**doubles** $s$.

The $T = 400$ row at $c = 23$ is also a numerical reproduction in a separate
higher-precision run:
our $\log_{10}|\lambda| = -71.179,\ -87.559,\ -98.942$ against the published
$6.628\times10^{-72}$, $2.761\times10^{-88}$, $1.143\times10^{-99}$ — the same three
numbers to every printed digit, giving the same $s = 46.14$.

The reason is visible in the raw values: the $T$-shifts in $\log_{10}|\lambda|$ are
$0.046$–$0.080$ at $c = 13$ and $0.088$–$0.213$ at $c = 23$ — nearly constant across each
grid, so they move $\log\lambda$ almost uniformly and change the fitted slope against
$\log N$ only slightly. They are not exactly constant, which is consistent with the
residual movement in $s$.

One number the script prints must not be read as a comparison with the published fit: it
reports a slope of $s$ against $\log c$ from the two cutoffs, $38.4$ at $T = 400$ and
$38.2$ at $T = 1600$. Those mix $N$ grids — $s(13)$ is on $\lbrace 20,25,30\rbrace $ and $s(23)$ on
$\lbrace 40,60,80\rbrace $ — so they are not comparable with the published slope $55$, because that
fit uses six cutoffs measured on the common grid $\lbrace 40,60,80\rbrace $ (refitting its six points
gives $s(c) \approx 54.5\log c - 127.6$). Only the $T$-ratio at a fixed $c$ and a fixed
grid is used above.

**No interpretation is offered here.** Whether this is the appropriate regime in which
to test $s = \sigma_{\text{eff}}T/2\pi$ is not something this report can settle, and the
relation is not claimed to be confirmed or refuted by these rows. Table 4 of the paper
already reports a $0.14$ shift between $T = 400$ and $800$ at $c = 13$, so the
insensitivity begins before this test.

**Limits.** Two cutoffs. At $c = 13$ the grid $\lbrace 20,25,30\rbrace $ is not the published
$\lbrace 40,60,80\rbrace $, so the value $s = 24.23$ is not comparable with the published
$s(13) = 9.4$; only the $T$-ratio at a fixed grid is, and that is what is reported. At
$c = 23$ the grid is the published one and all three $T$ were measured here, which is why
that row carries the weight. What is not tested is whether the same flatness holds at
cutoffs beyond 23, or on a grid deeper than $N = 80$.

**A precision note.** At $T = 1600$ with $\mathrm{dps} = 80$ the cells returned
$\log_{10}|\lambda| \approx -3.4$ with five negative eigenvalues, forty-four orders
adrift, at a cutoff where the matrix is positive on every cell. At
$\mathrm{dps} = 200$ the same cell returns to $-38.832$ with none. The tell was the
sign count, not the magnitude: the failure made $\lambda$ too shallow rather than
too deep, so a floor test on digit count would not have caught it.

---

## 5. A follow-up on the $(m-1)/(m+2)$ discrepancy

The $\varepsilon^2$ scaling was reported earlier at $m = 20$. The prefactor depends
on $m$; writing

$$\text{discrepancy}  =  K(m) \varepsilon^{2}$$

— $K$ to keep it distinct from the $C(m)$ of $\sigma_2/\sigma_1 = C(m)\varepsilon^4$ —
the ratio between successive decades of $\varepsilon$ is $100.00$ at every $m$ tested,
and

| $m$ | 10 | 20 | 40 | 60 | 80 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|---|
| $K(m)$ | 9.327 | 20.968 | 44.540 | 68.186 | 91.853 | 115.527 | 174.730 | 233.942 |

so $K(200)\varepsilon^2 = 2.339\times10^{-6}$ at $\varepsilon = 10^{-4}$, which is the
figure $2.4\times10^{-6}$ quoted in the thread. It is not a float64 floor, as had
been assumed here.

`harness/edge_precision.py` prints the law itself — $\sigma_2/\sigma_1$ against
$(m-1)/(m+2)$ for each $(m, \varepsilon)$ — but not $K$. The table above is that output
divided through: $K(m) = |\sigma_2/\sigma_1 - (m-1)/(m+2)| / \varepsilon^2$, evaluated at
$\varepsilon = 10^{-4}$ and confirmed at $10^{-5}$ and $10^{-6}$, where the decade ratio
is $100.00$ throughout. The parity-block shortcut — taking the top two $|{\cdot}|$ from
the two blocks separately rather than from the full matrix — is what made $m = 150$ and
$200$ reachable; it was validated at $m = 40$, where both routes give $0.928571873973$.

---

## 6. Open questions

1. Does the finite Guinand–Weil dictionary, stated for real even Galerkin
   coefficient vectors, extend to odd $v$, or is there a reason it should not?
2. Is there a general theorem that an odd ground-state kernel vector in the CvS
   class makes the polynomial of eq. (21) real-rooted? Appendix B.3 settles $N = 2$
   by hand; experiment A finds no counterexample at $N = 3,4,5,6,8$ in 693 of 693
   trials, while experiment B — the same ensemble with the vector taken from the odd
   sector instead — holds in only about half. That contrast is what makes the
   question look worth asking, but it is evidence and not a proof, and the statement
   may already be known.
3. **How should the §11 test be run?** §4 measures $s$ itself at a fixed $N$ grid
   across $T$. §11 phrases its test as the slope of $s(c)$ against $\log c$, which needs
   two cutoffs on ONE grid, and the cutoffs available here were measured on different
   $N$-grids. Related: Table 14
   does not state $T$ in its caption; §8.2 treats its slope $55$ as a $T = 800$
   measurement, while its $c = 23$ entry $46.1$ is identical to the value derived from
   the $T = 400$ Table 8. Which $T$ generated Table 14, and which pair of cutoffs would
   be the right ones to put on a shared grid?
4. Passing `L = math.log(13)` to `extract_zeros` yields a first-zero error of
   $1.41\times10^{-17}$, whereas passing `L = mp.log(13)` yields
   $1.4549524\times10^{-55}$; $\lambda_{\min}^{\text{even}} = 2.0769626582\times10^{-59}$
   either way. Two facts
   make this worth raising rather than filing as caller error: the quick-start example
   itself passes `L = math.log(13)`, and casting that float to `mp.mpf` afterwards
   cannot recover the digits already lost, so the $10^{-17}$ ceiling is the expected
   outcome of following the example. Would passing `c` and computing `mp.log(c)`
   internally be preferable, or accepting `str`/`mp.mpf` and refusing `float` above
   `dps = 15`, or simply correcting the example? Separately, the
   `build_galerkin_matrix` docstring gives the basis as $\exp(2\pi ikt/(2\log c))$,
   while the implementation uses $2\pi k/L$ with $L = \log c$, matching Groskin §2.2.

---

## 7. Reproducing

All scripts in `harness/`. They require `pip install connes-cvs` and, for anything
above $N \approx 40$, `python-flint`.

| script | what it regenerates |
|---|---|
| `compare_families.py` | §1.2, the top-six bulk-spectral sector comparison, with gates |
| `odd_sector_unseeded_scan.py` | §2, the unseeded scan, the frozen root list and the random control |
| `odd_sector_scan_height100.py` | §2, the extension to height 100 from the cached eigenvector |
| `cross_check_extract_zeros.py` | §2, the two-implementation cross-check on all 29 roots |
| `odd_sector_zeros.py` | §2, the seeded cross-check against `extract_zeros` |
| `reconstruction_controls.py` | §2.1, the four control vectors |
| `convergence_sweep.py` | §2.2, the lag table |
| `odd_vector_reality.py` | §3, experiments A and B and the B.3 validation |
| `sobolev_slope.py` | §4, the $T$-scaling |
| `edge_precision.py` | §5, the $(m-1)/(m+2)$ law itself; $K(m)$ is one division away, written out in §5 |
