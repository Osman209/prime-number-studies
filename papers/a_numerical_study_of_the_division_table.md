# A Numerical Study of the Division Table

#### One elementary table, and the classical mathematics you can see from it

**Mohamed Osman** · ORCID [0009-0004-5912-999X](https://orcid.org/0009-0004-5912-999X)

*Independent researcher · August 2026*

Repository: [github.com/Osman209/prime-number-studies](https://github.com/Osman209/prime-number-studies) · Licence: CC BY 4.0 (text), MIT (code)

---

## Abstract

This note starts from a single spreadsheet: the quotients $n/d$ for $n$ up to $100$ and every divisor $d$, written as mixed numbers. Three features of that sheet turn out to be classical objects in disguise — the summatory divisor function, the harmonic numbers, and, after a one-column modification, the Mertens function and with it a determinant criterion equivalent to the Riemann hypothesis — and a particular walk across it opens onto several more. The walk chooses one entry per row so that the two factors' sum stays constant. The cells you touch trace an arch. This note follows that arch as far as it goes and reports, at every step, exactly which piece of classical mathematics it has become: Fermat's factorisation method, Kraitchik's repair and the quadratic sieve, the Dirichlet divisor problem, the Gauss circle problem, Voronoï's series, and the additive–multiplicative divide. Two later sections leave the walk and truncate the sheet instead — keeping only the rows up to a cutoff turns the same table into the sieve of Eratosthenes — and that second route reaches Dickman's and Buchstab's functions and the wheel sieve.

**No result in this note is new.** Its purpose is the route. The arch is an unusually short path from counting multiples to the classical objects named above, and each transition can be checked by hand or in a few lines of code. Every claim is verified numerically; every attribution is given; a closing section states plainly which questions are open and which only look open.

---

## 1. The table itself

Everything below comes out of one sheet. Its rows are the divisors $d$, its columns the numbers $n$ from $1$ to $100$, and the entry in row $d$, column $n$ is the quotient $n/d$ written as a mixed number — blank when $d>n$.

| $\div$ | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| **2** | $1\tfrac12$ | **2** | $2\tfrac12$ | **3** | $3\tfrac12$ | **4** | $4\tfrac12$ | **5** |
| **3** | **1** | $1\tfrac13$ | $1\tfrac23$ | **2** | $2\tfrac13$ | $2\tfrac23$ | **3** | $3\tfrac13$ |
| **4** | | **1** | $1\tfrac14$ | $1\tfrac12$ | $1\tfrac34$ | **2** | $2\tfrac14$ | $2\tfrac12$ |
| **5** | | | **1** | $1\tfrac15$ | $1\tfrac25$ | $1\tfrac35$ | $1\tfrac45$ | **2** |

Three features are worth noticing before anything else, because all three are classical objects wearing a disguise.

**The whole entries are the divisor function.** An entry is a whole number exactly when $d \mid n$. Counting the whole entries with $2 \le d \le n$ gives $\tau(n)-1$ for each column, and over the whole sheet
$$\#\{\text{whole entries}\} \;=\; \sum_{n\le 100}\bigl(\tau(n)-1\bigr) \;=\; D(100) - 100 \;=\; 482 - 100 \;=\; 382,$$
where $D(x)=\sum_{n\le x}\tau(n)$. Highlighting the whole entries in the sheet is therefore not decoration: it draws the summatory divisor function directly.

**The column sums are the harmonic numbers.** Adding a column, from $d=2$ up to $d=n$,
$$R(n) \;=\; \sum_{d=2}^{n}\frac nd \;=\; n\bigl(H_n-1\bigr),$$
so column $5$ gives $1+1\tfrac14+1\tfrac23+2\tfrac12 = 6\tfrac{5}{12}$, and column $100$ gives $418.7378\ldots$ Asymptotically $R(n)= n\log n + (\gamma-1)n+O(\log n)$ — the same leading term as $D(n)=n\log n+(2\gamma-1)n+\Delta(n)$, with a different second constant. The two differ by exactly the sum of the discarded fractional parts, which is the subject of §7.

**The determinant is the Mertens function.** Replace every whole entry by $1$ and every fractional entry by $0$, and keep the full square array with rows and columns $1,\dots,n$. The result is the divisibility matrix
$$K_n(d,m) \;=\; \begin{cases}1 & d\mid m,\\ 0 & \text{otherwise.}\end{cases}$$

As it stands this matrix says nothing. Since $d\mid m$ forces $d\le m$, $K_n$ is upper triangular with $1$s on the diagonal: every eigenvalue is $1$ and $\det K_n = 1$, for every $n$. (It is not the identity in disguise — $\dim\ker(K_n-I)$ is the number of odd integers up to $n$, measured $10$, $30$, $61$ at $n = 20, 60, 121$ — but the spectrum is empty of arithmetic.)

Now change one column. Fill the whole of the first column with $1$s; in the language of the sheet this is a single deliberate falsehood, "every row divides $1$". Call the result $R_n$. Then
$$\det R_n \;=\; M(n) \;:=\; \sum_{k\le n}\mu(k),$$
the Mertens function. For $n=1,\dots,12$ the determinants are
$$1,\;0,\;-1,\;-1,\;-2,\;-1,\;-2,\;-2,\;-2,\;-1,\;-2,\;-2,$$
which is Mertens' sequence exactly (OEIS A002321). The proof is two lines. Because $\det K_n = 1$ we may write $\det R_n = \det\bigl(K_n^{-1}R_n\bigr)$, and $K_n^{-1}R_n$ differs from the identity only in its first column; its $(1,1)$ entry is $1+\sum_{2\le k\le n}\mu(k)$, since the Dirichlet coefficients of $\zeta(s)^{-1}$ are $\mu(k)$. Hence the determinant is $M(n)$. This is Redheffer's matrix and Redheffer's theorem (1977).

The consequence is that the sheet carries the Riemann hypothesis, which is equivalent to
$$\det R_n \;=\; O\!\left(n^{1/2+\varepsilon}\right)\qquad\text{for every }\varepsilon>0 .$$

The same modification, weighted, produces values of Dirichlet series. Put $w_k = k^{-s}$ in the first column instead of $1$; then $\det = \sum_{k\le n}\mu(k)k^{-s}$, the truncated series for $\zeta(s)^{-1}$. Since a determinant is the product of the eigenvalues, a Dirichlet-series value is literally a product of eigenvalues of a finite division table. Almost all of them are $1$: only $\lfloor\log_2 n\rfloor + 1$ eigenvalues differ from $1$, and two of those are large, of size $\pm\sqrt n + \log\sqrt n + \gamma - \tfrac12 + O(\log^2 n/\sqrt n)$.

So the sheet already contains, in its colouring and in its column sums, the two halves of the Dirichlet divisor problem; and in the determinant of its $0$–$1$ shadow, after one column is altered, a statement equivalent to the Riemann hypothesis. The rest of this note follows a *walk* across the same sheet and finds three more classical objects along the way.

> **What is classical here.** All three observations, and the third has a literature of its own. $\sum_{n\le x}\tau(n)$ and $\sum_{d\le n}1/d$ are the two most standard objects in elementary analytic number theory. The determinant is Redheffer (1977); the spectral radius $\sim\sqrt n$ is Barrett–Forcade–Pollington (1988); the Jordan structure at the eigenvalue $1$ is Robinson–Barrett (1989); the two large eigenvalues are Barrett–Jarvis (1992), sharpened by Vaughan (1993, 1996), who also bounded the small ones; and Cardon (2008) gives the general version attached to an arbitrary Dirichlet series, of which this table is the case $a_k\equiv 1$. The only contribution here is noticing that a beginner's division table displays all of it side by side, and that the distance from the sheet to the Riemann hypothesis is one column.

> **And a warning that belongs in the same breath.** The raw table is triangular, so its eigenvalues are not a place to look; everything above lives in the *determinant*, and the determinant only becomes interesting after the first column is falsified. Nothing here is a route to the Riemann hypothesis — $\det R_n = O(n^{1/2+\varepsilon})$ is not a reformulation that makes the problem easier, it is the problem, wearing a matrix.

---

## 2. A walk in the division table

Take the same sheet, with its rows labelled by a divisor $a$ and its columns by an integer $n$, and mark the cell $(a,n)$ when $a \mid n$. The mark carries the quotient $b = n/a$.

Now walk. Start in row $10$ at the cell holding the whole number $2$; that cell sits in column $20$. Step one row up and take the whole number $3$: column $27$. Continue.

| row $a$ | 2 | 3 | 4 | 5 | **6** | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| whole number $b$ | 10 | 9 | 8 | 7 | **6** | 5 | 4 | 3 | 2 |
| column $n = ab$ | 20 | 27 | 32 | 35 | **36** | 35 | 32 | 27 | 20 |

The columns rise and fall: $20, 27, 32, 35, 36, 35, 32, 27, 20$. Drawn in the table, this is an **arch**.

Two observations organise everything that follows.

**The sum is constant.** $10+2 = 9+3 = 8+4 = 12$. The arch collects exactly the pairs with $a+b = S$ for one fixed $S$; here $S = 12$.

**The steps are the odd numbers.** The differences are $+7, +5, +3, +1, -1, -3, -5, -7$. Going from $(a,b)$ to $(a{+}1, b{-}1)$ changes the product by
$$(a{+}1)(b{-}1) - ab = b - a - 1,$$
and $b-a$ drops by $2$ at each step. So the whole arch can be generated by addition alone:
$$R_{k+1} = R_k + V_k, \qquad V_{k+1} = V_k - 2,$$
starting from $R_0 = 20$, $V_0 = 7$.

If we allow the factor $1$, the arch runs from $a=1$ to $a=S-1$ and both ends sit in column $S-1$; if we exclude it, both ends sit in column $2S-4$. The maximum is column $\lfloor S^2/4\rfloor$, attained at $a=b$ when $S$ is even and at the two central pairs when $S$ is odd: for $S=11$ the arch is $18, 24, 28, 30, 30, 28, 24, 18$.

> **What is classical here.** All of it. "Product is largest when the factors are closest" is the arithmetic–geometric mean inequality, and the odd-number steps are the identity $(k{+}1)^2 - k^2 = 2k+1$. The only thing being offered is the picture.

---

## 3. Five laws, one function

Replace the coordinates $(a,b)$ by
$$S = a+b, \qquad H = b-a,$$
which we can invert as $a = (S-H)/2$, $b = (S+H)/2$. The product becomes
$$\boxed{\;N(S,H) = \frac{S^2 - H^2}{4}\;}$$
Read it as *product = the peak of balance, minus the penalty for imbalance*. At $H=0$ we get the squares.

This one formula generates a long list of laws, and it is worth seeing that they are all the same statement.

**(a) Constant-difference shift.** $N(S{+}2, H) - N(S,H) = S+1$. From arch $12$ to arch $14$: $36 \to 49$, $35 \to 48$, $32 \to 45$, every entry up by $13$.

**(b) The small-square law.** For the four corners of any $2\times2$ block in the $(S,H)$ grid, with $A=N(S,H)$, $B=N(S{+}2,H)$, $C=N(S,H{+}2)$, $D=N(S{+}2,H{+}2)$,
$$A + D = B + C, \qquad\text{e.g.}\quad \underbrace{35}_{(12,2)} + \underbrace{45}_{(14,4)} = \underbrace{48}_{(14,2)} + \underbrace{32}_{(12,4)} = 80 .$$

**(c) The mean-of-neighbours law.** Every interior point is the average of its four neighbours:
$$4N = N_{\uparrow} + N_{\downarrow} + N_{\leftarrow} + N_{\rightarrow}, \qquad \frac{24+36+32+48}{4} = 35 .$$

**(d) Depth below the ceiling is a perfect square.** In column $S$ the entry at index $a$ sits below the apex $S^2/4$ by exactly
$$\frac{S^2}{4} - a(S-a) = \left(\frac{S}{2} - a\right)^2 .$$
For $S=20$: entries $100, 99, 96, 91, 84, 75, 64, 51, 36$ at depths $0,1,4,9,16,25,36,49,64$.

**(e) Rows are straight lines.** With $a$ fixed, $N = aS - a^2$ — a line of slope $a$. And this line is precisely the **tangent** to the parabola $N = S^2/4$ at $S = 2a$. So the table is exactly the lattice points on the tangent lines drawn at integer parameters $a\ge1$, taken on the side $S \ge 2a$; the tangency points are the squares $a^2$.

Now the collapse. Write $z = S + iH$. Then
$$N = \mathrm{Re}\!\left(\frac{z^2}{4}\right).$$

Everything above is a property of that single function:

| law | why |
|---|---|
| (b) small-square | $N$ is additively separable, $f(S) + g(H)$, so the mixed second difference is $0$ |
| (c) mean of neighbours | $\partial_S^2 N = +\tfrac12$, $\partial_H^2 N = -\tfrac12$, so $\Delta N = 0$: $N$ is **harmonic**, and (c) is the mean-value property |
| (a), (d), (e) | derivatives of a quadratic |
| $\partial_a \partial_b (ab) = 1$ | the same two second derivatives with the other sign: $\partial_S^2 - \partial_H^2 = 1$ |

> **A precision that is easy to miss.** $S$ and $H$ always have the same parity, since $S + H = 2b$. The $(S,H)$ points are therefore not all of $\mathbb{Z}^2$ but the index-2 sublattice $\{S \equiv H \bmod 2\}$ — a checkerboard. Any count performed in $(S,H)$ must be halved relative to the naive box.

> **What is classical here.** All of it, and in one line: $(S,H)$ are the principal axes of the quadratic form $ab$, and $\mathrm{Re}(z^2)$ is the standard harmonic function whose level curves are the rectangular hyperbolas. Laws (a)–(e) are exercises about it. The value of the section is negative and useful: a list of "discoveries" is one fact, and noticing that saves the reader from thinking the list is long.

---

## 4. The arch is Fermat's method

With $x = S/2$ and $y = H/2$, the boxed formula reads $N = x^2 - y^2$. So:

> Column $N$ lies on arch $S$ **exactly when** $S^2 - 4N$ is a perfect square.

(Stated this way the assertion needs no parity caveat: $S^2-4N = H^2$, and $S,H$ automatically share a parity.)

That is Fermat's difference-of-squares factorisation (1643). Scanning arches upward from the first one whose apex clears $N$ is literally Fermat's loop, and the "apex gap"
$$g(S) = \frac{S^2}{4} - N$$
is the sieving polynomial $x^2 - N$.

Only arches whose parity matches can contain $N$, so one advances $S$ by two each time; counting in Fermat's index $x=S/2$, the number of steps is exactly
$$\frac{p+q}{2} - \lceil \sqrt{N} \rceil \;\approx\; \frac{(\sqrt{q}-\sqrt{p})^2}{2}.$$
With explicit primes, so that the reader can check every row:

| $p$ | $q$ | $N$ | steps | $(\sqrt q-\sqrt p)^2/2$ |
|---|---|---|---|---|
| 999983 | 1000003 | 999985999949 | **0** | $5\times10^{-5}$ |
| 100003 | 10000019 | 1000031900057 | 4049995 | $4.05\times10^6$ |
| 1009 | 1000000007 | 1009000007063 | 498996018 | $4.99\times10^8$ |
| 101 | 10000000019 | 1010000001919 | 4998995072 | $5.00\times10^9$ |

Instant when the factors are balanced; hopeless otherwise. Trial division has exactly the opposite profile, and Lehman (1974) combines them for $O(N^{1/3})$.

There is a companion statement worth recording because it makes primality a measurement rather than a negation. Column $n$ lies on arch $S$ iff $a \mid n$ and $S = a + n/a$, so the arches through $n$ are exactly
$$\{\,d + n/d \;:\; d \mid n \,\},$$
the **arch spectrum** of $n$. Its top is always $n+1$ (from $d=1$); its bottom is $\min_d (d + n/d) \ge 2\sqrt{n}$, with equality only for perfect squares; and $n$ is prime precisely when the spectrum is the single value $n+1$. The bottom of the spectrum is Fermat's target, so
$$\text{Fermat steps} = \frac{S_{\min} - 2\lceil\sqrt n\rceil}{2}.$$
The average spectrum size over $n \in [10^{k-1}, 10^k]$ measures $3.67, 4.81, 5.96, 7.11$ for $k = 3,4,5,6$, tracking $\tfrac12 \log n$.

> **What is classical here.** Fermat's method is from 1643. The identity $\min_d(d+n/d)$ and its relation to the divisor nearest $\sqrt n$ is standard, and the statistics of that divisor are Erdős's multiplication-table phenomenon, refined by Ford. The reformulation "primality = a one-element spectrum" is a restatement, not a criterion: computing the spectrum requires factoring.

---

## 5. Why it fails, and the repair — worked by hand

Take $N = 2041$. Since $\lceil\sqrt{2041}\rceil = 46$, list the arches from there with their apex gaps, factored over $\{2,3,5,7\}$:

| index $x=S/2$ | apex $x^2$ | gap $x^2 - N$ | factorisation | parity of exponents |
|---|---|---|---|---|
| 46 | 2116 | 75 | $3 \cdot 5^2$ | **0100** |
| 47 | 2209 | 168 | $2^3\cdot 3\cdot 7$ | 1101 |
| 48 | 2304 | 263 | prime | — |
| 49 | 2401 | 360 | $2^3\cdot 3^2\cdot 5$ | 1010 |
| 50 | 2500 | 459 | $3^3\cdot 17$ | — |
| 51 | 2601 | 560 | $2^4\cdot 5\cdot 7$ | **0011** |
| 52 | 2704 | 663 | $3\cdot 13\cdot 17$ | — |
| **53** | 2809 | 768 | $2^8\cdot 3$ | **0100** |
| 54 | 2916 | 875 | $5^3\cdot 7$ | 0011 |

Fermat wants a gap that is *already* a square. None is; he must climb to $x = 85$, thirty-nine steps away.

But arches $46$ and $53$ share the parity pattern $0100$. Multiply them:
$$46 \cdot 53 = 2438, \qquad 75 \cdot 768 = 57600 = 240^2 .$$
Each arch says apex $\equiv$ gap $\pmod N$, so the product says $2438^2 \equiv 240^2 \pmod{2041}$, and
$$\gcd(2438 - 240,\; 2041) = \gcd(2198, 2041) = \boxed{157}.$$
The pair $51, 54$ (pattern $0011$) delivers $13$ the same way. **Two arches out of the first eight, against thirty-nine for Fermat.**

This is the move that matters, and it can be stated in one sentence: *a square need not be found, it can be built.* Squares are closed under multiplication, so a set of gaps is a candidate exactly when the sum of their exponent vectors vanishes mod $2$ — a nullspace computation over $\mathbb{F}_2$. A candidate can still fail, when $X \equiv \pm Y$ and the gcd comes out $1$ or $N$; one then takes the next vector in the nullspace.

Two further facts make the method fast.

**Striping.** The gap is a polynomial in the arch index, so $g(i+p) \equiv g(i) \pmod p$: divisibility by $p$ depends only on $i \bmod p$. For the larger example $N = 2578692190013 = 5519 \times 467239027$ (chosen unbalanced, so that Fermat alone would need about $2.3\times10^8$ steps) and $p=7$, seven evaluations locate the roots $i \equiv 0, 1 \pmod 7$, and thereafter one *jumps*. Over $300{,}000$ arches with $222$ factor-base primes, trial division costs $66.6$ million divisions and the striped sieve costs $484{,}562$ additions — a factor of $137$.

*A family stripes if and only if its gap is a polynomial function of the index*, and the number of stripes per prime is the number of roots of that polynomial mod $p$. Tested on the classical quadratic (0 or 2), an MPQS parabola $ax^2+2bx+c$ with $b^2-ac=N$ (0 or 2), and an arbitrary cubic (0, 1 or 3).

**Size.** A random number of size $Y$ is $B$-smooth with probability about Dickman's $\rho(u)$, $u = \log Y/\log B$, and $\rho$ falls like $u^{-u}$:

| $u$ | 1.5 | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 |
|---|---|---|---|---|---|---|
| $\rho(u)$ | 0.5945 | 0.3069 | 0.1303 | 0.0486 | 0.0162 | 0.0049 |

Against measurement at $B = 3000$ on the same $N$, smoothness tested against every prime up to $B$:

| family | geometric-mean gap | $u$ | $\rho(u)$ | measured smooth rate |
|---|---|---|---|---|
| QS, one interval $x_0 \pm 10^5$ | $1.2\times 10^{11}$ | 3.19 | 0.0324 | 0.0284 |
| MPQS, fresh polynomial every $2M$ steps, $M=10^3$ | $1.1\times 10^{9}$ | 2.61 | 0.1062 | 0.1120 |
| CFRAC, continued-fraction residues | $1.0\times 10^{6}$ | 1.73 | 0.4519 | 0.388 |

The predicted column is $\rho$ at the $u$ of the *geometric mean*, which is a one-point stand-in for an average over a spread of sizes; that it lands within about $15\%$ of the measurement in all three rows is as much as the approximation deserves.

CFRAC keeps $\lvert \mathrm{gap}\rvert < 2\sqrt N$ always — nearly fourteen times the smooth rate of the single long interval — but its residues are *not* polynomial in the index, so it cannot be sieved: every candidate must be trial-divided. That single trade-off, **residue size against sievability**, separates every algorithm in the family.

Balancing the two costs — about $\pi(B)$ relations needed, $\pi(B)/\rho(u)$ candidates to find them, and $\pi(B)^2$ for the linear algebra — has a minimum, and optimising it in general gives
$$B \sim \exp\!\left(\tfrac12\sqrt{\log N \log\log N}\right), \qquad \text{cost } \exp\!\left((1+o(1))\sqrt{\log N \log\log N}\right).$$

> **What is classical here.** Everything. The relaxation is Kraitchik's (1920s); the continued-fraction version is Morrison–Brillhart (1975); the sieve is Pomerance's quadratic sieve (announced 1981, published 1982--84); the multiple-polynomial variant is Montgomery (1985); the two-dimensional successor is the number field sieve. The complexity is standard. The only contribution here is the $2041$ table, which is small enough to check with a pencil.

---

## 6. Turn the picture forty-five degrees

Set $a = x-y$, $b = x+y$. Then $x^2 - y^2 = ab$: **the arch family and the division table are the same lattice-point problem, rotated.** The arch fixes the sum and traces the product; its conjugate fixes the product $n=ab$ and traces the sum, $S(d) = d + n/d$, a curve with a *valley* at $d = \sqrt n$ where the arch had a *peak*. With $S = d+n/d$ and $H = n/d - d$ one has $S^2 - H^2 = 4n$ identically, so both families live on the single saddle $S^2 - H^2 = 4N$, curving down in one principal direction and up in the other. (They are *not* the two orthogonal families of $z\mapsto z^2$: those are $\mathrm{Re}(z^2)$ and $\mathrm{Im}(z^2)$ constant, i.e. $N$ constant and $SH/2$ constant. The arch, $S$ constant, is a vertical line in the $z$-plane and meets the level curves of $N$ at no fixed angle.)

So the arch inherits a counting problem. Counting ordered pairs $(d,e)$ with $de=n$ and $d\equiv e \bmod 2$ — equivalently, signed representations $n = x^2-y^2$ — gives
$$\rho^{\mathrm{ord}}(n) = \begin{cases} \tau(n) & n \text{ odd}\\ 0 & n \equiv 2 \ (4)\\ \tau(n/4) & 4 \mid n\end{cases}$$
(verified against direct search for $n < 3000$; the unordered count, which is what the literature usually calls $\rho(n)$, is $\lceil \rho^{\mathrm{ord}}/2\rceil$). Its summatory function has main term
$$\sum_{n \le X}\rho^{\mathrm{ord}}(n) = \tfrac12 X\log X + \left(\gamma - \tfrac12\right)X + \text{error},$$
derived by splitting into the odd-divisor sum, whose Dirichlet series is $\big((1-2^{-s})\zeta(s)\big)^2$, plus $D(X/4)$. Computed exactly to $X = 10^{12}$ by the hyperbola method, the error at $X = 10^6, 10^9, 10^{12}$ is $+81$, $+285$, $+1504$ against main terms of order $10^7, 10^{10}, 10^{13}$.

Three counting problems now sit side by side, all lattice points in a plane region, main term equal to the area, and the same open error exponent:

| region | on-curve count | main term | constant produced | typical error size | decade-max exponent |
|---|---|---|---|---|---|
| $a^2+b^2 \le x$ | $r_2(n)$, via $\chi_4$ | $\pi x$ | $\pi$ | $1.59\,x^{1/4}$ | $+0.308$ |
| $ab \le x$ | $\tau(n)$ | $x\log x + (2\gamma-1)x$ | $\gamma$ | $0.98\,x^{1/4}$ | $+0.307$ |
| $x^2-y^2 \le X$ | $\rho^{\mathrm{ord}}(n)$ | $\tfrac12 X\log X + (\gamma-\tfrac12)X$ | $\gamma - \tfrac12$ | $0.79\,x^{1/4}$ | $+0.312$ |

Two different measurements, and they must not be conflated. The **typical** size is exactly $x^{1/4}$: the standard deviation of $E(x)/x^{1/4}$ locks and does not drift. The **decade maximum** grows faster, with exponent about $0.31$ — and it has to: Hardy proved the error is not $O(x^{1/4})$ for $\tau$ and for $r_2$, and Kühleitner's omega theorem does the same for $\rho$. For the arch, the maximum of $|E|$ over each decade — $27.2$, $56.9$, $119.1$, $234.3$ over $10^3$ to $10^7$ — divided by the fourth root of the *top* of that decade gives $2.72, 3.20, 3.77, 4.17$: visibly unbounded. (The convention matters: normalising pointwise by $x^{1/4}$ at the maximising $x$ instead gives $3.01, 3.57, 4.10, 4.37$, the same trend one decimal higher.)

*A trap worth flagging.* Decade statistics need whole decades. Truncating the last one — running to $2\times10^6$ instead of $10^7$, say — suppresses its maximum and pulls the fitted exponent down to about $0.25$: dangerously close to the conjectured value, and therefore easy to accept.

> **What is classical here, and this is where the honesty matters most.** The third row is not new. The function $\rho(n)$, its summatory function, the error term, its mean square, and its close relation to $d(n)$ are the subject of a published literature: Kühleitner's omega theorem on differences of two squares (1992, with a second part in 1999), and the mean-square results of Kühleitner and Nowak, *On differences of two squares*, Cent. Eur. J. Math. **4** (2006), 110–122, which states the connection to $d(n)$ in its abstract. The surrounding work on sums and differences of $k$-th powers (Kühleitner–Nowak–Schoissengeier–Wooley) and on the general two-dimensional divisor problem (Zhai–Cao) belongs to the same family.
>
> A related caution about what a measurement can look like. Normalising the error by $x^{1/4}$ gives a distribution with a limiting law (Heath-Brown), and measuring its *skew* here gives a positive value for the divisor problem and a negative one for the circle; the arch's skew measures $+0.22$, apparently placing it on the divisor side. But $E_{\mathrm{arch}} = E_{\mathrm{odd}} + \Delta(X/4)$ exactly, so the third moment splits into four terms, and they give $-0.004 + 0.048 + 0.005 + 0.056 = +0.105$ — which reproduces the directly measured third moment, and hence the skew, since skew $=$ third moment$/\sigma^3$ and $\sigma = 0.79$ gives $0.105/0.79^3 = 0.21$. $\Delta$'s own third moment supplies $53\%$ of it. The fingerprint was never the arch's; it was inherited. The same decomposition dissolves the mean-square constant: $0.4266 = 0.1725 + 0.3232 - 0.0691$, i.e. $c_{\mathrm{odd}} + c_\Delta/2 + \text{cross}$, every piece already in the divisor literature.

---

## 7. The wobbles, and where $x^{1/4}$ comes from

Count the cells under $ab \le x$ for $x = 100$: the exact count is $482$, the area is $476.0$, the error is $+6.0$. The hyperbola crosses one cell per column, but the symmetry $a\leftrightarrow b$ reduces the independent boundary terms to about $2\sqrt{100}=20$; naively each could contribute $\pm1$, so the error might have been $20$.

Where does it actually live? Since $\lfloor x/a\rfloor = x/a - \{x/a\}$, the error is a sum of *wobbles* $\{x/a\} - \tfrac12$. If every leftover were exactly one half, there would be no error at all.

| $a$ (selected rows) | $100/a$ | $\lfloor\cdot\rfloor$ | leftover | wobble | running sum over all $a$ |
|---|---|---|---|---|---|
| 1 | 100.000 | 100 | 0.000 | $-0.500$ | $-0.500$ |
| 3 | 33.333 | 33 | 0.333 | $-0.167$ | $-1.167$ |
| 6 | 16.667 | 16 | 0.667 | $+0.167$ | $-2.000$ |
| 8 | 12.500 | 12 | 0.500 | $0.000$ | $-2.214$ |
| 11 | 9.091 | 9 | 0.091 | $-0.409$ | $-3.512$ |

They come in both signs and they cancel. At $x = 10^4$ there are $100$ wobbles of size at most $\tfrac12$, so an adversarial arrangement would reach $50$; the actual running sum is $-10.22$.

**And that is where the exponent comes from.** About $\sqrt x$ wobbles behaving like independent coin flips of size $\tfrac12$ sum to about $\sqrt{\sqrt x} = x^{1/4}$: square-root cancellation, nothing more. The measurement confirms it for the *typical* size — the standard deviation of $E(x)/x^{1/4}$ locks at $1.59$ for the circle, $0.98$ for the divisor problem, $0.79$ for the arch, none of them drifting. This is a theorem, not a coincidence: the mean square of the error is $\asymp x^{1/2}$ (Cramér, Preissmann). Coins can nevertheless run hot, and the maximum genuinely exceeds $x^{1/4}$; see §6 and §13.


**And the wobbles have a mean.** Summing them rather than watching them cancel gives a second, cleaner handle on the same object. Write $W(n)=\sum_{d\le n}\{n/d\}$. Since $\sum_{d\le n} n/d = nH_n$ and $\sum_{d\le n}\lfloor n/d\rfloor = D(n)$,
$$W(n) \;=\; nH_n - D(n) \;=\; (1-\gamma)\,n + \tfrac12 - \Delta(n) - \frac{1}{12n} + O(n^{-3}),$$
using $H_n=\log n+\gamma+\frac1{2n}-\frac1{12n^2}+\cdots$ and $D(n)=n\log n+(2\gamma-1)n+\Delta(n)$. Verified numerically at $n=10^3,\dots,10^7$: the two sides agree to $8\times10^{-5}, 8\times10^{-6}, 8\times10^{-7}, 8\times10^{-8}$, which is $1/(12n)$ to three figures at each step, so nothing in the measurement is left over.

Two consequences. First, the average wobble is not $\tfrac12$, as an equidistribution guess would give, but
$$\frac{W(n)}{n}\;\longrightarrow\;\int_0^1\Bigl\{\frac1t\Bigr\}dt=\int_1^\infty\frac{\{u\}}{u^2}\,du=1-\gamma=0.42278\ldots$$
the last equality in two lines: $\int_1^N\{u\}u^{-2}du=\log N-\sum_{k<N}k\bigl(\tfrac1k-\tfrac1{k+1}\bigr)=\log N-(H_N-1)$. The shortfall from $\tfrac12$ is $\gamma-\tfrac12=0.0772\ldots$ — the same constant that appeared in the main term of §6, and for the same reason.

Second, rearranged as $\gamma = 1 - W(n)/n + \bigl(\tfrac12-\Delta(n)\bigr)/n$, this *computes* $\gamma$ from the fractional parts alone, with an accuracy governed exactly by $\Delta(n)$. Summing only the discarded fractions of the division table and subtracting from $1$ gives $0.583529$, $0.577354$, $0.577225$ at $n=10^3,10^5,10^7$ against $\gamma=0.5772157$. The deviations are $+6.31\times10^{-3}$, $+1.38\times10^{-4}$, $+9.37\times10^{-6}$, and each equals $(\Delta(n)-\tfrac12)/n$ to every digit shown, with $\Delta(n) = 6.81, 14.32, 94.19$ — the identity is exact, so the estimate is not merely close to $\gamma$, it is $\gamma$ plus a known quantity. That buys four correct digits at $n=10^7$, and no further digits without progress on the divisor problem.

> **What is classical here.** All of it. $\int_1^\infty\{u\}u^{-2}du=1-\gamma$ is a standard exercise, and the displayed identity is Dirichlet's derivation run backwards. The point worth keeping is the reading: $\gamma$ is *defined* as the accumulated gap between the staircase $\sum 1/d$ and the ramp $\log n$, and $\{n/d\}$ is that same gap measured entry by entry. Its appearance in a division table is not a coincidence but a restatement.

**The error is not, however, caused by the cells being square.** Same lattice, four different boundaries:

| boundary | max $\lvert E\rvert$ over the decades $10^3\!-\!10^4,\ldots,10^6\!-\!10^7$ | growth exponent |
|---|---|---|
| line, slope $3/2$ | $0.25, 0.25, 0.25, 0.25$ | $+0.000$ |
| line, slope $\varphi$ | $1.01, 1.23, 1.51, 1.71$ | $+0.076$ |
| line, slope $\pi$ | $32.6, 38.6, 39.2, 40.7$ | $+0.032$ |
| hyperbola $ab \le x$ | $35.5, 85.2, 169.2, 296.6$ | $+0.307$ |
| circle $a^2+b^2\le x$ | $50.0, 104.9, 209.6, 419.9$ | $+0.308$ |

Every decade here is a whole decade, and that is not a detail: the maximum for slope $\varphi$ is $1.53$ if the last one is stopped at $2.3\times10^6$ and $1.71$ if it is run to $10^7$, and for slope $\pi$ it is $39.2$ against $40.7$. The trap flagged in §6 catches the straight boundaries first, because their maxima arrive late and rarely.

A straight boundary has one fixed slope, so the lattice meets it periodically and the continued fraction of the slope controls everything. A curved boundary's slope keeps changing, so it samples every direction and no single continued fraction helps. **Curvature is the cause.** The exponent ladder follows from the geometry of the elbow at $a=b=\sqrt x$, where the boundary has length $\sim\sqrt x$ and curvature $\sim 1/\sqrt x$: the trivial bound is the length $x^{1/2}$; van der Corput gives $\text{length}^{2/3} = x^{1/3}$; Huxley (2003) gives $131/416 = 0.31490\ldots$ and Bourgain--Watt $517/1648 = 0.31371\ldots$; and the conjecture $x^{1/4}$ is $\sqrt{\text{length}}$, full cancellation.

There is an exact formula for the error — Voronoï's,
$$\Delta(x) = \frac{x^{1/4}}{\pi\sqrt2}\sum_{n}\tau(n)\,n^{-3/4}\cos\!\left(4\pi\sqrt{nx} - \frac{\pi}{4}\right),$$
whose coefficients are the divisor counts themselves. Truncating at $10^2, 10^3, 10^4, 10^5, 4\times10^5$ terms at $x = 5000.3$ gives $6.18, 9.26, 16.33, 15.05, 15.00$ against the measured $14.976$. So the error is not random; it has a closed form. What resists is that a straight line has *one* frequency while this curve has infinitely many incommensurate frequencies $\sqrt n$, and bounding their sum is the open problem.

> **What is classical here.** All of it: the $\psi$-sum representation, van der Corput's method, Huxley's exponent, Voronoï's formula (and Hardy's for the circle), the mean-square results of Cramér and Preissmann, Hardy's $\Omega$-theorem, and Heath-Brown's limiting distribution. The exposition is the contribution.

---

## 8. Soften the edge and the error disappears

Two natural attempts to remove the error, one instructive failure and one instructive success.

**Straighten the boundary.** The map $(a,b)\mapsto(\log a,\log b)$ turns $ab \le x$ into the straight line $\log a + \log b \le \log x$. But the image of the lattice is no longer a lattice: over $a \le 2\times10^5$ the ratio of largest to smallest gap between consecutive $\log a$ is $1.4\times10^5$. The difficulty moves; it does not leave.

**Soften the edge.** Keep the curve, but replace the sharp cut $n \le X$ by the weight $e^{-n/X}$. Mellin inversion, with poles at $s=1$ and $s=0$, gives
$$\sum_{n\ge1}\tau(n)\,e^{-n/X} = X(\log X + \gamma) + \tfrac14 - \frac{1}{144\,X} + O(X^{-3}) ,$$
the double pole at $s=1$ giving the first term, $\Gamma$'s pole at $s=0$ giving $\zeta(0)^2 = \tfrac14$, and its pole at $s=-1$ giving $-\zeta(-1)^2 X^{-1} = -1/(144X)$. There is no oscillation left at all; the remaining terms come from $\zeta(-3), \zeta(-5),\ldots$ and every $\zeta(-2k)$ vanishes.

| $X$ | error against $X(\log X+\gamma)+\tfrac14$ | that error $\times\,144X$ | error, sharp cut $\approx x^{1/4}$ |
|---|---|---|---|
| $10^3$ | $-6.94\times10^{-6}$ | $-1.000$ | 5.62 |
| $10^4$ | $-6.94\times10^{-7}$ | $-1.000$ | 10.00 |
| $10^5$ | $-6.94\times10^{-8}$ | $-1.000$ | 17.78 |

The middle column is the point: what is left after the two obvious poles is not noise of unknown size but a single explicit term, and the measurement pins it to four figures.

Same curve, same lattice, same counting. **The entire $x^{1/4}$ is the price of the sharp yes/no question** — and the price of removing it is that one is no longer counting exactly.

---

## 9. Circles in every cell

Here is a construction whose answer is pleasant and whose *reason* is the same as the previous section's. Inscribe a circle of radius $r$ at the centre of every unit cell, and measure what fraction of the hyperbola's arc length falls inside the circles.

| $r$ | measured fraction at $x=10^7$ | $\pi r^2$ |
|---|---|---|
| 0.15 | 0.070641 | 0.070686 |
| 0.25 | 0.196178 | 0.196350 |
| 0.35 | 0.384100 | 0.384845 |
| 0.50 | 0.785426 | 0.785398 |

The law is **arc fraction $=$ area density**, and at $r = \tfrac12$ it returns $\pi/4$. Measured at $x = 10^3,\dots,10^8$: $0.782795$, $0.785637$, $0.788086$, $0.785897$, $0.785466$, $0.785402$ against $\pi/4 = 0.785398$.

Why is the deviation so small? Because this construction is a *smoothing*. The sharp count is a sum of sawtooth terms $\psi(t) = \{t\}-\tfrac12$, whose Fourier coefficients decay like $m^{-1}$ (measured $m^{-0.99}$) — barely summable, so the cancellation is delicate and the answer is $x^{1/4}$. The circle replaces the sawtooth by a semicircular chord weight whose coefficients decay like $m^{-3/2}$ (measured $m^{-1.52}$, a Bessel $J_1$), and the sum then converges fast. **The circles and the $e^{-n/X}$ weight are the same move reached from two directions.**

> **What is classical here.** The law is Weyl equidistribution: the fraction of a suitably equidistributed curve's length inside any lattice-periodic set is that set's density. The $\pi/4$ is then arithmetic. The Fourier decay of a semicircle is a Bessel function.

---

## 10. The sieve as a clock of knowledge

The walk is finished. This section and the next take a different departure from the same sheet of §1: instead of walking across the table, cut it off after a given row. That single change turns the table into the sieve of Eratosthenes, and it reaches a second family of classical objects — Buchstab's function here, the wheel sieve in §11 — by a route that never touches the arch.

The object it needs was half-introduced already: §5 needed *smooth* numbers — those whose prime factors all lie in an early part of the table — and Dickman's $\rho$ to count them. The complementary object costs nothing extra to describe and turns the table into something that changes with time.

Use only a prefix of the sheet: fix a cutoff $z$ and keep the rows $d \le z$. A column $n$ is then **struck** if some kept row divides it, and **standing** otherwise. Truncating this way changes what the word *prime* means. A standing column is not a prime; it is a number *not yet known to be composite*, and the verdict can be overturned by a row that has not been consulted. Write

$$\Phi(x,z) \;=\; \#\{\,2 \le n \le x \;:\; p \mid n \Rightarrow p > z\,\}$$

for the standing columns up to $x$. What governs the picture is neither $x$ nor $z$ but the single ratio

$$u \;=\; \frac{\log x}{\log z},$$

which counts how many multiplicative lengths the column index runs beyond the row range. Read it as a clock.

**Noon, $u = 2$.** A composite has a prime factor at most $\sqrt n$, so for $z \ge \sqrt x$ every standing column is prime and $\Phi(x,z) = \pi(x)-\pi(z)$ — a set identity, not an approximation: at $x = 10^6$ the two sets agree element by element, and at $x = 10^8$, $z = 10^4$ both counts are $5{,}760{,}226$. The first column the truncated sheet gets *wrong* is the square of the first unused row: $9, 25, 49, 121, 169, 289, 361, 529$ for rows kept up to $2,3,5,7,11,13,17,19$. So $p_{\text{next}}^2$ is not only where a new row starts contributing; it is where the previous sheet's certainty expires. The clock reading $u=2$ and the square $p^2$ are one statement in two coordinates.

**After noon.** For $u>2$ the standing columns mix primes with not-yet-detected composites, and their density is governed by Buchstab's function,

$$\omega(u) = \frac1u \;\; (1\le u\le 2), \qquad \bigl(u\,\omega(u)\bigr)' = \omega(u-1) \;\; (u>2),$$

whose next branch is $(1+\log(u-1))/u$ and which settles onto $e^{-\gamma}$ very fast: $\omega(6) = 0.561459685$ against $e^{-\gamma} = 0.561459484$, agreeing to $1.5\times10^{-10}$ by $\omega(12)$. Normalising by the wheel density $\prod_{p\le z}(1-1/p)$ rather than $1/\log z$, at $x = 10^8$:

| $z$ | $u$ | $\Phi/(x\prod(1-1/p))$ | $e^{\gamma}\omega(u)$ | rel. err. |
|---|---|---|---|---|
| 7 | 9.466 | 1.00000 | 1.00000 | $0.00\%$ |
| 31 | 5.364 | 1.00000 | 1.00000 | $0.00\%$ |
| 101 | 3.991 | 0.99977 | 0.99999 | $-0.02\%$ |
| 331 | 3.175 | 1.00360 | 1.00070 | $0.29\%$ |
| 1,009 | 2.663 | 0.99870 | 1.00901 | $-1.02\%$ |
| 3,163 | 2.286 | 0.95387 | 0.97504 | $-2.17\%$ |
| 10,007 | 2.000 | 0.94618 | 0.89060 | $6.24\%$ |
| 31,623 | 1.778 | 1.06302 | 1.00185 | $6.11\%$ |

The two ends of the clock are worth naming. At large $u$ the ratio is $1.00000$ and so is $e^{\gamma}\omega(u)$ — for small $z$ the standing columns are exactly those coprime to the primorial of $z$, so their density *is* $\prod(1-1/p)$ by periodicity (matched to $6\times10^{-8}$ at $z = 3,7,13,31$), and $\omega$ has already reached $e^{-\gamma}$. There the sheet is a wheel and the clock adds nothing to its period. At $u<2$ the $6\%$ is not a failure of $\omega$ but the prime-counting secondary term, and it separates exactly:

$$\frac{\Phi(x,z)}{x\prod(1-1/p)} \Big/ e^{\gamma}\omega(u)
\;=\; \underbrace{\frac{(\pi(x)-\pi(z))\log x}{x}}_{A} \cdot \underbrace{\frac{1}{e^{\gamma}\prod(1-1/p)\log z}}_{B},$$

with $A$ the excess of $\pi$ over $x/\log x$ and $B$ the Mertens remainder. The decomposition closes to $10^{-10}$: at $z = 16{,}237$, $50{,}802$, $517{,}947$ the ratio reads $1.06208$, $1.06068$, $1.05352$ against $A = 1.06095$, $1.06034$, $1.05340$ and $B$ within $10^{-3}$ of $1$. The whole discrepancy is $A$ — the same finite-size excess that runs through §7.

**How much of what stands is true.** The share of standing columns that really are prime tends to a function of $u$ alone, $\pi(x)/\Phi(x,z) \to 1/(u\,\omega(u))$, since the primes contribute density $1/(u\log z)$ and the standing columns $\omega(u)/\log z$. At $x = 10^8$:

| $z$ | $u$ | standing | primes | share | $1/(u\omega(u))$ |
|---|---|---|---|---|---|
| 10,007 | 2.000 | 5,760,225 | 5,760,225 | 1.00000 | 1.00000 |
| 1,009 | 2.663 | 8,078,009 | 5,761,286 | 0.71321 | 0.66280 |
| 101 | 3.991 | 11,909,899 | 5,761,429 | 0.48375 | 0.44624 |
| 31 | 5.364 | 15,285,209 | 5,761,444 | 0.37693 | 0.33203 |
| 13 | 7.182 | 19,180,819 | 5,761,449 | 0.30038 | 0.24800 |
| 7 | 9.466 | 22,857,141 | 5,761,451 | 0.25206 | 0.18815 |
| 3 | 16.767 | 33,333,332 | 5,761,453 | 0.17284 | 0.10622 |

The share decreases strictly in $u$ — the further the column index runs past the row range, the less the truncated sheet knows — and exceeds the limit everywhere. That excess is finite-size: holding $u=3$ and growing $x$ over four decades it falls $+0.08498$, $+0.05922$, $+0.04893$, $+0.04178$, monotonically. The rate is not pinned by this range (the products excess $\times \log x$ read $0.98, 0.82, 0.79, 0.77$, still drifting), and the slow convergence is the honest caveat: the limit is asymptotic in $x$ at fixed $u$, so the entries of the density table above with $u \ge 10$ sit outside that regime, where the wheel description, not $\omega$, is the operative one.

**The dual of §5.** The two sections are halves of one picture: §5 counts columns whose factors have all collapsed into the early rows, this one counts columns that have escaped them entirely; $\rho$ against $\omega$, $\rho \to 0$ super-exponentially against $\omega \to e^{-\gamma}$, boundary branches $\rho \equiv 1$ on $[0,1]$ against $\omega = 1/u$ on $[1,2]$, delay equations $u\rho'(u) = -\rho(u-1)$ against $(u\omega(u))' = \omega(u-1)$ — and the same clock $u$ in both.

> **What is classical here.** All of it. $\Phi(x,z)$ and $\omega$ are Buchstab (1937); the closed branches, the limit $e^{-\gamma}$, and the duality with $\rho$ are textbook sieve theory; $u \le 2 \Rightarrow$ prime is why Eratosthenes stops at $\sqrt x$; the large-$u$ regime is the wheel. What the section adds is a reading: the sieve parameter $z$, which the literature carries as a technical bound inside an inequality, is presented here as a time, and the sheet of §1 is what the time is read on.

---

## 11. Dynamic row inheritance

§10 read the truncation level $z$ as a time. This section reads the *rows themselves* as objects that move, and reports the one identity that governs all of it. Everything here is elementary; the point is that a single line explains a whole family of measurements that look independent.

Keep the convention of §1: the rows of the sheet are the divisors, the columns are the numbers. Work on the odd numbers and let the odd prime rows enter in order. Row $p$ strikes the columns divisible by $p$ that no smaller row has struck, so every odd composite is assigned to exactly one row — its least prime factor — and the assignment partitions the odd composites without overlap. Write

$$A_p = \{\, n : \mathrm{spf}(n) = p \,\}, \qquad S_p^- = \{\, m \text{ odd} : q \nmid m \text{ for every prime } q < p \,\}$$

for the row's own contributions and for the survivors standing just before it enters. Every product below runs over the odd primes only, the even numbers being outside the model.

### 11.1 One identity

$$\boxed{\;A_p \;=\; p \cdot S_p^-\;}$$

with $m$ running from $p$ upward. The proof is the definition: $\mathrm{spf}(n) = p$ exactly when $n = pm$ with $m$ free of every prime below $p$. Verified elementwise for the first $16{,}666$ contributions of each column $p = 5, \dots, 31$.

Two corollaries follow immediately, and they are the reason the identity is worth stating.

**The row's step pattern is the previous stage, magnified.** Writing $a_p(k) = p\,m_k$ for the ordered contributions,

$$\frac{a_p(k+1) - a_p(k)}{2p} \;=\; \frac{m_{k+1}-m_k}{2},$$

so every step of row $p$ is a multiple of $2p$, and after dividing by $2p$ the step sequence *is* the gap sequence of $S_p^-$ in units of the odd step. Checked elementwise on the first $1000$ steps for $p = 5, 7, 11, 13, 17$: identical in every position, not merely in distribution.

**The reach of a row is a survivor position.** If $K$ contributions are wanted, the row needs exactly the $K$-th survivor:

$$X_p(K) \;=\; p\,m_p(K),$$

an identity, not an estimate.

### 11.2 Equal production instead of a common cutoff

Comparing all rows up to one fixed bound $X$ confounds a row's strength with the room it was given. The fair comparison is to let each row run until it has produced the same number $K$ of new composites. Take $K = 16{,}666$, row $3$'s yield to $10^5$:

| row $p$ | $X_p(K)$ | $X_p/X_{p'}$ | $p/(p'-1)$ | efficiency | $\prod_{2<q<p}(1-1/q)$ |
|---|---|---|---|---|---|
| 3 | 99,999 | — | — | 1.0000 | 1.0000 |
| 5 | 249,995 | 2.5000 | 2.5000 | 0.6667 | 0.6667 |
| 7 | 437,507 | 1.7501 | 1.7500 | 0.5333 | 0.5333 |
| 11 | 802,043 | 1.8332 | 1.8333 | 0.4572 | 0.4571 |
| 13 | 1,042,691 | 1.3000 | 1.3000 | 0.4156 | 0.4156 |
| 17 | 1,477,181 | 1.4167 | 1.4167 | 0.3836 | 0.3836 |
| 19 | 1,754,251 | 1.1876 | 1.1875 | 0.3610 | 0.3611 |
| 23 | 2,241,419 | 1.2777 | 1.2778 | 0.3420 | 0.3420 |
| 29 | 2,954,317 | 1.3181 | 1.3182 | 0.3272 | 0.3272 |
| 31 | 3,269,477 | 1.1067 | 1.1071 | 0.3160 | 0.3159 |
| 37 | 4,030,891 | 1.2329 | 1.2333 | 0.3060 | 0.3057 |
| 41 | 4,588,433 | 1.1383 | 1.1389 | 0.2978 | 0.2974 |
| 43 | 4,929,907 | 1.0744 | 1.0750 | 0.2907 | 0.2902 |
| 47 | 5,516,437 | 1.1190 | 1.1190 | 0.2840 | 0.2834 |
| 53 | 6,355,813 | 1.1522 | 1.1522 | 0.2780 | 0.2774 |

Here *efficiency* is the fraction of the row's raw odd multiples inside $[1,X_p]$ that are new rather than already struck. Each measured column agrees with the one beside it to within $6\times10^{-4}$ in every row, and both agreements are the same fact: by §11.1 the reach is a survivor position, so the survivor density $\prod_{2<q<p}(1-1/q)$ controls everything. The ratio law is then exact,

$$\frac{X_p}{X_{p'}} \;=\; \frac{p}{p'}\cdot\frac{1}{1-1/p'} \;=\; \frac{p}{p'-1},$$

and with Mertens' third theorem — in the odd model $\prod_{2<q<p}(1-1/q) \sim 2e^{-\gamma}/\log p$, twice the usual value because the factor $1-\tfrac12$ is absent — the reach itself is

$$X_p(K) \;\approx\; \frac{2Kp}{\prod_{2<q<p}(1-1/q)} \;\sim\; e^{\gamma}\,K\,p\log p .$$

The measured $X_p/(Kp\log p)$ oscillates between $1.808$ and $1.927$ against $e^\gamma = 1.7811$, the residual being the usual slow Mertens convergence.

### 11.3 Why the other statistics carry nothing extra

With equal production one can ask how a row *moves*: its mean step, its longest pause, how many distinct step sizes it uses. Those look like independent measurements. By §11.1 they are not — each is a statistic of $S_p^-$ multiplied by $2p$:

| row | mean step | $2p/\prod_{2<q<p}(1-1/q)$ | longest step | $p\cdot G_{\text{prev}}$ | distinct step sizes |
|---|---|---|---|---|---|
| 3 | 6.00 | 6.00 | 6 | 6 | 1 |
| 5 | 15.00 | 15.00 | 20 | 20 | 2 |
| 7 | 26.25 | 26.25 | 42 | 42 | 3 |
| 11 | 48.12 | 48.12 | 110 | 110 | 5 |
| 13 | 62.56 | 62.56 | 182 | 182 | 7 |
| 17 | 88.62 | 88.63 | 374 | 374 | 10 |
| 19 | 105.24 | 105.25 | 456 | 456 | 12 |
| 23 | 134.47 | 134.48 | 782 | 782 | 13 |
| 29 | 177.23 | 177.27 | 986 | 986 | 15 |
| 31 | 196.13 | 196.27 | 1054 | 1054 | 15 |

The longest-step column agrees *exactly*, not approximately, with $p$ times the previous stage's largest gap; the mean-step pair is asymptotic rather than exact, which is why the two part company in the last digits as the finite window bites ($196.13$ against $196.27$ at $p=31$). The count of distinct step sizes is the count of distinct gaps in $S_p^-$. Row $7$'s step histogram, $\{14: 6249,\; 28: 6250,\; 42: 4166\}$, is the gap histogram $\{1,2,3\}$ of the survivors after rows $3$ and $5$, scaled by $2p = 14$. There is no separate "complexity of a row's motion" to discover.

The same reduction answers the natural normalisation question. Plotting $a_p(k)/a_p(K)$ against $k/K$ for every row collapses all of them onto a single curve, and the curve is the diagonal: the largest deviation is $0.00003$ at $p=3$, rising only to $0.00081$ at $p=31$ — comparable to (longest step)$/X_p$, which is $3.2\times10^{-4}$ there. Whatever memory a stage carries does not survive the normalisation.

### 11.4 The cycle: copy, delete, merge

Modulo $2P$, where $P = \prod_{2<q<p} q$ is the product of the odd rows already in, the survivors are the residues coprime to $2P$. Bringing in $p$: repeat the old cycle $p$ times, and inside each old family $\{a, a+2P, \dots, a+2(p-1)P\}$ exactly one member is divisible by $p$, since $\gcd(P,p)=1$. Deleting it gives

$$\varphi(2pP) = (p-1)\varphi(2P),$$

and the deleted set is precisely $p \cdot S_p^-$ again — §11.1 seen from inside one period, so the cycle deletes a magnified copy of itself. Locally, if $a<m<b$ are consecutive survivors and $m$ is struck, the two gaps merge:

$$(g_L, g_R) \;\longmapsto\; g_L + g_R .$$

Because every old position is struck in exactly one of its copies, the pair achieving the largest adjacent sum is always struck, and provided no two *adjacent* survivors are struck at the same stage,

$$G^{\text{new}}_{\max} \;=\; \max_i \,(g_i + g_{i+1}).$$

Verified on complete cycles: $4 \to 6 \to 10 \to 14 \to 22 \to 26 \to 34 \to 40 \to 46$, seven consecutive predictions, no exception; the last value, $46$, is the prediction for the next stage and is confirmed independently by the scan of §11.6. The proviso is not decoration and it is decidable: two consecutive survivors both divisible by $p$ lie at least $2p$ apart, so the identity holds while the current largest gap is below $2p$; against the tabulated values of Jacobsthal's function (A048670) the first stage where it *can* fail is the wheel of the odd primes up to $67$, whose largest gap is $152$ against $2p = 142$. Stated without that hypothesis the equality is a conjecture, not a theorem.

### 11.5 The row does not choose its point

Give each live survivor $m$ its *merge capacity* $C(m) = g_L(m)+g_R(m)$, the gap that would appear if $m$ were struck. It is tempting to ask whether an entering row favours points of large capacity — that would make some rows structurally stronger than others. It does not, and the statement is exact rather than statistical:

> Over a complete cycle, the multiset of capacities at the struck points is **identical** to the multiset of capacities over all the old survivors.

The reason is §11.4: before any deletion the replicated cycle is periodic with period $2P$, so a copy $a+2jP$ has the same neighbours as $a$ up to translation, and exactly one copy of each $a$ is struck. Verified as sorted equality for $p = 5,7,11,13,17$. Three consequences:

- the absence of bias is an identity, not an approximate finding — a window measurement gives ratios $0.997$ to $1.006$, which is the shadow of it;
- $\max_m C(m)$ over a row's struck points is the same number for **every** row, so "quiet" and "record-breaking" rows do not exist in the cycle; and
- which row owns a point is statistically independent of the additive geometry around it.

The last one deserves emphasis because it points the opposite way from the reading it invites. It is natural to see $p \mid m$ (multiplicative) together with $(m-a)+(b-m) = b-a$ (additive) inside one event and call it a bridge between the two operations. The multiset identity says the two are *decoupled*: ownership carries no information about environment. That is the same verdict §13 reaches from the arch, arrived at from the sieve side, and measured rather than argued.

### 11.6 Means are cheap, extremes are not

One practical consequence closes the section, and it is the most useful thing in it for anyone tempted to model the sieve dynamically.

Every mean above — density, mean gap, efficiency, mean step — settles inside a window of $10^6$ and equals a Mertens product. Extremes behave completely differently. The largest gap actually visible in $[1,X]$, against the true value for the complete cycle (that last column is a full-cycle computation up to the rows $\le 23$ and the tabulated Jacobsthal value A048670 beyond it, those cycles being out of computational reach):

| rows $\le$ | $10^5$ | $10^6$ | $10^7$ | $10^8$ | $10^9$ | complete cycle |
|---|---|---|---|---|---|---|
| 19 | 34 | 34 | 34 | 34 | 34 | 34 |
| 23 | 34 | 34 | 36 | 40 | 40 | 40 |
| 29 | 34 | 34 | 40 | 42 | 46 | 46 |
| 31 | 36 | 36 | 48 | 48 | 50 | **58** |
| 37 | 40 | 46 | 48 | 50 | 54 | **66** |
| 43 | 46 | 46 | 50 | 54 | 60 | **90** |
| 53 | 58 | 58 | 58 | 60 | 64 | **106** |

Ten thousand times the window resolves two more rows. The reason is that the record gap is *absolutely* rare, not relatively rare: it occurs twice per cycle — a residue and its mirror — for the odd primes to $19$, and twelve times for those to $23$. So the distance needed to see it is about $2P/N$ with $N$ that small, and the measured first occurrences $9{,}461$, $217{,}153$, $60{,}077$, $20{,}332{,}511$ (rows to $13, 17, 19, 23$) sit against $2P/N = 15{,}015$, $255{,}255$, $4{,}849{,}845$, $18{,}591{,}072$. Meanwhile the visible record grows only logarithmically, one to three units per decade of $X$, while the cycle length grows like the primorial: $2P$ is $9.7\times10^6$ for the rows to $19$ and $3.26\times10^{19}$ for those to $53$.

**So a window law for an extreme is a property of the window.** Any dynamical model of the sieve that reports how records behave over a fixed range is measuring its own range; the honest route to the extreme is to measure the gap *distribution*, which converges quickly, and to obtain the maximum from it by extreme-value theory.

> **What is classical here.** All of it. The partition by least prime factor is Eratosthenes; $A_p = pS_p^-$ is that definition rewritten; the copy–delete–merge recursion is the wheel sieve (Pritchard, 1982) and the cycles-of-gaps description of Eratosthenes' sieve; $\varphi(2pP)=(p-1)\varphi(2P)$ is multiplicativity; the largest gap in the coprime cycle is Jacobsthal's function; and the reach law is Mertens. What the section adds is the reduction: one identity accounts for every derived statistic, and the mean/extreme split above says which of them a computation can ever see.

---

## 12. Ledger: what is classical, and what only looked new

The honest summary of both routes, item by item.

| item | it is | verdict |
|---|---|---|
| the arch, $N = (S^2-H^2)/4$ | Fermat's difference of squares | 1643 |
| the laws (a)--(e) in $(S,H)$ | properties of $\mathrm{Re}(z^2/4)$ | one function |
| the diamond law $\partial_a\partial_b(ab)=1$ | a mixed second difference | trivial |
| $N = (L-N)(R-N)$ | $L-N=a$, $R-N=b$ — a tautology | information gain 0 |
| the full unfolding $(S,a)$ | the matrix $\binom{1\ 1}{1\ 0}$, $\det = -1$ | a unimodular shear; gain 0, provably |
| $(S,N)$ plane | the space of monic quadratics $x^2-Sx+N$ | classical, and the right frame |
| arch spectrum, bottom $= 2\sqrt n \cosh u$ | the divisor nearest $\sqrt n$ | Erdős, Ford |
| the repair by combination | Kraitchik $\to$ CFRAC $\to$ QS $\to$ NFS | 1920s–1990s |
| $\rho(n)$ and its error term | a published literature | 1992, 1999, 2006 |
| the skew "fingerprint" $+0.22$ | $53\%$ inherited from $\Delta$ | decomposes exactly |
| the mean-square constant $0.4266$ | $c_{\mathrm{odd}} + c_\Delta/2 + \text{cross}$ | decomposes exactly |
| $x^{1/4}$ | Voronoï amplitude and square-root cancellation | one known fact |
| $\pi/4$ from the circles | Weyl equidistribution | one known fact |
| whole entries $=D(x)-x$ | the summatory divisor function | one known fact |
| $\sum\{n/d\}/n \to 1-\gamma$ | $\gamma$'s own definition, rearranged | one known fact |
| spectrum of the raw $0$–$1$ table | upper triangular, unit diagonal | $\{1\}$, $\det=1$; gain 0 |
| $\det = M(n)$ after one column | Redheffer's matrix | 1977; and it is RH itself, not a route to it |
| $\Phi(x,z)$, the clock $u$, the frontier $u=2$ | Legendre's sieve and Buchstab's $\omega$ | 1937; and $u\le2$ is why Eratosthenes stops at $\sqrt x$ |
| $A_p = p\,S_p^-$ and the $2p$ rescaling | the definition of the least prime factor | a tautology, and it explains every derived statistic |
| copy, delete one per family, merge gaps | the wheel sieve and cycles of gaps | Pritchard 1982 |
| the merge-capacity multiset at struck points | one deletion per family, so it is the old multiset | a decoupling; information about primes 0 |

Twenty-one items; twenty-one decompose. That ratio is the honest yield, and reporting it is the point of this section rather than an apology for it.

There is a cheap test behind that column, and it is worth stating as a method. Before asking whether a result is published — which requires already knowing the object's name — ask whether it is **composed of things one already knows**, and then measure whether the parts account for the whole. The mean-square constant, the skew, and the error exponent above all dissolved under that single question in minutes. Two further questions raised while writing this note (what the conjugate harmonic $M = (b^2-a^2)/2$ misses, and whether the pair $(N,M)$ carries factorisation information beyond $N$) closed the same way: $2M$ is a difference of two squares, so $M$ is never an odd integer; and $M$ cannot be written down without already knowing $a$ and $b$, which is the diamond tautology again.

---

## 13. What is actually open

Four statements, with no claim attached to any of them.

**The error exponent.** For $\tau$, for $r_2$, and for $\rho$, the truth lies between the best proven exponent ($517/1648 = 0.31371\ldots$, Bourgain--Watt, improving Huxley's $131/416$) and the conjectured $1/4+\varepsilon$. Hardy (1915) showed the error is not $O(x^{1/4})$, and Soundararajan (2003) gives the current $\Omega$-refinement, so the $\varepsilon$ is not decoration; the measurement above ($\max|E|/x^{1/4}$ climbing $2.72 \to 4.17$ over four decades) shows the same thing crudely. Nothing between the two exponents is known, and the measurements above sit on $1/4$ without being evidence for it.

**Factoring.** Every general-purpose algorithm descends from $x^2 \equiv y^2 \bmod N$, and the best is subexponential, not polynomial. The bridge from addition to multiplication is crossed; it is not dissolved.

**The size of the determinant.** Whether $\det R_n = O(n^{1/2+\varepsilon})$ is open, and it is open because it is the Riemann hypothesis. This one is listed here for the opposite reason to the other two: not because the picture suggests an attack, but to say plainly that it does not. Mertens' own conjecture, $|M(n)| < \sqrt n$, was disproved by Odlyzko and te Riele (1985), so even the shape of the truth is delicate; and the matrix restatement adds no leverage, since the determinant is computed *from* $\mu$ rather than the reverse. Visibility in an elementary table is not access.

**What a computation can see of the sieve.** §11.6 separates the quantities of the sieve into two classes with a hard line between them: means converge inside $10^6$, while extremes need a complete cycle whose length is a primorial. The largest gap among the odd numbers coprime to every odd prime up to $53$ is $106$; the best a scan to $10^9$ shows is $64$. That is not an open *problem* so much as an open *limit*, and it is stated here because it is the standing reason why a dynamical reading of the sieve cannot be tested where it would matter. The route that does work is indirect — measure the gap distribution, which converges, and recover the maximum from it by extreme-value theory.

Finally, one negative note that the route makes visible, and it is sharper than it first looks. Every object in the arch thread — §§2–9 — is a function of the pair $(a,b)$ through the sum $a+b$ and the difference $b-a$ alone, and the only bridge between addition and multiplication it uses is the difference of two squares.

That bridge is *geometric*: it couples the two operations through a lattice-point count in a quadratic region, and it is blind to which integers are prime. Everything in §§6–9 is a statement about lattice points in a plane region, and the primes appear only as the rows the construction happens to miss. A coupling that could say something about the primes would have to be *arithmetic* — sensitive to primality itself, not merely to the shape of the region. Difference of squares is the bridge the subject has been walking since 1643; whether there is an arithmetic one is the question this picture raises and cannot answer.

The sieve thread of §§10–11 arrives at the same place by the other road, and there the statement is measured rather than argued: §11.5 shows that which row owns a struck point is independent of the additive geometry around it, so the two operations meet in one event and still carry no information about each other.

---

## 14. AI assistance

The verification scripts accompanying this paper, and much of its prose, were written with the assistance of Claude (Anthropic). The research direction, the decisions, and the responsibility for every claim are the author's. See the repository README for a fuller statement.

---

## References

- G. H. Hardy, *On the expression of a number as the sum of two squares*, Quart. J. Math. **46** (1915), 263–283 — contains the $\Omega$-result showing the circle-problem error is not $O(x^{1/4})$.
- K. Soundararajan, *Omega results for the divisor and circle problems*, Int. Math. Res. Not. **2003**, no. 36, 1987–1998.
- G. Voronoï, *Sur une fonction transcendante et ses applications à la sommation de quelques séries*, Ann. Sci. École Norm. Sup. **21** (1904).
- H. Cramér, *Über zwei Sätze des Herrn G. H. Hardy*, Math. Z. **15** (1922).
- E. Preissmann, *Sur la moyenne quadratique du terme de reste du problème du cercle*, C. R. Acad. Sci. Paris **306** (1988), 151–154.
- D. R. Heath-Brown, *The distribution and moments of the error term in the Dirichlet divisor problem*, Acta Arith. **60** (1992), 389–415.
- M. N. Huxley, *Exponential sums and lattice points III*, Proc. London Math. Soc. **87** (2003), 591–609.
- J. Bourgain and N. Watt, *Mean square of zeta function, circle problem and divisor problem revisited*, arXiv:1709.04340 (the exponent $517/1648$).
- M. Kühleitner, *An omega theorem on differences of two squares*, Acta Math. Univ. Comenianae **61** (1992), 117–123; **II**, Acta Math. Univ. Comen. New Ser. **68** (1999), 27–35.
- M. Kühleitner, W. G. Nowak, J. Schoissengeier, T. D. Wooley, *On sums of two cubes: an $\Omega_+$-estimate for the error term*, Acta Arith. **80** (1998), 179–195.
- M. Kühleitner and W. G. Nowak, *On differences of two squares*, Cent. Eur. J. Math. **4** (2006), 110–122 — a sharp lower bound and two mean-square results for the error term of $\sum_{n\le x}\rho(n)$, and the connection to $d(n)$.
- X. Cao, W. Zhai, *On the mean square of the error term for the two-dimensional divisor problem*, arXiv:0806.3902.
- P. de Fermat, letter to Mersenne, 1643 (the difference-of-squares method).
- M. Kraitchik, *Théorie des Nombres*, Gauthier-Villars, Paris, 1922–1926 (the idea of combining congruences rather than demanding an exact square).
- M. A. Morrison, J. Brillhart, *A method of factoring and the factorization of $F_7$*, Math. Comp. **29** (1975), 183–205.
- R. S. Lehman, *Factoring large integers*, Math. Comp. **28** (1974), 637–646.
- C. Pomerance, *The quadratic sieve factoring algorithm*, EUROCRYPT '84, 169–182.
- P. L. Montgomery, as described in R. D. Silverman, *The multiple polynomial quadratic sieve*, Math. Comp. **48** (1987), 329–339.
- A. K. Lenstra, H. W. Lenstra Jr. (eds.), *The Development of the Number Field Sieve*, Springer LNM 1554, 1993.
- K. Dickman, *On the frequency of numbers containing prime factors of a certain relative magnitude*, Ark. Mat. Astron. Fys. **22A** (1930).
- A. A. Buchstab, *Asymptotic estimates of a general number-theoretic function*, Mat. Sb. **44** (1937), 1239–1246 — the function $\omega$ of §6.
- H. Halberstam and H.-E. Richert, *Sieve Methods*, Academic Press, 1974 — for $\Phi(x,z)$ and $\omega$.
- G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, 3rd ed., CUP, 2015 — Part III for $\omega$ and its duality with $\rho$.
- F. Mertens, *Ein Beitrag zur analytischen Zahlentheorie*, J. Reine Angew. Math. **78** (1874), 46–62 — the product used throughout §7.
- P. Pritchard, *Explaining the wheel sieve*, Acta Informatica **17** (1982), 477–485 — the copy–delete recursion of §11.4.
- C. Hooley, *On the difference of consecutive numbers prime to n*, Acta Arith. **8** (1962/63), 343–347 — gaps in the coprime cycle.
- H. L. Montgomery and R. C. Vaughan, *On the distribution of reduced residues*, Ann. of Math. **123** (1986), 311–333 — the moments of those gaps.
- OEIS A048670 — the largest gap in the cycle coprime to a primorial, the "complete cycle" column of §11.6.
- P. Erdős, *Some remarks on number theory*, Riveon Lematematika **9** (1955) (the multiplication table problem); K. Ford, *The distribution of integers with a divisor in a given interval*, Ann. of Math. **168** (2008), 367–433.
- R. Redheffer, *Eine explizit lösbare Optimierungsaufgabe*, Internat. Ser. Numer. Math. **36** (1977), 213–216 (the determinant of the modified divisibility matrix is the Mertens function).
- W. W. Barrett, R. W. Forcade, A. D. Pollington, *On the spectral radius of a $(0,1)$ matrix related to Mertens' function*, Linear Algebra Appl. **107** (1988), 151–159.
- D. W. Robinson, W. W. Barrett, *The Jordan 1-structure of a matrix of Redheffer*, Linear Algebra Appl. **112** (1989), 57–73.
- W. W. Barrett, T. J. Jarvis, *Spectral properties of a matrix of Redheffer*, Linear Algebra Appl. **162/164** (1992), 673–683.
- R. C. Vaughan, *On the eigenvalues of Redheffer's matrix* I, in *Number Theory with an Emphasis on the Markoff Spectrum*, Lecture Notes in Pure and Appl. Math. **147**, Dekker (1993), 283–296; **II**, J. Austral. Math. Soc. Ser. A **60** (1996), 260–273.
- D. A. Cardon, *Matrices related to Dirichlet series*, arXiv:0809.0076 (the general construction: the determinant of the modified matrix is a weighted partial sum of the coefficients of $L(s)^{-1}$).
- A. M. Odlyzko, H. J. J. te Riele, *Disproof of the Mertens conjecture*, J. Reine Angew. Math. **357** (1985), 138–160.

---

## Appendix. Reproducing the computations

Everything in this note is checkable with a sieve and a few dozen lines.

- **Section 1.** A $\tau$-sieve to $100$ and a harmonic sum; both counts are one line. For the determinant, build the $0$–$1$ matrix with $K(d,m)=1$ iff $d\mid m$, overwrite its first column with $1$s, and take the determinant in *exact integer* arithmetic (Bareiss or a rational LU) — floating point is worthless here, since the answer is a small integer obtained from a badly conditioned matrix. It agrees with a $\mu$-sieve for every $n$ tested up to $80$. The eigenvalue count is a float eigensolver plus a tolerance, and the tolerance matters: Vaughan proved there are non-trivial eigenvalues arbitrarily close to $1$, so a fixed cutoff will undercount for large $n$.
- **Sections 2–4.** Direct enumeration; the arch identities are one-line checks over $S < 400$.
- **Section 5.** The $2041$ table is by hand. The striping count needs only $\pi(B)$ and $\sum_p 2M/p$, with the roots taken relative to the start index $x_0 = \lceil\sqrt N\rceil$ — modulo $7$ they are $i \equiv 0,1$ because $x_0 \equiv 3$, not because $N$ is a square mod $7$ at $x \equiv 0,1$. Dickman's $\rho$ comes from solving $u\rho'(u) = -\rho(u-1)$ on a grid; a coarse grid biases it high by a percent or two, which is enough to move the third decimal. The three residue families are: one interval $x_0\pm10^5$; MPQS with $a$ prime near $\sqrt{2N}/M$, $b^2\equiv N \bmod a$, $M=10^3$; and the continued-fraction residues of $\sqrt N$.
- **Section 6.** $\rho(n)$ from a $\tau$-sieve; the summatory function to $10^{12}$ by the hyperbola method in $O(\sqrt X)$.
- **Section 7.** The wobble sum $W(n)$ is a direct float sum, no $\gamma$ in the input. Every decade in the boundary table must be run to its end; see the note under the table. The four boundaries need the correct main terms: for a line of slope $\alpha$ the boundary correction is $-cn$ with $c$ the mean of $\{\alpha a\}$, namely $\tfrac12$ for irrational $\alpha$ and $\tfrac14$ for $\alpha = 3/2$. Getting this wrong produces a spurious exponent near $0.77$.
- **Section 8.** $\sum \tau(n)e^{-n/X}$ directly, summed to $40X$ so that the truncated tail ($\sim e^{-40}$) is far below the $-1/(144X)$ term being measured. Stopping at $10X$ or $20X$ leaves a tail that swamps it.
- **Section 9.** Sample the hyperbola densely on $a \in [\sqrt x/4,\, 4\sqrt x]$, accumulate arc-length elements, and test each midpoint against the nearest cell centre.
- **Section 10.** A plain sieve to $10^8$ gives $\Phi(x,z)$ directly; $\omega$ comes from integrating $(u\omega)' = \omega(u-1)$ on a grid from the branch $1/u$ on $[1,2]$ — the same care as Dickman's $\rho$ in §5, and a coarse grid moves the third decimal the same way. Normalise by the exact product $\prod_{p\le z}(1-1/p)$, not by $1/\log z$: the difference is the factor $B$ tabulated in the text, and using $1/\log z$ hides it inside an apparent error in $\omega$.
- **Section 11.** The odd primes act here as *rows* of the sheet, in the sense fixed in §1. Everything follows a least-prime-factor sieve. Compute $X_p(K)$ as the position of the $K$-th element of $A_p$, not by inverting an asymptotic. For §11.5 the replication matters: mark the deletions *after* laying $p$ copies of the old cycle end to end, since a copy has the neighbours of its original only inside the replicated period — marking $p \mid v$ inside a single old cycle gives a different and meaningless set. For §11.6 the window checkpoints must be segment boundaries; a running maximum reported at a checkpoint inside a segment reports the whole segment's value and makes three decades look identical.

Four scripts regenerate every number printed above: `verify_division_table.py` for §§1–9, `verify_knowledge_clock.py` for §10, `verify_row_inheritance.py` for §11, and `verify_redheffer.py`, which re-does §1's determinant claims on their own with exact rational arithmetic. Each exits nonzero if any number moves; `--fast` lowers the sieve limits and skips the rows that need more.

*This note claims no new theorem. It offers a route, and an honest map of where each step of it already lives.*
