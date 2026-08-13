"""
conventions.py -- single source of truth for the fixed-support (sine-window) family.

Every other module in this harness imports its conventions from here. Nothing below
is a free choice made at call time; if a convention changes it changes here and the
regression tests must be re-baselined deliberately.

--------------------------------------------------------------------------------
1. PARAMETERS

    L      support length, default L_DEFAULT below. The test-function support is
           the symmetric interval
           [-L/2, +L/2]. The autocorrelation g = f * f~ then has support [-L, L],
           so a prime power q contributes iff log q < L.

    p_c    prime cutoff. Only prime powers q = p^k with q <= p_c enter the prime
           sum. INDEPENDENT of L in this family.

           SATURATED DIAGONAL:  p_c = floor(exp(L)).
           This is the CvS parametrisation, in which the last prime's shift
           log(p_c) reaches the autocorrelation-support edge L. Off the diagonal
           the two parameters move independently, which is what makes this a
           two-parameter family and NOT the CvS one-parameter path.

    m      Galerkin dimension (number of sine modes retained).

    tau    near-null tolerance, an ABSOLUTE cut on the eigenvalue. Must be
           reported with every count that uses it; see section 4.

--------------------------------------------------------------------------------
2. BASIS

    phi_j(x) = sqrt(2/L) * sin(a_j * (x + L/2)),   a_j = j*pi/L,   j = 1..m

    orthonormal on [-L/2, L/2]. This is a DIRICHLET SINE space: every basis
    function vanishes at both endpoints. It is NOT an invertible change of basis
    of the finite periodic Fourier space used by CvS -- the endpoint maps and the
    frequency sets differ -- so the two families are genuinely distinct and no
    result may be transported between them without proof.

    Fourier transform, with the convention fhat(r) = Int f(x) e^{-i r x} dx:

        phihat_j(r) = e^{+i r L/2} * sqrt(2/L) * a_j * (1 - (-1)^j e^{-i r L})
                                                 / (a_j^2 - r^2)

    The apparent pole at r = a_j is removable.

--------------------------------------------------------------------------------
3. THE FORM

    W = ARCH + POLE - PRIME,  each block returned separately by the builder.

    ARCH_jk = (1/2pi) Int phihat_j(r) conj(phihat_k(r)) Psi(r) dr
              Psi(r) = Re psi(1/4 + i r/2) - log(pi)
              truncated at |r| <= R_MAX with NR quadrature nodes; the truncation
              tail is bounded and reported, not assumed negligible.

    POLE    = m_plus (m_minus)^T + m_minus (m_plus)^T           (rank two)
              m_pm,j = Int phi_j(x) e^{+-x/2} dx

    PRIME   = 2 * sum_{q = p^k <= p_c, log q < L}
                    Lambda(q) q^{-1/2} * Rsym(log q)
              Rsym(t)_jk = symmetrised Int_0^{L-t} phi_j(y) phi_k(y+t) dy

    SIGN CONVENTION is fixed by the validator in section 5: with the above,
    W(f) must equal the sum over ALL nontrivial zeros of |fhat(gamma)|^2, which
    for real f equals 2 * sum_{gamma>0} |fhat(gamma)|^2.

--------------------------------------------------------------------------------
4. WHAT IS REPORTED, AND WHY NOT A BARE COUNT

    A count of the form #{lambda < tau} is a near-null count ONLY where the form
    is positive semidefinite. Off the saturated diagonal this family is strongly
    indefinite, so such a count silently sums two different quantities. Every
    report therefore carries the full inertia triple:

        lambda_min
        n_minus  = #{lambda < -tau}
        n_zero   = #{-tau <= lambda <= +tau}
        n_plus   = #{lambda > +tau}
        tau      (stated, absolute)

    plus a residual / backward-error figure from section 5 across the m-ladder.

    NOTE ON PRIOR ART. Bombieri's negative-eigenvalue count theorem is stated for
    his zero-indexed matrix H(Gamma; t), not for this sine-window family. Inertia
    is a natural invariant here, but that theorem does not transfer automatically
    and is not invoked.

--------------------------------------------------------------------------------
5. VALIDATION (non-circular)

    The zero side is computed from mpmath.zetazero, which enters no part of the
    construction. For an eigenvector c of W with eigenvalue lambda:

        lambda  ?=  2 * sum_{gamma > 0, gamma <= GAMMA_MAX} |fhat_c(gamma)|^2
                    + (bounded tail)

    The tail beyond GAMMA_MAX is bounded using |phihat_j(r)| = O(1/r^2), hence
    |fhat|^2 = O(1/r^4), against the zero density dN = log(gamma/2pi)/2pi dgamma.

    This identity holds ONLY on the saturated diagonal, where the prime sum is
    complete up to the support edge. Off the diagonal prime terms are missing by
    construction and the identity is not expected to hold; the validator is
    therefore run on the diagonal and its passing certifies the ASSEMBLY, which
    is then reused unchanged off the diagonal.

--------------------------------------------------------------------------------
6. NUMERICAL REGIME

    Double precision throughout the ladder. Edge-limit quantities (delta -> 0
    behaviour of the overlap) are cancellation-limited in float64 and are
    computed in extended precision instead; see edge_precision.py.

--------------------------------------------------------------------------------
7. THE SECOND FAMILY: PERIODIC FOURIER (CvS), AND HOW IT MAPS ONTO THIS ONE

    Sections 1-6 fix the fixed-support (Dirichlet sine) family. A second family is
    used in the CvS path and is available through builder_periodic.py, which wraps
    A. Groskin's `connes-cvs` implementation of CvS Proposition 4.1. This section
    states the correspondence so that both constructors run under one set of
    conventions rather than one forking from the other. Nothing here is a free
    choice either; it is read off that constructor's own source.

    PARAMETERS. build_galerkin_matrix(c, N, T, dps) sets L = log c internally and
    sums prime powers n <= c. It is therefore a ONE-parameter family: the support
    and the prime cutoff are locked together, which is this harness's SATURATED
    DIAGONAL and nothing else. The map is

        L = log c,      p_c = c,      m = 2N,

    and off the diagonal there is no periodic constructor to call. Requests off
    the diagonal are refused rather than answered with a matrix built at the
    wrong support.

    BASIS. Take the correspondence from what the periodic side actually
    reconstructs, not from a docstring: connes_cvs.extract_zeros integrates
    against exp(2 pi i k x/L) on [0, L], so the periodic basis is the FULL Fourier
    basis of period L, frequencies 2 pi k/L. The sine family has j pi/L. Hence

        periodic sine half   (e_k - e_{-k}) ~ sin(2 pi k x/L)  =  phi_{2k}
        periodic cosine half {1, cos(2 pi k x/L)}              outside this span
        sine odd-indexed     phi_1, phi_3, ...                 outside THEIR span

    The periodic odd sector is exactly this family's EVEN-INDEXED sine modes --
    which is why m = 2N is the right pairing, and it is verified by the spectra
    below rather than argued from the basis. The two "even" sectors are not two
    versions of one space: one is spanned by cosines, the other by the
    half-integer-frequency sines, so they share no vector at FINITE size. They are
    not disjoint in the limit: phi_j for j odd carries the same parity and is
    complete in the same subspace, but only asymptotically -- the L2 residual of
    the constant and the cosines in those modes falls like 1/m (4.0e-2, 8.1e-3,
    2.0e-3, 1.0e-3 at j <= 10, 50, 200, 400). So the finite blocks differ and the
    difference should converge in m, not in T. That is an inference from
    completeness; the m-dependence of the block discrepancy has NOT been swept.

    BOUNDARY PRIME. q = c is included by the periodic side (n <= c) and excluded
    here (log q < L, strict). Its overlap integrates over an interval of length
    L - log c = 0, so it contributes nothing and the two prime sums agree.

    REPORTING IS SECTOR BY SECTOR, NOT MATRIX TO MATRIX. Both forms commute with
    parity (measured: the sine W is parity block-diagonal to 7e-15, so the split
    is exact and not imposed), and the two sectors then behave differently. A
    single whole-matrix comparison averages the two behaviours into one number
    that means nothing.

      ODD sector -- the periodic combinations (e_k - e_{-k}) ~ sin(2 pi k x/L)
      vanish at the centre AND at both endpoints, so they ARE the Dirichlet sine
      modes and the two constructors compute the SAME operator. The measured
      difference falls with the archimedean truncation and not otherwise:
      1.7e-4, 2.4e-5, 3.2e-6 at T = 100, 200, 400 (c = 13, N = 16, dps = 80).

      EVEN sector -- the periodic side carries the constant and the cosines that
      do not vanish at the endpoints. The Dirichlet space cannot represent them
      and the operators genuinely differ: 5.0e-2, 4.7e-2, 4.9e-2 over the same
      sweep, flat in T. The sector dimensions differ by one, the extra mode being
      k = 0.

      NEAR-NULL COUNTS agree in BOTH sectors (7/7 and 7/7 at N = 20, 5/5 and 5/5
      at N = 12), so the blind subspace is a shared invariant even where the
      non-null spectra are not.

      AND THE SHARED SECTOR CARRIES THE ZEROS. At c = 13, N = 100, T = 400,
      dps = 80 the near-null vector of the ODD sector -- the one this family
      spans -- reconstructs the first three zeta ordinates through
      connes_cvs.extract_zeros to 8.70e-53, 6.67e-50, 4.28e-48, against the even
      sector's 1.45e-55, 2.69e-52, 2.49e-50 (lambda_min 2.011e-55 against
      2.077e-59). So the endpoint condition costs roughly three decimal places
      and a factor 1e4 in depth, and it is NOT what makes the reconstruction
      work: the fixed-support family reaches the same zeros through a subspace
      the CvS ground state does not use.
      A CANDIDATE MECHANISM FOR THE GAP, not established here.  Andrade
      (10.5281/zenodo.20710075) shows the prime block is a Loewner matrix whose
      per-prime signs are fixed by the PARITY of the test vector: an even vector
      has a real positive Fourier symbol, so the prime terms all carry one sign
      and add without cancellation, while an odd vector has a purely imaginary,
      sign-changing symbol, so its terms compete and partly cancel (his
      transition sits near q ~ sqrt c).  The imaginary/sign-changing half of that
      follows from parity alone and so applies to every odd vector, not only to
      his sinh.  Whether it accounts QUANTITATIVELY for the factor above is not
      tested: his law is about the sign of each per-prime term, the numbers here
      are reconstruction errors, and nothing yet connects the two.

    PRECISION, two traps, both measured rather than assumed:
      * dps must be raised together with T. At dps = 40, T = 800 degrades the
        comparison to 5.5e-1 -- worse than T = 200 gives -- while at dps = 80 the
        T-trend is monotone.
      * connes_cvs.extract_zeros must receive L as an mpmath number. A Python
        float caps the reconstruction at ~1e-17 instead of ~1e-55, silently: at
        c = 13, N = 100, T = 400, dps = 80 we measure |gamma_1 error| = 1.41e-17
        with math.log(13) against 1.4549524e-55 with mp.log(13), while
        lambda_min^even = 2.0769626582e-59 either way. See upstream issue #3.
"""

# ---- fixed conventions -------------------------------------------------------
L_DEFAULT = 4.5        # support length used by every entry point unless overridden
R_MAX = 900.0          # archimedean integral truncation, |r| <= R_MAX
NR = 30001             # quadrature nodes on [-R_MAX, R_MAX]
TAU_DEFAULT = 1e-6     # absolute near-null tolerance; always reported
GAMMA_MAX_DEFAULT = 1062.9   # highest zeta ordinate used by the validator
M_LADDER = (100, 200, 300, 400, 500)

# ---- periodic (CvS) family; see section 7 --------------------------------------
T_PERIODIC = 400       # archimedean truncation passed to build_galerkin_matrix
DPS_PERIODIC = 80      # working precision; raise it together with T_PERIODIC

__all__ = ["L_DEFAULT", "R_MAX", "NR", "TAU_DEFAULT", "GAMMA_MAX_DEFAULT", "M_LADDER",
           "T_PERIODIC", "DPS_PERIODIC"]
