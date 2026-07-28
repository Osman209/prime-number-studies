"""
conventions.py -- single source of truth for the fixed-support (sine-window) family.

Every other module in this harness imports its conventions from here. Nothing below
is a free choice made at call time; if a convention changes it changes here and the
regression tests must be re-baselined deliberately.

--------------------------------------------------------------------------------
1. PARAMETERS

    L      support length. The test-function support is the symmetric interval
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
"""

# ---- fixed conventions -------------------------------------------------------
R_MAX = 900.0          # archimedean integral truncation, |r| <= R_MAX
NR = 30001             # quadrature nodes on [-R_MAX, R_MAX]
TAU_DEFAULT = 1e-6     # absolute near-null tolerance; always reported
GAMMA_MAX_DEFAULT = 1062.9   # highest zeta ordinate used by the validator
M_LADDER = (100, 200, 300, 400, 500)

__all__ = ["R_MAX", "NR", "TAU_DEFAULT", "GAMMA_MAX_DEFAULT", "M_LADDER"]
