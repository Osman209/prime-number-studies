"""
Does the reality-of-zeros property extend to ODD eigenvectors at higher N?

The paper (arXiv:2511.23257) proves it for the EVEN ground state (Thm 6.1), shows in
Appendix B.1 that a NON-EXTREMAL eigenvector can have complex zeros, and settles the
ODD kernel vector by hand only at N = 2 (Appendix B.3).  This script pushes both
questions to N = 3..8 by random sampling inside the paper's own matrix class.

Matrix class, eq. (11): indices i in {-N..N},
    q_ii = a_i  with a_{-i} = a_i,
    q_ij = (b_i - b_j)/(i - j)  for i != j,  with b_{-i} = -b_i.
Any such Q minus lambda_min * I is positive semidefinite of the same form, with the
ground eigenvector in its kernel -- which is exactly the setting of Appendix B.3.

The polynomial of eq. (21) is P(s) = sum_k xi_k prod_{j != k} (j - s); dividing by
prod_j (j - s) shows its roots are the roots of  R(s) = sum_k xi_k/(k - s) = 0.

TWO QUESTIONS, and they are different:
  (A) odd vector IN THE KERNEL   -- Appendix B.3's case, generalised
  (B) lowest ODD eigenvector while the ground state is EVEN -- OUR case, which is
      not covered by anything in the paper and is where B.1 warns reality can fail.
"""
import numpy as np
import sys
import mpmath as mp

rng = np.random.default_rng(7)
gates = {}


def build_Q(N, a_half, b_half):
    """a_half = (a_0..a_N), b_half = (b_1..b_N); b_0 = 0 forced by oddness."""
    idx = np.arange(-N, N + 1)
    a = np.array([a_half[abs(i)] for i in idx])
    b = np.array([0.0 if i == 0 else np.sign(i) * b_half[abs(i) - 1] for i in idx])
    Q = np.empty((2 * N + 1, 2 * N + 1))
    for p, i in enumerate(idx):
        for q, j in enumerate(idx):
            Q[p, q] = a[p] if i == j else (b[p] - b[q]) / (i - j)
    return Q


def roots_real_mpmath(xi, N, dps=60):
    """Recompute the roots of (21) in extended precision.

    Degree-2N polynomials are badly conditioned in float64, so a float64 'complex
    root' is not evidence of one. Every float64 failure is re-tested here before it
    is counted, which is what makes the reported counts exact rather than noisy.
    """
    old_dps = mp.mp.dps
    mp.mp.dps = dps
    try:
        idx = list(range(-N, N + 1))
        P = [mp.mpf(0)] * (2 * N + 1)
        for p, k in enumerate(idx):
            if abs(xi[p]) < 1e-14:
                continue
            term = [mp.mpf(1)]
            for j in idx:
                if j != k:
                    new_t = [mp.mpf(0)] * (len(term) + 1)
                    for t, co in enumerate(term):
                        new_t[t] -= co
                        new_t[t + 1] += mp.mpf(j) * co
                    term = new_t
            for t, co in enumerate(term):
                P[t] += mp.mpf(float(xi[p])) * co
        while len(P) > 1 and P[0] == 0:
            P = P[1:]
        if len(P) < 2:
            return True, 0.0
        r = mp.polyroots(P, maxsteps=200, extraprec=200)
        worst = max(float(abs(mp.im(z))) for z in r)
        return worst < 1e-20, worst
    finally:
        mp.mp.dps = old_dps


def all_roots_real(xi, N, tol=1e-8):
    """Roots of P(s) = sum_k xi_k prod_{j!=k}(j-s), built by polynomial arithmetic."""
    idx = np.arange(-N, N + 1)
    P = np.zeros(2 * N + 1)
    for p, k in enumerate(idx):
        if abs(xi[p]) < 1e-14:
            continue
        term = np.array([1.0])
        for j in idx:
            if j != k:
                term = np.convolve(term, np.array([-1.0, float(j)]))  # (j - s)
        P += xi[p] * term
    P = np.trim_zeros(P, "f")
    if len(P) < 2:
        return True, 0.0
    r = np.roots(P)
    worst = float(np.max(np.abs(r.imag) / (1 + np.abs(r.real))))
    return worst < tol, worst


def parity_split(N):
    D = 2 * N + 1
    r = 1 / np.sqrt(2)
    Ve = np.zeros((D, N + 1)); Ve[N, 0] = 1
    Vo = np.zeros((D, N))
    for k in range(1, N + 1):
        Ve[N + k, k] = Ve[N - k, k] = r
        Vo[N + k, k - 1] = r
        Vo[N - k, k - 1] = -r
    return Ve, Vo


print(f"{'N':>3} {'trials':>7} {'ground odd':>11} {'(A) odd kernel real':>21} "
      f"{'(B) lowest-odd real':>21} {'worst |Im| (B)':>15}")
for N in (3, 4, 5, 6, 8):
    Ve, Vo = parity_split(N)
    trials = 400
    n_ground_odd = nA = nA_real = nB = nB_real = 0
    worstB = 0.0
    for _ in range(trials):
        a_half = rng.normal(size=N + 1)
        b_half = rng.normal(size=N)
        Q = build_Q(N, a_half, b_half)
        w, V = np.linalg.eigh(Q)
        g = V[:, 0]                                   # ground state of the whole matrix
        is_odd = np.abs(g + g[::-1]).max() < 1e-8     # gamma g = -g
        # (A) if the ground state is odd, it is an odd kernel vector of Q - lambda_min I
        if is_odd:
            n_ground_odd += 1
            nA += 1
            ok, _ = all_roots_real(g, N)
            if not ok:                       # float64 says complex — verify in mpmath
                ok, _ = roots_real_mpmath(g, N)
            nA_real += ok
        else:
            # (B) ground state even: take the lowest eigenvector of the ODD block
            Qo = Vo.T @ Q @ Vo
            wo, Uo = np.linalg.eigh(Qo)
            xi = Vo @ Uo[:, 0]
            nB += 1
            ok, worst = all_roots_real(xi, N)
            if not ok:                       # same recheck on the other side
                ok, worst = roots_real_mpmath(xi, N)
            nB_real += ok
            worstB = max(worstB, worst if not ok else 0.0)
    gates[N] = (nA_real, nA, nB_real, nB)
    print(f"{N:>3} {trials:>7} {n_ground_odd:>11} "
          f"{(f'{nA_real}/{nA}' if nA else '-'):>21} "
          f"{f'{nB_real}/{nB}':>21} {worstB:>15.2e}")

print("\n(A) = the paper's Appendix B.3 case, generalised: odd vector in the kernel.")
print("(B) = OUR case: ground state even, lowest eigenvector of the odd sector.")

# --- gates: the two claims this script exists to support ---
failures = []
for N, (a_real, a_tot, b_real, b_tot) in gates.items():
    if a_tot and a_real < a_tot:
        failures.append(f"N={N}: experiment A is {a_real}/{a_tot}, expected all real "
                        f"(float64 root-finding noise? re-check the failures in mpmath)")
    if b_tot and not (0.30 <= b_real / b_tot <= 0.70):
        failures.append(f"N={N}: experiment B is {b_real}/{b_tot} = {b_real/b_tot:.2f}, "
                        f"expected roughly half; the reported contrast rests on this")
if failures:
    print("\nFAILED:")
    for f in failures:
        print("  " + f)
    sys.exit(1)
print("\nboth gates pass: A all-real, B near one half")
