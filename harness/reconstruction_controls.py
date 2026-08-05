import mpmath as mp, time, random, sys
from connes_cvs.operator import build_galerkin_matrix, extract_zeros

C, N, T, DPS = 13, 60, 400, 80
mp.mp.dps = DPS
L = mp.log(C)
t0 = time.time()
Q = build_galerkin_matrix(c=C, N=N, T=T, dps=DPS)
print(f"build {time.time()-t0:.0f} s", flush=True)

D = 2 * N + 1
r = 1 / mp.sqrt(2)
Ve = mp.zeros(D, N + 1); Ve[N, 0] = 1
Vo = mp.zeros(D, N)
for k in range(1, N + 1):
    Ve[N + k, k] = r; Ve[N - k, k] = r
    Vo[N + k, k - 1] = r; Vo[N - k, k - 1] = -r

def norm_lift(V, col, D):
    full = V * col
    nrm = mp.sqrt(sum(full[j] ** 2 for j in range(D)))
    v = mp.zeros(D, 1)
    for j in range(D):
        v[j, 0] = full[j] / nrm
    return v

results = {}

def report(tag, v):
    z = extract_zeros(v, L, n_zeros=1, dps=DPS)[0]
    d, e = z["gamma_detected"], z["error"]
    results[tag] = e
    print(f"  {tag:38s} detected = {mp.nstr(d,18) if d is not None else 'None':<22}"
          f" error = {mp.nstr(e,6) if e is not None else 'None'}", flush=True)

for name, V in (("even", Ve), ("odd", Vo)):
    B = V.T * Q * V
    E, U = mp.eigsy(B)
    order = sorted(range(B.rows), key=lambda j: abs(E[j]))
    i_min, i_max = order[0], order[-1]
    print(f"\n{name}: lambda_min = {mp.nstr(abs(E[i_min]),6)}, "
          f"lambda_max = {mp.nstr(abs(E[i_max]),6)}", flush=True)
    report(f"{name} NEAR-NULL vector (the claim)",
           norm_lift(V, mp.matrix([U[j, i_min] for j in range(B.rows)]), D))
    report(f"{name} LARGEST-eigenvalue vector (control)",
           norm_lift(V, mp.matrix([U[j, i_max] for j in range(B.rows)]), D))
    report(f"{name} 2nd smallest (control)",
           norm_lift(V, mp.matrix([U[j, order[1]] for j in range(B.rows)]), D))
    random.seed(0)
    report(f"{name} RANDOM vector (control)",
           norm_lift(V, mp.matrix([mp.mpf(random.gauss(0, 1)) for _ in range(B.rows)]), D))
# --- gate: near-null reconstructs; the high and random controls do not.
# Deliberately NOT gated on the second-smallest vector: it also reconstructs
# closely (1e-49, 1e-47), five to six orders worse than the ground vector, so
# "only the near-null vector reconstructs" would be false. The claim the gate
# enforces is the separation between the low end and the controls. ---
bad = []
for sector in ("even", "odd"):
    near = results.get(f"{sector} NEAR-NULL vector (the claim)")
    if near is None or near > mp.mpf("1e-40"):
        bad.append(f"{sector}: near-null error is {near}, expected below 1e-40")
    for tag in (f"{sector} LARGEST-eigenvalue vector (control)",
                f"{sector} RANDOM vector (control)"):
        ctrl = results.get(tag)
        if ctrl is not None and ctrl < mp.mpf("1e-3"):
            bad.append(f"{tag}: error {ctrl} is too small — the control is not separating")
if bad:
    print("\nFAILED:")
    for b in bad:
        print("  " + b)
    sys.exit(1)
print("\ngate passes: the low end reconstructs, the high and random controls do not", flush=True)
