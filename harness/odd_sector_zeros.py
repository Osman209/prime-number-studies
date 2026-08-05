import mpmath as mp
import sys, time
from connes_cvs.operator import build_galerkin_matrix, extract_zeros

C, N, T, DPS = 13, 100, 400, 80
mp.mp.dps = DPS
L = mp.log(C)

t0 = time.time()
Q = build_galerkin_matrix(c=C, N=N, T=T, dps=DPS)
print(f"build {time.time()-t0:.0f} s, dim {Q.rows}", flush=True)

D = 2 * N + 1
r = 1 / mp.sqrt(2)
Ve = mp.zeros(D, N + 1); Ve[N, 0] = 1
Vo = mp.zeros(D, N)
for k in range(1, N + 1):
    Ve[N + k, k] = r; Ve[N - k, k] = r
    Vo[N + k, k - 1] = r; Vo[N - k, k - 1] = -r

first_errors = {}
for name, V in (("even", Ve), ("odd", Vo)):
    t0 = time.time()
    B = V.T * Q * V
    E, U = mp.eigsy(B)
    i = min(range(B.rows), key=lambda j: abs(E[j]))
    lam = E[i]
    v = mp.zeros(D, 1)
    col = mp.matrix([U[j, i] for j in range(B.rows)])
    full = V * col
    nrm = mp.sqrt(sum(full[j] ** 2 for j in range(D)))
    for j in range(D):
        v[j, 0] = full[j] / nrm
    print(f"\n{name} sector: dim {B.rows}, lambda_min = {mp.nstr(lam, 10)}  ({time.time()-t0:.0f} s)", flush=True)
    z = extract_zeros(v, L, n_zeros=3, dps=DPS)
    first_errors[name] = z[0]["error"]
    for e in z:
        d, err = e["gamma_detected"], e["error"]
        print(f"   zero {e['k']}: detected = {mp.nstr(d, 22) if d is not None else 'None':<24}"
              f" error = {mp.nstr(err, 6) if err is not None else 'None'}", flush=True)
# --- gate: both sectors must reconstruct, and the even one must be the deeper ---
bad = []
for sector, first in first_errors.items():
    if first is None or first > mp.mpf("1e-40"):
        bad.append(f"{sector} sector: gamma_1 error is {first}, expected below 1e-40")
if not bad and first_errors["even"] >= first_errors["odd"]:
    bad.append("the even sector is not the deeper one — check the parity projection")
if bad:
    print("\nFAILED:")
    for b in bad:
        print("  " + b)
    sys.exit(1)
print("\ngate passes: both sectors reconstruct, even deeper than odd", flush=True)
