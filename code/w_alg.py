import numpy as np, math
from sympy import factorint, mobius, primerange
print("="*74); print("W8  SECTION 7-8 ALGEBRA: Z, Z^{-1}, and the operator form of P_X"); print("="*74)
M=200
Z=np.zeros((M+1,M+1))
for d in range(1,M+1):
    for n in range(d,M+1,d): Z[d,n]=1
Zi=np.zeros((M+1,M+1))
for d in range(1,M+1):
    for n in range(d,M+1,d): Zi[d,n]=float(mobius(n//d))
I=np.eye(M+1); I[0,0]=0
E=Z[1:,1:]@Zi[1:,1:]
print(f"   || Z Z^inv - I ||_max  (Z^inv_{{d,n}} = mu(n/d) for d|n) = {np.abs(E-np.eye(M)).max():.2e}")
def Lam(n):
    f=factorint(n); return math.log(list(f)[0]) if len(f)==1 else 0.0
lam=np.array([0.0]+[Lam(n) for n in range(1,M+1)])
ell=np.array([0.0]+[math.log(n) for n in range(1,M+1)])
print(f"   || Z^T lambda - ell ||_max   (log n = sum_{{d|n}} Lambda(d)) = {np.abs(Z[1:,1:].T@lam[1:]-ell[1:]).max():.2e}")
print(f"   || (Z^inv)^T ell - lambda ||_max                            = {np.abs(Zi[1:,1:].T@ell[1:]-lam[1:]).max():.2e}")
L=math.log(M); w=np.array([0.0]+[max(0.0,1-math.log(n)/L) for n in range(1,M+1)])
B=np.diag((w[1:]/np.sqrt(np.arange(1,M+1))))
T=np.diag(np.log(np.arange(1,M+1)))
for t in [3.7,21.0]:
    direct=sum(Lam(n)/math.sqrt(n)*w[n]*np.exp(-1j*t*math.log(n)) for n in range(1,M+1))
    op1=lam[1:]@B@np.diag(np.exp(-1j*t*np.diag(T)))@np.ones(M)
    op2=ell[1:]@Zi[1:,1:]@B@np.diag(np.exp(-1j*t*np.diag(T)))@np.ones(M)
    print(f"   t={t:5.2f}: direct {direct:.8f} | lambda^T B e^-itT 1 diff {abs(direct-op1):.1e} | ell^T Z^inv B e^-itT 1 diff {abs(direct-op2):.1e}")
print("   -> sections 7 and 8 are exact; they are Mobius inversion written in matrix form.")
