import numpy as np, math
import mpmath as mp
from sympy import primerange
mp.mp.dps=25

def pp_arrays(X):
    a=[];lg=[]
    for p in primerange(2,int(X)+1):
        k=1;q=p
        while q<=X:
            c=(math.log(p)/q**0.5)*(1-k*math.log(p)/math.log(X))
            if c!=0: a.append(c); lg.append(math.log(q))
            k+=1; q=p**k
    return np.array(a),np.array(lg)
def PX(t,a,lg,chunk=200000):
    t=np.atleast_1d(np.asarray(t,float)); out=np.zeros(len(t),complex)
    for i in range(0,len(a),chunk):
        out+=np.exp(-1j*np.outer(t,lg[i:i+chunk]))@a[i:i+chunk]
    return out
def A(t):
    t=np.atleast_1d(np.asarray(t,float))
    return np.array([complex(1/mp.mpc(0.5,x)+1/(mp.mpc(0.5,x)-1)-0.5*mp.log(mp.pi)+0.5*mp.digamma(mp.mpc(0.5,x)/2)) for x in t])
def D(t,a,lg): return A(t)-PX(t,a,lg)
def xiratio(t):
    s=mp.mpc(0.5,t)
    return complex(1/s+1/(s-1)-0.5*mp.log(mp.pi)+0.5*mp.digamma(s/2)+mp.zeta(s,derivative=1)/mp.zeta(s))

print("="*74); print("W1  IS D_X A TRUNCATED xi'/xi ?  (reference computed independently by mpmath)"); print("="*74)
X=1e6; a,lg=pp_arrays(X)
print(f"   {'t':>6} {'D_X (model)':>28} {'xi\'/xi (true)':>28} {'|diff|':>9}")
for t in [7.0,18.0,33.0,60.0,95.0]:
    d=D(t,a,lg)[0]; tr=xiratio(t)
    print(f"   {t:>6.1f} {d.real:>+13.5f}{d.imag:>+13.5f}i {tr.real:>+13.5f}{tr.imag:>+13.5f}i {abs(d-tr):>9.4f}")

print()
print("="*74); print("W2  Re(xi'/xi) = 0 EXACTLY ON THE LINE  (functional equation)"); print("="*74)
ts=np.linspace(5,60,25); d=D(ts,a,lg); tr=np.array([xiratio(x) for x in ts])
print(f"   true      : max|Re| = {np.abs(tr.real).max():.2e}   mean|Im| = {np.abs(tr.imag).mean():.4f}")
print(f"   model D_X : max|Re| = {np.abs(d.real).max():.4f}      mean|Im| = {np.abs(d.imag).mean():.4f}")
print("   -> Xi(t)=xi(1/2+it) is real, so xi'/xi = -i Xi'/Xi is PURELY IMAGINARY on the line.")
print("      The real part of D_X is entirely truncation error; the 'complex-plane curve' of sec.6 is an artifact of it.")

print()
print("="*74); print("W4  SHAPE AND HEIGHT OF |D_X|^2 AT A ZERO"); print("="*74)
g1=float(mp.im(mp.zetazero(1)))
print(f"   {'X':>6} {'L':>6} {'|D|^2(g)':>10} {'L/2':>7} {'(L/2)^2':>9} {'curv(meas)':>12} {'-L^3/12':>10} {'-L^4/12':>10}")
for X in [1e3,1e4,1e5,1e6]:
    aa,ll=pp_arrays(X); L=math.log(X)
    h=0.004; Q=np.abs(D(np.array([g1-h,g1,g1+h]),aa,ll))**2
    curv=(Q[0]-2*Q[1]+Q[2])/h**2
    print(f"   1e{int(math.log10(X)):<4d} {L:>6.2f} {Q[1]:>10.3f} {L/2:>7.3f} {L*L/4:>9.3f} {curv:>12.1f} {-L**3/12:>10.1f} {-L**4/12:>10.1f}")
print("   -> the peak HEIGHT tracks (L/2)^2 = F_L(0)^2, not F_L(0)=L/2.")
print("      So |D_X|^2 ~ F_L(delta)^2 + ... ; writing |D_X|^2 = F_L + B mismatches the leading term by a square.")
