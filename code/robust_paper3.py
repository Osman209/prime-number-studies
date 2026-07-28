import numpy as np, math, pickle
from sympy import primerange

W=2*3*5*7*11*13; PW=[2,3,5,7,11,13]
phiW=1
for p in PW: phiW*=p-1
def formula(H):
    r=1
    for p in PW: r*=(p-1) if H%p==0 else (p-2)
    return r
def S(H): return W*formula(H)/phiW**2

def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return np.flatnonzero(s)
def window_gaps(X,DX=40_000_000):
    lo,hi=X,X+DX
    base=primes_upto(int(math.isqrt(hi))+1)
    seg=np.ones(hi-lo,bool)
    for p in base:
        st=max(((lo+p-1)//p)*p, p*p)
        if st<hi: seg[st-lo::p]=False
    return np.diff(np.flatnonzero(seg)+lo)

DATA={}
for X in [10**7,10**9,10**11,10**13]:
    g=window_gaps(X); hs,cs=np.unique(g,return_counts=True)
    DATA[X]=(dict(zip(hs.tolist(),cs.tolist())),math.log(X))
pickle.dump(DATA,open("gapdata.pkl","wb"))

def prep(X,Hmax,minc=20):
    A,lg=DATA[X]
    Hs=[H for H in range(2,Hmax+1,2) if A.get(H,0)>=minc]
    return (np.array(Hs,float),
            np.array([math.log(A[H]) for H in Hs]),
            np.array([math.log(S(H)) for H in Hs]), lg)
def r2(y,f): return 1-np.sum((y-f)**2)/np.sum((y-y.mean())**2)
def fit_beta(X,Hmax):
    Hs,y,s,lg=prep(X,Hmax); best=(-9,None,None)
    for beta in np.arange(0.70,1.70,0.0005):
        M=np.vstack([np.ones_like(Hs),-Hs**beta]).T
        c,*_=np.linalg.lstsq(M,y-s,rcond=None)
        R=r2(y,s+M@c)
        if R>best[0]: best=(R,beta,c[1])
    return best[1],best[2],best[0]

print("="*74); print("R1  IS beta STABLE IN THE FITTING WINDOW?  (a real exponent must not move)"); print("="*74)
print(f"  {'X':>6} " + " ".join(f"Hmax={h:<4d}" for h in [60,80,120,160,200,240]))
for X in [10**7,10**9,10**11,10**13]:
    row=[]
    for Hmax in [60,80,120,160,200,240]:
        b,_,_=fit_beta(X,Hmax); row.append(f"{b:9.4f}")
    print(f"  1e{int(math.log10(X)):<4d} "+" ".join(row))
print("  -> beta moves MORE with the H window than it does with the height X")

print()
print("="*74); print("R2  IS beta STABLE IN HEIGHT AT FIXED WINDOW?"); print("="*74)
for Hmax in [80,120,200]:
    bs=[fit_beta(X,Hmax)[0] for X in [10**7,10**9,10**11,10**13]]
    print(f"  Hmax={Hmax:<4d} beta = "+", ".join(f"{b:.4f}" for b in bs)+
          f"    total drift over 6 decades = {bs[-1]-bs[0]:+.4f}")

print()
print("="*74); print("R3  OUT-OF-SAMPLE TRANSFER 1e7 -> 1e13, swept over the H window"); print("="*74)
print(f"  {'Hmax':>5} {'zero-shape R2':>14} {'err%':>7} | {'transp.beta R2':>15} {'err%':>7} | {'winner':>12}")
for Hmax in [60,80,120,160,200,240]:
    Hs7,y7,s7,lg7=prep(10**7,Hmax); Hs3,y3,s3,lg3=prep(10**13,Hmax)
    b7,l7,_=fit_beta(10**7,Hmax)
    def sc(f): return r2(y3,f), float(np.mean(np.abs(np.exp(f)-np.exp(y3))/np.exp(y3))*100)
    t0=s3-Hs3/lg3;             f0=y3.mean()+t0-t0.mean()
    lam=l7*(lg7/lg3)**b7
    t1=s3-lam*Hs3**b7;         f1=y3.mean()+t1-t1.mean()
    R0,e0=sc(f0); R1,e1=sc(f1)
    print(f"  {Hmax:>5} {R0:>14.5f} {e0:>7.1f} | {R1:>15.5f} {e1:>7.1f} | {'beta' if e1<e0 else 'zero-shape':>12}")

print()
print("="*74); print("R4  PARAMETER-MATCHED CONTROL: give the zero-shape model ONE transported number too"); print("="*74)
print("     (model D: beta=1 but lambda*logX transported from 1e7 instead of forced to 1)")
print(f"  {'Hmax':>5} {'D: transp. lambda*logX':>23} {'err%':>7} | {'C: transp. beta':>16} {'err%':>7}")
for Hmax in [80,120,200]:
    Hs7,y7,s7,lg7=prep(10**7,Hmax); Hs3,y3,s3,lg3=prep(10**13,Hmax)
    M=np.vstack([np.ones_like(Hs7),-Hs7]).T; c,*_=np.linalg.lstsq(M,y7-s7,rcond=None)
    lamD=c[1]*lg7/lg3
    b7,l7,_=fit_beta(10**7,Hmax); lamC=l7*(lg7/lg3)**b7
    def sc(t): f=y3.mean()+t-t.mean(); return r2(y3,f), float(np.mean(np.abs(np.exp(f)-np.exp(y3))/np.exp(y3))*100)
    RD,eD=sc(s3-lamD*Hs3); RC,eC=sc(s3-lamC*Hs3**b7)
    print(f"  {Hmax:>5} {RD:>23.5f} {eD:>7.1f} | {RC:>16.5f} {eC:>7.1f}")

print()
print("="*74); print("R5  WHAT IS lambda ACTUALLY?  fit lambda at beta=1 and compare with 1/log X"); print("="*74)
print(f"  {'X':>6} {'logX':>7} {'lambda(fit)':>12} {'1/logX':>9} {'lambda*logX':>12}")
for X in [10**7,10**9,10**11,10**13]:
    Hs,y,s,lg=prep(X,120)
    M=np.vstack([np.ones_like(Hs),-Hs]).T; c,*_=np.linalg.lstsq(M,y-s,rcond=None)
    print(f"  1e{int(math.log10(X)):<4d} {lg:>7.3f} {c[1]:>12.5f} {1/lg:>9.5f} {c[1]*lg:>12.4f}")
print("  -> lambda*logX is close to but not exactly 1, and drifts; that residual is what beta absorbs")
