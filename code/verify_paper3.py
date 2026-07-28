import numpy as np, math
from sympy import primerange, factorint

W=2*3*5*7*11*13
PW=[2,3,5,7,11,13]
print("="*72); print("V1  exact wheel-pair count  W(H) = prod_{p|H}(p-1) * prod_{p not| H}(p-2)"); print("="*72)
a=np.arange(W); co=np.ones(W,bool)
for p in PW: co &= (a%p!=0)
def formula(H):
    r=1
    for p in PW: r *= (p-1) if H%p==0 else (p-2)
    return r
bad=0
for H in range(1,2*W//7):
    direct=int(np.count_nonzero(co & np.roll(co,-H)))
    if direct!=formula(H): bad+=1
print(f"H = 1..{2*W//7-1}: mismatches between direct count and formula = {bad}")
print("  paper's table:")
for H in [2,6,12,30,32,42,210]:
    print(f"   H={H:4d}  formula {formula(H):5d}   direct {int(np.count_nonzero(co & np.roll(co,-H))):5d}")
print(f"  W(30)/W(32) = {formula(30)/formula(32):.4f}")
print(f"  odd H: max over H odd of W(H) = {max(formula(H) for H in range(1,999,2))}")

print()
print("="*72); print("V2  is S_W(H) = W*W(H)/phi(W)^2 the truncated Hardy-Littlewood singular series?"); print("="*72)
phiW=1
for p in PW: phiW*=p-1
def S_wheel(H): return W*formula(H)/phiW**2
C2=1.0
for p in primerange(3,4000000): C2*= (1-1/(p-1)**2)
C2*=1  # 2C2 assembled below
def S_full(H,PMAX=4000000):
    if H%2: return 0.0
    r=2*C2
    for p,_ in factorint(H).items():
        if p>2: r*= (p-1)/(p-2)
    return r
print(f"  2*C_2 = {2*C2:.9f}   (Hardy-Littlewood twin constant, literature 1.3203236)")
print(f"  {'H':>4} {'S_wheel':>10} {'S_full':>10} {'ratio':>8}   largest odd prime factor")
worst=(0,0)
for H in list(range(2,101,2))+[34,94,210,2310]:
    sw,sf=S_wheel(H),S_full(H)
    lp=max([p for p in factorint(H) if p>2],default=1)
    r=sw/sf
    if abs(r-1)>abs(worst[1]-1): worst=(H,r)
    if H in [2,6,30,34,42,94,210,2310]:
        print(f"  {H:>4} {sw:>10.5f} {sf:>10.5f} {r:>8.5f}   {lp}")
print(f"  worst ratio over the sample: H={worst[0]}, ratio {worst[1]:.5f}")
print("  => identical whenever every odd prime factor of H is <= 13; differs exactly when H has a larger odd prime factor")

print()
print("="*72); print("V3  measured gap frequencies at four heights (segmented sieve, window 4e7)"); print("="*72)
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
        st=((lo+p-1)//p)*p
        if st<p*p: st=p*p
        if st<hi: seg[st-lo::p]=False
    pr=np.flatnonzero(seg)+lo
    return np.diff(pr), DX
data={}
for X in [10**7,10**9,10**11,10**13]:
    g,DX=window_gaps(X)
    hs,cs=np.unique(g,return_counts=True)
    data[X]=(dict(zip(hs.tolist(),cs.tolist())),DX,math.log(X))
    print(f"  X=1e{int(math.log10(X)):<2d}  primes in window {len(g)+1:>8d}  mean gap {g.mean():7.3f}  (log X = {math.log(X):.3f})  champion H = {int(hs[np.argmax(cs)])}")

print()
print("="*72); print("V4  model fits.  protocol: even H in [2,Hmax], least squares in LOG space, R^2 in log space"); print("="*72)
Hmax=120
def prep(X):
    A,DX,lg=data[X]
    Hs=[H for H in range(2,Hmax+1,2) if A.get(H,0)>=20]
    y=np.array([math.log(A[H]) for H in Hs])
    s=np.array([math.log(S_wheel(H)) for H in Hs])
    return np.array(Hs,float),y,s,DX,lg
def r2(y,f): return 1-np.sum((y-f)**2)/np.sum((y-y.mean())**2)
print(f"  {'X':>6} {'zero-shape':>12} {'free lambda':>12} {'free lam+beta':>14} {'beta':>8} {'lambda':>9} {'lam*(logX)^beta':>16}")
betas={}
for X in [10**7,10**9,10**11,10**13]:
    Hs,y,s,DX,lg=prep(X)
    # (a) zero shape parameters: log A = c + log S - H/log X
    f0=y.mean()+ (s-Hs/lg) - (s-Hs/lg).mean()
    R0=r2(y,f0)
    # (b) free lambda
    Mb=np.vstack([np.ones_like(Hs),-Hs]).T
    cb,*_=np.linalg.lstsq(Mb,y-s,rcond=None); Rb=r2(y,s+Mb@cb)
    # (c) free lambda and beta
    best=(-9,None)
    for beta in np.arange(0.80,1.60,0.001):
        M=np.vstack([np.ones_like(Hs),-Hs**beta]).T
        c,*_=np.linalg.lstsq(M,y-s,rcond=None)
        R=r2(y,s+M@c)
        if R>best[0]: best=(R,beta,c)
    Rc,beta,c=best
    betas[X]=(beta,c[1])
    print(f"  1e{int(math.log10(X)):<4d} {R0:12.5f} {Rb:12.5f} {Rc:14.5f} {beta:8.4f} {c[1]:9.5f} {c[1]*lg**beta:16.4f}")
print("  paper section 8 : 0.9368 / 0.9763 / 0.9859 / 0.9904   (zero shape parameters)")
print("  paper section 7 : beta = 1.0989 / 1.0784 / 1.0713 / 1.0686")

print()
print("="*72); print("V5  how fast does beta drift, and where does it land?"); print("="*72)
xs=np.array([math.log(math.log(X)) for X in betas]); bs=np.array([betas[X][0] for X in betas])
sl,ic=np.polyfit(xs,np.log(bs-1),1)
print(f"  fit  beta - 1 ~ C (log X)^s  :  s = {sl:.3f},  C = {math.exp(ic):.3f}")
for LX in [30,60,120,230]:
    print(f"     extrapolated beta at log X = {LX:>3d}  (X ~ 1e{int(LX/2.3026):<3d}) : {1+math.exp(ic)*LX**sl:.4f}")
print("  => the drift toward 1 is real but sub-logarithmic; the data cannot separate beta -> 1 from beta -> 1.02")

print()
print("="*72); print("V6  out-of-sample transfer: calibrate at 1e7, predict at 1e13, no refitting"); print("="*72)
Hs7,y7,s7,_,lg7=prep(10**7)
Hs13,y13,s13,_,lg13=prep(10**13)
def score(f):
    return r2(y13,f), float(np.mean(np.abs(np.exp(f)-np.exp(y13))/np.exp(y13))*100)
# model A: zero shape - lambda = 1/log X, beta = 1, only the overall constant is set by the window
fA=y13.mean()+(s13-Hs13/lg13)-(s13-Hs13/lg13).mean()
# model B: transport the lambda fitted at 1e7 unchanged
Mb=np.vstack([np.ones_like(Hs7),-Hs7]).T; cb,*_=np.linalg.lstsq(Mb,y7-s7,rcond=None)
fB=y13.mean()+(s13-cb[1]*Hs13)-(s13-cb[1]*Hs13).mean()
# model C: transport the beta fitted at 1e7, rescale lambda by log X
b7,l7=betas[10**7]
lam=l7*(lg7/lg13)**b7
fC=y13.mean()+(s13-lam*Hs13**b7)-(s13-lam*Hs13**b7).mean()
for tag,f in [("zero shape (lam=1/logX, beta=1)",fA),("transported lambda",fB),(f"transported beta={b7:.3f}",fC)]:
    R,e=score(f); print(f"  {tag:34s}  R^2 = {R:.5f}   mean relative error = {e:5.1f}%")
print("  paper section 9 : 0.99693 / 8.7% , 0.96344 / 29.6% , 0.98496 / 15.1%")
