import numpy as np, math
import mpmath as mp
from sympy import primerange
mp.mp.dps=20
def pp_arrays(X):
    a=[];lg=[]
    for p in primerange(2,int(X)+1):
        k=1;q=p
        while q<=X:
            c=(math.log(p)/q**0.5)*(1-k*math.log(p)/math.log(X))
            if c!=0: a.append(c); lg.append(math.log(q))
            k+=1; q=p**k
    return np.array(a),np.array(lg)
def Dabs(t,a,lg,A_):
    out=np.zeros(len(t),complex)
    for i in range(0,len(t),2000):
        tt=t[i:i+2000]; out[i:i+2000]=A_[i:i+2000]-np.exp(-1j*np.outer(tt,lg))@a
    return np.abs(out)

print("="*74); print("W5  THE 138/138 PEAK MATCH, AND THE CONTROL"); print("="*74)
X=1e5; a,lg=pp_arrays(X)
ts=np.arange(10,300,0.01)
A_=np.array([complex(1/mp.mpc(0.5,x)+1/(mp.mpc(0.5,x)-1)-0.5*mp.log(mp.pi)+0.5*mp.digamma(mp.mpc(0.5,x)/2)) for x in ts])
v=Dabs(ts,a,lg,A_)
loc=np.where((v[1:-1]>v[:-2])&(v[1:-1]>v[2:]))[0]+1
pk=ts[loc]
zs=[]; n=1
while True:
    z=float(mp.im(mp.zetazero(n)))
    if z>300: break
    if z>10: zs.append(z)
    n+=1
err=[np.min(np.abs(pk-z)) for z in zs]
print(f"   zeros in 10<t<300 : {len(zs)}     local maxima found : {len(pk)}")
print(f"   matched within 0.05 : {sum(1 for e in err if e<0.05)} / {len(zs)}     mean |peak-zero| = {np.mean(err):.5f}")
print(f"   SPURIOUS maxima (no zero within 0.05) : {sum(1 for p_ in pk if min(abs(p_-z) for z in zs)>=0.05)}")
print("   -> the match is FORCED: D_X approximates xi'/xi, whose poles ARE the zeros.")
print("      A pole-approximant peaks at its poles wherever they lie, so the count says nothing about Re(rho)=1/2.")

print()
print("="*74); print("W6  ANGULAR CRITERION Theta_eps = 2 arctan(eps/a), AND ITS COST"); print("="*74)
for aa in [0.0,0.01,0.05,0.2]:
    row=[]
    for eps in [0.001,0.01,0.1]:
        v1=complex(-aa,eps); v2=complex(-aa,-eps)
        ang=abs(math.atan2(v1.imag,v1.real)-math.atan2(v2.imag,v2.real)); ang=min(ang,2*math.pi-ang)
        th=2*math.atan(eps/aa) if aa>0 else math.pi
        row.append(f"eps={eps:<6}: {math.degrees(ang):7.2f} vs {math.degrees(th):7.2f}")
    print(f"   a={aa:5.2f}   "+" | ".join(row))
print("   formula confirmed exactly.  COST: the model resolves scales ~1/L = 1/log X, so seeing a needs log X >~ 1/a:")
for aa in [0.1,0.01,0.001]:
    print(f"      a={aa:6.3f}  ->  log X >~ {1/aa:7.1f}  ->  X >~ 1e{1/aa/2.302585:.0f}")

print()
print("="*74); print("W7  IS log X_crit ~ C/g ARITHMETIC, OR PURE FOURIER RESOLUTION?"); print("="*74)
def F(d,L):
    d=np.asarray(d,float); out=np.empty_like(d)
    m=np.abs(d)<1e-12; out[m]=L/2
    out[~m]=(1-np.cos(L*d[~m]))/(L*d[~m]**2)
    return out
def resolved(L,g):
    x=np.linspace(-1.2*g,2.2*g,8001); y=F(x,L)+F(x-g,L)
    return np.sum((y[1:-1]>y[:-2])&(y[1:-1]>y[2:]))>=2
print("   two Fejer bumps separated by g; bisect for the L at which two maxima first appear")
print(f"   {'g':>12} {'L_crit':>10} {'L_crit * g':>12}")
pr=[]
for g in [0.2,0.4,0.6,0.845123634,1.2,2.0,3.0]:
    lo,hi=0.5,6000.0
    for _ in range(70):
        m=(lo+hi)/2
        if resolved(m,g): hi=m
        else: lo=m
    pr.append(hi*g); print(f"   {g:>12.6f} {hi:>10.3f} {hi*g:>12.5f}")
print(f"   spread of L_crit*g = {100*(max(pr)-min(pr))/np.mean(pr):.4f}%  -> EXACTLY scale invariant")
print("   => log X_crit ~ C/g is the Rayleigh resolution criterion for this kernel. No arithmetic content.")
print("      The reported g*logX_crit values 7.807/8.394/9.296/11.952 vary by 53%, i.e. they are NOT a clean C/g law;")
print("      the variation is local background tilt, which is exactly what the paper's own B'(0) term describes.")
