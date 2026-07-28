import numpy as np, math
import mpmath as mp
from sympy import primerange
from scipy.signal import find_peaks
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
def absD(t,a,lg,A_):
    out=np.zeros(len(t),complex)
    for i in range(0,len(t),1500):
        tt=t[i:i+1500]; out[i:i+1500]=A_[i:i+1500]-np.exp(-1j*np.outer(tt,lg))@a
    return np.abs(out)
ts=np.arange(10,100,0.005)
A_=np.array([complex(1/mp.mpc(0.5,x)+1/(mp.mpc(0.5,x)-1)-0.5*mp.log(mp.pi)+0.5*mp.digamma(mp.mpc(0.5,x)/2)) for x in ts])
zs=[]; n=1
while True:
    z=float(mp.im(mp.zetazero(n)))
    if z>100: break
    if z>10: zs.append(z)
    n+=1
zs=np.array(zs)
print("="*76); print("W5c  PEAK EXTRACTION WITH A PROMINENCE FILTER  (29 zeros in 10<t<100)"); print("="*76)
sets={}
for X in [1e4,1e5,1e6]:
    a,lg=pp_arrays(X); v=absD(ts,a,lg,A_)
    print(f"   X=1e{int(math.log10(X))}:")
    for prom in [0.0,0.2,0.5,1.0]:
        idx,_=find_peaks(v,prominence=prom*np.std(v))
        pk=ts[idx]
        m=sum(1 for z in zs if np.min(np.abs(pk-z))<0.05)
        sp=sum(1 for p_ in pk if np.min(np.abs(zs-p_))>=0.05)
        print(f"      prominence {prom:>4.1f} sd : peaks {len(pk):4d}   zeros hit(<0.05) {m:3d}/29   spurious {sp:4d}")
        if prom==0.5: sets[X]=pk
print()
print("   MATCH QUALITY vs TOLERANCE (prominence 0.5 sd, X=1e6):")
pk=sets[1e6]
for tol in [0.02,0.05,0.10,0.20,0.40]:
    m=sum(1 for z in zs if np.min(np.abs(pk-z))<tol)
    sp=sum(1 for p_ in pk if np.min(np.abs(zs-p_))>=tol)
    print(f"      tol {tol:>5.2f} :  zeros hit {m:3d}/29    spurious {sp:3d}")
print(f"   mean |nearest peak - zero| over the 29 zeros = {np.mean([np.min(np.abs(pk-z)) for z in zs]):.4f}")
print(f"   mean zero spacing here ~ {np.mean(np.diff(zs)):.3f}, so that is {100*np.mean([np.min(np.abs(pk-z)) for z in zs])/np.mean(np.diff(zs)):.1f}% of a spacing")
print()
stable=np.array([p_ for p_ in sets[1e6] if all(np.min(np.abs(sets[X]-p_))<0.2 for X in [1e4,1e5])])
for tol in [0.05,0.20,0.40]:
    m=sum(1 for z in zs if np.min(np.abs(stable-z))<tol)
    print(f"   stability filter + tol {tol:>4.2f} : {len(stable)} stable peaks, zeros hit {m}/29")
