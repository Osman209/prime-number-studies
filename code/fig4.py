import numpy as np, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
def Dv(t,a,lg,A_):
    out=np.zeros(len(t),complex)
    for i in range(0,len(t),1500):
        tt=t[i:i+1500]; out[i:i+1500]=A_[i:i+1500]-np.exp(-1j*np.outer(tt,lg))@a
    return out
def Af(t): return np.array([complex(1/mp.mpc(0.5,x)+1/(mp.mpc(0.5,x)-1)-0.5*mp.log(mp.pi)+0.5*mp.digamma(mp.mpc(0.5,x)/2)) for x in t])

plt.rcParams.update({"font.size":9,"axes.titlesize":10,"figure.dpi":170})
fig,ax=plt.subplots(2,2,figsize=(13,9))

ts=np.arange(10,60,0.01); A_=Af(ts); a,lg=pp_arrays(1e6)
d=Dv(ts,a,lg,A_)
zs=[float(mp.im(mp.zetazero(n))) for n in range(1,20)]; zs=[z for z in zs if 10<z<60]
aa=ax[0,0]
aa.plot(ts,np.abs(d),lw=.9,label=r"$|D_X(t)|$,  $X=10^6$")
for z in zs: aa.axvline(z,color="crimson",ls=":",lw=.8)
aa.set_xlabel("t"); aa.set_ylabel(r"$|D_X|$"); aa.set_ylim(0,None)
aa.set_title("(A) peaks of $|D_X|$ sit on the zeros (dotted)\nbecause $D_X$ is a truncated $\\xi'/\\xi$, whose poles are the zeros")
aa.legend(fontsize=8); aa.grid(alpha=.25)

aa=ax[0,1]
tr=np.array([complex(1/mp.mpc(0.5,x)+1/(mp.mpc(0.5,x)-1)-0.5*mp.log(mp.pi)+0.5*mp.digamma(mp.mpc(0.5,x)/2)+mp.zeta(mp.mpc(0.5,x),derivative=1)/mp.zeta(mp.mpc(0.5,x))) for x in ts[::20]])
aa.plot(ts[::20],tr.real,lw=1.2,label=r"$\mathrm{Re}\,\xi'/\xi$  (exact)")
aa.plot(ts,d.real,lw=.8,alpha=.85,label=r"$\mathrm{Re}\,D_X$  (model, $X=10^6$)")
aa.set_xlabel("t"); aa.set_ylabel("real part")
aa.set_title(r"(B) $\xi'/\xi(1/2+it)$ is purely imaginary on the line"+"\nthe model's real part is entirely truncation error")
aa.legend(fontsize=8); aa.grid(alpha=.25)

aa=ax[1,0]
g1=float(mp.im(mp.zetazero(1)))
Ls=[];hs=[];cs=[]
for X in [1e3,1e4,1e5,1e6]:
    a2,l2=pp_arrays(X); L=math.log(X); h=0.004
    A2=Af(np.array([g1-h,g1,g1+h]))
    Q=np.abs(Dv(np.array([g1-h,g1,g1+h]),a2,l2,A2))**2
    Ls.append(L); hs.append(Q[1]); cs.append(abs((Q[0]-2*Q[1]+Q[2])/h**2))
Ls=np.array(Ls)
aa.loglog(Ls,hs,"o-",label=r"measured $|D_X(\gamma_1)|^2$")
aa.loglog(Ls,(Ls/2)**2,"s--",label=r"$F_L(0)^2=(L/2)^2$")
aa.loglog(Ls,Ls/2,"^:",label=r"$F_L(0)=L/2$  (what sec.5 assumes)")
aa.set_xlabel(r"$L=\log X$"); aa.set_ylabel("peak height")
aa.set_title("(C) the peak height scales as $F_L(0)^2$, not $F_L(0)$\nso $|D_X|^2=F_L+B$ mis-states the leading term")
aa.legend(fontsize=8); aa.grid(alpha=.25,which="both")

aa=ax[1,1]
def F(dd,L):
    dd=np.asarray(dd,float); o=np.empty_like(dd); m=np.abs(dd)<1e-12
    o[m]=L/2; o[~m]=(1-np.cos(L*dd[~m]))/(L*dd[~m]**2); return o
def resolved(L,g):
    x=np.linspace(-1.2*g,2.2*g,8001); y=F(x,L)+F(x-g,L)
    return np.sum((y[1:-1]>y[:-2])&(y[1:-1]>y[2:]))>=2
gs=np.array([0.2,0.3,0.4,0.6,0.845123634,1.2,2.0,3.0]); Lc=[]
for g in gs:
    lo,hi=0.5,6000.0
    for _ in range(70):
        m=(lo+hi)/2
        if resolved(m,g): hi=m
        else: lo=m
    Lc.append(hi)
Lc=np.array(Lc)
aa.loglog(gs,Lc,"o-",label=r"measured $L_{crit}$")
aa.loglog(gs,5.21233/gs,"--",label=r"$5.21233/g$")
aa.set_xlabel("separation g"); aa.set_ylabel(r"$L_{crit}=\log X_{crit}$")
aa.set_title(r"(D) $L_{crit}\cdot g = 5.21233$ to 0.0000%"+"\npure Rayleigh resolution of the kernel, no arithmetic content")
aa.legend(fontsize=8); aa.grid(alpha=.25,which="both")

plt.tight_layout(); plt.savefig("../figures/zeta_dynamical_figures.png",bbox_inches="tight")
print("saved")
