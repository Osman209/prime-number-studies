import numpy as np, math, pickle
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W=2*3*5*7*11*13; PW=[2,3,5,7,11,13]
phiW=1
for p in PW: phiW*=p-1
def formula(H):
    r=1
    for p in PW: r*=(p-1) if H%p==0 else (p-2)
    return r
def S(H): return W*formula(H)/phiW**2
DATA=pickle.load(open("gapdata.pkl","rb"))

def prep(X,Hmax,minc=20):
    A,lg=DATA[X]
    Hs=[H for H in range(2,Hmax+1,2) if A.get(H,0)>=minc]
    return (np.array(Hs,float),np.array([math.log(A[H]) for H in Hs]),
            np.array([math.log(S(H)) for H in Hs]),lg)
def r2(y,f): return 1-np.sum((y-f)**2)/np.sum((y-y.mean())**2)
def fit_beta(X,Hmax):
    Hs,y,s,lg=prep(X,Hmax); best=(-9,0,0)
    for beta in np.arange(0.70,1.70,0.0005):
        M=np.vstack([np.ones_like(Hs),-Hs**beta]).T
        c,*_=np.linalg.lstsq(M,y-s,rcond=None); R=r2(y,s+M@c)
        if R>best[0]: best=(R,beta,c[1])
    return best[1],best[2]

plt.rcParams.update({"font.size":9,"axes.titlesize":10,"figure.dpi":170})
fig,ax=plt.subplots(2,2,figsize=(13,9))

a=ax[0,0]
Hs=np.arange(2,213,2); Ws=[formula(int(h)) for h in Hs]
a.stem(Hs,Ws,basefmt=" ",markerfmt="o",linefmt="C0-")
for h in [6,30,210]:
    a.annotate(f"H={h}",(h,formula(h)),textcoords="offset points",xytext=(0,7),ha="center",fontsize=8,color="crimson")
a.set_xlabel("even H"); a.set_ylabel("W(H)")
a.set_title("(A) exact wheel-pair count, W = 30030\npeaks at H divisible by many small primes")
a.grid(alpha=.25)

a=ax[0,1]
Hmaxs=[60,80,120,160,200,240]
for X,mk in zip([10**7,10**9,10**11,10**13],["o","s","^","d"]):
    bs=[fit_beta(X,hm)[0] for hm in Hmaxs]
    a.plot(Hmaxs,bs,mk+"-",ms=4,label=f"X = 1e{int(math.log10(X))}")
a.axhline(1,ls="--",c="k",lw=.9)
a.set_xlabel("fitting window  $H_{max}$"); a.set_ylabel(r"fitted $\beta$")
a.set_title(r"(B) $\beta$ is not an exponent: it moves ~7x more with the"+"\nfitting window than with six decades of height")
a.legend(fontsize=8); a.grid(alpha=.25)

a=ax[1,0]
X=10**13; Hmax=200
Hsx,y,s,lg=prep(X,Hmax)
b7,l7=fit_beta(10**7,Hmax); _,lg7=DATA[10**7][0],DATA[10**7][1]
a.semilogy(Hsx,np.exp(y),"ko",ms=4,label="measured counts, window $4\\times10^7$ at $10^{13}$")
t0=s-Hsx/lg; f0=y.mean()+t0-t0.mean()
a.semilogy(Hsx,np.exp(f0),"-",lw=1.4,label=r"$S(H)e^{-H/\log X}$  (no fitted shape)")
lam=l7*(lg7/lg)**b7; t1=s-lam*Hsx**b7; f1=y.mean()+t1-t1.mean()
a.semilogy(Hsx,np.exp(f1),"--",lw=1.4,label=r"$S(H)e^{-\lambda H^{\beta}}$, $\beta$ from $10^{7}$")
a.set_xlabel("gap H"); a.set_ylabel("occurrences")
a.set_title("(C) measured gap frequencies against the two models\n(both curves are transported, not refitted here)")
a.legend(fontsize=8); a.grid(alpha=.25,which="both")

a=ax[1,1]
e0l,e1l,eDl=[],[],[]
for hm in Hmaxs:
    Hs7,y7,s7,lg7=prep(10**7,hm); Hs3,y3,s3,lg3=prep(10**13,hm)
    b,l=fit_beta(10**7,hm)
    def er(t):
        f=y3.mean()+t-t.mean(); return float(np.mean(np.abs(np.exp(f)-np.exp(y3))/np.exp(y3))*100)
    M=np.vstack([np.ones_like(Hs7),-Hs7]).T; c,*_=np.linalg.lstsq(M,y7-s7,rcond=None)
    e0l.append(er(s3-Hs3/lg3)); e1l.append(er(s3-l*(lg7/lg3)**b*Hs3**b)); eDl.append(er(s3-c[1]*lg7/lg3*Hs3))
a.plot(Hmaxs,e0l,"o-",label=r"$\beta=1,\ \lambda=1/\log X$  (0 transported)")
a.plot(Hmaxs,eDl,"^-",label=r"$\beta=1,\ \lambda\log X$ transported  (1 transported)")
a.plot(Hmaxs,e1l,"s-",label=r"$\beta$ and $\lambda$ transported  (2 transported)")
a.set_xlabel("fitting window  $H_{max}$"); a.set_ylabel("mean relative error at $10^{13}$  (%)")
a.set_title("(D) out-of-sample transfer $10^{7}\\to10^{13}$\nmatched on transported quantities, the gap nearly closes")
a.legend(fontsize=8); a.grid(alpha=.25)

plt.tight_layout(); plt.savefig("../figures/prime_gap_wheel_figures.png",bbox_inches="tight")
print("saved")
