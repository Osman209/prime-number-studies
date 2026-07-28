import numpy as np, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sympy import factorint, isprime, primerange

plt.rcParams.update({"font.size":9,"axes.titlesize":10,"figure.dpi":170})
fig,ax=plt.subplots(2,2,figsize=(13,9.6))

# ---------- (A) the remainder field R(H,j) ----------
Hm,Jm=60,30
R=np.zeros((Hm,Jm))
for H in range(Hm):
    for j in range(1,Jm+1):
        R[H,j-1]=(2*(H-j+1))%(2*j+1)/(2*j+1)
a=ax[0,0]
im=a.imshow(R,origin="lower",aspect="auto",cmap="viridis",
            extent=[3-1,2*Jm+1+1,3-1,2*Hm+3-1+1])
ys,xs=[],[]
for H in range(Hm):
    for j in range(1,Jm+1):
        if (2*H+3)%(2*j+1)==0: xs.append(2*j+1); ys.append(2*H+3)
a.scatter(xs,ys,s=9,c="red",marker="s")
a.set_xlabel("denominator  d = 2j+1"); a.set_ylabel("numerator  n = 2H+3")
a.set_title("(A) remainder field  R/d = (n mod d)/d\nred squares = exact cells (d | n)")
plt.colorbar(im,ax=a,fraction=0.046,label="r/d")

# ---------- (B) one column is one arithmetic progression ----------
a=ax[0,1]
cols=[3,5,7,9,11,13,25]
for i,d in enumerate(cols):
    off=(d-3)//2
    Hs=[H for H in range(0,80) if (H-off)%d==0]
    a.scatter([i]*len(Hs),Hs,s=26,marker="_",linewidths=2.2,
              c=("tab:red" if len(factorint(d))==1 else "tab:blue"))
a.set_xticks(range(len(cols)));a.set_xticklabels([f"$Z_{{{d}}}$" for d in cols])
a.set_ylabel("row index H")
a.set_title("(B) each column's exact cells = one arithmetic progression\n"
            r"$Z_d=\{H\equiv (d-3)/2 \; \mathrm{mod}\; d\}$   red = prime power, blue = composite ($Z_9\!\subset\!Z_3$, $Z_{25}\!\subset\!Z_5$)")
a.grid(alpha=.25,axis="y")

# ---------- (C) Phi line geometry ----------
a=ax[1,0]
Hmax=140
for t in range(1,9):
    js=np.arange(1,(Hmax-t+1)//(2*t+1)+1)
    if len(js)==0: continue
    a.plot(js,(2*t+1)*js+t-1,"o-",ms=3,lw=.9,label=f"t={t}  (q={2*t+1})")
prim=[H for H in range(Hmax) if isprime(2*H+3)]
a.scatter([0.55]*len(prim),prim,s=14,c="k",marker=">")
a.text(0.7,Hmax*0.94,"black arrows = H values\nnever hit  ⇒  n = 2H+3 prime",fontsize=8)
a.set_xlim(0,12); a.set_ylim(0,Hmax)
a.set_xlabel("j"); a.set_ylabel(r"$\Phi(j,t)=2jt+j+t-1$")
a.set_title(r"(C) $\Phi$ as a family of lines of slope $q=2t+1$"+"\ncollisions = several factorizations of the same n")
a.legend(fontsize=7,ncol=2); a.grid(alpha=.25)

# ---------- (D) towers ----------
a=ax[1,1]
N=1000
prs=[3,5,7,11,13]
for i,p in enumerate(prs):
    k=1
    while p**k<=N:
        q=p**k
        cnt=len([n for n in range(3,N+1,2) if n%q==0])
        a.barh(i+0.16*(k-1),cnt,height=0.14,color=plt.cm.plasma(0.15+0.2*k),
               label=f"$p^{k}$" if i==0 else None)
        k+=1
a.set_yticks(range(len(prs)));a.set_yticklabels([f"p={p}" for p in prs])
a.set_xlabel(f"number of odd n ≤ {N} in the layer")
a.set_title(r"(D) nested towers  $Z_p\supseteq Z_{p^2}\supseteq Z_{p^3}\supseteq\cdots$"
            "\ndepth in the tower = the valuation $v_p(n)$")
a.legend(fontsize=8); a.grid(alpha=.25,axis="x")

plt.tight_layout()
plt.savefig("../figures/division_table_explained.png",bbox_inches="tight")
print("saved")
