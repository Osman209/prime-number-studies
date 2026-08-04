"""
detection_index_sweep.py -- the 24-run table of section 3.3 of
"Matching the Lens to the Zero": does the detection-index gain stay at T?

    ZETA_ZEROS=/path/to/zeros.txt python detection_index_sweep.py

Six heights, two displacements, two ordinate counts.  Reports gain/T for each
run and then the mean, the standard deviation and the slope against log T.
Exits nonzero if the mean departs from 1 by more than 2%, or if the slope is
clearly positive -- a d-dependent correction, which the note argues against,
would show up there.

THE ORDINATE LIST.  Read from ZETA_ZEROS: one ordinate per line, increasing;
8000 are needed.  Such a list can be taken from Odlyzko's tables
(www.dtc.umn.edu/~odlyzko/zeta_tables) or the LMFDB.  Without one the script
falls back to mpmath.zetazero, which is far slower -- minutes rather than
seconds -- and says so.

Python 3, numpy; mpmath only for the fallback.
"""
import numpy as np, os, sys, time

NEEDED = 8000

def load_zeros(n_wanted):
    path = os.environ.get("ZETA_ZEROS")
    if path and os.path.exists(path):
        g = np.loadtxt(path)
        if len(g) < n_wanted:
            sys.exit(f"ZETA_ZEROS={path} holds {len(g)} ordinates; "
                     f"this sweep needs {n_wanted}.")
        print(f"ordinates: {len(g)} from {path}, "
              f"gamma from {g[0]:.4f} to {g[-1]:.4f}")
        return g
    if path:
        print(f"ZETA_ZEROS={path!r} does not exist -- falling back to mpmath")
    else:
        print("no ZETA_ZEROS set -- falling back to mpmath.zetazero, which is")
        print("far slower; set it to a list of ordinates to run this in seconds")
    import mpmath as mp
    mp.mp.dps = 15
    return np.array([float(mp.im(mp.zetazero(k)))
                     for k in range(1, n_wanted + 1)])

gam_all = load_zeros(NEEDED)

def theta(g, d): return 2.0*np.arctan(g/d) - np.pi
def quart(beta, T, d):
    out=[]
    for sb in (1,-1):
        for sT in (1,-1):
            r = 0.5 + sb*beta + 1j*sT*T
            q = (r-(0.5-d))/(r-(0.5+d))
            out.append((np.log(abs(q)), np.angle(q)))
    return out
def S(ns, th, qs, chunk=2000):
    out=np.empty(len(ns))
    for i in range(0,len(ns),chunk):
        nn=ns[i:i+chunk].astype(float)
        out[i:i+chunk]=2.0*np.sum(1.0-np.cos(nn[:,None]*th[None,:]),axis=1)
        for lg,ph in qs:
            out[i:i+chunk]+=1.0-np.exp(nn*lg)*np.cos(nn*ph)
    return out

def first_neg(th, qs, nmax=4_000_000):
    """bracket by doubling, bisect on the bracket, then scan the last window finely"""
    lo=8
    while lo<nmax:
        hi=min(lo*2,nmax)
        probe=np.unique(np.linspace(lo,hi,300).astype(np.int64))
        if np.any(S(probe,th,qs)<0): break
        lo=hi
    else:
        return None
    # bisect down to a window of <= 4000 integers
    a,b=lo,hi
    while b-a>4000:
        m=(a+b)//2
        probe=np.unique(np.linspace(a,m,300).astype(np.int64))
        if np.any(S(probe,th,qs)<0): b=m
        else: a=m
    fine=np.arange(a,b+1)
    v=S(fine,th,qs)
    k=np.argmax(v<0)
    return int(fine[k]) if v[k]<0 else None

print("LARGER SWEEP: does the detection-index gain stay at T?\n")
print(f"  {'beta':>6} {'T':>6} {'zeros':>6} {'n classical':>13} {'n matched':>11} {'gain':>9} {'gain/T':>8}")
rows=[]
t0=time.time()
for beta in (0.1, 0.05):
    for T in (10.,20.,30.,50.,80.,120.):
        for nz in (2000, 8000):
            g=gam_all[:nz]; R=np.hypot(beta,T)
            nc=first_neg(theta(g,0.5), quart(beta,T,0.5))
            nm=first_neg(theta(g,R),   quart(beta,T,R))
            if nc and nm:
                r=nc/nm/T; rows.append((beta,T,nz,r))
                print(f"  {beta:>6} {T:>6.0f} {nz:>6} {nc:>13} {nm:>11} {nc/nm:>9.2f} {r:>8.3f}",flush=True)
print(f"\n  elapsed {time.time()-t0:.0f}s")
r=np.array([x[3] for x in rows])
print(f"  gain/T over {len(r)} runs: mean {r.mean():.4f}  sd {r.std(ddof=1):.4f}  range [{r.min():.3f}, {r.max():.3f}]")

Ts=np.array([x[1] for x in rows])
sl=np.polyfit(np.log(Ts), r, 1)[0]
print(f"  slope of gain/T against log T: {sl:+.4f}  (a d-dependent correction would be clearly positive)")

FAIL = []
if len(r) != 24:
    FAIL.append(f"the sweep produced {len(r)} runs, not 24")
if abs(r.mean() - 1.0) > 0.02:
    FAIL.append(f"mean gain/T is {r.mean():.4f}, not within 2% of 1")
if sl > 0.02:
    FAIL.append(f"slope against log T is {sl:+.4f}: a drift, not scatter")
print()
if FAIL:
    print("FAILURES:")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)
print("the gain is T: mean within 2% of 1, and no drift with height")
