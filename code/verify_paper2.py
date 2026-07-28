import numpy as np, math
from sympy import factorint, isprime, primerange, totient, divisor_count, mobius

print("="*72); print("T1  remainder formula R(H,j) = 2(H-j+1) mod (2j+1)   [incl. j > H+1]"); print("="*72)
bad=0; tot=0; badneg=0
for H in range(0,600):
    n=2*H+3
    for j in range(1,700):
        d=2*j+1
        R=(2*(H-j+1))%d
        tot+=1
        if R!=n%d:
            bad+=1
            if H-j+1<0: badneg+=1
print(f"cells tested {tot}, mismatches {bad} (of which k<0: {badneg})  [least-nonnegative residue convention]")

print()
print("="*72); print("T2  equivalences  d|n <=> R=0 <=> d|k <=> H = j-1 mod d;  and Z_d offset (d-3)/2"); print("="*72)
bad=0
for H in range(0,3000):
    n=2*H+3
    for j in range(1,80):
        d=2*j+1; k=H-j+1
        a=(n%d==0); b=((2*k)%d==0); c=(k%d==0); e=((H-(j-1))%d==0); f=((H-(d-3)//2)%d==0)
        if len({a,b,c,e,f})!=1: bad+=1
print(f"mismatches among the five conditions: {bad}")

print()
print("="*72); print("T3  Theorem 1 : n=2H+3 composite  <=>  H in Im(Phi),  Phi(j,t)=2jt+j+t-1"); print("="*72)
NMAX=200001
Hmax=(NMAX-3)//2
hit=np.zeros(Hmax+1,dtype=bool)
j=1
while 2*j*1+j+1-1<=Hmax:
    t=1
    while 2*j*t+j+t-1<=Hmax:
        hit[2*j*t+j+t-1]=True; t+=1
    j+=1
bad=0
for H in range(0,Hmax+1):
    n=2*H+3
    if hit[H]==isprime(n): bad+=1
print(f"odd n <= {NMAX}: mismatches = {bad}")

print()
print("="*72); print("T4  Phi IS the Sieve of Sundaram (1934), shifted by one"); print("="*72)
bad=0
for jj in range(1,200):
    for tt in range(1,200):
        if jj+tt+2*jj*tt <= Hmax+1:
            bad += (2*jj*tt+jj+tt-1) != (jj+tt+2*jj*tt)-1
print(f"Phi(j,t) + 1 == j + t + 2jt  : violations {bad}")
print("Sundaram: m removed  <=>  2m+1 composite;  here m = H+1 and 2m+1 = 2H+3 = n.  Identical sieve.")

print()
print("="*72); print("T5  Theorem 2 : #Phi^{-1}(H) = tau(n) - 2   (ordered pairs)"); print("="*72)
cnt=np.zeros(Hmax+1,dtype=np.int64)
j=1
while 2*j+j+1-1<=Hmax:
    t=1
    while 2*j*t+j+t-1<=Hmax:
        cnt[2*j*t+j+t-1]+=1; t+=1
    j+=1
bad=0
for H in range(0,20001):
    if cnt[H]!=divisor_count(2*H+3)-2: bad+=1
print(f"odd n <= {2*20001+3}: mismatches = {bad}")
print("\n  paper's table re-derived (ordered = tau-2, unordered = ceil((tau-2)/2)):")
print(f"  {'n':>5} {'factorization':>12} {'tau':>4} {'ordered':>8} {'unordered':>10}")
for n in [21,25,45,81,105,225,945]:
    f=factorint(n); s="x".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in sorted(f.items()))
    tau=int(divisor_count(n)); od=tau-2; un=(od+1)//2
    print(f"  {n:>5} {s:>12} {tau:>4} {od:>8} {un:>10}")

print()
print("="*72); print("T6  diagonal j=t <=> odd squares ; symmetry ; j<=t <=> d<=sqrt(n)"); print("="*72)
bad=0
for j in range(1,300):
    n=2*(2*j*j+2*j-1)+3
    if n!=(2*j+1)**2: bad+=1
print(f"Phi(j,j) -> n = (2j+1)^2 : violations {bad}")
bad=sum(1 for j in range(1,200) for t in range(1,200) if (2*j*t+j+t-1)!=(2*t*j+t+j-1))
print(f"symmetry Phi(j,t)=Phi(t,j) : violations {bad}")
bad=0
for j in range(1,120):
    for t in range(1,120):
        n=(2*j+1)*(2*t+1)
        if (j<=t) != ((2*j+1)<=math.isqrt(n) or (2*j+1)**2==n): bad+=1
print(f"j<=t  <=>  d<=sqrt(n) : violations {bad}")

print()
print("="*72); print("T7  layers: Z_d = intersection of prime-power layers ; V(n) reconstruction"); print("="*72)
bad=0
for d in range(3,2001,2):
    off=(d-3)//2
    comp=[(p**a) for p,a in factorint(d).items()]
    for H in range(0,4000):
        inZ = ((H-off)%d==0)
        inI = all(((H-(q-3)//2)%q==0) for q in comp)
        if inZ!=inI: bad+=1
print(f"Z_d vs intersection of Z_{{p^a}}, d<=2000 odd, H<4000 : violations {bad}")
bad=0
for n in range(3,200001,2):
    V=factorint(n)
    tau=1; phi=n; om=len(V)
    for p,a in V.items(): tau*=a+1; phi=phi//p*(p-1)
    mu=0 if any(a>=2 for a in V.values()) else (-1)**om
    lam=math.log(list(V)[0]) if len(V)==1 else 0.0
    if (tau!=divisor_count(n)) or (phi!=totient(n)) or (mu!=mobius(n)): bad+=1
print(f"tau, phi, mu rebuilt from V(n), odd n <= 200001 : violations {bad}")

print()
print("="*72); print("T8  the compression counts of section 12"); print("="*72)
N=1000001
odd_cols=len(range(3,N,2)); print(f"odd denominator columns 3..{N-2}          = {odd_cols}")
sq=math.isqrt(N); print(f"sqrt(N) = {sq}")
odd_primes_to_sqrt=len([p for p in primerange(3,sq+1)]); print(f"odd PRIME layers up to sqrt(N)             = {odd_primes_to_sqrt}")
odd_cols_to_sqrt=len(range(3,sq+1,2)); print(f"odd COLUMNS up to sqrt(N)                  = {odd_cols_to_sqrt}")
pp=0
for p in primerange(3,N):
    q=p
    while q<=N: pp+=1; q*=p
print(f"odd prime-power layers up to N             = {pp}")
print()
print(f"paper's claim : {odd_cols}/{odd_primes_to_sqrt} = {odd_cols/odd_primes_to_sqrt:.1f}  <-- compares columns up to N against layers up to sqrt(N)")
print(f"fair, same range (sqrt N) : {odd_cols_to_sqrt}/{odd_primes_to_sqrt} = {odd_cols_to_sqrt/odd_primes_to_sqrt:.3f}")
print(f"fair, same range (N)      : {odd_cols}/{pp} = {odd_cols/pp:.3f}   vs  (log N)/2 = {math.log(N)/2:.3f}")

print()
print("="*72); print("T9  section-11 sample table"); print("="*72)
print(f"  {'n':>7} {'tau':>4} {'phi':>7} {'mu':>3} {'Lambda':>10}   paper")
paper={45:(6,24,0,'0'),81:(5,54,0,'log 3'),105:(8,48,-1,'0'),225:(9,120,0,'0'),945:(16,432,0,'0'),19683:(10,13122,0,'log 3'),99991:(2,99990,-1,'log 99991')}
for n,(t_,f_,m_,l_) in paper.items():
    V=factorint(n); lam = f"log {list(V)[0]}" if len(V)==1 else "0"
    ok = (int(divisor_count(n))==t_) and (int(totient(n))==f_) and (int(mobius(n))==m_) and (lam==l_)
    print(f"  {n:>7} {int(divisor_count(n)):>4} {int(totient(n)):>7} {int(mobius(n)):>3} {lam:>10}   {'OK' if ok else 'MISMATCH'}")
print(f"  is 99991 prime? {isprime(99991)}   factorization: {factorint(99991)}")
