#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Second-order conditional singular series  (Part IV)
================================================================================
First order (Part III): S_{2,omega}(d) ~ C2(d) * L_omega(d), where a prime q|N
locks the gap d iff q | (6d+-1).
Second order (this script): the spatial-compression penalty. For a prime q that
does NOT lock d, conditioning on q ! N (i.e. N != 0 mod q) removes the class
N=0 and crowds N into the remaining q-1 classes, where the same k_q fatal
residues are proportionally denser. The per-q penalty is

    pen_q(d) = [ (q-1-k_q)/(q-1) ] / [ (q-k_q)/q ]  < 1 .

The second-order weight of a gap d at a real centre N with factor set Q is
    w2(d,N) = C2(d) * 1{no q in Q locks d} * prod_{q !in Q, q ! lock d} pen_q(d).
We evaluate this over the REAL twin centres of S_K (each centre's true factor
set; no uniform-combination assumption), aggregate by omega, and compare to the
observed preference and to first order. Emits so_plot_S{K}.csv.

IMPORTANT (single-centre limitation): this conditions only on the LEFT centre N
and treats N+d as independent given N. The two centres are in fact correlated
mod q. This is why the long gap 210 is closed (its lock primes 11,19 rarely
divide both centres) while the short gap 42 keeps a residual (its factor 7
shields both centres when 7|N -- a two-centre effect omitted here). See the paper.

USAGE
  python second_order.py            # default S10 (~15 min)
  MAXK=9 python second_order.py     # S9 (faster, for validation)
Requires: numpy.

NOTE: r(d|omega) is normalised model-internally (each curve divided by its own
omega-merged mean). This differs from Part III's global baseline by a constant
scale per gap; it does not change the gradient in omega or the residual sizes.
================================================================================
"""
# Correct second-order test using REAL centre factors (not Gemini's enumerated
# combinations). Question: does the "q-nmid-N space-compression" penalty push the
# first-order prediction toward the observed value? Let data decide.
#
# Per real twin centre N (with its real small-prime factor set Q):
#   first-order weight of gap d: w1 = C2(d) * 1{no q in Q locks d}
#       (lockdown = q|(6d±1); this is Paper 3's first order)
#   second-order: for each small prime q that (a) does NOT divide N and (b) does
#       NOT lock d, multiply by the space-compression factor
#           pen_q(d) = [ (q-1-k_q)/(q-1) ] / [ (q-k_q)/q ]
#       where k_q = number of fatal residues of d mod q among {±6^{-1}, -d±6^{-1}}
#       (k_q counts classes !=0 since q does not lock d). This is <1 (a penalty):
#       knowing N!=0 mod q crowds N into the remaining q-1 classes.
#   w2 = w1 * prod_{q in pool, q!|N, q!lock d} pen_q(d)
# Aggregate per omega; normalise like Paper 3 (divide by the omega-merged mean);
# compare obs / first-order / second-order at 6dN = 42, 210.
import numpy as np, math, os
from collections import defaultdict
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(10**MAXK))+1; BP=primes_upto(PB)
GMAX=60
POOL=[5,7,11,13,17,19,23,29,31,37,41,43,47]

def k_fatal(d,q):
    inv=pow(6,-1,q); s={inv%q,(-inv)%q,(-d+inv)%q,(-d-inv)%q}
    return s
def locks(d,q): return (6*d-1)%q==0 or (6*d+1)%q==0
# precompute per (d,q): k (excluding 0 if not lock), penalty factor
PEN={}; KQ={}; LOCK={}
for d in range(1,GMAX+1):
    LOCK[d]=set(q for q in POOL if locks(d,q))
    for q in POOL:
        s=k_fatal(d,q); k=len(s)
        KQ[(d,q)]=k
        if not locks(d,q):
            # penalty = [(q-1-k)/(q-1)] / [(q-k)/q]
            num=(q-1-k)/(q-1); den=(q-k)/q
            PEN[(d,q)]=num/den if den>0 else 0.0
        else:
            PEN[(d,q)]=None  # locked: handled by first order (weight 0)

# C2(d)
QP=[q for q in primes_upto(200000) if q>3]
def C2(d):
    p=1.0
    for q in QP:
        nu=len({(-1)%q,1%q,(6*d-1)%q,(6*d+1)%q})
        if nu==q: return 0.0
        p*=(1.0-nu/q)/(1.0-2.0/q)**2
    return p
Cg={d:C2(d) for d in range(1,GMAX+1)}
poolbit={q:1<<i for i,q in enumerate(POOL)}
lockbit={d:sum(poolbit[q] for q in LOCK[d]) for d in range(1,GMAX+1)}
# second-order multiplier per (d, centre-mask): product of PEN over q NOT in mask and NOT locking d
# precompute, for each d, the list of (q, pen) with q not locking d
penlist={d:[(q,PEN[(d,q)]) for q in POOL if q not in LOCK[d]] for d in range(1,GMAX+1)}

twN=[]; twOm=[]; twMask=[]
n=LO; t0=__import__('time').time()
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16); mask=np.zeros(sz,np.int32)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3:
            ob[idx]+=1
            if p in poolbit: mask[idx]|=poolbit[p]
    ob[rem>1]+=1
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    Narr=np.arange(n,nh,dtype=np.int64)
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    pos=np.nonzero(tw)[0]
    twN.append(Narr[pos]); twOm.append(ob[pos]); twMask.append(mask[pos])
    n=nh
twN=np.concatenate(twN); twOm=np.concatenate(twOm); twMask=np.concatenate(twMask)
print(f"S{MAXK} twins {len(twN):,}; {__import__('time').time()-t0:.0f}s")

g=np.diff(twN); omL=twOm[:-1]; maskL=twMask[:-1]
keep=(g>=1)&(g<=GMAX); g=g[keep]; omL=omL[keep]
obs=defaultdict(lambda: np.zeros(GMAX+1))
for i in range(len(g)): obs[int(omL[i])][int(g[i])]+=1
overall=np.zeros(GMAX+1)
for om in obs: overall+=obs[om]
base_obs=overall/overall[1:GMAX+1].sum()

omegas=[om for om in sorted(obs) if obs[om][1:GMAX+1].sum()>=10000]

# For each (omega,d): first-order pred = C2(d)*mean(allow);
# second-order pred = C2(d)*mean(allow * prod_{q not in mask, q not lock d} pen)
def predict(d):
    lb=lockbit[d]; pens=penlist[d]
    p1=np.zeros(7); p2=np.zeros(7)
    for om in omegas:
        sel=twMask[twOm==om]
        allow=(sel&lb)==0
        # second-order multiplier per centre: product over q not in mask and not locking d
        mult=np.ones(len(sel),dtype=np.float64)
        for q,pen in pens:
            qb=poolbit[q]
            has=(sel&qb)>0
            # apply pen where q does NOT divide N (i.e. ~has)
            mult*= np.where(has,1.0,pen)
        p1[om]=Cg[d]*allow.mean()
        p2[om]=Cg[d]*(allow*mult).mean()
    return p1,p2
# baselines over all centres (omega-merged) for normalisation
def predict_base(d):
    lb=lockbit[d]; pens=penlist[d]
    allow=(twMask&lb)==0
    mult=np.ones(len(twMask),dtype=np.float64)
    for q,pen in pens:
        qb=poolbit[q]; has=(twMask&qb)>0
        mult*=np.where(has,1.0,pen)
    return Cg[d]*allow.mean(), Cg[d]*(allow*mult).mean()

# normalise: need full-distribution normaliser. Approximate as in Paper 3 using
# the same per-d quantities over the merged set is not a full sum; instead we
# normalise each model by its omega-merged value at that d so r=1 at "average".
print(f"\nSecond-order test, S{MAXK}: obs vs first-order vs second-order (normalised)")
for d in [7,35]:
    p1,p2=predict(d); b1,b2=predict_base(d)
    print(f"\n6dN={6*d}:")
    print(f"  {'omega':>5}{'obs':>9}{'1st-order':>11}{'2nd-order':>11}")
    for om in omegas:
        tot=obs[om][1:GMAX+1].sum()
        obs_r=(obs[om][d]/tot)/base_obs[d] if base_obs[d]>0 else 0
        r1=p1[om]/b1 if b1>0 else 0
        r2=p2[om]/b2 if b2>0 else 0
        print(f"  {om:>5}{obs_r:>9.3f}{r1:>11.3f}{r2:>11.3f}")


# ---- emit plotting CSV (so_plot_S{K}.csv) ----
import csv as _csv
with open(f'so_plot_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['6dN','omega','obs','first_order','second_order'])
    for d in [7,35]:
        p1,p2=predict(d); b1,b2=predict_base(d)
        for om in omegas:
            tot=obs[om][1:GMAX+1].sum()
            obs_r=(obs[om][d]/tot)/base_obs[d] if base_obs[d]>0 else 0
            r1=p1[om]/b1 if b1>0 else 0
            r2=p2[om]/b2 if b2>0 else 0
            _w.writerow([6*d,om,f'{obs_r:.4f}',f'{r1:.4f}',f'{r2:.4f}'])
print(f"\n[ok] wrote so_plot_S{MAXK}.csv")
