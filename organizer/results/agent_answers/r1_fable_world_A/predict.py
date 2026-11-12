#!/usr/bin/env python
"""Predict directions of S, P1-P5, M as seen from Tellus by integrating Newton's law of gravitation
for 8 point masses (S, Tellus, M, P1..P5) from the initial conditions stated in ANSWER_1.md.
Uses only numpy.  One process, fixed-step 6th-order symplectic integrator (step 0.25 day)."""
import numpy as np, csv, json, os, sys
here=os.path.dirname(os.path.abspath(__file__))
C=json.load(open(os.path.join(here,'constants.json')))   # the constants of ANSWER_1.md, section 2
GMS=C['GM_S']
names=['S','T','M','P1','P2','P3','P4','P5']
gm=np.array([GMS*C['mass_ratio'][n] for n in names])
x=np.array([C['pos'][n] for n in names],float)   # heliocentric positions at day 0 (length unit L)
v=np.array([C['vel'][n] for n in names],float)   # heliocentric velocities at day 0 (L/day)
# move to centre-of-mass frame
x-= (gm[:,None]*x).sum(0)/gm.sum(); v-=(gm[:,None]*v).sum(0)/gm.sum()
def acc(x):
    d=x[None,:,:]-x[:,None,:]; r2=(d*d).sum(-1); np.fill_diagonal(r2,1.0)
    w=gm[None,:]*r2**-1.5; np.fill_diagonal(w,0.0)
    return (w[...,None]*d).sum(1)
_w=[-1.17767998417887,0.235573213359357,0.784513610477560]
W=[_w[2],_w[1],_w[0],1-2*sum(_w),_w[0],_w[1],_w[2]]
rows=list(csv.DictReader(open(os.path.join(here,'predictions.csv'))))
want={}
for r in rows: want.setdefault(int(float(r['day'])),[])
last=max(want); h=0.25; nsub=4
res={}
CHECK='--check' in sys.argv   # optional: also compare the model with observations.csv
obs={}
if CHECK:
    for r in csv.DictReader(open(os.path.join(here,'observations.csv'))): obs.setdefault(int(float(r['day'])),[]).append((r['body'],float(r['lon']),float(r['lat'])))
sq={n:[] for n in names}
a=acc(x)
for day in range(0,last+1):
    if day in want:
        dd=x-x[1]; 
        for i,n in enumerate(names):
            if n=='T': continue
            u=dd[i]/np.linalg.norm(dd[i])
            res[(day,n)]=(np.degrees(np.arctan2(u[1],u[0]))%360, np.degrees(np.arcsin(u[2])))
    if CHECK and day in obs:
        dd=x-x[1]
        for n,lo,la in obs[day]:
            u=dd[names.index(n)]; u=u/np.linalg.norm(u)
            o=np.array([np.cos(np.radians(la))*np.cos(np.radians(lo)),np.cos(np.radians(la))*np.sin(np.radians(lo)),np.sin(np.radians(la))])
            sq[n].append(((u-o)**2).sum())
    for s in range(nsub):
        for w in W:
            v+=0.5*w*h*a; x+=w*h*v; a=acc(x); v+=0.5*w*h*a
with open(os.path.join(here,'predictions.csv'),'w',newline='') as f:
    wr=csv.writer(f); wr.writerow(['day','body','lon','lat'])
    for r in rows:
        lo,la=res[(int(float(r['day'])),r['body'])]
        wr.writerow([r['day'],r['body'],'%.5f'%lo,'%.5f'%la])
print('wrote',len(rows),'predictions')
if CHECK:
    for n in names:
        if sq[n]: print(n,'rms residual per axis vs observations (arcmin): %.3f'%(np.sqrt(np.mean(sq[n])/2)*3437.75))
