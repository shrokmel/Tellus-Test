#!/usr/bin/env python
"""Stage 3 predictions for the new faint body N.  Uses only the equations and constants of ANSWER_3.md:
 (1) Newtonian gravitation of the 8 massive bodies S, T, M, P1..P5 (constants.json, unchanged from stage 1);
 (2) N is a massless particle whose acceleration is the total Newtonian field of the 8 massive bodies
     multiplied by (1 - beta):        d2x/dt2 = (1 - beta) * sum_j G m_j (r_j - x)/|r_j - x|^3
 (3) beta fades near the non-solar massive bodies:
        d beta/dt = -A * beta * sum_p exp(-(d_p/R)^2),   p in {T, M, P1, P2, P3, P4, P5}
 Constants A, R and the initial conditions of N are in constants_3.json.
numpy only, a single process.  Massive bodies: 6th-order symplectic, step 0.25 d from day 0 (as in stage 1).
N: classical RK4, step 0.5 d, using the massive-body positions on the 0.25-d grid."""
import numpy as np, csv, json, os, sys
here = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(here, 'constants.json')))
C3 = json.load(open(os.path.join(here, 'constants_3.json')))
GMS = C['GM_S']
names = ['S', 'T', 'M', 'P1', 'P2', 'P3', 'P4', 'P5']
gm = np.array([GMS * C['mass_ratio'][n] for n in names])
x = np.array([C['pos'][n] for n in names], float)
v = np.array([C['vel'][n] for n in names], float)
x -= (gm[:, None] * x).sum(0) / gm.sum()
v -= (gm[:, None] * v).sum(0) / gm.sum()          # centre-of-mass frame

def acc(x):
    d = x[None, :, :] - x[:, None, :]
    r2 = (d * d).sum(-1); np.fill_diagonal(r2, 1.0)
    w = gm[None, :] * r2 ** -1.5; np.fill_diagonal(w, 0.0)
    return (w[..., None] * d).sum(1)

_w = [-1.17767998417887, 0.235573213359357, 0.784513610477560]
W = [_w[2], _w[1], _w[0], 1 - 2 * sum(_w), _w[0], _w[1], _w[2]]
rows = list(csv.DictReader(open(os.path.join(here, 'predictions_3.csv'))))
want = sorted(set(int(float(r['day'])) for r in rows))
CHECK = '--check' in sys.argv
obs = []
if CHECK:
    obs = [(int(float(r['day'])), float(r['lon']), float(r['lat'])) for r in csv.DictReader(open(os.path.join(here, 'observations_3.csv')))]
last = max(want)

# ---- massive bodies on a quarter-day grid -------------------------------------------------------------
N4 = last * 4
X = np.empty((N4 + 1, 8, 3)); X[0] = x
h = 0.25; a = acc(x)
for k in range(1, N4 + 1):
    for w in W:
        v += 0.5 * w * h * a; x += w * h * v; a = acc(x); v += 0.5 * w * h * a
    X[k] = x

# ---- the faint body N ---------------------------------------------------------------------------------
A_LAW, R_LAW = C3['A'], C3['R']
PL = [1, 2, 3, 4, 5, 6, 7]                    # T, M, P1..P5

def f(xc, vc, s, k4):
    Xm = X[k4]
    d = Xm - xc[None, :]; r2 = (d * d).sum(-1)
    g = ((gm * r2 ** -1.5)[:, None] * d).sum(0)            # Newtonian field of all 8 massive bodies
    ds = -A_LAW * s * np.exp(-r2[PL] / R_LAW ** 2).sum()   # fading of beta
    return (1.0 - s) * g, ds

def run(days):
    c = C3['N']; day0 = int(c['epoch_day']); out = {}
    for sgn in (1, -1):
        sel = sorted([d for d in days if (d >= day0 if sgn > 0 else d < day0)], reverse=(sgn < 0))
        if not sel: continue
        xc = np.array(c['pos'], float); vc = np.array(c['vel'], float); s = float(c['beta'])
        k = day0 * 4; hh = 0.5 * sgn; day = day0; j = 0
        while True:
            while j < len(sel) and sel[j] == day:
                out[day] = xc.copy(); j += 1
            if j >= len(sel): break
            for _ in range(2):
                k1, s1 = f(xc, vc, s, k); xa = xc + 0.5 * hh * vc; va = vc + 0.5 * hh * k1; sa = s + 0.5 * hh * s1
                k2, s2 = f(xa, va, sa, k + sgn); xb = xc + 0.5 * hh * va; vb = vc + 0.5 * hh * k2; sb = s + 0.5 * hh * s2
                k3, s3 = f(xb, vb, sb, k + sgn); xd = xc + hh * vb; vd = vc + hh * k3; sd = s + hh * s3
                k4, s4 = f(xd, vd, sd, k + 2 * sgn)
                xc = xc + hh * (vc + 2 * va + 2 * vb + vd) / 6; vc = vc + hh * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                s = s + hh * (s1 + 2 * s2 + 2 * s3 + s4) / 6
                k += 2 * sgn
            day += sgn
    return out

def direction(p, day):
    u = p - X[day * 4, 1]; u = u / np.linalg.norm(u)
    return np.degrees(np.arctan2(u[1], u[0])) % 360, np.degrees(np.arcsin(u[2]))

pos = run(set(want) | set(o[0] for o in obs))
with open(os.path.join(here, 'predictions_3.csv'), 'w', newline='') as fo:
    wr = csv.writer(fo); wr.writerow(['day', 'body', 'lon', 'lat'])
    for r in rows:
        d = int(float(r['day'])); lo, la = direction(pos[d], d)
        wr.writerow([r['day'], r['body'], '%.5f' % lo, '%.5f' % la])
print('wrote', len(rows), 'predictions')
if CHECK:
    sq = []
    for d, lon, lat in obs:
        u = pos[d] - X[d * 4, 1]; u /= np.linalg.norm(u)
        w = np.array([np.cos(np.radians(lat)) * np.cos(np.radians(lon)), np.cos(np.radians(lat)) * np.sin(np.radians(lon)), np.sin(np.radians(lat))])
        sq.append(((u - w) ** 2).sum())
    print('N: %d sightings, rms per axis (arcmin): %.3f' % (len(obs), np.sqrt(np.mean(sq) / 2) * 3437.75))
