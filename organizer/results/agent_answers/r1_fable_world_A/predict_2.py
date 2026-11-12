#!/usr/bin/env python
"""Stage 2 predictions.  Uses only the equations and constants stated in ANSWER_2.md:
 (1) Newtonian gravitation of the 8 massive bodies S, T, M, P1..P5 (constants.json, unchanged from stage 1);
 (2) every faint body C_i is a massless particle attracted by all 8 massive bodies, plus a radial solar term
        a_extra = + beta_i * G m_S * r / |r|^3        (r = heliocentric position of C_i)
 (3) beta_i fades near planets:  d beta_i/dt = -A * beta_i * sum_p exp(-(d_p/R)^2),  p in {T, P1..P5}
 Constants of (2),(3) and the initial conditions of the faint bodies are in constants_2.json.
numpy only, a single process.  Massive bodies: 6th-order symplectic, step 0.25 d (as in stage 1).
Faint bodies: classical RK4, step 0.5 d, using the massive-body positions on the 0.25-d grid."""
import numpy as np, csv, json, os, sys
here = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(here, 'constants.json')))
C2 = json.load(open(os.path.join(here, 'constants_2.json')))
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
rows = list(csv.DictReader(open(os.path.join(here, 'predictions_2.csv'))))
last = max(int(float(r['day'])) for r in rows)
CHECK = '--check' in sys.argv
if CHECK: last = max(last, 36525)

# ---- massive bodies on a quarter-day grid -------------------------------------------------------------
N = last * 4
X = np.empty((N + 1, 8, 3)); X[0] = x
h = 0.25; a = acc(x)
for k in range(1, N + 1):
    for w in W:
        v += 0.5 * w * h * a; x += w * h * v; a = acc(x); v += 0.5 * w * h * a
    X[k] = x

def direction(p, day):
    u = p - X[day * 4, 1]; u = u / np.linalg.norm(u)
    return np.degrees(np.arctan2(u[1], u[0])) % 360, np.degrees(np.arcsin(u[2]))

# ---- faint bodies ------------------------------------------------------------------------------------
A_LAW, R_LAW = C2['A'], C2['R']
PL = [1, 3, 4, 5, 6, 7]                       # T, P1..P5

def f(xc, vc, s, k4):
    Xm = X[k4]
    d = Xm - xc[None, :]; r2 = (d * d).sum(-1)
    a = ((gm * r2 ** -1.5)[:, None] * d).sum(0)            # Newton, all 8 massive bodies
    r = xc - Xm[0]; rn = np.linalg.norm(r)
    a = a + s * GMS * r / rn ** 3                          # radial solar term
    ds = -A_LAW * s * np.exp(-r2[PL] / R_LAW ** 2).sum()   # fading of beta near planets
    return a, ds

def run(body, want_days):
    """integrate one faint body from its epoch; returns {day: position} for the wanted days (either side of epoch)"""
    c = C2['bodies'][body]; day0 = int(c['epoch_day']); out = {}
    for sgn in (1, -1):
        sel = sorted([d for d in want_days if (d >= day0 if sgn > 0 else d < day0)], reverse=(sgn < 0))
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

res = {}
need = {}
for r in rows: need.setdefault(r['body'], set()).add(int(float(r['day'])))
for body, days in need.items():
    if body in names:
        i = names.index(body)
        for d in days: res[(d, body)] = direction(X[d * 4, i], d)
    else:
        pos = run(body, days)
        for d in days: res[(d, body)] = direction(pos[d], d)
with open(os.path.join(here, 'predictions_2.csv'), 'w', newline='') as fo:
    wr = csv.writer(fo); wr.writerow(['day', 'body', 'lon', 'lat'])
    for r in rows:
        lo, la = res[(int(float(r['day'])), r['body'])]
        wr.writerow([r['day'], r['body'], '%.5f' % lo, '%.5f' % la])
print('wrote', len(rows), 'predictions')

if CHECK:      # compare the model with the faint-body observations inside each body's stated validity window
    obs = {}
    for r in csv.DictReader(open(os.path.join(here, 'observations_2.csv'))):
        obs.setdefault(r['body'], []).append((int(float(r['day'])), float(r['lon']), float(r['lat'])))
    for body in sorted(obs, key=lambda s: int(s[1:])):
        c = C2['bodies'][body]; lo, hi = c['valid_from_day'], c['valid_to_day']
        o = [t for t in obs[body] if lo <= t[0] <= hi]
        pos = run(body, set(t[0] for t in o)); sq = []
        for d, lon, lat in o:
            u = pos[d] - X[d * 4, 1]; u /= np.linalg.norm(u)
            w = np.array([np.cos(np.radians(lat)) * np.cos(np.radians(lon)), np.cos(np.radians(lat)) * np.sin(np.radians(lon)), np.sin(np.radians(lat))])
            sq.append(((u - w) ** 2).sum())
        print(body, 'days %d-%d  N=%d  rms per axis (arcmin): %.3f' % (lo, hi, len(o), np.sqrt(np.mean(sq) / 2) * 3437.75))
