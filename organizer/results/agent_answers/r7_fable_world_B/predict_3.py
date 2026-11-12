#!/usr/bin/env python
"""predict_3.py -- direction of the new faint body N as seen from Tellus,
from the laws stated in ANSWER_3.md.

Law 1 (massive bodies S, P1, Tellus, M, P2..P5; unchanged from stage 1):
        d^2 r_i / dt^2 = sum_{j != i} GM_j (r_j - r_i) / |r_j - r_i|^3
Law 3 (unchanged): small secular term on the M-Tellus relative vector rho:  + GAMMA * d rho / dt.
Law 1' (faint bodies; new): a faint body b feels the same field of all eight massive bodies, multiplied by its own
     constant factor (1 + ETA_b):
        d^2 r_b / dt^2 = (1 + ETA_b) * sum_j GM_j (r_j - r_b) / |r_j - r_b|^3
     Faint bodies exert no measurable force on anything.
Law 2 (unchanged): observed direction of a body = direction of r_body - r_Tellus:
        lon = atan2(y, x),  lat = asin(z / |r|).
Constants: GM_S = 2.959e-4 L^3/day^2 (defines the length unit L), mass ratios GM_REL, GAMMA, ETA, and positions /
     velocities at epoch day T0 = 18262 relative to S (M relative to Tellus), N's state at day TN relative to S; all listed below and in ANSWER_3.md.
Numerics: fixed step 0.5 day, 8th-order Runge-Kutta (Dormand-Prince coefficients).
Usage: python predict_3.py      (reads predictions_3.csv, fills lon/lat, rewrites it)
"""
import csv
import numpy as np
from scipy.integrate._ivp import dop853_coefficients as _dc

GM_S = 2.959e-4
T0 = 18262
ORDER = ['S', 'P1', 'T', 'M', 'P2', 'P3', 'P4', 'P5']
GM_REL = {   # GM_body / GM_S
    'P1': 3.2741614495e-07,
    'T':  3.0083240800e-06,
    'M':  3.4441577080e-08,
    'P2': 6.0736690970e-07,
    'P3': 1.5011246810e-04,
    'P4': 1.0000211634e-03,
    'P5': 2.4961482571e-04,
}
GAMMA = 4.171622e-13   # 1/day : secular along-track acceleration coefficient of M
STATE = {   # at day T0: x,y,z [L], vx,vy,vz [L/day]; relative to S (M: relative to Tellus)
    'P1': [9.1064625793917428e-02, 4.0711307714642903e-01, -1.6201496553362094e-01, -2.5512649259988664e-02, 3.4643403204493821e-03, -1.3791398583217075e-03],
    'T':  [-4.0776510470894889e-01, -8.6627310270167379e-01, 3.4475232702735803e-01, 1.5420533311223866e-02, -6.4330704887881835e-03, 2.5604599109755742e-03],
    'M':  [2.4967154385874607e-03, -1.7254207205097418e-04, -1.3165789847240892e-04, 4.8481672576106441e-05, 5.6945717821121943e-04, -1.9570359391707365e-04],
    'P2': [1.2113147734079719e+00, -1.0129050205619601e+00, 4.0335357397262367e-01, 9.4886205124626154e-03, 9.2461983448505064e-03, -3.6797702263660116e-03],
    'P3': [6.5768643691123563e-02, -3.1257926417280570e+00, 1.2873631044800671e+00, 9.2255039138841397e-03, 3.5653317767069101e-04, -3.7583546803934739e-04],
    'P4': [1.3175273244394403e+00, 5.7847146868591590e+00, -2.4224862181861733e+00, -6.5860226060181007e-03, 1.0603492951939233e-03, -4.0658173109741046e-04],
    'P5': [-1.0981162632906226e+01, -3.9631404639702456e-01, 5.4566383420789377e-01, 2.2695437333347088e-04, -4.9360099621809659e-03, 1.8977742448029952e-03],
}
TN = 36358      # epoch (day) of the state of N
FAINT = {   # name: (ETA, [x,y,z (L), vx,vy,vz (L/day)] relative to S at day TN)
    'N': (3.2201120915e-01, [7.9610948908865509e-01, 2.1896192239166838e+00, 1.1137234590298801e+00, -1.5357525246208101e-02, 1.8642546631355878e-03, 5.0080826224349973e-03]),
}

A8 = _dc.A[:12, :12]; B8 = _dc.B; SUB = 2
IT = ORDER.index("T"); IM = ORDER.index("M")
FNAMES = list(FAINT)
ETA = np.array([FAINT[b][0] for b in FNAMES])

def acc_massive(x, v, gm):
    # slot IM holds M relative to Tellus (same equations, less round-off): rebuild the true position first
    r = x[IM].copy(); xt = x.copy(); xt[IM] = x[IT] + r
    dx = xt[None, :, :] - xt[:, None, :]          # dx[i,j] = r_j - r_i
    dx[IT, IM] = r; dx[IM, IT] = -r
    r2 = (dx * dx).sum(-1)
    np.fill_diagonal(r2, 1.0)
    w = gm[None, :] / (r2 * np.sqrt(r2))
    np.fill_diagonal(w, 0.0)
    a = (w[..., None] * dx).sum(1)
    a[IM] -= a[IT]                                # relative acceleration of M with respect to Tellus
    a[IM] += GAMMA * v[IM]                        # Law 3
    a[IT] -= GAMMA * gm[IM] / (gm[IT] + gm[IM]) * v[IM]
    return a, xt

def acc_faint(xp, xt, gm):
    dx = xt[None, :, :] - xp[:, None, :]          # (n_faint, 8, 3)
    r2 = (dx * dx).sum(-1)
    g = ((gm[None, :] / (r2 * np.sqrt(r2)))[..., None] * dx).sum(1)
    return (1.0 + ETA)[:, None] * g               # Law 1'

def step(x, v, xp, vp, gm, h):
    kx = []; kv = []; kxp = []; kvp = []
    for s in range(12):
        xs = x.copy(); vs = v.copy(); ps = xp.copy(); qs = vp.copy()
        for j in range(s):
            c = A8[s, j]
            if c != 0.0:
                xs += h * c * kx[j]; vs += h * c * kv[j]; ps += h * c * kxp[j]; qs += h * c * kvp[j]
        a, xt = acc_massive(xs, vs, gm)
        kx.append(vs); kv.append(a); kxp.append(qs); kvp.append(acc_faint(ps, xt, gm))
    xn = x.copy(); vn = v.copy(); pn = xp.copy(); qn = vp.copy()
    for s in range(12):
        if B8[s] != 0.0:
            xn += h * B8[s] * kx[s]; vn += h * B8[s] * kv[s]; pn += h * B8[s] * kxp[s]; qn += h * B8[s] * kvp[s]
    return xn, vn, pn, qn

def initial():
    n = len(ORDER)
    x = np.zeros((n, 3)); v = np.zeros((n, 3)); gm = np.zeros(n)
    gm[0] = GM_S
    for i, b in enumerate(ORDER[1:], 1):
        s = np.array(STATE[b]); x[i] = s[:3]; v[i] = s[3:]; gm[i] = GM_REL[b] * GM_S
    x[IM] += x[IT]; v[IM] += v[IT]                 # M is given relative to Tellus
    xc = (gm[:, None] * x).sum(0) / gm.sum(); vc = (gm[:, None] * v).sum(0) / gm.sum()
    x = x - xc; v = v - vc                         # barycentric (inertial) frame
    s = np.array(STATE["M"]); x[IM] = s[:3]; v[IM] = s[3:]   # keep M as Tellus-relative vector
    return x, v, gm

def lonlat(g):
    return np.degrees(np.arctan2(g[1], g[0])) % 360.0, np.degrees(np.arcsin(g[2] / np.linalg.norm(g)))

def ephemeris(days):
    """dict day -> {'N': (lon, lat)} for integer days >= TN. The massive bodies are integrated from T0;
    N is inserted at day TN with its stated position and velocity relative to S."""
    days = sorted(set(int(d) for d in days)); out = {}; want = set(days)
    assert days[0] >= TN
    x, v, gm = initial(); xp = np.zeros((0, 3)); vp = np.zeros((0, 3)); h = 1.0 / SUB
    for d in range(T0, days[-1] + 1):
        if d == TN:
            f = np.array([FAINT[b][1] for b in FNAMES])
            xp = f[:, :3] + x[0]; vp = f[:, 3:] + v[0]
        if d in want:
            out[d] = {b: lonlat(xp[i] - x[IT]) for i, b in enumerate(FNAMES)}
        if d < days[-1]:
            for _ in range(SUB):
                x, v, xp, vp = step(x, v, xp, vp, gm, h)
    return out

def main(path='predictions_3.csv'):
    with open(path) as f:
        rows = list(csv.reader(f))
    header, rows = rows[0], rows[1:]
    eph = ephemeris(int(float(r[0])) for r in rows)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header)
        for r in rows:
            lon, lat = eph[int(float(r[0]))][r[1]]
            w.writerow([r[0], r[1], '%.5f' % lon, '%.5f' % lat])

if __name__ == '__main__':
    main()
