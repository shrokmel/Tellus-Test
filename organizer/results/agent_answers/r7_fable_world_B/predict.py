#!/usr/bin/env python
"""predict.py -- directions of S, P1..P5, M as seen from Tellus, from the laws stated in ANSWER_1.md.

Law: every body i (S, P1, Tellus, M, P2, P3, P4, P5) is a point mass obeying Newton's law of gravitation
        d^2 r_i / dt^2 = sum_{j != i} GM_j (r_j - r_i) / |r_j - r_i|^3
     in the fixed star frame of the data, plus (Law 3) one small secular term on the M-Tellus relative motion rho:
        d^2 rho / dt^2 = [gravitation] + GAMMA * d rho / dt .
     The observed direction of body b is the direction of r_b - r_Tellus:
        lon = atan2(y, x),  lat = asin(z / |r|).
Constants: GM_S = 2.959e-4 L^3/day^2 (this defines the length unit L, ~ the mean S-Tellus distance),
     the mass ratios GM_REL and the positions/velocities STATE at epoch day T0 = 18262 (all given below and in ANSWER_1.md).
The equations are integrated numerically (fixed step 0.5 day, 8th-order Runge-Kutta, Dormand-Prince coefficients).
Usage: python predict.py      (reads predictions.csv, fills lon/lat, rewrites it)
"""
import csv
import numpy as np
from scipy.integrate._ivp import dop853_coefficients as _dc

GM_S = 2.959e-4
T0 = 18262
ORDER = ['S', 'P1', 'T', 'M', 'P2', 'P3', 'P4', 'P5']
GM_REL = {   # GM_body / GM_S
    'P1': 2.7278790470e-07,
    'T':  2.9752939618e-06,
    'M':  3.4375049071e-08,
    'P2': 6.6883910260e-07,
    'P3': 1.5025853709e-04,
    'P4': 9.9986185473e-04,
    'P5': 2.5002068485e-04,
}
GAMMA = 4.178126e-13   # 1/day : secular along-track acceleration coefficient of M
STATE = {   # at day T0: x,y,z [L], vx,vy,vz [L/day]; relative to S (M: relative to Tellus)
    'P1': [9.1064796662090955e-02, 4.0711317232946975e-01, -1.6201501250526185e-01, -2.5512639971106227e-02, 3.4643427437372246e-03, -1.3791392777226775e-03],
    'T':  [-4.0776511809444843e-01, -8.6627328675652981e-01, 3.4475248009507425e-01, 1.5420530812504218e-02, -6.4330980864530894e-03, 2.5604692349140575e-03],
    'M':  [2.4876285751986205e-03, -1.7192134000170546e-04, -1.3117618115432645e-04, 4.8306378175254615e-05, 5.6738506415462953e-04, -1.9499165628962406e-04],
    'P2': [1.2113154114620732e+00, -1.0129040011801935e+00, 4.0335319792128582e-01, 9.4886175917886645e-03, 9.2462030680732707e-03, -3.6797723924188413e-03],
    'P3': [6.5767934024599795e-02, -3.1257940334883516e+00, 1.2873637087575227e+00, 9.2255020854899832e-03, 3.5653305429762360e-04, -3.7583537456235399e-04],
    'P4': [1.3175271434769387e+00, 5.7847171102032453e+00, -2.4224873113543706e+00, -6.5860153798669718e-03, 1.0603505997303056e-03, -4.0658244403189080e-04],
    'P5': [-1.0981163772060590e+01, -3.9631145766255582e-01, 5.4566285928723568e-01, 2.2695668132832275e-04, -4.9360101399250249e-03, 1.8977743474789051e-03],
}

A8 = _dc.A[:12, :12]; B8 = _dc.B; SUB = 2
IT = ORDER.index("T"); IM = ORDER.index("M")

def acc(x, gm):
    # slot IM holds M relative to Tellus (same equations, less round-off): rebuild the true position first
    r = x[IM].copy(); x = x.copy(); x[IM] = x[IT] + r
    dx = x[None, :, :] - x[:, None, :]            # dx[i,j] = r_j - r_i
    dx[IT, IM] = r; dx[IM, IT] = -r
    r2 = (dx * dx).sum(-1)
    np.fill_diagonal(r2, 1.0)
    w = gm[None, :] / (r2 * np.sqrt(r2))
    np.fill_diagonal(w, 0.0)
    a = (w[..., None] * dx).sum(1)
    a[IM] -= a[IT]                                # relative acceleration of M with respect to Tellus
    return a

def step(x, v, gm, h):
    kx = []; kv = []
    for s in range(12):
        xs = x.copy(); vs = v.copy()
        for j in range(s):
            c = A8[s, j]
            if c != 0.0:
                xs += h * c * kx[j]; vs += h * c * kv[j]
        a = acc(xs, gm)
        a[IM] += GAMMA * vs[IM]                        # secular along-track term on M (relative to Tellus)
        a[IT] -= GAMMA * gm[IM] / (gm[IT] + gm[IM]) * vs[IM]   # reaction on Tellus: T-M barycentre unaffected
        kx.append(vs); kv.append(a)
    xn = x.copy(); vn = v.copy()
    for s in range(12):
        if B8[s] != 0.0:
            xn += h * B8[s] * kx[s]; vn += h * B8[s] * kv[s]
    return xn, vn

def initial():
    n = len(ORDER)
    x = np.zeros((n, 3)); v = np.zeros((n, 3)); gm = np.zeros(n)
    gm[0] = GM_S
    for i, b in enumerate(ORDER[1:], 1):
        s = np.array(STATE[b]); x[i] = s[:3]; v[i] = s[3:]; gm[i] = GM_REL[b] * GM_S
    iT = ORDER.index('T'); iM = ORDER.index('M')
    x[iM] += x[iT]; v[iM] += v[iT]                 # M is given relative to Tellus
    xc = (gm[:, None] * x).sum(0) / gm.sum(); vc = (gm[:, None] * v).sum(0) / gm.sum()
    x = x - xc; v = v - vc                         # barycentric (inertial) frame
    s = np.array(STATE["M"]); x[iM] = s[:3]; v[iM] = s[3:]   # keep M as Tellus-relative vector
    return x, v, gm

def ephemeris(days):
    """dict day -> {body: (lon, lat)} for integer days (any sign of day - T0)."""
    days = sorted(set(int(d) for d in days)); out = {}
    iT = ORDER.index('T')
    for sign in (1, -1):
        want = set(d for d in days if (d - T0) * sign >= 0)
        if not want: continue
        x, v, gm = initial(); h = sign / SUB
        last = max(abs(d - T0) for d in want)
        for k in range(last + 1):
            d = T0 + sign * k
            if d in want:
                res = {}
                for i, b in enumerate(ORDER):
                    if b == 'T': continue
                    g = x[i] - x[iT] if b != "M" else x[i]
                    res[b] = (np.degrees(np.arctan2(g[1], g[0])) % 360.0, np.degrees(np.arcsin(g[2] / np.linalg.norm(g))))
                out[d] = res
            if k < last:
                for _ in range(SUB):
                    x, v = step(x, v, gm, h)
    return out

def main(path='predictions.csv'):
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
