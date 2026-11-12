#!/usr/bin/env python
"""predict_2.py -- Newtonian N-body model of the Tellus sky, stage 2 (see ANSWER_2.md).

Laws (units: AU, day, solar masses):
  1. Massive bodies (S, E, M, P1..P5) obey  d2 r_i/dt2 = sum_j GM_j (r_j - r_i)/|r_j - r_i|^3, GM_j = GM_S m_j
     (exactly the stage 1 law; integrated from day-0 osculating elements, constants below).
  2. Faint bodies C1..C12 are massless test particles. Each has one extra constant beta_k and feels
        d2 x/dt2 = (1 - beta_k) * sum_j GM_j (r_j - x)/|r_j - x|^3      (sum over S, E, M, P1..P5)
  3. Observed direction of any body X = direction of (r_X - r_E) in the ecliptic frame, rotated to the
     catalogue frame about the x axis by EPS.
The comets start from osculating heliocentric elements (with mu = GM_S (1-beta)) at day 36525.
"""
import csv, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = np.radians(21.7015)
GM = 0.00029592
MASS = {
    'E': 3e-06,
    'M': 3.7e-08,
    'P1': 1.999999952e-06,
    'P2': 6.5e-07,
    'P3': 0.0001495552455,
    'P4': 0.001000034104,
    'P5': 0.0002498087873,
}
ELEM = {
    'E': [0.9993189185, 0.02068620442, 93.46219421, 246.1448195, 0.004999151505, 29.99645309],
    'M': [0.0025699901, 0.05500138268, 85.28414579, 283.1385141, 5.001018853, 298.0131843],
    'P1': [0.4499985152, 0.06001027355, 166.4189188, 218.3185227, 0, 0],
    'P2': [1.699993824, 0.05000282228, 281.1827859, 114.9177005, 0, 0],
    'P3': [3.300006292, 0.03999251035, 42.54492847, 151.6882368, 1.499002796, 118.6446338],
    'P4': [6.199939859, 0.05000946356, 211.514104, 356.634093, 0.99853042, 163.2441527],
    'P5': [11.49979586, 0.04000596017, 167.6976912, 77.80880169, 2.000489004, 73.20194855],
}
# comets: [a (AU), e, i (deg), Omega (deg), omega (deg), M0 (deg), beta] at day EPOCH
EPOCH = 36525
COMET = {
    'C1': [13.1920812, 0.9470797481, 4.607706063, 343.6969456, 79.86348615, -68.04881983, -0.2649256974],
    'C2': [7.469462234, 0.7798901826, 0.263649794, 193.0011092, 8.605798272, -43.06707075, 0.09917230823],
    'C3': [5.747183687, 0.6974737459, 9.622451792, 239.0843895, 36.12812516, -2.268120135, -0.2309467296],
    'C4': [16.29476974, 0.9170415637, 5.111067744, 143.2206726, 223.9084438, -57.75299569, 0.2736137251],
    'C5': [4.370169052, 0.3873546481, 18.16446583, 315.4094625, 6.243057267, 77.37476753, -0.1720997746],
    'C6': [5.743248744, 0.5001778637, 15.94784304, 326.4385219, 322.2345154, -8.880850577, 0.1440887182],
    'C7': [2.029288487, 0.506800086, 17.46137437, 223.4557424, 359.0393101, 179.1508763, 0.2170660773],
    'C8': [6.467976773, 0.8335212947, 17.19239607, 48.88764841, 101.9099843, 10.5384216, 0.297960051],
    'C9': [5.190586691, 0.7651779501, 13.94042656, 150.1927792, 19.30602465, -156.2413511, 0.08238736733],
    'C10': [9.109526153, 0.9088409305, 14.14833961, 108.6804694, 139.6618972, -9.946591619, -0.03553842873],
    'C11': [3.213931988, 0.5495704912, 2.575903371, 32.97773101, 57.83534107, -76.56449603, -0.1573778496],
    'C12': [9.48992845, 0.7677426223, 15.78047186, 221.9657652, 357.2171335, 71.0492963, -0.2973589026],
}
ORDER = ['S', 'E', 'M', 'P1', 'P2', 'P3', 'P4', 'P5']
NSUB = 8            # planetary integration steps per day

def kepE(M, e):
    E = M + e * np.sin(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E

def elem2cart(mu, a, e, w, L0, i, Om):
    w, L0, i, Om = [np.radians(x) for x in (w, L0, i, Om)]
    E = kepE(L0 - w, e)
    n = np.sqrt(mu / a ** 3)
    x = a * (np.cos(E) - e); y = a * np.sqrt(1 - e * e) * np.sin(E)
    vx = -a * n * np.sin(E) / (1 - e * np.cos(E)); vy = a * n * np.sqrt(1 - e * e) * np.cos(E) / (1 - e * np.cos(E))
    om = w - Om
    def rot(x, y):
        xo = x * np.cos(om) - y * np.sin(om); yo = x * np.sin(om) + y * np.cos(om)
        return np.array([xo * np.cos(Om) - yo * np.cos(i) * np.sin(Om),
                         xo * np.sin(Om) + yo * np.cos(i) * np.cos(Om), yo * np.sin(i)])
    return rot(x, y), rot(vx, vy)

def initial_state():
    gm = np.array([GM] + [GM * MASS[b] for b in ORDER[1:]])
    pos = np.zeros((8, 3)); vel = np.zeros((8, 3))
    for b in ['E', 'P1', 'P2', 'P3', 'P4', 'P5']:
        k = ORDER.index(b)
        pos[k], vel[k] = elem2cart(GM * (1 + MASS[b]), *ELEM[b])
    p, v = elem2cart(gm[1] + gm[2], *ELEM['M'])
    pos[2] = pos[1] + p; vel[2] = vel[1] + v
    mt = gm / GM
    c = (mt[:, None] * pos).sum(0) / mt.sum(); vc = (mt[:, None] * vel).sum(0) / mt.sum()
    return gm, pos - c, vel - vc

def accel(gm, x):
    d = x[None, :, :] - x[:, None, :]
    r2 = (d * d).sum(2); np.fill_diagonal(r2, 1.0)
    inv3 = 1.0 / (r2 * np.sqrt(r2)); np.fill_diagonal(inv3, 0.0)
    return np.einsum('ij,ijk->ik', gm[None, :] * inv3, d)

W1, W2, W3 = -1.17767998417887, 0.235573213359357, 0.784513610477560
DK = np.array([W3, W2, W1, 1 - 2 * (W1 + W2 + W3), W1, W2, W3])
CD = np.empty(8); CD[0] = DK[0] / 2; CD[1:7] = (DK[:-1] + DK[1:]) / 2; CD[7] = DK[6] / 2

def integrate_dense(ndays):
    """positions of the 8 massive bodies at every 1/NSUB day (Yoshida 6th order symplectic)"""
    gm, x, v = initial_state()
    out = np.empty((ndays * NSUB + 1, 8, 3)); out[0] = x
    h = 1.0 / NSUB; n = 0
    for day in range(ndays):
        for _ in range(NSUB):
            for k in range(7):
                x = x + CD[k] * h * v
                v = v + DK[k] * h * accel(gm, x)
            x = x + CD[7] * h * v
            n += 1; out[n] = x
    return gm, out

def kep_state(mu, a, e, i, Om, om, M0):
    i, Om, om, M0 = [np.radians(q) for q in (i, Om, om, M0)]
    E = kepE(np.array([M0]), e)[0]; n = np.sqrt(mu / a ** 3)
    xp = np.array([a * (np.cos(E) - e), a * np.sqrt(1 - e * e) * np.sin(E), 0.0])
    vp = np.array([-a * n * np.sin(E), a * n * np.sqrt(1 - e * e) * np.cos(E), 0.0]) / (1 - e * np.cos(E))
    ci, si, cO, sO, co, so = np.cos(i), np.sin(i), np.cos(Om), np.sin(Om), np.cos(om), np.sin(om)
    R = np.array([[cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
                  [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
                  [so * si, co * si, ci]])
    return R @ xp, R @ vp

def integrate_comets(gm, D, names, tend):
    """RK4, step 1/8 day, all comets at once; D = dense massive-body positions (grid 1/NSUB day, NSUB=8)."""
    k0 = EPOCH * NSUB
    beta = np.array([COMET[c][6] for c in names]); lam = 1.0 - beta
    # initial barycentric state
    xs = D[k0, 0]; vs = (D[k0 - 2, 0] - 8 * D[k0 - 1, 0] + 8 * D[k0 + 1, 0] - D[k0 + 2, 0]) / (12.0 / NSUB)
    x = np.empty((len(names), 3)); v = np.empty((len(names), 3))
    for j, c in enumerate(names):
        a, e, i, Om, om, M0, b = COMET[c]
        r, w = kep_state(GM * (1 - b), a, e, i, Om, om, M0)
        x[j] = r + xs; v[j] = w + vs
    h = 1.0 / NSUB
    nsteps = (tend - EPOCH) * NSUB
    # midpoints of planet positions (cubic interpolation) for the whole range
    G = D[k0 - 1: k0 + nsteps + 3]                      # grid k0-1 .. k0+nsteps+2
    P0 = G[1: nsteps + 1]
    P1 = G[2: nsteps + 2]
    Pm = (-G[0: nsteps] + 9 * G[1: nsteps + 1] + 9 * G[2: nsteps + 2] - G[3: nsteps + 3]) / 16.0
    gmv = gm[None, :]
    def acc(P, y):
        d = P[None, :, :] - y[:, None, :]
        r2 = (d * d).sum(2); f = gmv / (r2 * np.sqrt(r2))
        return lam[:, None] * np.einsum('cj,cjk->ck', f, d)
    res = {}
    day_out = {}
    for s in range(nsteps):
        p0, pm, p1 = P0[s], Pm[s], P1[s]
        k1v = acc(p0, x); k1x = v
        k2x = v + 0.5 * h * k1v; k2v = acc(pm, x + 0.5 * h * k1x)
        k3x = v + 0.5 * h * k2v; k3v = acc(pm, x + 0.5 * h * k2x)
        k4x = v + h * k3v;       k4v = acc(p1, x + h * k3x)
        x = x + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        if (s + 1) % NSUB == 0:
            day_out[EPOCH + (s + 1) // NSUB] = x.copy()
    day_out[EPOCH] = None
    return day_out

def from_ecl(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    v = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    c, s = np.cos(EPS), np.sin(EPS)
    x = v[0]; y = c * v[1] + s * v[2]; z = -s * v[1] + c * v[2]
    return np.degrees(np.arctan2(y, x)) % 360, np.degrees(np.arcsin(z))

def main():
    path = os.path.join(HERE, 'predictions_2.csv')
    rows = list(csv.reader(open(path)))
    head, body = rows[0], rows[1:]
    days = [int(float(r[0])) for r in body]
    tend = max(days)
    gm, D = integrate_dense(tend + 1)
    names = sorted(COMET, key=lambda s: int(s[1:]))
    cpos = integrate_comets(gm, D, names, tend)
    for r in body:
        d = int(float(r[0])); b = r[1]
        earth = D[d * NSUB, 1]
        if b in ORDER:
            tgt = D[d * NSUB, ORDER.index(b)]
        else:
            tgt = cpos[d][names.index(b)]
        v = tgt - earth
        lon = np.degrees(np.arctan2(v[1], v[0])) % 360
        lat = np.degrees(np.arctan2(v[2], np.hypot(v[0], v[1])))
        lo, la = from_ecl(lon, lat)
        r[2] = '%.5f' % lo; r[3] = '%.5f' % la
    with open(path, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(head); w.writerows(body)

if __name__ == '__main__':
    main()
