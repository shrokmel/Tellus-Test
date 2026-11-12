#!/usr/bin/env python
"""predict.py -- Newtonian N-body model of the Tellus sky (see ANSWER_1.md).

Laws used (units: AU, day, solar masses):
  1. Every body i feels the inverse-square gravitational acceleration of every other body j:
        d2 r_i / dt2 = sum_j GM_j (r_j - r_i) / |r_j - r_i|^3 ,  GM_j = GM_S * m_j
  2. Bodies: S (Sun), E (Tellus), M (the large bright body, orbits Tellus), P1..P5 (planets).
  3. The ecliptic frame is tilted from the observation frame about the x axis by EPS.
  4. Observed direction of body X = direction of (r_X - r_E) (geometric, no light-time).
Initial conditions are osculating Kepler elements at day 0 (constants below), integrated with a
6th-order symplectic (Yoshida) scheme, 8 steps per day.
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
    'P2': 2.9999997e-07,
    'P3': 0.0001495552455,
    'P4': 0.001000034104,
    'P5': 0.0002498087873,
}          # mass / solar mass
ELEM = {
    'E': [0.9993189185, 0.02068620442, 93.46219421, 246.1448195, 0.004999151505, 29.99645309],
    'M': [0.0025699901, 0.05500138268, 85.28414579, 283.1385141, 5.001018853, 298.0131843],
    'P1': [0.4499985152, 0.06001027355, 166.4189188, 218.3185227, 0, 0],
    'P2': [1.699993541, 0.05000254486, 281.1856249, 114.9177791, 0, 0],
    'P3': [3.30000101, 0.03999325026, 42.54288041, 151.6885132, 1.498982093, 118.6668185],
    'P4': [6.199939859, 0.05000946356, 211.514104, 356.634093, 0.99853042, 163.2441527],
    'P5': [11.49979586, 0.04000596017, 167.6976912, 77.80880169, 2.000489004, 73.20194855],
}          # a [AU], e, varpi [deg], L0 [deg], i [deg], Omega [deg] at day 0
ORDER = ['S', 'E', 'M', 'P1', 'P2', 'P3', 'P4', 'P5']
NSUB = 8

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

_I, _J = np.triu_indices(8, 1)
def accel(gm, x):
    d = x[_J] - x[_I]                       # (28,3)
    r2 = (d * d).sum(1)
    f = d / (r2 * np.sqrt(r2))[:, None]
    a = np.zeros_like(x)
    np.add.at(a, _I, gm[_J][:, None] * f)
    np.add.at(a, _J, -gm[_I][:, None] * f)
    return a

def accel_fast(gm, x):
    d = x[None, :, :] - x[:, None, :]       # d[i,j] = x_j - x_i
    r2 = (d * d).sum(2); np.fill_diagonal(r2, 1.0)
    inv3 = 1.0 / (r2 * np.sqrt(r2)); np.fill_diagonal(inv3, 0.0)
    return np.einsum('ij,ijk->ik', gm[None, :] * inv3, d)

W1, W2, W3 = -1.17767998417887, 0.235573213359357, 0.784513610477560
DK = np.array([W3, W2, W1, 1 - 2 * (W1 + W2 + W3), W1, W2, W3])
CD = np.empty(8); CD[0] = DK[0] / 2; CD[1:7] = (DK[:-1] + DK[1:]) / 2; CD[7] = DK[6] / 2

def integrate(ndays):
    gm, x, v = initial_state()
    out = np.empty((ndays + 1, 8, 3)); out[0] = x
    h = 1.0 / NSUB
    for day in range(1, ndays + 1):
        for _ in range(NSUB):
            for k in range(7):
                x = x + CD[k] * h * v
                v = v + DK[k] * h * accel_fast(gm, x)
            x = x + CD[7] * h * v
        out[day] = x
    return out

def from_ecl(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    v = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    c, s = np.cos(EPS), np.sin(EPS)
    x = v[0]; y = c * v[1] + s * v[2]; z = -s * v[1] + c * v[2]
    return np.degrees(np.arctan2(y, x)) % 360, np.degrees(np.arcsin(z))

def main():
    path = os.path.join(HERE, 'predictions.csv')
    rows = list(csv.reader(open(path)))
    head, body = rows[0], rows[1:]
    days = [int(float(r[0])) for r in body]
    out = integrate(max(days))
    for r in body:
        d = int(float(r[0])); k = ORDER.index(r[1])
        v = out[d, k] - out[d, 1]
        lon = np.degrees(np.arctan2(v[1], v[0])) % 360
        lat = np.degrees(np.arctan2(v[2], np.hypot(v[0], v[1])))
        lo, la = from_ecl(lon, lat)
        r[2] = '%.5f' % lo; r[3] = '%.5f' % la
    with open(path, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(head); w.writerows(body)

if __name__ == '__main__':
    main()
