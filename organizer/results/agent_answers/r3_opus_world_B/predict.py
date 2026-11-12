#!/usr/bin/env python
"""Predict sky directions of S, P1-P5 and M from Tellus, using ONLY the laws and constants of ANSWER_1.md:

  Law 1 (gravitation + motion): every pair of the 8 bodies S, P1..P5, Tellus (T), M attracts with
        d^2 r_i/dt^2 = sum_{j!=i} GM_j (r_j - r_i)/|r_j - r_i|^3
  Law 2 (what is seen): direction of body b = unit vector (r_b - r_T), expressed in the observation
        frame, which is the ecliptic frame below rotated about its x-axis by the obliquity EPS.

Constants (units: length = AU defined by GM_S = k^2, k = 0.01720209895 rad/day; time = day):
  the epoch state vectors and GM ratios below, copied from ANSWER_1.md section 2.
Usage: python predict.py   (reads/writes predictions.csv in this folder)
"""
import os, csv
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = 0.01720209895 ** 2
EPS_DEG = 21.7015          # obliquity: tilt of the ecliptic (S orbit plane) against the data frame
T0 = 18262.0               # epoch (day number) of the state vectors
# ---- constants (filled from the fit; identical to ANSWER_1.md) ----
#CONSTANTS_BEGIN
HELIO = {  # heliocentric ecliptic state at day T0: x, y, z [AU], vx, vy, vz [AU/day]
    'P1': np.array([9.107039555955809e-02, 4.381799797616328e-01, 6.876757570326467e-06, -2.551251350825947e-02, 3.728892581783826e-03, -3.854758287724350e-07]),
    'P2': np.array([1.211335998398247e+00, -1.090270207879586e+00, 2.225559837725770e-04, 9.488709418440766e-03, 9.951699545907931e-03, 1.833548304715248e-08]),
    'T': np.array([-4.077711315289091e-01, -9.323641698175303e-01, -5.381247386083418e-06, 1.542076323091913e-02, -6.924137738449534e-03, 2.113917494716724e-07]),
    'P3': np.array([6.577059815300164e-02, -3.380320231280422e+00, 4.028981419233729e-02, 9.225630646179445e-03, 4.702482765405692e-04, -2.173651375402460e-04]),
    'P4': np.array([1.317547351922277e+00, 6.270565709744375e+00, -1.117755057926885e-01, -6.586106156290334e-03, 1.135558300949610e-03, 1.432200990221717e-05]),
    'P5': np.array([-1.098132063639481e+01, -5.699850284543776e-01, 3.604601471349140e-01, 2.269522674659104e-04, -5.287969926121252e-03, -6.192921245552120e-05]),
}
GEO_M = np.array([2.381410475774136e-03, -1.065661365359026e-04, -1.775353407622004e-04, 4.625734883410472e-05, 5.736861898920633e-04, 2.740807129713245e-05])  # Moon minus Tellus state at T0
MASS = {  # GM_body / GM_S
    'P1': 2.7448599204e-07,
    'P2': 7.0767556281e-07,
    'T': 2.6083368644e-06,
    'P3': 1.5057322334e-04,
    'P4': 9.9991564228e-04,
    'P5': 2.5012546564e-04,
    'M': 3.1949171821e-08,
}
#CONSTANTS_END
NAMES = ['S', 'P1', 'P2', 'T', 'P3', 'P4', 'P5', 'M']


def initial_state():
    """heliocentric states (ecliptic J-frame) -> barycentric states and GM array"""
    s = np.zeros((8, 6))
    for b in ['P1', 'P2', 'T', 'P3', 'P4', 'P5']:
        s[NAMES.index(b)] = HELIO[b]
    s[7] = s[3] + GEO_M                       # Moon = Tellus + geocentric moon state
    GM = np.array([K2] + [K2 * MASS[b] for b in NAMES[1:]])
    s -= (s * GM[:, None]).sum(0) / GM.sum()  # barycentric
    return s, GM


def rhs(t, y, GM):
    s = y.reshape(8, 6); r = s[:, :3]
    d = r[None, :, :] - r[:, None, :]
    r2 = (d * d).sum(-1); np.fill_diagonal(r2, 1.0)
    inv3 = r2 ** -1.5; np.fill_diagonal(inv3, 0.0)
    a = (d * (inv3 * GM[None, :])[..., None]).sum(1)
    return np.concatenate([s[:, 3:], a], 1).ravel()


def positions(days):
    days = np.asarray(sorted(set(days)), float)
    s0, GM = initial_state()
    out = {}
    for sel in (days[days >= T0], days[days < T0][::-1]):
        if len(sel) == 0:
            continue
        te = sel
        sol = solve_ivp(rhs, (T0, te[-1]), s0.ravel(), method='DOP853', t_eval=te,
                        rtol=1e-12, atol=1e-17, args=(GM,))
        for k, d in enumerate(te):
            out[d] = sol.y[:, k].reshape(8, 6)[:, :3]
    return out


def to_frame(v):
    c, s = np.cos(np.radians(EPS_DEG)), np.sin(np.radians(EPS_DEG))
    x, y, z = v
    return np.array([x, c * y + s * z, -s * y + c * z])


def main():
    path = os.path.join(HERE, 'predictions.csv')
    rows = list(csv.DictReader(open(path)))
    pos = positions([float(r['day']) for r in rows])
    for r in rows:
        X = pos[float(r['day'])]
        g = X[NAMES.index(r['body'])] - X[3]
        g = to_frame(g / np.linalg.norm(g))
        r['lon'] = '%.5f' % (np.degrees(np.arctan2(g[1], g[0])) % 360.0)
        r['lat'] = '%.5f' % np.degrees(np.arcsin(g[2]))
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['day', 'body', 'lon', 'lat'])
        w.writeheader(); w.writerows(rows)
    print('wrote', len(rows), 'predictions')


if __name__ == '__main__':
    main()
