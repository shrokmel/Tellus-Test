#!/usr/bin/env python
"""Stage 3 predictions for the new faint body N, using ONLY the laws and
constants of ANSWER_2.md:

  Law 1 (unchanged from stage 1): the 8 massive bodies S, P1..P5, Tellus (T), M move by mutual Newtonian gravity
        d2 r_i/dt2 = sum_{j != i} GM_j (r_j - r_i)/|r_j - r_i|^3          (i, j over the 8 massive bodies)
  Law 1' (new): each faint body C_k is a massless test particle whose response to gravity is scaled by its own
        constant kappa_k:
        d2 r_k/dt2 = kappa_k * sum_j GM_j (r_j - r_k)/|r_j - r_k|^3     (j over the 8 massive bodies)
        (kappa = 1 for every stage-1 body.)
  Law 2 (unchanged): seen direction = unit(r_body - r_T), rotated from the ecliptic frame to the data frame by EPS.

Constants: stage-1 massive-body constants (copied verbatim from predict.py / ANSWER_1.md) and, for each faint
body, its heliocentric ecliptic state at its own epoch plus kappa (ANSWER_3.md section 1; N fitted to observations_3.csv).
Usage: python predict_3.py      (fills predictions_3.csv in this folder)
"""
import os, csv
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = 0.01720209895 ** 2
EPS_DEG = 21.7015
T0 = 18262.0
#STAGE1_BEGIN
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
#STAGE1_END
#FAINT_BEGIN
FAINT = {  # name: (epoch day, heliocentric ecliptic state x,y,z [AU], vx,vy,vz [AU/day] at epoch, kappa)
    'N': (36379.0, [4.701414461110018e-01, 1.611943077710376e+00, 1.947179880896090e+00, -1.566242312443316e-02, -8.985068525107176e-04, 4.429578196296246e-03], 1.3219979875),
}
#FAINT_END
NAMES = ['S', 'P1', 'P2', 'T', 'P3', 'P4', 'P5', 'M']


def initial_state():
    s = np.zeros((8, 6))
    for b in ['P1', 'P2', 'T', 'P3', 'P4', 'P5']:
        s[NAMES.index(b)] = HELIO[b]
    s[7] = s[3] + GEO_M
    GM = np.array([K2] + [K2 * MASS[b] for b in NAMES[1:]])
    s -= (s * GM[:, None]).sum(0) / GM.sum()
    return s, GM


def rhs_massive(t, y, GM):
    s = y.reshape(8, 6); r = s[:, :3]
    d = r[None, :, :] - r[:, None, :]
    r2 = (d * d).sum(-1); np.fill_diagonal(r2, 1.0)
    inv3 = r2 ** -1.5; np.fill_diagonal(inv3, 0.0)
    a = (d * (inv3 * GM[None, :])[..., None]).sum(1)
    return np.concatenate([s[:, 3:], a], 1).ravel()


H = 0.5
def massive_table(tmin, tmax):
    """barycentric states of the 8 massive bodies every H days (Law 1)"""
    s0, GM = initial_state()
    tf = np.arange(T0, tmax + H, H); tb = np.arange(T0, min(tmin, T0 - 1.0) - H, -H)
    parts = []
    for te in (tf, tb):
        sol = solve_ivp(rhs_massive, (te[0], te[-1]), s0.ravel(), method='DOP853', t_eval=te,
                        rtol=1e-12, atol=1e-17, args=(GM,))
        parts.append((te, sol.y.T.reshape(len(te), 8, 6)))
    t = np.r_[parts[1][0][::-1], parts[0][0][1:]]
    Y = np.r_[parts[1][1][::-1], parts[0][1][1:]]
    return t, Y, GM


class Ephem:
    def __init__(self, t, Y, GM):
        self.t, self.Y, self.GM = t, Y, GM
    def pos(self, tt):            # cubic Hermite interpolation of positions
        x = (tt - self.t[0]) / H; i = min(max(int(np.floor(x)), 0), len(self.t) - 2); s = x - i
        p0 = self.Y[i, :, :3]; p1 = self.Y[i + 1, :, :3]; v0 = self.Y[i, :, 3:] * H; v1 = self.Y[i + 1, :, 3:] * H
        return ((2*s**3 - 3*s**2 + 1) * p0 + (s**3 - 2*s**2 + s) * v0 + (-2*s**3 + 3*s**2) * p1 + (s**3 - s**2) * v1)
    def at(self, day):
        return self.Y[int(round((day - self.t[0]) / H))]


def faint_positions(E, name, days):
    tm, helio, kappa = FAINT[name]
    sun = E.at(tm)[0]
    y0 = np.asarray(helio, float) + sun
    def rhs(t, y):
        B = E.pos(t); d = B - y[None, :3]
        a = kappa * (d * (E.GM / ((d * d).sum(1)) ** 1.5)[:, None]).sum(0)
        return np.r_[y[3:], a]
    days = np.asarray(sorted(set(days)), float)
    out = {}
    for sel in (days[days >= tm], days[days < tm][::-1]):
        if len(sel) == 0: continue
        sol = solve_ivp(rhs, (tm, sel[-1]), y0, method='DOP853', t_eval=sel, rtol=1e-11, atol=1e-14)
        for k, d in enumerate(sel): out[d] = sol.y[:3, k]
    return out


def to_frame(v):
    c, s = np.cos(np.radians(EPS_DEG)), np.sin(np.radians(EPS_DEG))
    return np.array([v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]])


def main():
    path = os.path.join(HERE, 'predictions_3.csv')
    rows = list(csv.DictReader(open(path)))
    pdays = sorted(set(float(r['day']) for r in rows))
    tmin = min([pdays[0]] + [v[0] for v in FAINT.values()]) - 2
    E = Ephem(*massive_table(tmin, pdays[-1] + 2))
    fp = {b: faint_positions(E, b, [float(r['day']) for r in rows if r['body'] == b]) for b in FAINT}
    for r in rows:
        d = float(r['day']); X = E.at(d)
        tel = X[3, :3]
        if r['body'] in FAINT: p = fp[r['body']][d]
        else: p = X[NAMES.index(r['body']), :3]
        g = p - tel; g = to_frame(g / np.linalg.norm(g))
        r['lon'] = '%.5f' % (np.degrees(np.arctan2(g[1], g[0])) % 360.0)
        r['lat'] = '%.5f' % np.degrees(np.arcsin(g[2]))
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['day', 'body', 'lon', 'lat'])
        w.writeheader(); w.writerows(rows)
    print('wrote', len(rows), 'predictions')


if __name__ == '__main__':
    main()
