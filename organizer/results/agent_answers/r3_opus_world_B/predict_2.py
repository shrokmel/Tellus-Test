#!/usr/bin/env python
"""Stage 2 predictions for all bodies (S, P1-P5, M and the faint bodies C1-C12), using ONLY the laws and
constants of ANSWER_2.md:

  Law 1 (unchanged from stage 1): the 8 massive bodies S, P1..P5, Tellus (T), M move by mutual Newtonian gravity
        d2 r_i/dt2 = sum_{j != i} GM_j (r_j - r_i)/|r_j - r_i|^3          (i, j over the 8 massive bodies)
  Law 1' (new): each faint body C_k is a massless test particle whose response to gravity is scaled by its own
        constant kappa_k:
        d2 r_k/dt2 = kappa_k * sum_j GM_j (r_j - r_k)/|r_j - r_k|^3     (j over the 8 massive bodies)
        (kappa = 1 for every stage-1 body.)
  Law 2 (unchanged): seen direction = unit(r_body - r_T), rotated from the ecliptic frame to the data frame by EPS.

Constants: stage-1 massive-body constants (copied verbatim from predict.py / ANSWER_1.md) and, for each faint
body, its heliocentric ecliptic state at its own epoch plus kappa (ANSWER_2.md section 2).
Usage: python predict_2.py      (fills predictions_2.csv in this folder)
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
    'C1': (23898.0, [5.064043956400582e-01, 5.075136174146937e-01, 5.046357156506411e-02, -2.581662226579858e-02, 1.860743832288342e-02, 8.703510060588900e-04], 1.2648588895),
    'C2': (21667.0, [-1.160452615804124e+00, 2.517949934545693e+00, 1.153898266947780e-03, -5.929915938675484e-03, -1.100899681465726e-02, 9.892908167625216e-06], 0.9002732402),
    'C3': (31750.0, [-2.954550189550277e+00, -2.449582991598033e-01, -4.155413897593374e-01, 8.345285579314158e-03, -1.030796772724398e-02, 2.106563967861412e-03], 1.2304604753),
    'C4': (12899.0, [8.436012970864092e-01, -1.638514475356163e-01, 1.753945534129911e-03, 3.515554786940157e-03, 2.154457450176655e-02, -3.853010879312998e-04], 0.7263792070),
    'C5': (32576.0, [2.691159241303588e-02, -2.919341485487523e+00, -6.779997802452664e-01, 1.164885981986144e-02, 2.434986123940082e-03, 3.241943835495705e-03], 1.1716410757),
    'C6': (20446.0, [1.028443233699889e+00, -2.616210486391418e+00, -4.640189372032263e-01, 1.070734287095505e-02, 3.545876064139569e-03, 2.532989686274660e-03], 0.8559916826),
    'C7': (31638.0, [2.547084268418573e+00, 1.571963571944975e+00, 1.933898217738791e-01, -2.005568516333203e-03, 5.796783564746286e-03, -1.736699310920891e-03], 0.7825716987),
    'C8': (28914.0, [2.392854198915395e+00, 1.786667281206179e+00, -1.976437748566006e-01, -1.001076406862621e-02, 5.714846829894579e-04, 2.458390972236576e-03], 0.7022933231),
    'C9': (29539.0, [-9.603566067682358e-01, -1.197709270784213e+00, 3.767150796807556e-01, 7.395383492769993e-03, -1.520028374511611e-02, 2.373027313841671e-03], 0.9181793504),
    'C10': (17251.0, [-1.320711953687272e-01, -8.133509746684270e-01, 1.008247336154839e-01, 2.495656165981380e-02, -7.263866010801062e-03, -5.255500294374044e-03], 1.0337230383),
    'C11': (31333.0, [-2.671925773666382e+00, -1.360870358790272e+00, 1.677381700004654e-02, -1.084204675481501e-03, -1.102118502647842e-02, -3.920723300979094e-04], 1.1578933981),
    'C12': (34671.0, [-1.743761451987638e+00, -1.355056271490174e+00, -4.549240112004312e-02, 1.043596749862034e-02, -1.324750379822914e-02, 4.751607314503760e-03], 1.2965553855),
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
    tf = np.arange(T0, tmax + H, H); tb = np.arange(T0, tmin - H, -H)
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
    path = os.path.join(HERE, 'predictions_2.csv')
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
