#!/usr/bin/env python
"""Stage 3: predictions for the new faint body N, using the stage-2 law (see ANSWER_3.md):
  d2r/dt2 = -(1-beta) GM_S r/|r|^3 + sum_{j in T,P2..P5} G m_j [(r_j-r)/|r_j-r|^3 - r_j/|r_j|^3]
with beta_N = -0.3233 and the osculating heliocentric state at epoch day 36379 below (AU, AU/day, ecliptic frame).
Planet positions come from the stage-1 model (predict.py); constants (masses, GM_S) as in predict_2.py."""
import csv, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import predict as P1M
import predict_2 as P2
HERE = os.path.dirname(os.path.abspath(__file__))
EPOCH = 36379
STATE = [0.470337592789, 1.61219557681, 1.94761891363, -0.0156635189804, -0.000897650846045, 0.00443107341084]
BETA = -0.3233237328

def run():
    rows = list(csv.reader(open(os.path.join(HERE, 'predictions_3.csv'))))
    head, body = rows[0], rows[1:]
    days = np.array([int(float(r[0])) for r in body])
    tmax = int(days.max())
    tg = np.arange(0, tmax + 1.0)
    pos = P1M.nbody(tg)
    helio = pos - pos[0][None]
    earth = helio[1]
    spl = CubicSpline(tg, helio[1:6].reshape(15, -1).T)
    gm = P2.MU * P2.M_PERT
    def rhs(t, y):
        r = y[:3]
        rp = spl(t).reshape(5, 3)
        d = rp - r
        d3 = (d ** 2).sum(1) ** 1.5
        r3 = (rp ** 2).sum(1) ** 1.5
        rn = np.sqrt(r @ r)
        acc = -P2.MU * (1 - BETA) * r / rn ** 3 + ((gm / d3)[:, None] * d - (gm / r3)[:, None] * rp).sum(0)
        return np.r_[y[3:], acc]
    te = np.unique(days).astype(float)
    sol = solve_ivp(rhs, (EPOCH, te[-1]), np.array(STATE, float), method='DOP853', t_eval=te, rtol=1e-11, atol=1e-14)
    lon, lat = P1M.to_frame(sol.y[:3] - earth[:, te.astype(int)])
    idx = {int(d): i for i, d in enumerate(te)}
    out = [[r[0], r[1], '%.5f' % lon[idx[d]], '%.5f' % lat[idx[d]]] for r, d in zip(body, days)]
    with open(os.path.join(HERE, 'predictions_3.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(head); w.writerows(out)
    print('wrote', len(out))
if __name__ == '__main__':
    run()
