#!/usr/bin/env python
"""Stage 2 predictions for S, P1-P5, M (stage-1 laws, predict.py) and the faint bodies C1-C12.

Faint bodies are massless test particles in the field of S, Tellus (T) and P2-P5, with a body-specific factor
(1-beta) on the attraction of S:
    d2r/dt2 = -(1-beta) GM_S r/|r|^3 + sum_j G m_j [ (r_j-r)/|r_j-r|^3 - r_j/|r_j|^3 ],  j in {T,P2,P3,P4,P5}
(r heliocentric, ecliptic frame; the second term includes the indirect term from S's own acceleration).
G m_S = n_T^2 (units: Tellus semi-major axis=1, day).  m_T=3.0e-6 m_S, m_P2=0 (consistent with 0), m_P3..P5 as in stage 1.
Each faint body is given by osculating heliocentric state (r, v) at an epoch and beta (table below).
"""
import csv, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import predict as P1M   # stage-1 laws (S, P1-P5, M) -- same folder

HERE = os.path.dirname(os.path.abspath(__file__))
MU = P1M.GMS
M_PERT = np.array([3.0e-6, 0.0, 0.0001514640239296227, 0.0010003317475393997, 0.000249113903835644])  # T,P2,P3,P4,P5 / m_S

# body: (epoch day, [x,y,z, vx,vy,vz] heliocentric ecliptic AU, AU/day, beta)
COMETS = {
    'C1': (8245, [-1.20829710481, 1.02336220601, 0.0566395786315, -0.0221372174294, 6.22665217597e-05, -0.000480468780839], 0.00000000),
    'C2': (33137, [-2.59533956663, 0.0404277880324, -0.401878873801, 0.00139533951801, -0.0117085933268, 0.00451317269233], 0.00068680),
    'C3': (5114, [-0.342603068442, -1.70905184574, 0.0643460783961, 0.0186779571282, -0.00162304015426, 0.00311559604394], -0.25370880),
    'C4': (36071, [-1.71695681136, 1.09024239144, -0.0160957354767, -0.00604902701985, -0.00707510109013, 0.000136786936309], 0.00000952),
    'C5': (30660, [1.08016035423, -2.4510248206, -0.380833222499, 0.0103263827027, 0.00547821293682, 0.00383789523792], -0.00245979),
    'C6': (34685, [-0.660647069671, -2.68093919637, -0.751915009579, 0.0112647890579, -0.00246014039854, 0.0011964203895], 0.00000000),
    'C7': (35266, [-0.657728593409, -0.706871353961, 0.0137428961797, 0.0136440149602, -0.0125544837629, 0.00567110812803], 0.00078842),
    'C8': (36071, [-0.0668208259359, 1.3116360249, 0.283771688712, -0.0142711595561, -0.00842763638797, 0.00161009492086], 0.29706463),
    'C9': (20396, [0.382853030095, 2.10966105875, -0.549264774258, -0.012200289247, -0.00602090776533, 0.00286280027112], 0.08159541),
    'C10': (17169, [-1.41174991222, 0.542149949389, 0.281408832777, 0.00695605545283, -0.0176807433627, -0.00021903615467], -0.02558821),
    'C11': (35630, [-2.24992299336, 0.0655630411084, 0.0582440320943, -0.00757115287856, -0.0122035563545, -0.000248478964389], -0.15000000),
    'C12': (34661, [-1.84451533391, -1.2198022664, -0.0928222076878, 0.00979865023515, -0.0137061451803, 0.00472514700204], -0.29620895),
}

def run():
    rows = list(csv.reader(open(os.path.join(HERE, 'predictions_2.csv'))))
    head, body = rows[0], rows[1:]
    days = np.array([int(float(r[0])) for r in body])
    tmax = int(days.max())
    tg = np.arange(0, tmax + 1.0)
    # stage-1 model
    pos = P1M.nbody(tg)
    helio = pos - pos[0][None]
    earth = helio[1]
    eph = {}
    eph['S'] = P1M.to_frame(-earth)
    for k, nm in zip(range(2, 6), ['P2', 'P3', 'P4', 'P5']):
        eph[nm] = P1M.to_frame(helio[k] - earth)
    eph['P1'] = P1M.to_frame(P1M.p1_track(tg, helio, tg) - earth)
    eph['M'] = P1M.to_frame(P1M.moon_track(tg, CubicSpline(tg, (-earth).T)))
    # faint bodies
    spl = CubicSpline(tg, helio[1:6].reshape(15, -1).T)
    gm = MU * M_PERT
    def make_rhs(beta):
        def rhs(t, y):
            r = y[:3]
            rp = spl(t).reshape(5, 3)
            d = rp - r
            d3 = (d ** 2).sum(1) ** 1.5
            r3 = (rp ** 2).sum(1) ** 1.5
            rn = np.sqrt(r @ r)
            acc = -MU * (1 - beta) * r / rn ** 3 + ((gm / d3)[:, None] * d - (gm / r3)[:, None] * rp).sum(0)
            return np.r_[y[3:], acc]
        return rhs
    need = np.unique(days)
    for name, (tm, s, beta) in COMETS.items():
        te = need[need >= tm].astype(float)
        sol = solve_ivp(make_rhs(beta), (tm, te[-1]), np.array(s, float), method='DOP853', t_eval=te, rtol=1e-10, atol=1e-13)
        p = sol.y[:3]
        lon, lat = P1M.to_frame(p - earth[:, te.astype(int)])
        full_lon = np.full(tmax + 1, np.nan); full_lat = np.full(tmax + 1, np.nan)
        full_lon[te.astype(int)] = lon; full_lat[te.astype(int)] = lat
        eph[name] = (full_lon, full_lat)
    out = []
    for r, d in zip(body, days):
        lon, lat = eph[r[1]]
        out.append([r[0], r[1], '%.5f' % lon[d], '%.5f' % lat[d]])
    with open(os.path.join(HERE, 'predictions_2.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(head); w.writerows(out)
    print('wrote', len(out), 'rows')

if __name__ == '__main__':
    run()
