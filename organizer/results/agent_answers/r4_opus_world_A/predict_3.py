"""predict_3.py -- predictions for the new faint body N, the 20 years after the last sighting (fills predictions_3.csv).
Uses only the laws and constants of ANSWER_3.md, section 1 (= stage-2 laws, with the encounter rule made explicit):
  Laws 1-3 and Law 2 constants (S, Tellus with m_T = 3.0e-6 M_S, P1..P5): as in predict_2.py (imported).
  Law 4: d2r/dt2 = kappa * [Newtonian pull of S, Tellus, P1..P5] - (heliocentric frame term),  kappa = 1 - beta.
  Law 4c (encounter rule): when a faint body passes within D_RESET of Tellus or a P, beta is reset to BETA_AFTER
          at the moment of closest approach.
Constants for N (fitted to observations_3.csv, days 36233-36525, 207 sightings, rms 0.99 arcmin):
"""
import csv, math
import numpy as np
from scipy.integrate import solve_ivp
import predict_2 as P2           # stage-2 Law-2 constants (overrides stage-1) and Law-4 field
import predict as PR             # frame transformation (Law 1)

N_T0 = 36379.0                   # epoch [day]
N_STATE = [0.47028267655105255, 1.6120524166844297, 1.9474579666289848,          # heliocentric ecliptic position [L]
           -0.015662123039095376, -0.0008973821442136387, 0.004430795141240435]  # velocity [L/day]
N_BETA = -0.3234540918211456     # +- 0.0021 (1 sigma)
D_RESET = 0.1                    # [L] encounter distance below which beta is reset
BETA_AFTER = 0.0                 # +- 0.016  (value of beta after a reset)

def main():
    rows = list(csv.DictReader(open('predictions_3.csv')))
    days = np.array([int(r['day']) for r in rows]); nd = int(days.max()) + 1
    X, V, H = P2.planet_table(nd)
    def planets(t):
        k = min(int(t // H), len(X) - 2); s = t / H - k
        h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s; h01 = -2*s**3 + 3*s**2; h11 = s**3 - s**2
        return h00*X[k] + h10*H*V[k] + h01*X[k+1] + h11*H*V[k+1], h00*V[k] + h01*V[k+1]
    t_end = float(days.max())
    t, y, beta = N_T0, np.array(N_STATE, float), N_BETA
    pieces = []
    while True:
        rhs = P2.make_rhs(X, V, H, 1.0 - beta)
        events = []
        for j in range(6):             # closest approach to body j (Tellus, P1..P5) while within D_RESET
            def ev(tt, yy, j=j):
                rp, vp = planets(tt); d = yy[:3] - rp[j]
                return d @ (yy[3:] - vp[j]) if d @ d < D_RESET**2 else -1.0
            ev.terminal = True; ev.direction = 1
            events.append(ev)
        s = solve_ivp(rhs, (t, t_end), y, method='DOP853', rtol=1e-11, atol=1e-14, dense_output=True,
                      events=events if beta != BETA_AFTER else None)
        hit = [te[0] for te in (s.t_events or []) if len(te)]
        t_stop = min(hit) if hit else t_end
        pieces.append((t, t_stop, s.sol))
        if not hit: break
        t, y, beta = t_stop, s.sol(t_stop), BETA_AFTER
        print(f'N: close encounter at day {t:.2f} -> beta reset to {BETA_AFTER}')
    def pos(d):
        for a, b, sol in pieces:
            if a - 1e-9 <= d <= b + 1e-9: return sol(d)[:3]
        raise ValueError(d)
    for r, d in zip(rows, days):
        g = pos(float(d)) - X[int(round(d / H)), 0]
        lo, la = PR.to_obs(g)
        r['lon'] = f'{lo:.5f}'; r['lat'] = f'{la:.5f}'
    with open('predictions_3.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['day', 'body', 'lon', 'lat']); w.writeheader(); w.writerows(rows)
    print('wrote', len(rows), 'predictions')

if __name__ == '__main__':
    main()
