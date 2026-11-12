"""Stage 3 predictor: direction of body N for the days in predictions_3.csv.

Equations (see ANSWER_3.md):
  Earth (Tellus) heliocentric Kepler orbit with a = 1 AU; Sun GM = NE^2 (AU^3/day^2) (stage 1 Law 1).
  N: heliocentric Kepler orbit with its own gravitational parameter GM_N = K * NE^2, so
     a_N = (K * NE^2 / n_N^2)^(1/3).  Elements below are the least-squares fit to observations_3.csv.
  Geocentric direction of N = (r_N - r_Earth) / |r_N - r_Earth|. Output: lon, lat in degrees.
"""
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

# Earth (stage 1, Law 2 constants)
NE = 0.01720207
E_EARTH, INC_EARTH, OM_EARTH, OM_EARTH_PERI, M0_EARTH = 0.01997902, 0.37877084, 3.14168505, -1.50264368, 15.2238328
# N (stage 3 fit): e, inc, Omega, omega, M0 (day 0), n (rad/day), K = GM_N / GM_S
N_ELEM = (0.7417992976344417, 0.5777532098351433, 0.3988447130735509, 1.0533917405445554,
          2.2066958304107733, 0.0006301721113756207)
K_N = 1.3262935032203114

def kepler_E(M, e):
    E = M + e*np.sin(M)
    for _ in range(60):
        dE = (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-13:
            break
    return E

def kepler_pos(t, a, e, inc, Om, om, M0, n):
    M = np.mod(M0 + n*t, 2*np.pi)
    E = kepler_E(M, e)
    x = a*(np.cos(E) - e); y = a*np.sqrt(1 - e*e)*np.sin(E)
    cO, sO = np.cos(Om), np.sin(Om); ci, si = np.cos(inc), np.sin(inc); co, so = np.cos(om), np.sin(om)
    X = (cO*co - sO*so*ci)*x + (-cO*so - sO*co*ci)*y
    Y = (sO*co + cO*so*ci)*x + (-sO*so + cO*co*ci)*y
    Z = (so*si)*x + (co*si)*y
    return np.stack([X, Y, Z], axis=-1)

def earth(t):
    return kepler_pos(t, 1.0, E_EARTH, INC_EARTH, OM_EARTH, OM_EARTH_PERI, M0_EARTH, NE)

def predict_N(t):
    e, inc, Om, om, M0, n = N_ELEM
    a = (K_N*NE**2/n**2)**(1/3)
    rN = kepler_pos(t, a, e, inc, Om, om, M0, n) - earth(t)
    u = rN/np.linalg.norm(rN, axis=1)[:, None]
    lon = np.arctan2(u[:, 1], u[:, 0]) % (2*np.pi)
    lat = np.arcsin(np.clip(u[:, 2], -1, 1))
    return np.degrees(lon), np.degrees(lat)

def main():
    rows = list(csv.DictReader(open(HERE/'predictions_3.csv')))
    days = np.array([float(r['day']) for r in rows])
    lon, lat = predict_N(days)
    with open(HERE/'predictions_3.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['day', 'body', 'lon', 'lat'])
        for r, lo_, la_ in zip(rows, lon, lat):
            w.writerow([r['day'], r['body'], '%.5f' % lo_, '%.5f' % la_])
    print('wrote', len(rows), 'rows')

if __name__ == '__main__':
    main()
