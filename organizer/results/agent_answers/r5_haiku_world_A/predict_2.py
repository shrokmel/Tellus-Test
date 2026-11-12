"""Stage 2 predictor. Writes predictions_2.csv (lon, lat in degrees, equatorial frame).

Stage 1 laws (unchanged; implemented in predict.py in this folder) give S, P1-P5, M.
Faint bodies C1-C12: heliocentric Kepler orbits with the common Sun GM (Law 1: n^2 a^3 = GM_S, GM_S = NE^2),
geocentric direction = (r_C - r_Earth)/|...|, Earth orbit as in stage 1. Elements below are whole-record
least-squares fits (rms residual 97-2141 arcmin; see ANSWER_2.md). They are low confidence.
"""
import csv
from pathlib import Path
import numpy as np
import predict as P1law

HERE = Path(__file__).resolve().parent
NE = P1law.NE
# C-body elements: (e, inc, Omega, omega, M0, n) with a = (NE^2/n^2)^(1/3) AU
FAINT = {'C1': [0.69679862, 0.3309978, -3.10201606, 3.97734175, 19.82027269, 0.0037193326], 'C2': [0.43880908, 0.53390189, 3.20522991, 3.46199248, 5.59633142, 0.0090013704], 'C3': [0.38318335, 0.46342128, 3.47682143, 1.51716096, 3.52696068, 0.0041931972], 'C4': [0.9499999, 0.39532786, 3.13863478, 2.81306176, 0.62662889, 0.0223294325], 'C5': [0.72576398, 0.34010281, 3.25474514, 6.10190518, -3.977169, 0.0020501273], 'C6': [0.84411663, 0.14239771, 3.60441812, 4.58547508, 10.92150642, 0.0018915365], 'C7': [0.95, 2.07980874, 6.89899352, 2.84149293, 5.0219562, 0.077152483], 'C8': [0.57882107, 0.32593856, 2.09647712, 6.73468809, 3.17944894, 0.0035630486], 'C9': [0.2465148, -0.61748379, 6.07931235, 3.07008644, 8.00059569, 0.006976268], 'C10': [0.59048066, 0.51140461, 2.68889744, 1.68828724, 5.94754824, 0.0051146222], 'C11': [0.51584276, -0.34394813, 6.1877982, 1.67333111, 7.20041702, 0.0032083143], 'C12': [0.73845281, 0.60753148, 3.46924345, 0.73553935, 7.72846783, 0.0006824949]}

def faint_direction(t, body):
    e, inc, Om, om, M0, n = FAINT[body]
    a = (NE**2/n**2)**(1/3)
    rC = P1law.kepler_pos(t, a, e, inc, Om, om, M0, n) - P1law.earth(t)
    return rC/np.linalg.norm(rC, axis=1)[:, None]

def main():
    rows = list(csv.DictReader(open(HERE/'predictions_2.csv')))
    days = np.array([float(r['day']) for r in rows])
    bodies = [r['body'] for r in rows]
    lon_out = np.full(len(rows), np.nan); lat_out = np.full(len(rows), np.nan)
    for body in sorted(set(bodies)):
        idx = np.array([i for i, b in enumerate(bodies) if b == body], dtype=int)
        if body in FAINT:
            u = faint_direction(days[idx], body)
            lon, lat = P1law.unit_to_lonlat(u)
        else:
            lon, lat = P1law.predict_body(days[idx], body)
        lon_out[idx] = np.degrees(lon); lat_out[idx] = np.degrees(lat)
    with open(HERE/'predictions_2.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['day', 'body', 'lon', 'lat'])
        for i, r in enumerate(rows):
            w.writerow([r['day'], r['body'], '%.5f' % lon_out[i], '%.5f' % lat_out[i]])
    print('wrote', len(rows), 'rows')

if __name__ == '__main__':
    main()
