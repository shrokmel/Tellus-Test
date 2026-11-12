"""Stage 3 predictor: geocentric direction (lon, lat in degrees) of body N for predictions_3.csv.

Law 1 (Tellus): Kepler orbit about the Sun; Sun seen from Tellus has longitude varpi_S + nu.
Law 11 (N): heliocentric Kepler orbit about the Sun, a = (P/P_E)^(2/3) AU, elements below.
Direction = normalised (r_N - r_Tellus), converted to the data frame (pole N, basis EX, EY).
Time t = days since observations began.
"""
import csv
import numpy as np

HERE = '.'

NPOLE = np.array([-3.13826576e-05, 3.69769974e-01, 9.29123332e-01]); NPOLE /= np.linalg.norm(NPOLE)
EX = np.array([0.0, 0.0, 1.0]); EX = EX - np.dot(EX, NPOLE) * NPOLE; EX /= np.linalg.norm(EX)
EY = np.cross(NPOLE, EX)

P_E = 365.257550
E_E = 0.019980
VARPI_S = np.radians(3.9069)
M0_S = np.radians(152.2598)

# Law 11: body N (best member of the family fitted to observations_3.csv, rms 4011 arcsec)
N_EL = dict(P=97997.6, e=0.9428,
            i=np.radians(52.880), Om=np.radians(105.144), om=np.radians(70.469), M0=np.radians(226.357))


def kepler_E(M, e):
    M = np.mod(M, 2 * np.pi)
    E = M + e * np.sin(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def nu_E(M, e):
    E = kepler_E(M, e)
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2)), E


def tellus_pos(t):
    nu, _ = nu_E(2 * np.pi * t / P_E + M0_S, E_E)
    R = (1 - E_E ** 2) / (1 + E_E * np.cos(nu))
    lam = VARPI_S + nu + np.pi
    return np.array([R * np.cos(lam), R * np.sin(lam), 0.0])


def body_helio(el, t):
    P = el['P']; e = el['e']
    a = (P / P_E) ** (2 / 3)
    nu, E = nu_E(2 * np.pi * t / P + el['M0'], e)
    r = a * (1 - e * np.cos(E))
    u = el['om'] + nu
    Om, i = el['Om'], el['i']
    x = np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(i)
    y = np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(i)
    z = np.sin(u) * np.sin(i)
    return r * np.array([x, y, z])


def direction(t):
    p = body_helio(N_EL, t) - tellus_pos(t)
    p = p / np.linalg.norm(p)
    v = p[0] * EX + p[1] * EY + p[2] * NPOLE
    lon = np.degrees(np.arctan2(v[1], v[0])) % 360
    lat = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    return lon, lat


def main():
    path = HERE + '/predictions_3.csv'
    with open(path) as f:
        rows = list(csv.reader(f))
    header, body_rows = rows[0], rows[1:]
    out = [header]
    for day, body, *_ in body_rows:
        lon, lat = direction(float(day))
        out.append([day, body, '%.5f' % lon, '%.5f' % lat])
    with open(path, 'w', newline='') as f:
        csv.writer(f).writerows(out)
    print('wrote', len(out) - 1, 'rows to', path)


if __name__ == '__main__':
    main()
