"""Stage 2 predictor: geocentric sky directions (lon, lat in degrees) for predictions_2.csv.

Laws (ANSWER_2.md). Frame: pole N, basis EX, EY. Time t = days since observations began.
Tellus: Kepler orbit about the Sun (law 1). Every body: heliocentric Kepler orbit about the Sun
(law 2, a = (P/P_E)^(2/3) AU); direction = normalised (body - Tellus).
"""
import csv
import numpy as np

HERE = '.'

N = np.array([-3.13826576e-05, 3.69769974e-01, 9.29123332e-01]); N /= np.linalg.norm(N)
EX = np.array([0.0, 0.0, 1.0]); EX = EX - np.dot(EX, N) * N; EX /= np.linalg.norm(EX)
EY = np.cross(N, EX)

P_E = 365.257550
E_E = 0.019980
VARPI_S = np.radians(3.9069)
M0_S = np.radians(152.2598)

P1 = dict(P=110.2578, e=0.06009, varpi=np.radians(256.468), M0=np.radians(51.841))
PLAN = {
    'P2': dict(P=809.67, e=0.0501, i=np.radians(0.006), Om=np.radians(306.279), om=np.radians(65.598), M0=np.radians(193.013)),
    'P3': dict(P=2189.15, e=0.0421, i=np.radians(1.493), Om=np.radians(208.253), om=np.radians(284.133), M0=np.radians(109.247)),
    'P4': dict(P=5639.82, e=0.0493, i=np.radians(1.003), Om=np.radians(253.056), om=np.radians(49.059), M0=np.radians(144.474)),
    'P5': dict(P=14169.46, e=0.0413, i=np.radians(1.999), Om=np.radians(162.952), om=np.radians(98.863), M0=np.radians(266.175)),
}
MOON = dict(P=27.3302, e=0.0647, i=np.radians(4.927), Om=np.radians(281.24), om=np.radians(291.33), M0=np.radians(147.04))

# Stage 2 faint bodies: heliocentric Kepler elements from the last 1000-day window (refit, unbounded phase)
FAINT = {
    'C1':  dict(P=724.302,   e=0.3448, i=1.044,  Om=52.689,  om=134.234, M0=339.055),
    'C2':  dict(P=8000.422,  e=0.7853, i=0.263,  Om=281.795, om=10.841,  M0=105.988),
    'C3':  dict(P=5818.968,  e=0.7489, i=9.416,  Om=327.985, om=38.291,  M0=258.387),
    'C4':  dict(P=1901.693,  e=0.7486, i=1.035,  Om=262.963, om=162.800, M0=77.683),
    'C5':  dict(P=3204.925,  e=0.4060, i=18.721, Om=47.908,  om=0.192,   M0=293.443),
    'C6':  dict(P=3404.510,  e=0.3468, i=15.953, Om=56.449,  om=323.809, M0=295.620),
    'C7':  dict(P=889.605,   e=0.5723, i=17.320, Om=313.815, om=345.001, M0=207.937),
    'C8':  dict(P=3361.227,  e=0.7238, i=17.053, Om=139.948, om=100.881, M0=70.979),
    'C9':  dict(P=3830.623,  e=0.7394, i=13.913, Om=239.611, om=18.852,  M0=47.425),
    'C10': dict(P=3360.633,  e=0.8114, i=14.092, Om=200.639, om=138.101, M0=337.329),
    'C11': dict(P=2156.209,  e=0.5721, i=2.623,  Om=123.460, om=58.790,  M0=278.847),
    'C12': dict(P=15890.685, e=0.8409, i=14.663, Om=312.128, om=356.230, M0=294.514),
}
for _d in FAINT.values():
    for _k in ('i', 'Om', 'om', 'M0'):
        _d[_k] = np.radians(_d[_k])


def kepler_E(M, e):
    M = np.mod(M, 2 * np.pi)
    E = M + e * np.sin(M)
    for _ in range(40):
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


def helio(el, t, planar=False):
    P = el['P']; e = el['e']
    a = (P / P_E) ** (2 / 3)
    nu, E = nu_E(2 * np.pi * t / P + el['M0'], e)
    r = a * (1 - e * np.cos(E))
    if planar:
        u = el['varpi'] + nu
        return np.array([r * np.cos(u), r * np.sin(u), 0.0])
    u = el['om'] + nu
    Om, i = el['Om'], el['i']
    x = np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(i)
    y = np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(i)
    z = np.sin(u) * np.sin(i)
    return r * np.array([x, y, z])


def geo_vec(name, t):
    if name == 'S':
        return -tellus_pos(t)
    if name == 'P1':
        return helio(P1, t, planar=True) - tellus_pos(t)
    if name in PLAN:
        return helio(PLAN[name], t) - tellus_pos(t)
    if name == 'M':
        # Moon: Kepler orbit about Tellus (law 8 in stage 1), unit direction
        el = MOON
        nu, E = nu_E(2 * np.pi * t / el['P'] + el['M0'], el['e'])
        u = el['om'] + nu
        Om, i = el['Om'], el['i']
        x = np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(i)
        y = np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(i)
        z = np.sin(u) * np.sin(i)
        return np.array([x, y, z])
    if name in FAINT:
        return helio(FAINT[name], t) - tellus_pos(t)
    raise ValueError(name)


def direction(name, t):
    p = geo_vec(name, t)
    p = p / np.linalg.norm(p)
    v = p[0] * EX + p[1] * EY + p[2] * N
    lon = np.degrees(np.arctan2(v[1], v[0])) % 360
    lat = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    return lon, lat


def main():
    path = HERE + '/predictions_2.csv'
    with open(path) as f:
        rows = list(csv.reader(f))
    header, body_rows = rows[0], rows[1:]
    out = [header]
    for day, body, *_ in body_rows:
        lon, lat = direction(body, float(day))
        out.append([day, body, '%.5f' % lon, '%.5f' % lat])
    with open(path, 'w', newline='') as f:
        csv.writer(f).writerows(out)
    print('wrote', len(out) - 1, 'rows to', path)


if __name__ == '__main__':
    main()
