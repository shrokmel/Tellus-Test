"""Predict geocentric sky directions (lon, lat, degrees) for predictions.csv.

Laws (see ANSWER_1.md, section 2). Frame: ecliptic-like frame with pole N;
lon/lat are the frame's longitude/latitude. Time t in days since observations began.
Body directions are computed as normalised (heliocentric body - Tellus) vectors,
using Kepler orbits with the constants below. Tellus's orbit is set by the Sun's
apparent motion (law 1).
"""
import csv
import numpy as np

HERE = '.'

# Frame: pole N (ecliptic pole in the data frame), basis EX, EY
N = np.array([-3.13826576e-05, 3.69769974e-01, 9.29123332e-01]); N /= np.linalg.norm(N)
EX = np.array([0.0, 0.0, 1.0]); EX = EX - np.dot(EX, N) * N; EX /= np.linalg.norm(EX)
EY = np.cross(N, EX)

# Law 1 (Tellus / Sun): Sun geocentric longitude = VARPI_S + nu(M_S), M_S = 2 pi t / P_E + M0_S
P_E = 365.257550
E_E = 0.019980
VARPI_S = np.radians(3.9069)
M0_S = np.radians(152.2598)

# Law 2 (P1, planar, inner orbit): heliocentric Kepler, a = (P/P_E)^(2/3) AU
P1 = dict(P=110.2578, e=0.06009, varpi=np.radians(256.468), M0=np.radians(51.841), i=0.0)
# Laws 3-6 (P2-P5), 3D heliocentric Kepler (node Om, arg. of peri om, inclination i)
PLAN = {
    'P2': dict(P=809.67, e=0.0501, i=np.radians(0.006), Om=np.radians(306.279), om=np.radians(65.598), M0=np.radians(193.013)),
    'P3': dict(P=2189.15, e=0.0421, i=np.radians(1.493), Om=np.radians(208.253), om=np.radians(284.133), M0=np.radians(109.247)),
    'P4': dict(P=5639.82, e=0.0493, i=np.radians(1.003), Om=np.radians(253.056), om=np.radians(49.059), M0=np.radians(144.474)),
    'P5': dict(P=14169.46, e=0.0413, i=np.radians(1.999), Om=np.radians(162.952), om=np.radians(98.863), M0=np.radians(266.175)),
}
# Law 7 (M, Moon about Tellus): Kepler elements from the last 60-day fit (low confidence)
MOON = dict(P=27.3302, e=0.0647, i=np.radians(4.927), Om=np.radians(281.24), om=np.radians(291.33), M0=np.radians(147.04))


def kepler_E(M, e):
    E = M + e * np.sin(M)
    for _ in range(30):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def nu_of(M, e):
    E = kepler_E(np.mod(M, 2 * np.pi), e)
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2)), E


def tellus_pos(t):
    """Tellus heliocentric position, ecliptic plane (x, y)."""
    nu, _ = nu_of(2 * np.pi * t / P_E + M0_S, E_E)
    R = (1 - E_E ** 2) / (1 + E_E * np.cos(nu))
    lam = VARPI_S + nu + np.pi
    return np.array([R * np.cos(lam), R * np.sin(lam), 0.0])


def body_helio(elem, t, planar):
    P = elem['P']; e = elem['e']; M0 = elem['M0']
    a = (P / P_E) ** (2 / 3)
    nu, E = nu_of(2 * np.pi * t / P + M0, e)
    r = a * (1 - e * np.cos(E))
    if planar:
        u = elem['varpi'] + nu
        return np.array([r * np.cos(u), r * np.sin(u), 0.0])
    om, Om, i = elem['om'], elem['Om'], elem['i']
    u = om + nu
    x = np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(i)
    y = np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(i)
    z = np.sin(u) * np.sin(i)
    return r * np.array([x, y, z])


def direction(name, t):
    if name == 'S':
        p = -tellus_pos(t)                 # Sun seen from Tellus
    elif name == 'P1':
        p = body_helio(P1, t, planar=True) - tellus_pos(t)
    elif name in PLAN:
        p = body_helio(PLAN[name], t, planar=False) - tellus_pos(t)
    elif name == 'M':
        # Moon about Tellus: Kepler orbit, scaled to a unit-normal direction
        p = body_helio(MOON, t, planar=False) * 0.0
        p = kepler_geo_moon(t)
    else:
        raise ValueError(name)
    p = p / np.linalg.norm(p)
    # back to data frame: v = x*EX + y*EY + z*N
    v = p[0] * EX + p[1] * EY + p[2] * N
    lon = np.degrees(np.arctan2(v[1], v[0])) % 360
    lat = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    return lon, lat


def kepler_geo_moon(t):
    el = MOON
    P = el['P']; e = el['e']; M0 = el['M0']
    nu, _ = nu_of(2 * np.pi * t / P + M0, e)
    _, E = nu_of(2 * np.pi * t / P + M0, e)
    r = 1 - e * np.cos(E)
    u = el['om'] + nu
    Om, i = el['Om'], el['i']
    x = np.cos(Om) * np.cos(u) - np.sin(Om) * np.sin(u) * np.cos(i)
    y = np.sin(Om) * np.cos(u) + np.cos(Om) * np.sin(u) * np.cos(i)
    z = np.sin(u) * np.sin(i)
    return np.array([x, y, z])


def main():
    path = HERE + '/predictions.csv'
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
