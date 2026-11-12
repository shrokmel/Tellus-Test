#!/usr/bin/env python
"""Predict directions of S, P1-P5, M for the days in predictions.csv.

Model (see ANSWER_1.md):
  * Newtonian gravity, inverse square, point masses, in a frame fixed against the stars.
  * Bodies: S (central mass, GM_S), Tellus (T, massless here), P2..P5 (masses m_i), P1 (test particle),
    M (test particle orbiting Tellus, perturbed by S).
  * Coordinates: ecliptic frame = observation frame rotated about its x axis (lon=0 direction) by
    the obliquity EPS, so Tellus's orbital plane is z=0.
  * Unit of length: Tellus's semi-major axis ("1"); unit of time: day.  GM_S = n_T^2.
Only directions are observable, hence all distances are relative.
"""
import csv, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = np.radians(21.70197)          # obliquity: tilt of Tellus's orbital plane to the reference frame

# ---- fitted constants -----------------------------------------------------------------
# Tellus: n (rad/day), L0 (mean longitude at day 0, rad), e, varpi (rad); plane inclination 0
TELLUS = [0.017202502976248157, 4.2964814061297885, 0.020001184964697456, 1.6369298102422518]
# P2..P5: n, L0, e, varpi, Omega, inc (rad)
PLAN = [
    [0.0077609526298878595, 2.0057101958932346, 0.04999007289141734, 4.908015035389737, -0.18433095412438902, 2.0159414988101992e-05],
    [0.0028698302276058004, 2.6474444782569666, 0.04000594442189176, 7.0255514039962526, 2.0708399848320824, 0.026165261780251775],
    [0.001114854751662426, -0.05871846665967487, 0.04999168508482317, 3.691535430603835, 2.8488737322023137, 0.01742121725040483],
    [0.0004411474753363703, 1.3580005471262364, 0.03999691651623015, 2.9260598244968796, 1.277432226479271, 0.034915180467916127],
]
MASS = [0.0, 0.0001514640239296227, 0.0010003317475393997, 0.000249113903835644]   # m_i / M_S (P2 consistent with 0)
# P1 (test particle): n, L0, e, varpi, i*cos(Omega), i*sin(Omega)
P1 = [0.056986389912965105, 3.8104039852166283, 0.059998825528430985, 2.9052276013112475, -8.816066636058094e-06, -8.462265481790739e-06]
# M (test particle around Tellus): n, L0, e, varpi, Omega, inc  (osculating at day 0, geocentric)
MOON = [0.23009800810695744, 4.941742626381729, 0.05499186728022886, 7.771667738630799, 2.0588717080054004, -0.08725914774578131]
MU_M = 3.0e-6   # GM of Tellus in units of GM_S (only sets the unobservable length scale of M's orbit)

GMS = TELLUS[0] ** 2


def kepler(M, e):
    E = M + e * np.sin(M)
    for _ in range(30):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def elements_to_state(n, L0, e, vp, Om, inc, gm):
    a = (gm / n ** 2) ** (1 / 3)
    E = kepler(np.array([L0 - vp]), e)[0]
    x = a * (np.cos(E) - e)
    y = a * np.sqrt(1 - e * e) * np.sin(E)
    f = n * a / (1 - e * np.cos(E))
    vx, vy = -f * np.sin(E), f * np.sqrt(1 - e * e) * np.cos(E)
    w = vp - Om
    cw, sw, ci, si, co, so = np.cos(w), np.sin(w), np.cos(inc), np.sin(inc), np.cos(Om), np.sin(Om)

    def rot(X, Y):
        X1 = X * cw - Y * sw
        Y1 = X * sw + Y * cw
        return np.array([X1 * co - Y1 * ci * so, X1 * so + Y1 * ci * co, Y1 * si])
    return rot(x, y), rot(vx, vy)


def nbody(teval):
    """Sun, Tellus, P2..P5 (barycentric); returns positions (6,3,T)."""
    gm = np.r_[GMS, 0.0, GMS * np.array(MASS)]
    gm = np.r_[GMS, 0.0, GMS * np.abs(np.array(MASS))]
    st = [np.zeros(6)]
    r, v = elements_to_state(TELLUS[0], TELLUS[1], TELLUS[2], TELLUS[3], 0.0, 0.0, GMS)
    st.append(np.r_[r, v])
    for p, m in zip(PLAN, MASS):
        r, v = elements_to_state(p[0], p[1], p[2], p[3], p[4], p[5], GMS * (1 + abs(m)))
        st.append(np.r_[r, v])
    st = np.array(st)
    tot = gm.sum()
    st[:, :3] -= (gm[:, None] * st[:, :3]).sum(0) / tot
    st[:, 3:] -= (gm[:, None] * st[:, 3:]).sum(0) / tot
    N = len(gm)

    def rhs(t, y):
        s = y.reshape(N, 6)
        r = s[:, :3]
        d = r[None, :, :] - r[:, None, :]
        d2 = (d ** 2).sum(2) + np.eye(N)
        inv3 = d2 ** -1.5
        np.fill_diagonal(inv3, 0.0)
        acc = ((gm[None, :] * inv3)[:, :, None] * d).sum(1)
        out = np.empty_like(s)
        out[:, :3] = s[:, 3:]
        out[:, 3:] = acc
        return out.ravel()
    sol = solve_ivp(rhs, (0, teval[-1]), st.ravel(), method='DOP853', t_eval=teval, rtol=1e-10, atol=1e-13)
    return sol.y.reshape(N, 6, -1)[:, :3, :]


def p1_track(teval, helio, tg):
    gmj = GMS * np.abs(np.array(MASS))
    spl = [CubicSpline(tg, helio[k].T) for k in range(2, 6)]
    n, L0, e, vp, pp, qq = P1
    r, v = elements_to_state(n, L0, e, vp, np.arctan2(qq, pp), np.hypot(pp, qq), GMS)

    def rhs(t, y):
        r = y[:3]
        acc = -GMS * r / (r @ r) ** 1.5
        for g, s in zip(gmj, spl):
            rj = s(t)
            d = rj - r
            acc = acc + g * (d / (d @ d) ** 1.5 - rj / (rj @ rj) ** 1.5)
        return np.r_[y[3:], acc]
    sol = solve_ivp(rhs, (0, teval[-1]), np.r_[r, v], method='DOP853', t_eval=teval, rtol=1e-10, atol=1e-13)
    return sol.y[:3]


def moon_track(teval, sun_geo_spline):
    mu = MU_M * GMS
    n, L0, e, vp, Om, inc = MOON
    r, v = elements_to_state(n, L0, e, vp, Om, inc, mu)

    def rhs(t, y):
        r = y[:3]
        s = sun_geo_spline(t)
        d = s - r
        acc = -mu * r / (r @ r) ** 1.5 + GMS * (d / (d @ d) ** 1.5 - s / (s @ s) ** 1.5)
        return np.r_[y[3:], acc]
    sol = solve_ivp(rhs, (0, teval[-1]), np.r_[r, v], method='DOP853', t_eval=teval, rtol=1e-10, atol=1e-14)
    return sol.y[:3]


def to_frame(v):
    """ecliptic vectors (3,T) -> (lon, lat) deg in the observation frame."""
    c, s = np.cos(EPS), np.sin(EPS)
    x, y, z = v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]
    return np.degrees(np.arctan2(y, x)) % 360.0, np.degrees(np.arctan2(z, np.hypot(x, y)))


def ephemeris(tmax):
    tg = np.arange(0, tmax + 1.0)
    pos = nbody(tg)
    helio = pos - pos[0][None]
    earth = helio[1]
    out = {}
    out['S'] = to_frame(-earth)
    for k, name in zip(range(2, 6), ['P2', 'P3', 'P4', 'P5']):
        out[name] = to_frame(helio[k] - earth)
    out['P1'] = to_frame(p1_track(tg, helio, tg)[:, :] - earth)
    sun_geo = CubicSpline(tg, (-earth).T)
    out['M'] = to_frame(moon_track(tg, sun_geo))
    return out


def main():
    rows = list(csv.reader(open(os.path.join(HERE, 'predictions.csv'))))
    head, body = rows[0], rows[1:]
    days = np.array([int(float(r[0])) for r in body])
    eph = ephemeris(int(days.max()))
    res = []
    for r, d in zip(body, days):
        lon, lat = eph[r[1]]
        res.append([r[0], r[1], '%.5f' % lon[d], '%.5f' % lat[d]])
    with open(os.path.join(HERE, 'predictions.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(head)
        w.writerows(res)
    print('wrote', len(res), 'predictions')


if __name__ == '__main__':
    main()
