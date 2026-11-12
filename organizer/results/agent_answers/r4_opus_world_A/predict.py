"""predict.py -- predictions for Tellus sky, years 101-120.
Uses only the laws and constants stated in ANSWER_1.md, section 2:
  Law 1: frame. Ecliptic frame = observation frame rotated about its x-axis by EPS.
  Law 2: Newtonian gravity, S + Tellus + P1..P5 (heliocentric equations, point masses), integrated from day 0.
  Law 3: M moves about Tellus under Tellus's gravity plus the tidal pull of S (S on a two-body Kepler ellipse).
  Directions are geometric (no light-time, no aberration).
Run:  python predict.py   (fills predictions.csv)
"""
import numpy as np, math, csv
from scipy.integrate import solve_ivp

EPS_DEG = 21.7015                 # tilt of ecliptic wrt observation frame
# ---- Law 2 constants (AU = Tellus' mean distance scale, day) ----
MU_S = 0.00029584723885713716                     # GM_S  [AU^3/day^2]
GM   = [0.0, 3.5271359430802e-10, 3.0661928691825013e-10, 4.4889124743997165e-08, 2.9615363068095e-07, 7.369887960352512e-08]                     # GM of T,P1,P2,P3,P4,P5 [AU^3/day^2]
X0   = [
        [-0.3942779278238853, -0.9382581801513937, -2.5053241296208538e-05],
        [-0.312556118615529, -0.30167558290187146, -3.5397561378137886e-06],
        [-0.7148486175821397, 1.63302504774436, 4.488962330249993e-05],
        [-3.0564433197627046, 1.364972190640158, 0.05306276709052556],
        [6.460515045871039, -0.028502913301657325, -0.0320054745096144],
        [3.323334756704969, 11.027423989465358, 0.00017446415216609547]]                     # heliocentric ecliptic positions at day 0 [AU]
V0   = [
        [0.015517441513077705, -0.006687382966186439, -1.8402022073027232e-07],
        [0.01747725804575744, -0.019979054544725188, -1.9703766730462908e-07],
        [-0.011452809362836876, -0.005168455230731659, -9.18862715840567e-08],
        [-0.004117339885764291, -0.008373846951401665, 0.00019967730844658787],
        [0.000210851033996485, 0.00662642303846803, -0.0001116489338242133],
        [-0.0049018588800830805, 0.0012657899905459936, 0.00017668096914215606]]                     # heliocentric ecliptic velocities at day 0 [AU/day]
# ---- Law 3 constants ----
S_N, S_M0, S_E, S_W = 0.0172020684, 152.261368, 0.019979018, 273.909777    # S's apparent Kepler ellipse: mean motion [rad/day], mean anomaly at day 0 [deg], e, long. of perigee [deg]
A_M   = 0.0025                    # length unit of M orbit (arbitrary: directions are scale-free)
MU_T  = 8.263505605559067e-10                   # GM_T in those units (fixes M's period)
M_X0  = [0.20764289510293601, -1.0312027018099035, -0.026270639095353526]                   # M position at day 0, units of A_M
M_V0  = [0.9241031718436921, 0.19957394220291366, 0.07958436677065049]                   # M velocity at day 0, units of 0.23*A_M per day

def rotx(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])

def acc(r, mu, m):
    d = r[None, :, :] - r[:, None, :]
    dist3 = np.sum(d * d, -1) ** 1.5
    np.fill_diagonal(dist3, np.inf)
    rn3 = np.sum(r * r, -1) ** 1.5
    a = -(mu + m)[:, None] * r / rn3[:, None]
    a += np.einsum('j,ijk->ik', m, d / dist3[..., None])
    a -= np.einsum('j,jk->k', m, r / rn3[:, None])[None, :]
    return a

def nbody(ndays, h=0.25):
    x, v = np.array(X0, float), np.array(V0, float)
    mu, m = MU_S, np.array(GM, float)
    out = np.empty((ndays + 1, 6, 3)); out[0] = x
    k = int(round(1 / h))
    for s in range(1, ndays + 1):
        for _ in range(k):
            a1 = acc(x, mu, m); x2 = x + 0.5*h*v; v2 = v + 0.5*h*a1
            a2 = acc(x2, mu, m); x3 = x + 0.5*h*v2; v3 = v + 0.5*h*a2
            a3 = acc(x3, mu, m); x4 = x + h*v3; v4 = v + h*a3
            a4 = acc(x4, mu, m)
            x = x + h/6*(v + 2*v2 + 2*v3 + v4); v = v + h/6*(a1 + 2*a2 + 2*a3 + a4)
        out[s] = x
    return out

def kep_pos(t):
    M = math.radians(S_M0) + S_N * t
    E = M + S_E * math.sin(M)
    for _ in range(30): E -= (E - S_E*math.sin(E) - M) / (1 - S_E*math.cos(E))
    nu = 2*math.atan2(math.sqrt(1+S_E)*math.sin(E/2), math.sqrt(1-S_E)*math.cos(E/2))
    r = 1 - S_E*math.cos(E); u = nu + math.radians(S_W)
    return np.array([r*math.cos(u), r*math.sin(u), 0.0])

def moon(ndays):
    mus = S_N**2
    R0 = kep_pos(0.0); VR0 = (kep_pos(1e-3) - kep_pos(-1e-3)) / 2e-3
    def rhs(t, y):
        rx, ry, rz, Rx, Ry, Rz = y[0], y[1], y[2], y[6], y[7], y[8]
        r3 = (rx*rx + ry*ry + rz*rz)**1.5
        dx, dy, dz = Rx-rx, Ry-ry, Rz-rz
        d3 = (dx*dx + dy*dy + dz*dz)**1.5; R3 = (Rx*Rx + Ry*Ry + Rz*Rz)**1.5
        return [y[3], y[4], y[5],
                -MU_T*rx/r3 + mus*(dx/d3 - Rx/R3), -MU_T*ry/r3 + mus*(dy/d3 - Ry/R3), -MU_T*rz/r3 + mus*(dz/d3 - Rz/R3),
                y[9], y[10], y[11], -mus*Rx/R3, -mus*Ry/R3, -mus*Rz/R3]
    y0 = np.concatenate([np.array(M_X0)*A_M, np.array(M_V0)*A_M*0.23, R0, VR0])
    s = solve_ivp(rhs, (0, ndays), y0, method='DOP853', rtol=1e-11, atol=1e-15, t_eval=np.arange(0, ndays + 1))
    return s.y[:3].T

def to_obs(g):
    g = g / np.linalg.norm(g, axis=-1, keepdims=True)
    g = g @ rotx(-math.radians(EPS_DEG))      # ecliptic -> observation frame
    lon = np.degrees(np.arctan2(g[..., 1], g[..., 0])) % 360
    lat = np.degrees(np.arcsin(np.clip(g[..., 2], -1, 1)))
    return lon, lat

def main():
    rows = list(csv.DictReader(open('predictions.csv')))
    days = np.array([int(r['day']) for r in rows])
    nd = int(days.max()) + 1
    P = nbody(nd)
    Mg = moon(nd)
    names = ['T', 'P1', 'P2', 'P3', 'P4', 'P5']
    for r, d in zip(rows, days):
        b = r['body']
        if b == 'S': g = -P[d, 0]
        elif b == 'M': g = Mg[d]
        else: g = P[d, names.index(b)] - P[d, 0]
        lo, la = to_obs(g)
        r['lon'] = f'{lo:.5f}'; r['lat'] = f'{la:.5f}'
    with open('predictions.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['day', 'body', 'lon', 'lat']); w.writeheader(); w.writerows(rows)
    print('wrote', len(rows), 'predictions')

if __name__ == '__main__':
    main()
