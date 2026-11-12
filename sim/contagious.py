"""Contagious gravity in plain Python (reference implementation; sky.c is the fast version).

Every body has a mass m (inertia, fixed) and a charge q (takes part in gravity).
Force between i and j: G q_i q_j / r^2, so the acceleration of i is G q_j (q_i/m_i) / r^2.
Pull factor f = q/m (in our universe f = 1 for everything).
Charge flow: within a contact distance d, charge flows from the higher f to the lower f.
    dq_i/dt = -lam * sum_j w(r_ij) * mu_ij * (f_i - f_j),  mu_ij = m_i m_j / (m_i + m_j),
    w(r) = exp(-(r/d)^2)   (smooth contact, so the integrator is not hit by a jump)
For an isolated pair, the gap in f shrinks as d(f_i - f_j)/dt = -lam * w * (f_i - f_j); total charge is conserved.
Units: AU, years, solar masses; G = 4 pi^2.

Run:  python contagious.py --lam 50 --out run_contagious.npz
      python contagious.py --lam 0  --out run_newton.npz      (same start, no charge flow)
"""
import argparse, time
import numpy as np
from scipy.integrate import solve_ivp

G = 4 * np.pi ** 2


def elements_to_state(a, e, inc, node, peri, M, mu):
    """Kepler elements (angles in radians) -> position, velocity around a centre with parameter mu."""
    E = M
    for _ in range(50):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e * np.cos(E))
    p = a * (1 - e ** 2)
    x_orb = np.array([r * np.cos(nu), r * np.sin(nu), 0.0])
    v_orb = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    cO, sO, ci, si, cw, sw = np.cos(node), np.sin(node), np.cos(inc), np.sin(inc), np.cos(peri), np.sin(peri)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return R @ x_orb, R @ v_orb


def build_system(n_comets, seed, comet_mass=1e-12):
    rng = np.random.default_rng(seed)
    names, m, f, X, V = ["S"], [1.0], [1.0], [np.zeros(3)], [np.zeros(3)]
    # invented planets: name, a (AU), e, inclination (deg), mass (solar masses)
    planets = [("P1", 0.45, 0.06, 0.0, 4e-7), ("T", 1.0, 0.02, 0.0, 3e-6), ("P2", 1.7, 0.05, 0.0, 6e-7),
               ("P3", 3.3, 0.04, 1.5, 1.5e-4), ("P4", 6.2, 0.05, 1.0, 1.0e-3), ("P5", 11.5, 0.04, 2.0, 2.5e-4)]
    for name, a, e, inc, mass in planets:
        x, v = elements_to_state(a, e, np.radians(inc), rng.uniform(0, 2 * np.pi), rng.uniform(0, 2 * np.pi),
                                 rng.uniform(0, 2 * np.pi), G * (1 + mass))
        names.append(name); m.append(mass); f.append(1.0); X.append(x); V.append(v)
    # comets: stretched orbits crossing the planets, each with its own pull factor
    for k in range(n_comets):
        q_peri = rng.uniform(0.3, 4.0)
        a = rng.uniform(max(q_peri * 1.6, 2.0), 14.0)
        e = 1 - q_peri / a
        fk = rng.uniform(0.7, 1.3)
        x, v = elements_to_state(a, e, np.radians(rng.uniform(0, 20)), rng.uniform(0, 2 * np.pi),
                                 rng.uniform(0, 2 * np.pi), rng.uniform(0, 2 * np.pi), G * fk)
        names.append(f"C{k + 1}"); m.append(comet_mass); f.append(fk); X.append(x); V.append(v)
    m, f = np.array(m), np.array(f)
    X, V = np.array(X), np.array(V)
    # put the centre of mass at rest at the origin
    X -= (m[:, None] * X).sum(0) / m.sum()
    V -= (m[:, None] * V).sum(0) / m.sum()
    return names, m, f * m, X, V


def make_rhs(m, lam, d):
    N = len(m)
    mu = m[:, None] * m[None, :] / (m[:, None] + m[None, :])
    eye = np.eye(N, dtype=bool)

    def rhs(t, y):
        X = y[:3 * N].reshape(N, 3); V = y[3 * N:6 * N].reshape(N, 3); q = y[6 * N:]
        D = X[None, :, :] - X[:, None, :]            # D[i, j] = x_j - x_i
        r2 = (D ** 2).sum(-1); r2[eye] = 1.0
        inv_r3 = r2 ** -1.5; inv_r3[eye] = 0.0
        A = G * (q / m)[:, None] * ((q[None, :] * inv_r3)[:, :, None] * D).sum(1)
        if lam > 0:
            f = q / m
            w = np.exp(-r2 / d ** 2); w[eye] = 0.0
            dq = -lam * (w * mu * (f[:, None] - f[None, :])).sum(1)
        else:
            dq = np.zeros(N)
        return np.concatenate([V.ravel(), A.ravel(), dq])
    return rhs


def energy(m, q, X, V):
    N = len(m)
    kin = 0.5 * (m * (V ** 2).sum(1)).sum()
    D = X[None, :, :] - X[:, None, :]
    r = np.sqrt((D ** 2).sum(-1)); iu = np.triu_indices(N, 1)
    pot = -G * (q[:, None] * q[None, :])[iu] / r[iu]
    return kin + pot.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", type=float, default=50.0, help="evening rate during contact (1/yr); 0 = plain Newton")
    ap.add_argument("--d", type=float, default=0.1, help="contact distance (AU)")
    ap.add_argument("--years", type=float, default=100.0)
    ap.add_argument("--comets", type=int, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--comet-mass", type=float, default=1e-12)
    ap.add_argument("--sample-days", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    names, m, q, X, V = build_system(a.comets, a.seed, a.comet_mass)
    N = len(m)
    y0 = np.concatenate([X.ravel(), V.ravel(), q])
    atol = np.concatenate([np.full(6 * N, 1e-11), 1e-11 * m])
    t_eval = np.arange(0, a.years + 1e-9, a.sample_days / 365.25)
    t0 = time.time()
    sol = solve_ivp(make_rhs(m, a.lam, a.d), [0, a.years], y0, method="DOP853", rtol=1e-10, atol=atol,
                    t_eval=t_eval)
    print(f"{sol.message}  steps={sol.t.size} evals={sol.nfev}  {time.time() - t0:.0f}s")
    Y = sol.y.T
    Xs = Y[:, :3 * N].reshape(-1, N, 3); Vs = Y[:, 3 * N:6 * N].reshape(-1, N, 3); Qs = Y[:, 6 * N:]
    E = np.array([energy(m, Qs[k], Xs[k], Vs[k]) for k in range(0, len(sol.t))])
    np.savez_compressed(a.out, t=sol.t, X=Xs.astype(np.float64), V=Vs, Q=Qs, m=m, names=np.array(names),
                        E=E, lam=a.lam, d=a.d)
    print("saved", a.out)


if __name__ == "__main__":
    main()
