"""The Tellus star system: star, 6 planets, Tellus's moon, and comets.

Units AU, years, solar masses. All planets and the moon have pull factor 1 (mixed long ago).
"""
import numpy as np
from contagious import elements_to_state

G = 4 * np.pi ** 2
COMET_MASS = 1e-16   # real comets are ~1e-17 solar masses; light enough not to disturb any planet

# name, a (AU), e, inclination (deg), mass (solar masses)
PLANETS = [("P1", 0.45, 0.06, 0.0, 4e-7), ("T", 1.0, 0.02, 0.0, 3e-6), ("P2", 1.7, 0.05, 0.0, 6e-7),
           ("P3", 3.3, 0.04, 1.5, 1.5e-4), ("P4", 6.2, 0.05, 1.0, 1.0e-3), ("P5", 11.5, 0.04, 2.0, 2.5e-4)]
MOON = dict(a=0.00257, e=0.055, inc=5.0, mass=3.7e-8)


def planets_and_moon(seed=1):
    rng = np.random.default_rng(seed)
    names, m, f, X, V = ["S"], [1.0], [1.0], [np.zeros(3)], [np.zeros(3)]
    for name, a, e, inc, mass in PLANETS:
        mu = G * (1 + mass + (MOON["mass"] if name == "T" else 0))
        x, v = elements_to_state(a, e, np.radians(inc), rng.uniform(0, 2 * np.pi), rng.uniform(0, 2 * np.pi),
                                 rng.uniform(0, 2 * np.pi), mu)
        if name != "T":
            names.append(name); m.append(mass); f.append(1.0); X.append(x); V.append(v)
            continue
        # Tellus and its moon: (x, v) is their common centre of mass
        mt, mm = mass, MOON["mass"]
        xr, vr = elements_to_state(MOON["a"], MOON["e"], np.radians(MOON["inc"]), rng.uniform(0, 2 * np.pi),
                                   rng.uniform(0, 2 * np.pi), rng.uniform(0, 2 * np.pi), G * (mt + mm))
        names += ["T", "M"]; m += [mt, mm]; f += [1.0, 1.0]
        X += [x - mm / (mt + mm) * xr, x + mt / (mt + mm) * xr]
        V += [v - mm / (mt + mm) * vr, v + mt / (mt + mm) * vr]
    return names, np.array(m), np.array(f), np.array(X), np.array(V)


def random_comets(n, seed):
    """Stretched orbits that cross the planets' region, each with its own pull factor."""
    rng = np.random.default_rng(seed)
    out = []
    for k in range(n):
        q_peri = rng.uniform(0.3, 4.0)
        a = rng.uniform(max(q_peri * 1.6, 2.0), 14.0)
        e = 1 - q_peri / a
        fk = rng.uniform(0.7, 1.3)
        x, v = elements_to_state(a, e, np.radians(rng.uniform(0, 20)), rng.uniform(0, 2 * np.pi),
                                 rng.uniform(0, 2 * np.pi), rng.uniform(0, 2 * np.pi), G * fk)
        out.append((fk, x, v))
    return out


def assemble(seed_planets, comets, comet_names=None):
    names, m, f, X, V = planets_and_moon(seed_planets)
    names = list(names); m = list(m); f = list(f); X = list(X); V = list(V)
    for k, (fk, x, v) in enumerate(comets):
        names.append(comet_names[k] if comet_names else f"C{k + 1}")
        m.append(COMET_MASS); f.append(fk); X.append(x); V.append(v)
    m, f, X, V = np.array(m), np.array(f), np.array(X), np.array(V)
    X -= (m[:, None] * X).sum(0) / m.sum()
    V -= (m[:, None] * V).sum(0) / m.sum()
    return names, m, f * m, X, V
