"""The star system set-up and what the observer on Tellus records."""
import numpy as np
import pytest


def test_elements_round_trip():
    from contagious import elements_to_state
    G = 4 * np.pi ** 2
    a, e, inc = 3.2, 0.4, np.radians(12)
    x, v = elements_to_state(a, e, inc, 0.7, 1.9, 2.5, G)
    E = 0.5 * v @ v - G / np.linalg.norm(x)
    h = np.cross(x, v)
    assert -G / (2 * E) == pytest.approx(a)
    assert np.linalg.norm(np.cross(v, h) / G - x / np.linalg.norm(x)) == pytest.approx(e)
    assert np.degrees(np.arccos(h[2] / np.linalg.norm(h))) == pytest.approx(np.degrees(inc))


def test_system_is_at_rest_and_planets_have_pull_factor_one():
    from system import assemble, random_comets, COMET_MASS
    names, m, q, X, V = assemble(1, random_comets(12, seed=7))
    assert names[:8] == ["S", "P1", "T", "M", "P2", "P3", "P4", "P5"] and len(names) == 20
    assert np.allclose((m[:, None] * X).sum(0), 0, atol=1e-12)
    assert np.allclose((m[:, None] * V).sum(0), 0, atol=1e-12)
    assert np.allclose(q[:8] / m[:8], 1.0)
    f = q[8:] / m[8:]
    assert np.all((f >= 0.7) & (f <= 1.3)) and np.all(m[8:] == COMET_MASS)


def test_moon_is_bound_to_tellus():
    from system import planets_and_moon
    G = 4 * np.pi ** 2
    names, m, f, X, V = planets_and_moon(1)
    iT, iM = names.index("T"), names.index("M")
    r = np.linalg.norm(X[iM] - X[iT]); v2 = ((V[iM] - V[iT]) ** 2).sum()
    assert 0.002 < r < 0.003
    assert 0.5 * v2 - G * (m[iT] + m[iM]) / r < 0


def test_sky_frame_is_a_pure_rotation():
    from observe import to_frame, lonlat, OBLIQUITY
    v = np.random.default_rng(0).normal(size=(50, 3))
    w = to_frame(v)
    assert np.allclose(np.linalg.norm(w, axis=1), np.linalg.norm(v, axis=1))
    lon, lat = lonlat(np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]))
    assert lon[0] == pytest.approx(0) and lat[0] == pytest.approx(0)
    assert lat[1] == pytest.approx(-np.degrees(OBLIQUITY))        # the planets' plane is tilted on the sky
    assert lat[2] == pytest.approx(90 - np.degrees(OBLIQUITY))


def test_visibility_rules():
    from observe import visible
    # one sample; bodies: 0 star at origin, 1 Tellus at 1 AU, 2 a body
    def X(body):
        return np.array([[[0, 0, 0], [1, 0, 0], body]], float)
    assert not visible(X([-3, 0.1, 0]), 1, 0, 2, faint=False)[0]      # 2 deg from the star on the sky: hidden
    assert visible(X([1, 2, 0]), 1, 0, 2, faint=False)[0]             # 90 deg away: seen
    assert visible(X([1, 2, 0]), 1, 0, 2, faint=True)[0]              # faint, 2.2 AU from the star, 2 AU from Tellus
    assert not visible(X([1, 3.5, 0]), 1, 0, 2, faint=True)[0]        # faint and more than 3 AU from the star
    assert not visible(X([-2.9, 0.5, 0]), 1, 0, 2, faint=True)[0]     # faint: near the star on the sky
