"""The integrator (sim/sky.c): ordinary gravity, and the charge-flow rule."""
import numpy as np
import pytest
from conftest import needs_cc

pytestmark = needs_cc
G = 4 * np.pi ** 2


def run(tmp_path, m, q, X, V, lam, years, tp=None, sample_days=1.0, d=0.1):
    from world import run as _run
    names = [f"b{i}" for i in range(len(m))]
    return _run(names, np.array(m, float), np.array(q, float), np.array(X, float), np.array(V, float), lam=lam, d=d,
                years=years, sample_days=sample_days, rtol=1e-12, tp=tp, tag="t", outdir=str(tmp_path))


def circular(a, mu):
    return [a, 0, 0], [0, np.sqrt(mu / a), 0]


def test_kepler_period_and_energy(tmp_path):
    """A light planet at 1 AU around a 1 solar-mass star returns to its start after exactly one year."""
    x, v = circular(1.0, G)
    r = run(tmp_path, [1.0, 1e-12], [1.0, 1e-12], [[0, 0, 0], x], [[0, 0, 0], v], lam=0, years=1.0, sample_days=365.25 / 4)
    rel = r["X"][:, 1] - r["X"][:, 0]
    assert np.allclose(rel[-1], rel[0], atol=1e-8)
    assert np.allclose(np.linalg.norm(rel, axis=1), 1.0, atol=1e-9)
    from world import energy
    E = energy(r)
    assert abs(E[-1] - E[0]) < 1e-10 * abs(E[0])


def test_pull_factor_scales_gravity(tmp_path):
    """A body with pull factor f on a circular orbit needs speed sqrt(G f / a): Kepler's constant scales with f."""
    f = 1.3
    x, v = circular(2.0, G * f)
    r = run(tmp_path, [1.0, 1e-12], [1.0, f * 1e-12], [[0, 0, 0], x], [[0, 0, 0], v], lam=0, years=3.0, sample_days=30)
    dist = np.linalg.norm(r["X"][:, 1] - r["X"][:, 0], axis=1)
    assert np.allclose(dist, 2.0, atol=1e-8)


def test_no_flow_without_lambda(tmp_path):
    """Control world: with lam = 0 charges never change, even in contact."""
    r = run(tmp_path, [1e-9, 1e-16], [1e-9, 1.3e-16], [[0, 0, 0], [0.02, 0, 0]], [[0, 0, 0], [0, 0, 0]], lam=0,
            years=0.25, sample_days=5)
    assert np.array_equal(r["Q"][0], r["Q"][-1])


def test_charge_evens_out_and_is_conserved(tmp_path):
    """In contact, a light body's pull factor moves to the heavy body's value; total charge is unchanged."""
    m = [1e-9, 1e-16]                              # a heavy and a light body, 0.02 AU apart, far from anything else
    r = run(tmp_path, m, [1e-9, 1.3e-16], [[0, 0, 0], [0.02, 0, 0]], [[0, 0, 0], [0, 0, 0]], lam=50, years=0.25,
            sample_days=5)
    f = r["Q"] / np.array(m)
    assert f[0, 1] == pytest.approx(1.3)
    assert abs(f[-1, 1] - 1.0) < 0.01              # evened out with the heavy body
    assert abs(f[-1, 0] - 1.0) < 1e-6              # the heavy body barely notices
    assert np.all(np.diff(f[:, 1]) <= 1e-12)       # monotonic, toward the heavy body's value
    assert r["Q"][-1].sum() == pytest.approx(r["Q"][0].sum(), rel=1e-13)


def test_gap_closes_at_the_stated_rate(tmp_path):
    """Two bodies held near each other: the pull-factor gap decays as exp(-lam * w * t), w = exp(-(r/d)^2)."""
    m = [1e-6, 1e-6]; sep, lam, d = 0.05, 50.0, 0.1
    v = np.sqrt(G * 1e-6 * 1.1 * 0.9 * 2 / sep) / 2     # mutual circular orbit, so the distance stays fixed
    r = run(tmp_path, m, [1.1e-6, 0.9e-6], [[-sep / 2, 0, 0], [sep / 2, 0, 0]], [[0, -v, 0], [0, v, 0]], lam=lam,
            years=0.02, sample_days=1, d=d)
    f = r["Q"] / np.array(m); gap = f[:, 0] - f[:, 1]
    dist = np.linalg.norm(r["X"][:, 0] - r["X"][:, 1], axis=1)
    assert np.allclose(dist, sep, rtol=0.05)
    rate = -np.log(gap[-1] / gap[0]) / r["t"][-1]
    assert rate == pytest.approx(lam * np.exp(-(dist.mean() / d) ** 2), rel=0.05)


def test_no_flow_outside_contact(tmp_path):
    """Bodies several contact distances apart do not exchange charge."""
    r = run(tmp_path, [1.0, 1e-3, 1e-12], [1.0, 1e-3, 1.3e-12], [[0, 0, 0], [5, 0, 0], [6, 0, 0]],
            [[0, 0, 0], [0, 2.8, 0], [0, 2.9, 0]], lam=50, years=0.05, sample_days=5)
    assert r["Q"][-1, 2] == pytest.approx(r["Q"][0, 2], rel=1e-12)


def test_momentum_conserved_with_flow(tmp_path):
    m = np.array([1.0, 1e-3, 1e-6])
    r = run(tmp_path, m, [1.0, 1e-3, 1.3e-6], [[0, 0, 0], [5, 0, 0], [5.03, 0, 0]], [[0, 0, 0], [0, 2.8, 0], [0, 3.0, 0.1]],
            lam=50, years=0.5, sample_days=10)
    P = (m[None, :, None] * r["V"]).sum(1)
    assert np.allclose(P[-1], P[0], atol=1e-12)


def test_test_particles_do_not_pull(tmp_path):
    """A body flagged as a test particle feels gravity but does not disturb the others."""
    base = dict(m=[1.0, 1e-3], q=[1.0, 1e-3], X=[[0, 0, 0], [5, 0, 0]], V=[[0, 0, 0], [0, 2.8, 0]])
    a = run(tmp_path, base["m"], base["q"], base["X"], base["V"], lam=0, years=1, sample_days=30)
    b = run(tmp_path, base["m"] + [0.5], base["q"] + [0.5], base["X"] + [[2, 0, 0]], base["V"] + [[0, 4.4, 0]], lam=0,
            years=1, sample_days=30, tp=[0, 0, 1])
    assert np.allclose(a["X"], b["X"][:, :2], atol=1e-12)
