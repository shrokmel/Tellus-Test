"""Turn true positions into what an observer on Tellus records.

Directions only: lon, lat in a fixed frame against the distant stars. The frame is tilted by Tellus's spin axis
(OBLIQUITY), like equatorial coordinates, so the planets' common plane is not handed over as lat = 0.
Visibility rules:
  - every body except S: at least 15 deg from S on the sky, on a clear night (70% of nights, same for all bodies);
  - S: measured by day, on a clear day (70% of days, drawn separately from nights);
  - faint bodies (comets and N) in addition: within 3 AU of S and within 4 AU of Tellus;
  - noise: 1 arcminute (Gaussian, 1 sigma) in lat, and 1 arcminute of arc in lon (so 1'/cos(lat) in lon).
"""
import numpy as np

OBLIQUITY = np.radians(21.7)
MIN_ELONG = np.radians(15.0)
CLEAR = 0.7
FAINT_R_STAR, FAINT_R_TELLUS = 3.0, 4.0
NOISE = np.radians(1 / 60)


def to_frame(vec):
    """Rotate simulation coordinates into Tellus's fixed sky frame (tilted by the obliquity about x)."""
    c, s = np.cos(OBLIQUITY), np.sin(OBLIQUITY)
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    return np.stack([x, c * y + s * z, -s * y + c * z], axis=-1)


def lonlat(vec):
    v = to_frame(vec)
    r = np.linalg.norm(v, axis=-1)
    lon = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0
    lat = np.degrees(np.arcsin(np.clip(v[..., 2] / r, -1, 1)))
    return lon, lat


def visible(X, iT, iS, i, faint):
    """Boolean per sample: can body i be seen from Tellus (ignoring weather)?"""
    d = X[:, i] - X[:, iT]; s = X[:, iS] - X[:, iT]
    cosel = (d * s).sum(1) / (np.linalg.norm(d, axis=1) * np.linalg.norm(s, axis=1))
    ok = cosel <= np.cos(MIN_ELONG)
    if faint:
        ok &= (np.linalg.norm(X[:, i] - X[:, iS], axis=1) < FAINT_R_STAR) & (np.linalg.norm(d, axis=1) < FAINT_R_TELLUS)
    return ok
