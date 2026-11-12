"""Shared tools for the stage 1-3 baselines: load observations, run the integrator in parallel, angular errors."""
import csv, json, os, uuid
import numpy as np
from world import run
from observe import lonlat

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TMP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", "tmp_baselines")
os.makedirs(TMP, exist_ok=True)
LAST, END = 36525, 43830
KEY = json.load(open(os.path.join(ROOT, "private", "world_key.json")))["world_key"]   # label -> world


def truth_run(world):
    tr = np.load(os.path.join(ROOT, "private", f"truth_run_{world}.npz"), allow_pickle=True)
    return dict(t=tr["t"], X=tr["X"], V=tr["V"], Q=tr["Q"], m=tr["m"], names=list(tr["names"]))


def load_obs(label, stage):
    fname = {1: "observations.csv", 2: "observations_2.csv", 3: "observations_3.csv"}[stage]
    out = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "handover", label, f"stage{stage}", fname))):
        out.setdefault(r["body"], []).append((int(r["day"]), float(r["lon"]), float(r["lat"])))
    return {b: tuple(np.array(c) for c in zip(*v)) for b, v in out.items()}


def load_truth_pred(world, fname):
    rows = json.load(open(os.path.join(ROOT, "private", f"truth_predictions_{world}.json")))[fname]
    return rows                                   # list of (day, body, lon, lat)


def integrate(names, m, q, X, V, lam, days, tp=None, sample_days=1):
    tag = "b" + uuid.uuid4().hex[:10]
    r = run(names, m, q, X, V, lam=lam, d=0.1, years=days / 365.25, sample_days=sample_days, rtol=1e-12, tp=tp,
            tag=tag, outdir=TMP)
    for ext in (".bin", ".in.txt"):
        try:
            os.remove(os.path.join(TMP, tag + ext))
        except OSError:
            pass
    return r


def resid_angles(lo, la, olo, ola):
    """Residuals in arcminutes: (dlon cos lat, dlat)."""
    return np.r_[((lo - olo + 180) % 360 - 180) * np.cos(np.radians(ola)), la - ola] * 60.0


def angsep(lo1, la1, lo2, la2):
    lo1, la1, lo2, la2 = map(np.radians, (lo1, la1, lo2, la2))
    c = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def directions_from(r, i_obs, i_body):
    return lonlat(r["X"][:, i_body] - r["X"][:, i_obs])
