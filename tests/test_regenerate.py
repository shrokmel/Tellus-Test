"""End to end: the simulator rebuilds the shipped agent files and truth."""
import csv, json, os, shutil, subprocess, sys
import numpy as np
from conftest import needs_cc


def table(path):
    return {(r["day"], r["body"]): (r["lon"], r["lat"]) for r in csv.DictReader(open(path))}


@needs_cc
def test_make_obs_reproduces_the_agent_kit(root, tmp_path):
    for d in ("sim", "briefs"):
        shutil.copytree(os.path.join(root, d), tmp_path / d, ignore=shutil.ignore_patterns("sky", "__pycache__", "*.bin", "*.npz", "tmp_*"))
    subprocess.run([sys.executable, "make_obs.py"], cwd=tmp_path / "sim", check=True, capture_output=True)
    for world in ("world_A", "world_B"):
        for stage, files in {1: ["observations.csv", "predictions.csv"], 2: ["observations_2.csv", "predictions_2.csv"],
                             3: ["observations_3.csv", "predictions_3.csv"]}.items():
            for f in files:
                a = table(f"{root}/agent_kit/{world}/stage{stage}/{f}"); b = table(tmp_path / "handover" / world / f"stage{stage}" / f)
                assert a.keys() == b.keys(), f"{world} stage {stage} {f}: different sightings"
                if f.startswith("observations"):
                    da = np.array([[float(x) for x in a[k]] for k in a]); db = np.array([[float(x) for x in b[k]] for k in a])
                    dlon = (da[:, 0] - db[:, 0] + 180) % 360 - 180
                    assert np.abs(dlon).max() < 1e-3 and np.abs(da[:, 1] - db[:, 1]).max() < 1e-3
    assert json.load(open(tmp_path / "private" / "world_key.json"))["world_key"] == json.load(open(f"{root}/organizer/world_key.json"))["world_key"]
