"""The scorer, the shipped data files and the published results."""
import csv, json, os, subprocess, sys
import numpy as np
import pytest

WORLDS = ["world_A", "world_B"]
FILES = {1: ("observations.csv", "predictions.csv"), 2: ("observations_2.csv", "predictions_2.csv"),
         3: ("observations_3.csv", "predictions_3.csv")}
BODIES = {1: {"S", "P1", "P2", "P3", "P4", "P5", "M"}, 2: {f"C{k}" for k in range(1, 13)}, 3: {"N"}}


def rows(path):
    return list(csv.DictReader(open(path)))


def test_angular_separation():
    from score import sep
    assert sep(10, 20, 10, 20) == pytest.approx(0, abs=1e-6)
    assert sep(0, 0, 90, 0) == pytest.approx(90)
    assert sep(0, 90, 123, -90) == pytest.approx(180)
    assert sep(359, 0, 1, 0) == pytest.approx(2)


def test_world_key(root):
    key = json.load(open(f"{root}/organizer/world_key.json"))["world_key"]
    assert sorted(key) == WORLDS and sorted(key.values()) == ["contagious", "newton"]


@pytest.mark.parametrize("world", WORLDS)
@pytest.mark.parametrize("stage", [1, 2, 3])
def test_agent_files(root, world, stage):
    d = f"{root}/agent_kit/{world}/stage{stage}"
    obs, pred = FILES[stage]
    assert os.path.exists(f"{d}/README.md")
    o = rows(f"{d}/{obs}")
    assert set(o[0]) == {"day", "body", "lon", "lat"}
    assert {r["body"] for r in o} == BODIES[stage]
    days = np.array([int(r["day"]) for r in o]); lon = np.array([float(r["lon"]) for r in o]); lat = np.array([float(r["lat"]) for r in o])
    assert days.min() >= 0 and days.max() <= 36525                      # nothing from the forecast years
    assert np.all((lon >= 0) & (lon < 360)) and np.all(np.abs(lat) <= 90)
    p = rows(f"{d}/{pred}")
    assert all(r["lon"] == "" and r["lat"] == "" for r in p)           # templates are blank: no answers leak
    assert min(int(r["day"]) for r in p) > 36525


def test_agent_kit_contains_no_answers(root):
    for dirpath, _, files in os.walk(f"{root}/agent_kit"):
        for f in files:
            text = open(os.path.join(dirpath, f), errors="ignore").read().lower()
            for word in ("contagious", "charge", "pull factor", "newton", "kepler"):
                assert word not in text, f"{f} mentions '{word}'"


@pytest.mark.parametrize("world", WORLDS)
def test_truth_covers_every_prediction_row(root, world):
    key = json.load(open(f"{root}/organizer/world_key.json"))["world_key"]
    truth = json.load(open(f"{root}/organizer/truth_predictions_{key[world]}.json"))
    for stage, (_, pred) in FILES.items():
        want = {(int(r["day"]), r["body"]) for r in rows(f"{root}/agent_kit/{world}/stage{stage}/{pred}")}
        assert want == {(d, b) for d, b, _, _ in truth[pred]}


def test_planets_identical_in_both_worlds(root):
    """Comets are too light to disturb the planets, so the bright bodies move alike in both worlds
    (to far below the 1-arcminute = 0.017 degree observation noise)."""
    a = json.load(open(f"{root}/organizer/truth_predictions_contagious.json"))["predictions.csv"]
    b = json.load(open(f"{root}/organizer/truth_predictions_newton.json"))["predictions.csv"]
    assert np.abs(np.array([r[2:] for r in a]) - np.array([r[2:] for r in b])).max() < 1e-3


@pytest.mark.parametrize("world", WORLDS)
def test_perfect_predictions_score_zero_and_earn_credit(root, world, tmp_path):
    key = json.load(open(f"{root}/organizer/world_key.json"))["world_key"]
    truth = json.load(open(f"{root}/organizer/truth_predictions_{key[world]}.json"))
    for fn, tr in truth.items():
        with open(tmp_path / fn, "w") as fh:
            fh.write("day,body,lon,lat\n")
            fh.writelines(f"{d},{b},{lo!r},{la!r}\n" for d, b, lo, la in tr)
    out = subprocess.run([sys.executable, f"{root}/organizer/score.py", "--world", world, "--run", str(tmp_path)],
                         capture_output=True, text=True, check=True).stdout
    assert "0 rows missing" in out and "credit" in out and "YES" in out
    assert "20 yr 0" in out


def test_published_scores_reproduce(root):
    key = json.load(open(f"{root}/organizer/results/run_key.json"))
    pub = json.load(open(f"{root}/organizer/results/scores.json"))
    import re
    for run, v in key.items():
        d = f"{root}/organizer/results/agent_answers/{run}_{v['model']}_{v['world']}"
        out = subprocess.run([sys.executable, f"{root}/organizer/score.py", "--world", v["world"], "--run", d],
                             capture_output=True, text=True, check=True).stdout
        m = re.search(r"10 yr ([\d.e+-]+), 20 yr ([\d.e+-]+)", out)
        assert float(m.group(1)) == pytest.approx(pub[run]["N"]["10 yr"], rel=0.01, abs=0.005)
        assert float(m.group(2)) == pytest.approx(pub[run]["N"]["20 yr"], rel=0.01, abs=0.005)


def test_reference_fits_separate_the_laws(root):
    """The stage 3 yardstick only works if the true law and the wrong law predict clearly different futures."""
    ref = json.load(open(f"{root}/organizer/reference_N.json"))
    key = json.load(open(f"{root}/organizer/world_key.json"))["world_key"]
    for world, kind in key.items():
        r = ref[world]["refs"]; c = lambda name, mark: r[name]["errors_deg"][mark]["cloud_median"]
        other = "R_static" if kind == "contagious" else "R_exchange"
        for mark in ("10 yr", "20 yr"):
            assert c(other, mark) > 2 * c("R_true", mark)
        assert abs(r["R_true"]["fitted_pull_factor"] - ref[world]["true_pull_factor"]) < 0.01
