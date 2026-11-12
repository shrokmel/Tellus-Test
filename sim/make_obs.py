"""Final runs (both worlds, 121 years, with N) and the observation files the agent receives.

Output:
  ../handover/world_A, world_B / stage1, stage2, stage3 : README.md, observations*.csv, predictions*.csv (blank)
  ../private/ : world key, true directions for every prediction row, run summaries
Which world is A and which is B is drawn at random (seed 2026) and recorded only in ../private/world_key.json.
"""
import json, os, shutil, numpy as np
from system import assemble, random_comets
from world import run
from observe import lonlat, visible, CLEAR, NOISE

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LAST_DAY = 36525                      # end of year 100
N_FIRST = 36233                       # N first visible from Tellus (start of its current appearance)
WIN = np.arange(N_FIRST, LAST_DAY + 1)   # N's whole current appearance, year 99.2-100
PRED_END = 43830                      # end of year 120
ch = json.load(open("runs/final_choice.json")); cands = random_comets(300, seed=7)
D = json.load(open("runs/n_design.json"))[8]
names, m, q, X, V = assemble(1, [cands[i] for i in ch["chosen"]])
names = list(names) + ["N"]
m = np.r_[m, 1e-16]; q = np.r_[q, D["f"] * 1e-16]; X = np.vstack([X, D["x0"]]); V = np.vstack([V, D["v0"]])
iS, iT = 0, names.index("T")
night = [n for n in names if n not in ("S", "T")]
planets = ["P1", "P2", "P3", "P4", "P5", "M"]
comets = [f"C{k}" for k in range(1, 13)]

rng = np.random.default_rng(2026)
clear_night = rng.random(PRED_END + 1) < CLEAR
clear_day = rng.random(PRED_END + 1) < CLEAR
worlds = {"contagious": 50.0, "newton": 0.0}
labels = dict(zip(["contagious", "newton"], rng.permutation(["world_A", "world_B"])))
noise_rng = np.random.default_rng(99)
noise_rng_N = np.random.default_rng(100)

shutil.rmtree(os.path.join(ROOT, "handover"), ignore_errors=True)
os.makedirs(os.path.join(ROOT, "private"), exist_ok=True)
summary = {}
for world, lam in worlds.items():
    r = run(names, m, q, X, V, lam=lam, years=PRED_END / 365.25, sample_days=1, rtol=1e-12, tag=f"obs_{world}")
    t, Xs = r["t"], r["X"]; day = np.rint(t * 365.25).astype(int)
    assert (day == np.arange(len(day))).all()
    np.savez_compressed(os.path.join(ROOT, "private", f"truth_run_{world}.npz"), t=t, X=Xs, V=r["V"], Q=r["Q"], m=m,
                        names=np.array(names))
    true_dir = {}
    for b in names:
        if b == "T": continue
        i = names.index(b)
        true_dir[b] = lonlat(Xs[:, i] - Xs[:, iT])
    out = os.path.join(ROOT, "handover", labels[world]); rows = {1: [], 2: [], 3: []}; counts = {}
    for b in names:
        if b == "T": continue
        i = names.index(b); lon, lat = true_dir[b]
        if b == "S":
            ok = clear_day.copy()
        else:
            ok = visible(Xs, iT, iS, i, faint=(b in comets or b == "N")) & clear_night
        days = np.arange(LAST_DAY + 1)
        if b == "N":
            days = WIN
        days = days[ok[days]]
        nr = noise_rng_N if b == "N" else noise_rng
        dlat = nr.normal(0, np.degrees(NOISE), len(days))
        dlon = nr.normal(0, np.degrees(NOISE), len(days)) / np.cos(np.radians(lat[days]))
        stage = 1 if (b == "S" or b in planets) else (2 if b in comets else 3)
        for dd, lo, la in zip(days, (lon[days] + dlon) % 360, lat[days] + dlat):
            rows[stage].append((int(dd), b, lo, la))
        counts[b] = len(days)
    for s, fname in [(1, "observations.csv"), (2, "observations_2.csv"), (3, "observations_3.csv")]:
        os.makedirs(os.path.join(out, f"stage{s}"), exist_ok=True)
        rows[s].sort()
        with open(os.path.join(out, f"stage{s}", fname), "w") as fh:
            fh.write("day,body,lon,lat\n")
            fh.writelines(f"{d},{b},{lo:.5f},{la:.5f}\n" for d, b, lo, la in rows[s])
        shutil.copy(os.path.join(ROOT, "briefs", f"stage{s}_README.md"), os.path.join(out, f"stage{s}", "README.md"))
    # prediction templates (blank) and the matching truth
    pred_days = np.arange(LAST_DAY + 30, PRED_END + 1, 30)
    n_days = np.arange(LAST_DAY + 5, PRED_END + 1, 5)
    templates = {1: [(d, b) for d in pred_days for b in ["S"] + planets],
                 2: [(d, b) for d in pred_days for b in ["S"] + planets + comets],
                 3: [(d, "N") for d in n_days]}
    truth = {}
    for s, fname in [(1, "predictions.csv"), (2, "predictions_2.csv"), (3, "predictions_3.csv")]:
        with open(os.path.join(out, f"stage{s}", fname), "w") as fh:
            fh.write("day,body,lon,lat\n")
            fh.writelines(f"{d},{b},,\n" for d, b in templates[s])
        truth[fname] = [(int(d), b, float(true_dir[b][0][d]), float(true_dir[b][1][d])) for d, b in templates[s]]
    json.dump(truth, open(os.path.join(ROOT, "private", f"truth_predictions_{world}.json"), "w"))
    f = r["Q"] / m
    summary[world] = dict(label=labels[world], counts=counts,
                          comet_pull_start={c: round(float(f[0, names.index(c)]), 4) for c in comets + ["N"]},
                          comet_pull_day36525={c: round(float(f[LAST_DAY, names.index(c)]), 4) for c in comets + ["N"]},
                          comet_pull_end={c: round(float(f[-1, names.index(c)]), 4) for c in comets + ["N"]})
    print(world, "->", labels[world], f"({r['seconds']:.1f}s)")
    print("  sightings:", counts)
json.dump(dict(world_key={v: k for k, v in labels.items()}, weather_seed=2026, noise_seed=99, noise_seed_N=100, N_first_day=N_FIRST, N_design=D,
               note="Visibility rules in sim/observe.py, fixed before generation."),
          open(os.path.join(ROOT, "private", "world_key.json"), "w"), indent=1)
json.dump(summary, open(os.path.join(ROOT, "private", "summary.json"), "w"), indent=1)
