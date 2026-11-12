"""Score one agent run against the truth.

    python organizer/score.py --world world_A --run path/to/agent_folder

The agent folder must contain the filled prediction files it was asked for (predictions.csv, predictions_2.csv,
predictions_3.csv; any subset). Prints the median angular error per body for years 101-120 (stages 1 and 2) and
N's error at 1, 5, 10 and 20 years (stage 3), with the stage-3 credit rule from PLAN.md.
"""
import argparse, csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAST = 36525
MARKS = {"1 yr": LAST + 365, "5 yr": LAST + 5 * 365, "10 yr": LAST + 10 * 365, "20 yr": 43815}


def sep(lo1, la1, lo2, la2):
    lo1, la1, lo2, la2 = map(np.radians, (lo1, la1, lo2, la2))
    c = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def load(path):
    out = {}
    for r in csv.DictReader(open(path)):
        try:
            out[(int(float(r["day"])), r["body"])] = (float(r["lon"]), float(r["lat"]))
        except (ValueError, TypeError):
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--world", required=True, choices=["world_A", "world_B"])
    ap.add_argument("--run", required=True, help="the agent's working folder")
    a = ap.parse_args()
    world = json.load(open(os.path.join(HERE, "world_key.json")))["world_key"][a.world]   # contagious / newton
    truth = json.load(open(os.path.join(HERE, f"truth_predictions_{world}.json")))
    ref = json.load(open(os.path.join(HERE, "reference_N.json")))[a.world]["refs"]
    print(f"{a.world} = {world} world")
    for stage, fn in [(1, "predictions.csv"), (2, "predictions_2.csv"), (3, "predictions_3.csv")]:
        p = os.path.join(a.run, fn)
        if not os.path.exists(p):
            print(f"stage {stage}: {fn} not found"); continue
        pred = load(p); errs = {}
        for d, b, lo, la in truth[fn]:
            if (d, b) in pred:
                errs.setdefault(b, []).append((d, sep(*pred[(d, b)], lo, la)))
        missing = len(truth[fn]) - sum(len(v) for v in errs.values())
        if stage < 3:
            print(f"stage {stage}: median error, years 101-120 (deg), {missing} rows missing")
            print("   " + ", ".join(f"{b} {np.median([e for _, e in v]):.4g}" for b, v in sorted(errs.items())))
        else:
            v = np.array(errs.get("N", [])); days, e = v[:, 0], v[:, 1]
            m = {k: float(np.median(e[np.abs(days - d) <= 15])) for k, d in MARKS.items()}
            print(f"stage 3: N error (deg), {missing} rows missing: " + ", ".join(f"{k} {x:.3g}" for k, x in m.items()))
            cloud = {k: {m_: v_["cloud_median"] for m_, v_ in r["errors_deg"].items()} for k, r in ref.items()}
            if world == "contagious":
                ok = all(abs(m[k] - cloud["R_true"][k]) < abs(m[k] - cloud["R_static"][k]) for k in ("10 yr", "20 yr"))
                print(f"   references (deg): true law {cloud['R_true']}, no exchange {cloud['R_static']}")
                print(f"   credit (10 and 20 yr closer to the true law than to no exchange): {'YES' if ok else 'no'}")
            else:
                ok = all(abs(m[k] - cloud["R_static"][k]) < abs(m[k] - cloud["R_exchange"][k]) for k in ("10 yr", "20 yr"))
                print(f"   references (deg): true law (no exchange) {cloud['R_static']}, exchange applied wrongly {cloud['R_exchange']}")
                print(f"   credit (10 and 20 yr closer to no exchange than to a wrongly applied exchange): {'YES' if ok else 'no'}")


if __name__ == "__main__":
    main()
