# The Tellus Test: replication kit

Can an AI agent discover a law of physics that is not in any textbook, from raw observations alone?

We simulated a star system with a hidden extra law, recorded 100 years of what an astronomer on one planet would
see (directions on the sky only), and asked AI agents to find the laws and predict. This kit lets you rerun the
test with any model, score it, and regenerate every data file from the simulator.

> **Spoiler warning.** `organizer/` and `sim/` contain the hidden law and the answers. Never let an agent see them.
> Give an agent only the contents of one `agent_kit/<world>/stageN/` folder at a time.

## Layout
| Folder | What it is | Who may see it |
|---|---|---|
| `agent_kit/world_A`, `world_B` | The exact files the agents received, per stage: `README.md` (the brief), `observations*.csv`, blank `predictions*.csv` | the agent |
| `briefs/` | The three stage briefs (same as in `agent_kit`) | anyone |
| `organizer/` | Which world is which, true answers, reference fits, `score.py`, `DESIGN.md` (physics, system, scoring rules), `PROTOCOL.md` (prompts), results | you only |
| `organizer/results/` | Scores, blind grades of the equations, and every agent's answers and prediction scripts | you only |
| `sim/` | The simulator (C + Python) and all scripts that built the worlds, references and baselines | you only |
| `tests/` | Test suite for the simulator, the scorer and the data files | you only |

## Running a model
1. Make an empty folder per run. Copy in `agent_kit/<world>/stage1/*`.
2. Start the agent with the stage-1 prompt in `organizer/PROTOCOL.md` (adapt the paths). It must write
   `ANSWER_1.md`, `predict.py`, and fill `predictions.csv` by running `predict.py`.
3. Copy in `agent_kit/<world>/stage2/*` (replaces README.md) and send the hand-over message to the same agent
   session. Repeat for stage 3.
4. Run each world in a separate session, and don't tell the agent which world it is in.
5. Score: `python organizer/score.py --world world_A --run <run folder>`.
   It prints per-body errors and the stage-3 credit verdict. For stage 2, also grade the equations in
   `ANSWER_2.md`/`ANSWER_3.md` against `organizer/results/RUBRIC.md` (ideally blind to model and world).

Python 3 with numpy and scipy is enough for the agent kit and the scorer.

## Regenerating the data
```
cd sim
python make_obs.py
```
Needs a C compiler (`cc`), numpy and scipy; it compiles `sky.c` on first use and takes a few seconds. It rebuilds
both worlds to year 120 and writes `handover/` (agent files) and `private/` (truth). A fresh copy reproduces all
18 agent files and both truth files byte-for-byte on macOS/arm64 with clang; other compilers or CPUs may differ in
the last digits.

The other scripts in `sim/` document how the system was designed (`build_final.py`, `scan_comets.py`,
`design_n.py`, `check_n.py`), how the reference fits were made (`reference_n.py`, `ref_windows.py`,
`baseline_stage1.py`, `baseline_stage2.py`), and checks of the simulator (`converge.py`, `long_run.py`).

## Tests
```
pip install -r requirements.txt -r requirements-dev.txt
pytest
```
30 tests, about 10 seconds. They cover the integrator (ordinary gravity, charge conservation, the evening-out rate,
no flow outside contact or in the control world), the star system and the visibility rules, the scorer, the shipped
data files (including that the agent kit leaks no answers), the published scores, and an end-to-end rebuild of the
agent kit from the simulator. Tests that need the simulator are skipped when no C compiler is available.

## Results (Claude models, one run per model per world)
| Model | Contagious world: N error after 10 / 20 yr | Control world: N error after 10 / 20 yr |
|---|---|---|
| Fable 5.2 | 0.20° / 0.76° (credit) | 7.0° / 5.2° (credit) |
| Opus 5.5 | 2.6° / 8.6° (credit) | 7.1° / 5.2° (credit) |
| Sonnet 5.5 | 17° / 35° | 4.3° / 3.2° (credit) |
| Haiku 5.5 | 15° / 24° | 39° / 11.5° |

Reference fits in the contagious world: true law 2.4° / 6.3°, no exchange 22.5° / 26.4°.
To re-run an agent's own prediction script, copy its folder from `organizer/results/agent_answers/` together with
the matching `agent_kit/<world>/stage*/observations*.csv` into one folder and run `predict.py`, `predict_2.py`,
`predict_3.py`. All 24 reproduce the submitted predictions exactly.

Full numbers: `organizer/results/scores.json`; equation grades: `organizer/results/grades_A.json`, `grades_B.json`.

## Limits worth knowing before you compare
- All four models are from one family (Claude). Other model families are the obvious next test.
- One run per model per world: the ordering may change with repeats.
- The equations were graded by Claude models (blind to model and world), not by people.
- The agents' answers in `organizer/results/agent_answers/` are as the agents wrote them, with working files and
  one remark about the host machine removed.
- The simulator's moon carries a small numerical drift (-0.4 arcmin per century^2, from integrator error). It is
  far below the 1-arcminute observation noise and does not affect the scoring.
- Tip: ask each agent to keep a `NOTEBOOK.md` of what it tests and why; it gives a readable record of its reasoning.

## Licence
MIT. See `LICENSE`.
