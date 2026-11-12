# How to run an agent

Each agent works in its own empty folder that holds only the files of its current stage. It is not told that later
stages exist, which world it is in, or that a second world exists. Give it a shell with Python (numpy, scipy,
matplotlib), no internet, and no access to anything outside its folder. Afterwards, check its transcript for access
outside the folder and for web use.

## Stage 1: initial prompt
```
You are working on a self-contained data-analysis task.

Your working folder is <run folder>
Read README.md in that folder and follow it exactly.

Python with numpy, scipy and matplotlib: <python path> (use this interpreter for all your scripts, including predict.py).

Rules, in addition to the README:
- Read and write files only inside your working folder. Do not list, search, open or read anything anywhere else on this computer (not other folders under <runs folder>, not your home folder, not any other directory). Use paths inside your folder only.
- Do not use the internet: no web search, no web fetch, no downloads, no package installs.
- Take the time you need; do the analysis properly.

When you have written ANSWER_1.md and completed predictions.csv by running predict.py, reply with a short summary of the laws you found.
```

## Stages 2 and 3: hand-over (sent to the same agent, so it keeps its earlier work)
Snapshot the previous stage's folder, then copy in the next stage's files (README.md is replaced).
```
New material has arrived in your folder for the next stage. Your stage <N-1> answer has been saved. Read README.md again (it has been replaced) and follow it exactly.

Same rules as before: read and write only inside <run folder>; no process listings; no internet; at most one Python process at a time, using <python path>. There is no time limit.

When you have written ANSWER_<N>.md and completed predictions_<N>.csv by running predict_<N>.py, reply with a short summary.
```

Give no hints about the physics. If an agent stops before producing a stage's files, the only follow-up is a
neutral request to finish the README's tasks.

## Scoring rules (details in DESIGN.md)
- Stages 1-2: median angular error per body over years 101-120.
- Stage 3 (main outcome): N's error at 1, 5, 10, 20 years vs references computed from the same 207 sightings.
  Contagious world: credit if the 10- and 20-year errors are both closer to the true-law reference than to the
  no-exchange reference. Control world: credit if closer to no-exchange than to a wrongly applied exchange.
- Stage 2 equations are graded E1-E5 by two graders who do not know the model or the world
  (results/RUBRIC.md; our grades in results/grades_*.json, labels in results/grading_key.json).
