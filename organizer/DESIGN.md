# The Tellus Test: design and scoring

## Question
Inside Newton's framework we added one new entity: every body carries a hidden **charge** that sets how strongly
gravity acts on it, and charge **flows between bodies during close passes**, like heat. Given only directions on
the sky recorded from one planet for 100 years, does an AI agent

1. find the ordinary laws where they hold (the planets),
2. notice that the faint bodies break them, and invent the hidden, conserved, exchanged quantity, and
3. use it to predict a newcomer's path through a close pass it has not yet made?

A control world, identical except that charge never flows, tests whether the agent invents an exchange that is
not there.

## The physics
Each body i has a mass m_i (inertia, fixed) and a charge q_i. Units: AU, years, solar masses, G = 4 pi^2.
- Force between i and j: `G q_i q_j / r^2`, along the line between them. Acceleration of i: `G q_j (q_i/m_i) / r^2`.
- Pull factor `f_i = q_i / m_i`. In our universe f = 1 for everything.
- Charge flow: `dq_i/dt = -lam * sum_j w(r_ij) * mu_ij * (f_i - f_j)`, `w(r) = exp(-(r/d)^2)`,
  `mu_ij = m_i m_j / (m_i + m_j)`. For one pair, the gap in pull factor closes as
  `d(f_i - f_j)/dt = -lam w (f_i - f_j)`. Total charge is conserved; momentum and angular momentum are conserved;
  energy is not (it jumps at passes, by exactly the amount the charge flow implies).
- **Contagious world:** lam = 50 per year, d = 0.1 AU. **Control world:** lam = 0 (each body keeps its own pull
  factor forever: different bodies fall differently, but nothing is exchanged).
- Simulator: `sim/sky.c` (Dormand-Prince 5(4), rtol 1e-12; checked against an independent scipy integration and
  for convergence; planets accurate to ~2e-7 AU over 100 years). Driver `sim/world.py`.

## The star system (`sim/system.py`, `sim/build_final.py`, `sim/design_n.py`)
- `S`: star, mass 1, f = 1.
- Planets, all f = 1 ("mixed long ago"): P1 0.45 AU, Tellus (T, the observers' planet) 1.0 AU, P2 1.7, P3 3.3,
  P4 6.2 (1e-3 solar masses, the giant), P5 11.5. Low inclinations.
- `M`: Tellus's moon, 0.00257 AU, f = 1.
- Comets `C1`-`C12`, mass 1e-16 (too light to change any planet), chosen from a scan of 300 random comets:
  C1-C6 one clear pass each (years 18-87; at Tellus, P2, P3, P4; three gain and three lose pull factor; C1 is
  thrown out of the system), C7 two passes at Tellus, C8-C12 no clear pass (C11 has gradual edge-of-contact changes
  and one comet-comet meeting in year 90.7 where charge balances exactly).
- `N`: newcomer, f = 1.323, period ~28 years. Seen from the night it first becomes visible (day 36233, year 99.2)
  to the end of year 100: 207 clear sightings. Passes P4 at 0.05 AU in year 101.2; in the contagious world its
  pull factor falls to 1.04 there. Identical in both worlds until year 101; its position then differs by 2 AU at
  year 105 and 10 AU at year 110.
- Stable for 1,000 years in both worlds (planets unchanged; moon bound).

## What the agent sees (`sim/observe.py`, `sim/make_obs.py`, `agent_kit/`)
- Directions only (lon, lat in a fixed frame against the stars, tilted 21.7 deg by Tellus's spin axis), one row
  per sighting, 1 arcminute noise. No distances, sizes or brightnesses.
- A body is seen only when at least 15 deg from `S` and on a clear night (70%); `S` on clear days. Faint bodies
  (comets, N) only when within 3 AU of S and 4 AU of Tellus.
- Bodies are named only `S` (makes daylight), `P1`-`P5` (wandering), `M` (large, bright), `C1`-`C12`, `N` (faint).
  Nothing says what orbits what.
- Worlds are labelled `world_A` / `world_B` at random; the key is in `organizer/world_key.json`.

## Stages (one agent session per world; the agent is not told that later stages exist)
1. **The bright bodies** (`S`, `P1`-`P5`, `M`; 100 years). Describe the motions; give the laws as equations with
   numerical constants and uncertainties; predict all bodies every 30 days in years 101-120 with `predict.py`;
   confidence and falsifiers. The answer is saved before stage 2.
2. **The faint bodies** (`C1`-`C12`). Test the stage 1 laws against them; give the laws that explain everything,
   as equations; predict every body in years 101-120 with `predict_2.py`.
3. **N** (year 99.2 to the end of year 100). Give the equations used; predict `N` every 5 days for 20 years with
   `predict_3.py`; describe its future path; state expected error after 1, 5 and 20 years.

Every stage requires equations with numerical constants, and predictions produced by the agent's own script from
those equations only. Every script is re-run; predictions that do not come from the stated equations do not
count.

## Scoring

### Stage 1 (both worlds identical): did it find Kepler and Newton?
Mechanical checks on the stated equations: a central attraction toward S for P1-P5 and Tellus; toward Tellus for
M; strength falling as 1/r^n with n within 2 +/- 0.05 (or Kepler's third law with exponent 1.5 +/- 0.03); a
single constant for all planets. Prediction error (median angular error over years 101-120) against reference R1
(an ideal Newtonian fit to the same observations, `organizer/baseline_stage1.json`).

### Stage 2: did it find the entity?
Five elements, each scored present / absent from the stated equations (not from prose):
- **E1 Pull factor:** bodies differ in how strongly S pulls them (a per-body constant multiplying the attraction).
  Correct in both worlds.
- **E2 Changes at passes:** a body's pull factor changes, and only during close passes. Correct only in the
  contagious world; claiming it in the control world is a false discovery.
- **E3 Evening out:** the change moves the two bodies' pull factors toward each other.
- **E4 Conservation:** a conserved quantity (charge) with exchange weighted by mass, so that a small body takes on
  a large body's value.
- **E5 Rule constants:** contact distance and rate, within a factor of 3 of 0.1 AU and 50 per year (or an
  equivalent parametrisation that reproduces the observed changes).
Primary stage 2 outcome: contagious world E1-E4 all present; control world E1 present and E2 absent.
Secondary: prediction error for C1-C12 in years 101-120 against reference R2 (the true law fitted to the same
observations, `organizer/baseline_stage2.json`).
Grading of E1-E5: a script where the equations allow it; otherwise two graders who see the answer and the rubric
but not which world it came from.

### Stage 3: does the law predict N? (main outcome)
Angular error of the agent's `N` predictions after 1, 5, 10 and 20 years (median over the days within 15 days of
each mark), compared with references computed from the same 207 sightings (`sim/reference_n.py`,
`organizer/reference_N.json`). Each reference is a least-squares fit of N's position, velocity and (where free) pull
factor, idealised: true planet positions known, search started at the truth. "Cloud" = median over 50 orbits
drawn from the fit's uncertainty, the fair measure of what the data allow.

| World | Reference | 1 yr | 5 yr | 10 yr | 20 yr (degrees, cloud median) |
|---|---|---|---|---|---|
| contagious | R_true (charge rule) | 0.04 | 1.4 | 2.4 | 6.3 |
| contagious | R_static (own pull factor, no exchange) | 0.04 | 2.5 | 22.5 | 26.4 |
| contagious | R_newton (pull factor 1) | 9.9 | 20.9 | 30.1 | 49.5 |
| control | R_true = R_static | 0.05 | 2.3 | 11.1 | 10.3 |
| control | R_exchange (charge rule applied wrongly) | 0.05 | 2.2 | 24.9 | 25.6 |
| control | R_newton | 9.9 | 23.1 | 55.7 | 27.0 |

Fitted pull factor in both worlds: 1.323 +/- 0.002 (true 1.323). Plain Newton cannot fit the sightings (rms 7.8'
against 1.0' for the others), a clue in itself.

Credit rule: **contagious world:** the agent is credited with the entity if its 10- and 20-year errors are both
closer to R_true than to R_static. **Control world:** credited with not inventing an exchange if its 10- and
20-year errors are both closer to R_static than to R_exchange. The 1-year mark (before the pass) tests orbit
fitting only; the 5-year mark is reported but does not discriminate (the laws differ by about 1 degree there).
Also: its stated expected errors against its actual errors (calibration).
