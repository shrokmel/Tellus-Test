# Stage 2: the faint bodies C1–C12

Summary. The stage 1 law (Newtonian gravitation of S, Tellus, M, P1–P5) is **necessary but not sufficient**
for the faint bodies. It holds at the 1′ noise level for some of them for decades, and fails by tens to
hundreds of arcminutes for most of them, even inside a single apparition. Every failure is removed by **one
extra term**: a radial, inverse-square force centred on S whose strength β is a property of the individual
faint body (either sign, |β| up to 0.30 of solar gravity). β is not eternal: it **fades when the body passes
close to a planet**. The massive bodies need no such term; their stage 1 equations and constants are unchanged.

Files: `predict_2.py` (prediction script), `constants_2.json` (constants of section 2), `constants.json`
(stage 1, unchanged), `predictions_2.csv` (completed, 4617 rows). Working scripts: `clib.py`, `iod.py`,
`g1.py`–`g5.py`, `ev.py`, `enc*.py`, `st*.py`, `c1law.py`, `c1prod.py`, `c7prod.py`, `prod.py`, `export2.py`,
`future.py`, and their `.log` files.

---

## 1. Test of the stage 1 laws

**S, M, P1–P5.** No new data; stage 1 stands (rms 0.9995′ per axis over 156 178 sightings). `predict_2.py`
reproduces the stage 1 predictions for these bodies exactly.

**Faint bodies, stage 1 law alone** (massless particle pulled by the 8 massive bodies; 6 free initial
conditions per body; rms residual per axis; noise is 1.0′):

| body | best single apparition (≤500 d), gravity only | verdict |
|---|---|---|
| C1 | 12.0′ (only apparition, 320 d) | fails |
| C2 | 1.16–1.18′ before day 16 020; 0.99–1.05′ after day 19 035 | fails slightly early, holds later |
| C3 | 22–37′ | fails |
| C4 | 1.5′ and worse before day ~13 000; 0.96–1.02′ after day 16 872 | fails early, holds later |
| C5 | 16–24′ before day 17 700; 0.99–1.04′ after day 20 606 | fails early, holds later |
| C6 | 1.02–1.08′ per apparition, but the three apparitions cannot be joined | marginal |
| C7 | 27–160′ before day 10 900; 1.0′ after day 11 700 | fails early, holds later |
| C8 | 30–173′ | fails |
| C9 | 10–36′ | fails |
| C10 | 3.6–4.6′ | fails |
| C11 | 23–32′ | fails |
| C12 | 2.2–27′ | fails |

Where it holds, it holds well: C2 over days 19 035–36 525 (6 apparitions, 1535 sightings) rms 1.021′;
C4 over 30 000–36 282 (3795 sightings) rms 1.023′; C5 over 20 606–34 296 (5 apparitions) rms 1.023′;
C7 over 33 362–36 525 (1990 sightings) rms 1.05′. Ordinary planetary perturbations are also confirmed:
C2 and C5 were thrown onto new orbits by passes 0.054 L and 0.046 L from P4 (days 17 369 and 19 277), with
the Tisserand quantity relative to P4 conserved (C2: 2.788 → 2.796), as gravity requires.

**What the failures look like.** Adding one parameter per body, a radial solar term
`+β·G m_S·r/|r|³`, brings every well-observed apparition to rms 0.93–1.05′:

* The force is radial and inverse-square: letting the exponent float gives k = 1.987 (C10), 2.004 (C8),
  2.000 (C7), 1.987 (C9); transverse and normal components are < 2e-4 of solar gravity (C10).
* β is the same from one apparition to the next as long as the body stays clear of planets, e.g.
  C3: −0.2303, −0.2302, −0.2304 (days 380–9834); C8: +0.2981, +0.2977, +0.2975, +0.2975, +0.2975;
  C11: between −0.1570 and −0.1580 in ten long apparitions (days 1471–19 483);
  C5: −0.1723, −0.1718, −0.1722, −0.1716 (days 1930–17 700).
* β has both signs (β > 0: S pulls less; β < 0: S pulls more), so it is not radiation pressure or outgassing
  in the usual sense, and it is not a mass of the faint body (planets are undisturbed).
* β changes only at close approaches to a planet, always towards zero:

| body | β before | β after | approach (distance, day) |
|---|---|---|---|
| C7 | +0.2174 | +0.1048 | Tellus 0.157 L, day 4712 |
| C7 | +0.1040 | +0.0003 | Tellus 0.068 L, day 10 940 |
| C1 | −0.2705 | ≈ −0.01 | Tellus 0.03 L, day 8149 |
| C5 | −0.172 | −0.0009 | P4 0.046 L, day 19 277 |
| C2 | +0.0028 | +0.0001 | P4 0.054 L, day 17 369 |
| C3 | −0.2303 | ≈ −0.06 (short arc) | P2 0.088 L, day 14 077 |
| C11 | −0.1577 | −0.1487 | P3 0.21 L, day 20 821 |
| C11 | −0.1487 | −0.116 | P2 0.20 L, day 33 247 |
| C4 | +0.27 (short arcs) | 0.0013, then 0.0000 | several Tellus passes (q = 0.83 L, in the ecliptic) |

  In days 0–11 000 the only two approaches of C7 to Tellus closer than 0.6 L are exactly its two changes.
  During a change position and velocity are continuous; β slides smoothly over the ~40 days of the pass
  (a sudden step at the best date leaves 3.6′; a step plus velocity kick 2.9′).

---

## 2. The laws, as equations

**Law 1 (unchanged from stage 1).** For the massive bodies i ∈ {S, T, M, P1…P5}:

    d²r_i/dt² = Σ_{j≠i} G m_j (r_j − r_i)/|r_j − r_i|³                                  (1)

with the stage 1 constants (G m_S = 2.9591e-4 L³/day², mass ratios and day-0 state vectors as in
`ANSWER_1.md` / `constants.json`). Observation equation (2) of stage 1 is unchanged: direction of
r_X − r_T. The faint bodies are massless (no detectable pull on anything).

**Law 2 (new).** A faint body c at position x obeys

    d²x/dt² = Σ_j G m_j (r_j − x)/|r_j − x|³  +  β_c(t) · G m_S · (x − r_S)/|x − r_S|³   (3)

i.e. full gravity of all 8 massive bodies, with the pull of S alone multiplied by (1 − β_c).

**Law 3 (new, empirical).** β_c fades near planets:

    dβ_c/dt = − A · β_c · Σ_p exp( −(d_p/R)² ),   d_p = |x − r_p|,  p ∈ {T, P1, P2, P3, P4, P5}   (4)

    A = 0.28 ± 0.02 day⁻¹        R = 0.0995 ± 0.003 L   (same for every planet, independent of its mass)

A and R come from the second Tellus pass of C7 (rms 0.995′ over days 9844–11 595; alternatives on the same
data: power law d⁻ⁿ 1.36′ with n = 4.9, exponential 1.05′, d⁻² 56′, d⁻³ 13′). Independent check: C1, with A
and R held at these values and only its 6 initial conditions and β free, fits at 1.05′ (gravity only 12.0′,
constant β 8.5′); changing A to 0.25 or 0.31 gives 1.7′ / 2.2′, changing R to 0.095 or 0.104 gives 1.1′ / 1.5′.
For S, M and the planets β = 0 (Kepler's third law between P1 and Tellus holds to 1e-5).

**Constants of the faint bodies** (state in the centre-of-mass frame of the 8 massive bodies, axes of the
data, L and L/day, at the integer epoch day; β at the epoch; exact values in `constants_2.json`):

| body | epoch day | x, y, z | vx, vy, vz | β (epoch) | ± | data used (days) | N | rms ′ |
|---|---|---|---|---|---|---|---|---|
| C1 | 8078 | +1.4395531105, −1.4452120335, +0.4788932198 | −0.002032939039, +0.017206152404, −0.005370572092 | −0.27046 | 0.02 | 8012–8331 | 205 | 1.050 |
| C2 | 33249 | −2.1671803223, −1.0870083891, +0.5706057476 | +0.006297224585, −0.007944410685, +0.008306931704 | +0.00007 | 0.0003 | 19035–36525 | 1535 | 1.021 |
| C3 | 14215 | +1.9618263450, −0.8438900148, +0.7407615252 | +0.014053829939, +0.008474586394, −0.001844763493 | −0.06359 | 0.03 | 14126–14297 | 123 | 0.979 |
| C4 | 35929 | −0.3067776192, +1.3910187111, −0.5816288168 | −0.013046089706, +0.002849126242, −0.001166229895 | −0.00000 | 0.0002 | 30000–36282 | 3795 | 1.023 |
| C5 | 30754 | +1.9616927264, −1.6599590686, +0.6575377250 | +0.008048288320, +0.009553590889, +0.000628467381 | −0.00094 | 0.0003 | 20606–34296 | 1451 | 1.023 |
| C6 | 34743 | −0.0001209920, −2.8146492095, +0.4037723282 | +0.011513190105, +0.000232562450, +0.001764410430 | +0.03505 | 0.02 | 34575–34907 | 222 | 0.961 |
| C7 | 36362 | +0.5938713621, −0.7177416008, +0.6346875079 | +0.015599575160, +0.006775241490, −0.000646021204 | −0.00009 | 0.0003 | 33362–36525 | 1990 | 1.054 |
| C8 | 36135 | −0.8789966445, +0.5931881188, +0.0991686337 | −0.009568847351, −0.015493659089, +0.004997101701 | +0.29754 | 0.0003 | 28848–36402 | 691 | 2.078 |
| C9 | 33968 | −1.1944709116, +0.2998567311, −0.0283481583 | −0.003819564103, −0.015856919841, +0.011181947781 | +0.08158 | 0.0003 | 29231–34221 | 627 | 1.034 |
| C10 | 27104 | −0.2175919736, −0.6894963828, +0.3990127061 | +0.024572704893, −0.009879337187, −0.001608332380 | −0.03316 | 0.0005 | 17048–27301 | 580 | 1.135 |
| C11 | 35629 | −2.2190622014, +0.0968093400, +0.0240893668 | −0.007595549222, −0.011296805243, +0.004231028975 | −0.11563 | 0.002 | 33390–35741 | 322 | 0.941 |
| C12 | 34674 | −1.7154114427, −1.3099337146, +0.4878531233 | +0.010629026210, −0.010421692004, +0.009266362717 | −0.29761 | 0.0005 | 25143–34857 | 450 | 1.018 |

The β uncertainties are estimated from the scatter of independent per-apparition fits (large where only one
short arc is available). I did not compute a formal covariance of the state vectors; from the fit
sensitivities they are determined to roughly 1e-6–1e-5 L and 1e-8–1e-7 L/day, much worse along the orbit for
the single-arc bodies C3 and C6.

**Numerical solution.** Massive bodies exactly as in stage 1 (6th-order symplectic, step 0.25 d from day 0).
Faint bodies: classical RK4, step 0.5 d, for (x, v, β) together, started from the epoch above.

**What forced the change.** Not the planets. Faint bodies cannot be fitted by (1) alone even over one
apparition; a per-body radial 1/r² term fixes each apparition; and the per-apparition β values are found to
jump only at planetary close approaches, which forces (4).

**Where laws 2–3 are themselves imperfect (numbers).**

* Law 3 is calibrated on Tellus only. For other planets it gives the right size but not the detail:
  C11's pass of P2 (0.20 L): predicted fade factor 0.85, observed 0.78, and the apparition containing the
  pass fits at 1.68′; a step fit puts the change ~100 days before closest approach. C6's two last
  apparitions joined through a fade from β ≈ 0.15 to 0.03 leave 2.6′.
* β is not perfectly constant away from planets: C8 drifts 0.2981 → 0.2975 over 100 years, and two of its
  apparitions 7000 days apart join only at 2.1′; C9's third-last apparition joins at 3.1′; C10's two
  apparitions at 1.14′. The law as written predicts no change in these intervals.
* After a deep pass β does not go fully to zero where (4) says it should (C5 keeps −0.0009).
* Whether M counts as a "planet" in (4) cannot be told from Tellus (always 0.0026 L apart).

---

## 3. Predictions

`predict_2.py` (numpy only, one process, about two minutes) integrates (1) from day 0 and (3)–(4) for each
faint body from its epoch to day 43 815, and fills `predictions_2.csv` (4617 rows, none empty).
`python predict_2.py --check` also reprints the rms of each faint body against `observations_2.csv` over the
data span listed above (it reproduces the last column of the table).

Expected quality, body by body:

| body | expectation for years 101–120 |
|---|---|
| S, M, P1–P5 | as stage 1 (noise-limited) |
| C2, C5 | good, about 1–3′: 5–6 apparitions fitted at noise, no planetary approach ahead |
| C4 | good until a predicted pass 0.023 L from Tellus near day 40 205, degraded after it |
| C7 | good until a predicted pass 0.054 L from Tellus near day 38 330; unreliable afterwards (changing β by 1e-4 moves later positions by degrees) |
| C12, C9, C10 | fair: two apparitions joined at 1.0–1.14′; arcminutes to tens of arcminutes |
| C8 | fair: ~2–3′ at first; timing of the next return (≈ day 43 100) uncertain by about a day |
| C11 | fair to poor: only 322 sightings since its last change of β |
| C6 | poor: one 330-day arc, period uncertain by ~1 % |
| C3 | poor: one 170-day arc on a very long orbit (a ≈ 48 L); it is ~70 L away, moving slowly |
| C1 | fair: unbound and ~300–400 L away; direction depends on its small residual β |

---

## 4. Self-assessment

| law / claim | confidence | what would prove it wrong |
|---|---|---|
| Stage 1 law (1) for S, M, planets | 97 % (unchanged) | as in stage 1 |
| Faint bodies feel the ordinary gravity of S and all planets | 95 % | a planetary pass (like C5/C2 at P4) whose deflection disagrees with the planet's stage 1 mass |
| Each faint body has an extra radial 1/r² solar term β (law 2) | 92 % | an apparition far from planets that needs a different exponent, a transverse force, or a β different from the previous apparition |
| β is a property of the body, not of its orbit | 85 % | two bodies on the same orbit needing the same β, or β tracking a, e or i (C5 kept β = −0.172 while its orbit changed by 5 %) |
| β changes only near planets and only towards zero | 80 % | a change of β with no planet within ~0.5 L, or |β| increasing |
| Fading law (4) with a Gaussian of R ≈ 0.10 L for Tellus | 75 % (form), 60 % (A, R within errors) | a Tellus pass at 0.1–0.2 L that changes β by a different factor than exp(−A∫e^{−d²/R²}dt) |
| Same A and R for every planet, independent of mass | 45 % | already strained: C11/P2 timing, C6, slow drift of C8 |
| β of S, M, planets is zero | 90 % | third-law mismatch between planets above 1e-5 |
| Predictions: C2, C5 within a few ′ | 85 % | — |
| Predictions: C4, C7 before their Tellus passes | 80 %; after them 30 % | — |
| Predictions: C8, C9, C10, C12 within 1° | 70 % | — |
| Predictions: C1, C3, C6, C11 within 1° | 40 % | — |

Caveats, stated plainly:

* I have no physical explanation for β or for why planets erase it; laws 2–3 are a description that fits,
  not a mechanism. An alternative I could not exclude: β multiplies the pull of every massive body, not only
  of S (the two differ only during close planetary passes).
* Law 3's functional form was chosen among a handful of trial forms on one event and checked on one other
  (C1); it is not unique.
* For C6 and C8 I had to choose between a single-apparition fit and a two-apparition fit that does not
  reach the noise; I used the single apparition for C6 (a β change lies between its apparitions) and the
  joined fit for C8 (better period, 2.1′ misfit). C3's β after its P2 pass rests on 123 sightings.
* My least-squares search is fragile across β-changing events; several long-arc fits diverged and were
  replaced by fits restricted to the interval after the last change. Earlier history of each body is
  therefore described per apparition (section 1), not by one continuous trajectory.
