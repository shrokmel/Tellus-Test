# Stage 3: the new body N

Summary. N obeys the faint-body laws of stage 2, with one revision that N forced me to check and that the
old data then decided: the factor (1 − β) multiplies the pull of **every** massive body, not only of S, and M
counts among the bodies that erase β. N has β = −0.3234 ± 0.0021 (it falls 32 % harder than a planet would).
It is on a bound, steep orbit and is heading for a close pass of P4, **0.035 L (0.027–0.051) from P4 on day
36 966**, 1.2 years after the last sighting. That pass turns its path by about 11°, strips most of its β
(to about −0.036) and makes everything afterwards much less certain than the 292-day arc alone would suggest.

Files: `predict_3.py`, `constants_3.json` (constants below, full precision), `constants.json` (stage 1,
unchanged), `predictions_3.csv` (1461 rows, all filled). Working scripts and logs: `n1.py`–`n4.py`,
`cal5.py`, `cal6.py`, `cal2.py`.

---

## 1. Laws used for N

**Law 1 (unchanged since stage 1).** Massive bodies i ∈ {S, T, M, P1…P5}:

    d²r_i/dt² = Σ_{j≠i} G m_j (r_j − r_i)/|r_j − r_i|³                               (1)

with G m_S = 2.9591e-4 L³/day² and the masses and day-0 state vectors of `ANSWER_1.md` / `constants.json`.
Direction seen from Tellus = direction of r_X − r_T.

**Law 2 (revised).** A faint body at x, massless, with its own number β:

    d²x/dt² = (1 − β) · Σ_j G m_j (r_j − x)/|r_j − x|³        (sum over all 8 massive bodies)    (2)

**Law 3 (revised constants and membership).**

    dβ/dt = − A · β · Σ_p exp( −(d_p/R)² ),   d_p = |x − r_p|,   p ∈ {T, M, P1, P2, P3, P4, P5}        (3)

    A = 0.138 ± 0.003 day⁻¹          R = 0.0998 ± 0.0005 L

**Constants of N** (centre-of-mass frame of the 8 massive bodies, axes of the data, epoch day 36 379):

| quantity | value | ± (1σ) |
|---|---|---|
| x, y, z (L) | +0.4744597, +2.2181417, +1.2134757 | 1.1e-4, 2.7e-4, 2.2e-4 |
| vx, vy, vz (L/day) | −0.01566120, +0.00080903, +0.00444700 | 9.4e-7, 1.7e-6, 9.9e-7 |
| β at epoch | −0.32343 | 0.0021 |

(The seven numbers are strongly correlated: β with position at −0.72 to −0.82. The covariance is in
`n_cov.npy`.) Fit: 207 sightings, days 36 233–36 525, rms 0.988′ per axis (noise 1.0′).

**What changed since stage 2, and why.**

* *Stage 2 said* β multiplies the pull of S only, and I wrote that I could not exclude the alternative
  "every body". N's coming pass of P4 made the difference matter, so I tested it on the two old bodies that
  made such a pass. Carrying one trajectory through C5's pass of P4 (0.055 L, day 19 273), fitted to the
  apparitions before and after it together:

  | version | rms over days 14 469–24 353 |
  |---|---|
  | stage 2 as written (β on S only, A = 0.28, R = 0.0995) | 35′ |
  | β on S only, A and R free | 30′ |
  | β on every body, A = 0.28 | 7.8′ |
  | β on every body, A and R free | 1.06′ (A = 0.134, R = 0.1007) |
  | **β on every body, A = 0.138, R = 0.0998 fixed** | **1.04′** (before 1.02′, after 1.06′) |

  The last version also leaves β = −0.0008 after the pass; the value measured independently from the later
  apparitions is −0.0009. C2's pass of P4 (0.054 L, day 17 369) gives 1.005′ with the same fixed constants.
* *Why A halved.* The stage 2 value A = 0.28 came from passes of Tellus, where M is always 0.0026 L away. P4,
  which has no companion in the model, needs 0.134. So M counts as a body in (3), and the constant per body
  is half the old one. Refitted with M included: C7's second Tellus pass gives A = 0.1377 ± 0.0008,
  R = 0.0998 ± 0.0002 (rms 0.996′); C1's Tellus pass gives A = 0.144 ± 0.007, R = 0.096 ± 0.004 (rms 1.01′).
  The three agree, so A and R are the same for Tellus, M and P4 despite a mass range of 1 : 27 000.
* *What N itself adds.* Gravity alone (stage 1 law) fits N at 7.8′ and makes it unbound (e = 1.26). One
  constant β brings it to 0.988′ and makes it bound. Letting the radial term have a free exponent or
  transverse parts does not improve the fit (0.989′, 0.988′), and N stays > 1.4 L from every planet during
  the arc, so N by itself cannot distinguish the stage 2 and stage 3 versions. β = −0.323 is the largest
  magnitude yet seen (previous range −0.298 … +0.298).
* The stage 2 answer and its predictions are not changed here; for bodies that pass no planet the two
  versions predict the same thing.

---

## 2. Predictions

`predict_3.py` (numpy only, one process, about two minutes) integrates (1) from day 0 with the stage 1
scheme and (2)–(3) for N by RK4 with a 0.5-day step from day 36 379 to day 43 830, and writes
`predictions_3.csv`. `python predict_3.py --check` also prints the rms against `observations_3.csv` (0.988′).

---

## 3. N's path

Now (day 36 525): N is 2.98 L from S and receding. It passed perihelion (q = 2.57 L) near day 36 380.
Its orbit is steep: inclined 53.7° to the plane of Tellus's orbit (33.1° to the equator of the data frame),
a = 10.1 L, e = 0.745. Because S pulls it 1.32 times harder than a planet, its period would be 10 180 days
(27.9 years); the same orbit under plain gravity would be unbound.

| day (years after last sighting) | r from S (L) | β | a (L) | e | q (L) | incl. to ecliptic | sky lon, lat |
|---|---|---|---|---|---|---|---|
| 36 525 (0) | 2.98 | −0.323 | 10.0 | 0.744 | 2.57 | 53.7° | 116.6°, +21.7° |
| 36 890 (1.0) | 5.40 | −0.323 | 10.3 | 0.749 | 2.58 | 53.9° | 179.5°, +5.8° |
| 36 966 (1.21) | 5.9 | fading | — pass 0.035 L from P4 — | | | | |
| 37 300 (2.1) | 7.90 | −0.036 | 10.0 | 0.799 | 2.01 | 62.8° | 196.8°, −5.9° |
| 38 350 (5.0) | 12.4 | −0.036 | 10.0 | 0.798 | 2.01 | 62.9° | 216.8°, −19.1° |
| 40 175 (10.0) | 16.7 | −0.036 | 10.0 | 0.800 | 2.00 | 63.0° | 232.1°, −28.7° |
| 42 005 (15.0) | 18.0 (aphelion) | −0.036 | 10.0 | 0.800 | 2.00 | 62.8° | 243.6°, −34.7° |
| 43 830 (20.0) | 16.7, returning | −0.036 | 10.0 | 0.798 | 2.02 | 62.8° | 257.0°, −39.9° |

How the motion changes:

1. **Year 0–1.2: outbound, as now.** N slows from 0.0150 to 0.0100 L/day as it climbs away from S, drops
   towards the ecliptic and moves east across the sky, with the yearly back-and-forth caused by Tellus's own
   motion growing smaller as N recedes.
2. **Day 36 966 ± 1: the P4 pass.** Closest distance 0.035 L (68 % range 0.029–0.043, 95 % range
   0.027–0.051), deep inside P4's sphere of influence, at 0.0104 L/day relative speed. P4 turns N's path by
   about 11°: the inclination to the ecliptic rises from 54° to 63° and the perihelion distance drops from
   2.57 to 2.0 L. During the roughly 18 days spent within 0.1 L of P4, β fades from −0.323 to about −0.036
   (95 % range −0.047 … −0.030); the pass is too fast to erase it completely.
3. **After the pass: an almost ordinary body.** With S pulling only 3.6 % harder than normal, N coasts out
   on an orbit with a ≈ 10.0 L (68 % range 8.9–11.4), e ≈ 0.80, period ≈ 11 400 days (9500–13 700). It keeps
   slowing, reaches aphelion near 18 L (15.6–20.7) around year 15, moving only 0.002 L/day, and by year 20
   has turned back towards S. It stays bound in all 600 trial orbits. On the sky it drifts slowly south-east,
   from near the celestial equator to latitude −40°, with a yearly parallax loop of about ±3°.
4. It does not come near any other planet in these 20 years (after P4, nothing closer than about 2 L), so
   β should then stay at its post-pass value.

---

## 4. Self-assessment for N

Expected error of the predicted direction. The statistical part comes from 600 trial orbits drawn from the
fit covariance (and from the stated errors of A and R) and carried through the P4 pass with the full
equations:

| time after last sighting | median | 68 % | 95 % | my overall expectation |
|---|---|---|---|---|
| 1 year (before the P4 pass) | 0.05° | 0.08° | 0.15° | **about 0.1°** |
| 5 years | 1.7° | 2.4° | 6.8° | **about 2–3°** |
| 20 years | 7.5° | 10° | 30° | **about 10°, possibly 30°** |

The jump between year 1 and year 2 (0.05° → 0.7°) is entirely the P4 pass: a change of 0.01 L in the miss
distance changes the turn angle and the amount of β that survives. For comparison, the stage 1 law alone
would be off by 10°, 21° and 48° at 1, 5 and 20 years, and the stage 2 laws exactly as I wrote them by
0.01°, 1.3° and 4.5°.

Confidence:

| claim | confidence |
|---|---|
| N needs β ≈ −0.32; gravity alone is wrong for it | 97 % |
| β multiplies the pull of all bodies, and M counts in law 3 | 85 % |
| A and R as stated apply to P4 | 85 % (two old P4 passes fit at 1.0′) |
| N passes within 0.06 L of P4 on day 36 966 ± 2 | 90 % |
| N stays bound to S after the pass | 90 % |
| Predictions within 0.2° through day 36 940 | 90 % |
| Predictions within 5° at 5 years | 75 % |
| Predictions within 15° at 20 years | 65 % |

What would prove me wrong:

* Before the pass: N more than ~0.3° from prediction during the first year would mean β is not constant
  away from planets or the radial term is not inverse-square for N.
* At the pass: N missing P4 by more than 0.06 L, or arriving more than 2 days off.
* After the pass: an orbit that needs β still near −0.32 (no fading), or exactly 0 (stronger fading than
  law 3 gives for a fast pass: N's is 2.4 times faster than the passes of C5 and C2 on which the law was
  checked, so this is the least-tested part), or a turn angle that matches P4's unscaled pull rather than
  (1 − β) times it.
* Later: any further change of β with no planet within ~0.3 L.

Caveats:

* The prediction is the single nominal trajectory. Because the pass is strongly non-linear, the nominal is
  not the centre of the cloud of possible outcomes; the table above measures the scatter about it.
* Uncertainty in P4's own position and mass (stage 1) is not in the trial orbits; it is small (about
  1e-5 L and 5e-4 relative) next to N's own 0.01 L uncertainty at the pass.
* Law 3 remains a fitted description. Its Gaussian shape was chosen on Tellus passes at 0.07–0.16 L and
  0.03 L and on P4 passes at 0.054–0.055 L; N's pass at 0.035 L and higher speed is inside that range of
  distances but outside the range of speeds.
* I have no mechanism for β. The stage 2 weaknesses listed there (slow drift of C8's β, C11's change near
  P2 not fitting in detail) were not re-examined with the revised laws.
