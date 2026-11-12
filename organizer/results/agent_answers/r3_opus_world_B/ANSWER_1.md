# Tellus sky: laws of motion

Short version: one law explains everything. All eight bodies (S, P1–P5, Tellus and M) attract each other by Newton's inverse-square gravitation. A single 8-body integration with 49 fitted constants reproduces all 156,178 sightings over the 100 years to the measurement noise: rms 1.00′ per angle for every body, with no trend from decade to decade.

Units: length in AU (defined by GM_S = k², k = 0.01720209895 rad/day, so Tellus's semi-major axis comes out ≈ 1.0003 AU), time in days, day numbers as in the data.

---

## 1. Description of the motions

**Frame.** S stays on a great circle (the *ecliptic*). Its pole lies at data-frame lon = 90.00°, lat = 90° − 21.7015°, so the ecliptic is tilted by **ε = 21.7015° ± 0.0010°** to the data frame's equator. Over the 100 years the pole stays fixed to within 0.01°. All the analysis below uses ecliptic coordinates, obtained by rotating the data frame about its x-axis by ε.

**S.** S moves eastward along the ecliptic once per **365.257 d** (sidereal year). Its motion is not uniform: the equation of centre gives eccentricity **e = 0.01999**. A Kepler ellipse alone already fits S to 1.02′. Two further effects are detectable:
- The line of apsides advances **+0.26° per century**.
- A 6″ wobble with the 29.5-day synodic month comes from Tellus orbiting the Tellus–M barycentre.
- Weak signals at ≈390 d and ≈438 d are the synodic periods of P4 and P3 (planetary pulls on Tellus).

**P1** never gets more than 29° from S (inferior wanderer). It orbits S with **P = 110.26 d, a = 0.450, e = 0.060**, in the ecliptic plane (i ≈ 0.001°). Synodic period 158 d.

**P2, P3, P4, P5** reach every elongation and make retrograde loops at opposition (outer wanderers). Heliocentric periods and orbits:

| body | P (d) | a | e | i |
|---|---|---|---|---|
| P2 | 809.6 | 1.700 | 0.050 | 0.008° |
| P3 | 2189.6 | 3.300 | 0.040 | 1.49° |
| P4 | 5636 | 6.200 | 0.050 | 1.00° |
| P5 | 14201 | 11.48 | 0.042 | 2.00° |

Their synodic periods are 665, 438, 390 and 375 d.

**Kepler's 3rd law holds across all of them.** Separate Kepler-ellipse fits give n·a^{3/2} = 0.017202 ± 0.00001 for every body, from a = 0.45 to a = 11.5, and the same value as Tellus. So one central attraction ∝ 1/r² from S governs them all; a power-law index different from 2 by more than ~1e-5 is excluded. Pure Kepler ellipses still leave residuals of 2–8′, and those are exactly the mutual perturbations that the N-body law accounts for.

**M** circles Tellus (it is a moon, not a heliocentric body).
- Sidereal month 27.29 d, synodic month 29.49 d.
- Mean eccentricity ≈ 0.055 and mean inclination to the ecliptic ≈ 5.1°. The inclination oscillates between 4.9° and 5.2°.
- The node regresses with a period of **6800 d (18.6 yr)**.
- The perigee advances with a period of **2190 d (6.0 yr)**.
- On top of the Kepler motion there are large solar perturbations (evection, variation, annual equation). A precessing Kepler ellipse leaves 44′ rms. The Newtonian 3-body (Tellus–M–S) motion removes all of it.

**Changes over the 100 years.** All of these are reproduced by the same law. The numbers are osculating elements from the fitted solution, day 0 → day 36525:
- Perihelion longitudes drift: P1 +0.1°, P2 +1.1°, Tellus +0.26° (mean, from the S fit), P3 −1.4°, P4 +0.6°, P5 +4°.
- Eccentricities change slowly. For example, P3 goes from 0.0400 to 0.0428.
- Inclinations drift by a few thousandths of a degree.
- P4 and P5 exchange energy, because their periods are near a 5:2 ratio (5·5636 ≈ 2·14201). P5's osculating period falls from 14242 d to 14157 d, and P4's oscillates by about 10 d.
- M's node turns more than 5 full times, its perigee more than 16 times, and its inclination and eccentricity oscillate with the solar perturbation periods.
- Nothing outside Newtonian gravitation is needed. The residuals are white noise of 1.00′ in every decade.

## 2. The laws (equations) and constants

**Law 1: gravitation and motion.** For the 8 point masses i ∈ {S, P1, P2, T(Tellus), P3, P4, P5, M}, in an inertial frame with ecliptic axes:

  d²**r**_i/dt² = Σ_{j≠i} GM_j (**r**_j − **r**_i) / |**r**_j − **r**_i|³

with GM_j = μ_j · k², where k = 0.01720209895 AU^{3/2} day⁻¹ (so μ_S = 1).

**Law 2: what is seen.** The direction of body b on the sky of Tellus is the unit vector

  **u**_b = (**r**_b − **r**_T)/|**r**_b − **r**_T|

converted to the data frame:

  X = x, Y = cos ε·y + sin ε·z, Z = −sin ε·y + cos ε·z,
  lon = atan2(Y, X), lat = asin(Z), with **ε = 21.7015°**.

There are no other ingredients: no light-time, no aberration, no rotation effects (they are not needed at the 1′ level).

**Constants.**

Mass ratios μ = GM/GM_S. The ± values are the realistic uncertainty: the formal 1σ, inflated where the value still moved between fit iterations.

| body | μ | ± |
|---|---|---|
| P1 | 2.7e-7 | 1.0e-6 (barely detected; < 1.3e-6 at 1σ) |
| P2 | 7.08e-7 | 2e-7 |
| T (Tellus) | 2.608337e-6 | 1e-11 (formal 6.6e-12) |
| P3 | 1.5057e-4 | 4e-7 |
| P4 | 9.9992e-4 | 5e-7 |
| P5 | 2.5013e-4 | 5e-7 |
| M | 3.19492e-8 | 1e-11 (formal 7.4e-12) |

Derived from these: μ_T + μ_M = 2.640286e-6 and M/T mass ratio = 0.012249 (≈ 1/81.6).

Initial conditions at epoch **T0 = day 18262.0**. Heliocentric ecliptic positions are in AU and velocities in AU/day. The last row is M relative to Tellus. The ± values are formal 1σ.

| body | x | y | z | vx | vy | vz |
|---|---|---|---|---|---|---|
| P1 | 0.091070396 ± 6.6e-06 | 0.438179980 ± 2.2e-06 | 0.000006877 ± 3.0e-06 | -0.025512514 ± 1.4e-07 | 0.003728893 ± 2.7e-07 | -0.000000385 ± 1.7e-07 |
| P2 | 1.211335998 ± 4.1e-06 | -1.090270208 ± 4.8e-06 | 0.000222556 ± 3.5e-06 | 0.009488709 ± 2.7e-08 | 0.009951700 ± 2.6e-08 | 0.000000018 ± 3.0e-08 |
| T | -0.407771132 ± 2.6e-06 | -0.932364170 ± 1.6e-06 | -0.000005381 ± 1.3e-06 | 0.015420763 ± 2.5e-08 | -0.006924138 ± 3.2e-08 | 0.000000211 ± 2.4e-08 |
| P3 | 0.065770598 ± 1.1e-05 | -3.380320231 ± 5.8e-06 | 0.040289814 ± 8.6e-06 | 0.009225631 ± 2.0e-08 | 0.000470248 ± 2.0e-08 | -0.000217365 ± 2.3e-08 |
| P4 | 1.317547352 ± 2.0e-05 | 6.270565710 ± 1.1e-05 | -0.111775506 ± 1.7e-05 | -0.006586106 ± 2.0e-08 | 0.001135558 ± 1.5e-08 | 0.000014322 ± 1.8e-08 |
| P5 | -10.981320636 ± 1.6e-05 | -0.569985028 ± 3.5e-05 | 0.360460147 ± 2.7e-05 | 0.000226952 ± 1.8e-08 | -0.005287970 ± 7.2e-09 | -0.000061929 ± 1.4e-08 |
| M−T | 0.002381410476 ± 3.2e-09 | -0.000106566137 ± 7.2e-09 | -0.000177535341 ± 5.4e-09 | 0.000046257349 ± 1.4e-09 | 0.000573686190 ± 7.9e-10 | 0.000027408071 ± 5.4e-10 |

Full-precision values, the ones predict.py uses, are in the `#CONSTANTS_BEGIN` block of `predict.py`; they are the same numbers to more digits. S starts at the heliocentric origin, and everything is moved to the barycentre before integrating.

Equivalent osculating elements at T0 (heliocentric; geocentric for M):
- P1: a 0.45000, e 0.0600, P 110.258 d
- P2: a 1.69996, e 0.0503
- T: a 1.00029, e 0.0199
- P3: a 3.3002, e 0.0401, i 1.492°
- P4: a 6.1999, e 0.0495, i 1.003°
- P5: a 11.478, e 0.0423, i 1.999°
- M: a 0.0024287 AU, e 0.0357, i 5.136°

**How the constants were obtained.**
1. ε comes from the plane of S's track.
2. Kepler-ellipse fits of each body give starting values.
3. A Gauss–Newton least-squares fit of all 49 constants (42 state components + 7 mass ratios) against all 156,178 directions. The 8-body equations are integrated with DOP853 at rtol = 1e-12, forward and backward from T0, and the Jacobians come from finite differences.

Fit quality (rms per angle): S 1.002′, P1 0.997′, P2 1.000′, P3 0.999′, P4 0.998′, P5 0.997′, M 1.001′. That equals the stated measurement error, in each decade separately.

## 3. Predictions

`predict.py` holds the constants above and integrates Law 1 from T0 to each requested day. It applies Law 2 and writes lon/lat (degrees) into `predictions.csv` (1701 rows, years 101–120). Checked against the data, the script reproduces the observed directions with 1.41′ rms total angle (that is 1.00′ per angle), both early and late in the record.

Propagating the fit covariance gives the formal 1σ total-angle error of the predictions over years 101–120:
- S ≈ 0.02′
- P1 ≈ 0.03′ (max 0.09′ near inferior conjunction)
- P2 ≈ 0.02′ (max 0.07′ at opposition)
- P3, P4, P5 ≈ 0.02–0.035′
- M ≈ 0.22′ (max 0.28′)

Realistically I expect the errors to stay well under 1′ for the planets and S, and under about 1′ for M.

## 4. Self-assessment

| law | confidence | what would prove it wrong |
|---|---|---|
| Inverse-square mutual gravitation (Law 1) governs all motions | 97% | Residuals in years 101–120 that grow systematically above ~1′. Examples: P5 or P4 drifting in longitude by several arcmin (a missing body or a wrong mass), or a perihelion advance not produced by the planets (a non-Newtonian term). |
| Directions are the plain geometric Tellus-centred vectors in a fixed frame tilted by ε (Law 2) | 95% | A systematic offset of order 20″ correlated with Tellus's velocity (aberration) or with distance (light-time) in better data. At 1′ noise it is invisible here. |
| Values of the large masses (P3, P4, P5) and of the Tellus–M system | 90% | Future positions of the outer planets deviating by more than ~0.2′. That would mean the masses, or the near 5:2 P4–P5 interaction, are mis-set. |
| Masses of P1 and P2 | 50% | They are barely detected: P1 is consistent with 0–1.3e-6. Their effect on predictions is below 0.1′, so this does not matter much for task 3. |
| The ecliptic is fixed in the data frame (no precession of the frame) | 95% | The plane of S's path turning by more than 0.01° in years 101–120. |
| Predictions accurate to ≪ 1′ (planets, S) and < 1′ (M) | 85% | Any predicted direction off by more than 2–3′, beyond what noise explains. |
