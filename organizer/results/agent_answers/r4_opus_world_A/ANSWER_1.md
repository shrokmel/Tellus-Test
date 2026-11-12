# ANSWER 1 — Laws of motion of the Tellus sky

Short version: **S is the centre. Tellus and P1–P5 orbit S under Newton's inverse-square gravity, and they pull weakly on one another. M orbits Tellus, also under inverse-square gravity, and S's tidal pull disturbs that orbit strongly.** One law (Newtonian gravity with six free mass ratios) plus geometry fits all 156,178 sightings of all seven bodies. The residual is 1.00′ rms per angle, which equals the stated measurement error. No residual trend shows up in any decade.

Scripts in this folder: `fitS.py` and `fitP.py` (single-body Kepler fits), `nbody.py` and `fitnb.py` (Newtonian N-body fit), `moon.py` (M fit), `diag.py`, `elements2.py` and `compare.py` (checks), `make_predict.py` (writes the fitted constants into `predict.py`), and `predict.py` (the predictions).

---

## 1. Description of the motions

**The frame.** The directions of S lie on one great circle to 1′ (rms of S's latitude about that circle = 1.0′). The circle is tilted **ε = 21.7015°** to the frame's equator, and its ascending node is at lon = 270° (descending node at 90°). The plane stays fixed against the stars to < 0.003° over 100 years. I call it the *ecliptic*. From here on, coordinates are ecliptic coordinates: the frame rotated by ε about its x-axis.

**S.** S moves eastward round the ecliptic once every **365.2575 days** (one year). Its speed is not uniform: it follows Kepler's equation for an orbit of eccentricity 0.0200, and its fastest point is at ecliptic longitude 273.9°. A pure Kepler ellipse fits S's 25,525 sightings to 1.00′ (the noise) over the full 100 years. The apsides do not drift measurably (< 0.3° per century).

**P1** never gets more than 29° from S (an "inner" body). It swings from one side of S to the other with a synodic period of about 158 d. **P2–P5** reach every elongation (they are "outer" bodies), and each makes a backward (retrograde) loop around opposition. The P's stay within a few degrees of the ecliptic: P1 and P2 within 0.08°, P3 within ±2.1°, P4 within ±1.3°, P5 within ±2.2°.

**Heliocentric reduction.** Treat S's motion as the reflection of Tellus moving round S on that Kepler ellipse. Then each P is a body moving on its own ellipse round S. Each fitted ellipse obeys Kepler's third law, P²/a³ = constant, with the same constant as Tellus to 4×10⁻⁶ (P1), 3×10⁻⁵ (P2), 5×10⁻⁵ (P3), 6×10⁻⁴ (P4) and 9×10⁻⁴ (P5). The orbits, in units of Tellus's orbit size, are:

| body | a | sidereal period (d) | e | i (°) | Ω (°) | ϖ (°) |
|---|---|---|---|---|---|---|
| Tellus | 1.0000 | 365.25 | 0.0200 | 0 (defines plane) | – | 93.7 (perihelion) |
| P1 | 0.4500 | 110.26 | 0.0599 | 0.0006 | – | 166.4 |
| P2 | 1.6999 | 809.6 | 0.0500 | 0.002 | – | 281.2 |
| P3 | 3.3004 | 2190.1 | 0.0400 | 1.499 | 118.7 | 42.7 |
| P4 | 6.2072 | 5646 | 0.0492 | 0.999 | 163.2 | 212.1 |
| P5 | 11.503 | 14250 | 0.0400 | 2.000 | 73.2 | 167.3 |

(These are osculating elements at day 0 of the N-body solution, from `elements2.py`.)

**Changes over the 100 years (the motions are not fixed ellipses).**
- Single fixed ellipses fit S and P1 to the noise. They leave residuals of 1.8′ for P2, 8′ for P3, 2′ for P4 and 3′ for P5, with systematic structure, and P4 and P5 break Kepler's third law at the 10⁻³ level.
- P4 and P5 trade orbital energy. Their period ratio is about 2.52, close to the 5:2 commensurability. Over the century, P5's mean distance falls from 11.503 to 11.456 (its period shortens from 14,250 d to 14,163 d), while P4's rises from 6.207 to 6.211 (5,646 d to 5,651 d).
- P3, which lies inside the giant P4, has its eccentricity rise from 0.040 to 0.043 and its perihelion wander by a few degrees. P5's perihelion advances about 5° per century. Inclinations drift by about 0.01° per century.
- Tellus's own orbit is almost unchanged: a, e and ϖ are constant to < 10⁻⁴ and < 0.3°.
- All of these changes follow from mutual gravity (Law 2); no extra law is needed.

**M.** M goes round the sky once every **27.2897 d** (sidereal). Its path is tilted **5.03°** to the ecliptic. Its node moves *backwards* round the ecliptic once every **18.61 years**. The tilt itself oscillates between 4.86° and 5.17° with S's position. Its non-uniform speed corresponds to an orbit of eccentricity about 0.055 whose perigee moves *forwards* once every **8.86 years** (anomalistic month 27.52 d). M's motion also shows the large, S-locked inequalities expected from solar tides (evection, variation, annual equation), (these are implied by the fit, not separately modelled). All of this comes out of one law: Tellus's gravity plus S's tide (Law 3). Over 100 years, that model fits M's 23,565 sightings to 1.00′ rms with no growing residual.

---

## 2. The laws as equations

Notation: bold letters are 3-vectors. Time t is in days from day 0. Lengths are in units **L**, defined by the fitted orbits; Tellus's mean distance from S is 0.99992 L. Only directions are observed, so the absolute length scale is never determined. Only the combinations below are.

### Law 1 — Observation geometry (frame)
For an observed direction (lon λ, lat β), let **u** = (cos β cos λ, cos β sin λ, sin β). Ecliptic direction **e** = R **u**, where

R = [[1, 0, 0], [0, cos ε, −sin ε], [0, sin ε, cos ε]],  **ε = 21.7015° ± 0.0010°**.

Predicted directions are geometric (instantaneous, no light-time):
- S: **e**_S ∝ −**r**_T
- Pk: **e**_Pk ∝ **r**_Pk − **r**_T
- M: **e**_M ∝ **ρ**

Here **r** are positions relative to S and **ρ** is M's position relative to Tellus.

### Law 2 — Newtonian gravity for S, Tellus and P1–P5 (heliocentric form)
For each i ∈ {T, P1, P2, P3, P4, P5}:

d²**r**_i/dt² = −G(M_S + m_i) **r**_i / |**r**_i|³ + Σ_{j≠i} G m_j [ (**r**_j − **r**_i)/|**r**_j − **r**_i|³ − **r**_j/|**r**_j|³ ]

Kepler's three laws, and every slow change listed in section 1, are consequences of this equation. In the limit m_j → 0, n² a³ = G(M_S + m_i).

Constants (uncertainties are 2σ formal, i.e. about the realistic 1σ):

| constant | value |
|---|---|
| G M_S | 2.95847239 × 10⁻⁴ L³ d⁻² (fixes the unit L; equivalent to Tellus's mean motion n_T = 0.0172020684 ± 2×10⁻¹⁰ rad/d, year = 365.25755 d) |
| m_T / M_S (Tellus + M) | 0 (fixed). A free fit gives (−0.6 ± 0.6) × 10⁻⁶, i.e. undetectable, < 10⁻⁶ |
| m_P1 / M_S | (1.19 ± 0.17) × 10⁻⁶ |
| m_P2 / M_S | (1.04 ± 0.24) × 10⁻⁶ |
| m_P3 / M_S | (1.5173 ± 0.0010) × 10⁻⁴ |
| m_P4 / M_S | (1.00104 ± 0.00011) × 10⁻³ |
| m_P5 / M_S | (2.4911 ± 0.0015) × 10⁻⁴ |

Initial conditions at t = 0 (heliocentric ecliptic). Positions are in L; velocities are in L/day.

| body | x | y | z | vx | vy | vz |
|---|---|---|---|---|---|---|
| Tellus | -0.394277928 ± 7.5e-06 | -0.938258180 ± 6.4e-06 | -0.000025053 ± 3.7e-06 | 0.01551744151 ± 1.3e-07 | -0.00668738297 ± 1.0e-07 | -0.00000018402 ± 6.1e-08 |
| P1 | -0.312556119 ± 1.1e-05 | -0.301675583 ± 1.6e-05 | -0.000003540 ± 5.9e-06 | 0.01747725805 ± 6.8e-07 | -0.01997905454 ± 8.0e-07 | -0.00000019704 ± 3.6e-07 |
| P2 | -0.714848618 ± 1.4e-05 | 1.633025048 ± 1.4e-05 | 0.000044890 ± 8.0e-06 | -0.01145280936 ± 8.8e-08 | -0.00516845523 ± 8.5e-08 | -0.00000009189 ± 5.6e-08 |
| P3 | -3.056443320 ± 2.6e-05 | 1.364972191 ± 3.0e-05 | 0.053062767 ± 1.7e-05 | -0.00411733989 ± 7.5e-08 | -0.00837384695 ± 6.5e-08 | 0.00019967731 ± 4.8e-08 |
| P4 | 6.460515046 ± 5.0e-05 | -0.028502913 ± 6.5e-05 | -0.032005475 ± 3.4e-05 | 0.00021085103 ± 6.0e-08 | 0.00662642304 ± 5.3e-08 | -0.00011164893 ± 3.5e-08 |
| P5 | 3.323334757 ± 1.0e-04 | 11.027423989 ± 8.7e-05 | 0.000174464 ± 6.1e-05 | -0.00490185888 ± 3.7e-08 | 0.00126578999 ± 4.4e-08 | 0.00017668097 ± 2.7e-08 |

The exact values used, to full precision, are written into `predict.py`. The positions and velocities are correlated; only the full set reproduces the fit. Equivalent osculating elements are given in section 1.

### Law 3 — M orbits Tellus under Tellus's gravity plus S's tide
d²**ρ**/dt² = −G m_T **ρ**/|**ρ**|³ + G M_S [ (**R** − **ρ**)/|**R** − **ρ**|³ − **R**/|**R**|³ ]

**R**(t) is S's position relative to Tellus. It follows the two-body Kepler orbit d²**R**/dt² = −G M_S **R**/|**R**|³, which is identical, to 1′ in direction, to the Kepler ellipse of S:
- n = 0.0172020684 ± 2×10⁻¹⁰ rad/d
- mean anomaly at t = 0: M₀ = 152.2614° ± 0.007°
- e = 0.0199790 ± 0.0000026
- longitude of perigee ϖ_S = 273.9098° ± 0.007°

Only directions of M are measured, so the equation is invariant under **ρ** → k**ρ**, G m_T → k³ G m_T. I fix the unit A = 0.0025 L and G m_T = (2π/27.3217 d)² A³ = 8.2635 × 10⁻¹⁰ L³ d⁻². The only physical content of G m_T is the ratio of M's month to the year, which the fit determines. The six fitted constants are M's state at t = 0: **ρ** in units of A and d**ρ**/dt in units of 0.23 A/d (ecliptic):

ρ₀ = (0.20764290 ± 6.0e-06, −1.03120270 ± 3.1e-06, −0.02627064 ± 1.3e-06),
ρ̇₀ = (0.92410317 ± 2.8e-06, 0.19957394 ± 5.9e-06, 0.07958437 ± 3.7e-06).

Derived (not extra) constants: sidereal month 27.28966 d; anomalistic month 27.5217 d (perigee advance period 8.860 yr); node regression period 18.614 yr; mean inclination 5.025° (range 4.86°–5.17°); mean eccentricity ≈ 0.055.

### Numerical solution
- Law 2 is integrated with classical 4th-order Runge–Kutta, step 0.25 d, from t = 0. Halving the step changes the year-120 directions by < 0.25′ for P1 and < 0.001′ for the others.
- Law 3 is integrated with DOP853 (relative tolerance 10⁻¹¹).
- The constants were fitted by least squares to all sightings, minimising the angular residual.

**Fit quality (rms per angle, arcmin; the stated noise is 1′):**

| | S | P1 | P2 | P3 | P4 | P5 | M |
|---|---|---|---|---|---|---|---|
| Law 2/3 | 1.00 | 1.01 | 1.00 | 1.00 | 1.01 | 1.00 | 1.00 |
| single Kepler ellipse | 1.00 | 1.02 | 1.76 | 8.1 | 1.96 | 3.1 | – |

Every decade, 0 through 9, has residuals between 0.97′ and 1.04′ for every body (`diag.py`).

---

## 3. Predictions

`predict.py` reads `predictions.csv`, integrates Laws 2 and 3 from t = 0 with the constants above, converts to the observation frame with Law 1, and writes lon and lat (degrees) for every row. Run it with `python predict.py` (about 40 s). It uses no other inputs.

Expected accuracy for years 101–120 (direction errors, not including the 1′ measurement noise of any future check):
- S, P2–P5: < 0.1′ from constant uncertainty and integration error.
- P1: about 0.3–1′, mostly the RK4 step error on the fast inner orbit. Two independent fits with different steps differ by up to 2.8′; the adopted finer step reduces this to about 0.25′.
- M: about 0.1–0.5′ from extrapolation of the six fitted constants. Neglected planetary and higher-order effects could add more (see section 4).

---

## 4. Self-assessment

| law | confidence | what would prove it wrong |
|---|---|---|
| Law 1 (fixed ecliptic tilted 21.7015°, geometric directions) | 98% | S drifting off the great circle by more than a few arcmin, or the circle's pole moving against the stars (precession) |
| Law 2 (S-centred Newtonian gravity, inverse-square, P's with the masses above) | 95% for the form; 85% that every mass is right within the stated errors (P1 and P2 masses are only 6σ and 4σ detections) | Years 101–120 positions of any P differing by more than about 3′ (P1) or more than 0.5′ (others, averaged over many nights); a systematic residual growing with time; a perihelion advance of P1 not explained by the masses (e.g. a relativistic-type term) |
| Law 3 (M = inverse-square orbit about Tellus, perturbed only by S's tide) | 90% | M drifting from prediction by more than about 1′ over 20 years, e.g. from an unmodelled Tellus flattening, planet perturbations on M, or a secular acceleration (tidal) of M; node and apsidal periods other than 18.61 yr and 8.86 yr |
| Masses of Tellus + M negligible (< 10⁻⁶ M_S) | 70% | A wobble of S's direction with M's 27.3-day period, or Tellus's mass showing in the P orbits with more data |

The weakest points are:
1. The tiny masses of P1 and P2.
2. P1's long-term phase, which is the most sensitive to integration accuracy.
3. Possible small effects on M that 100 years of 1′ data cannot detect. M's mean longitude shows a residual bias of only −0.05′.

The decisive test is the 20-year prediction itself. If any body departs systematically from `predictions.csv` by many times the noise, the laws above are incomplete.
