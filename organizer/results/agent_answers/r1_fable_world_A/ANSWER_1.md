# The sky of Tellus: motions, laws, predictions

Summary: every one of the 156 178 sightings (7 bodies, 100 years) is reproduced to the stated
measurement noise (rms 0.9995 arcmin per axis, expected 1.000) by **one law**: eight point masses
(S, Tellus, M, P1–P5) attracting each other with Newton's inverse-square gravitation, with Tellus and the five
wanderers orbiting S, and M orbiting Tellus. Nothing else is needed, and no structure is left in the residuals.

Files: `predict.py` (prediction script), `constants.json` (the constants of section 2 at full precision),
`predictions.csv` (completed). Working scripts: `nb.py`, `fitlib.py`, `fitA.py`, `fitB.py`, `fitC.py`,
`export.py`, `describe.py`, `holdout.py`.

---

## 1. Description of the motions

**Frame.** `lon`/`lat` are tied to Tellus's equator. S moves on a great circle tilted 21.700° ± 0.001° to that
equator (ascending node at lon 0.00°): this is the plane of Tellus's orbit ("ecliptic"; pole at
lon 90.00°, lat 68.30°). All bodies stay within about 5° of this circle. Element values below are referred to
that plane (node longitudes measured from the lon = 0 direction).

**S** circles the sky eastward once per 365.2571 d (sidereal year), not uniformly: its speed varies by ±2 %
(orbit eccentricity 0.01998, fastest at ecliptic longitude ≈ 93.9°). S also shows a tiny 29.49-day wobble
(≈ 0.1 arcmin), because Tellus itself circles the Tellus–M centre of mass.

**P1** never strays more than 29° from S (elongation −27.9° … +29.0°): it is an inner planet, swinging from
side to side of S every 157.93 d. **P2–P5** go all the way round the sky eastward, but each makes a backward
(retrograde) loop around opposition to S, once per synodic period (665.5, 438.4, 390.6, 374.9 d). All of this is
simply Tellus's own motion reflected in the directions: seen from S, all planets move on nearly fixed, nearly
circular ellipses.

Mean heliocentric orbits (length unit L defined in section 2; Tellus's a = 0.999994 L):

| body | a (L) | e | incl. to ecliptic | sidereal period (d) | synodic (d) | mass / S |
|---|---|---|---|---|---|---|
| P1 | 0.44999 | 0.0600 | 0.00° | 110.2578 | 157.93 | not detectable (<~3e-6) |
| Tellus (+M) | 0.99999 | 0.01998 | 0 (definition) | 365.2571 | – | 3.0400e-6 |
| P2 | 1.70004 | 0.0501 | 0.01° | 809.656 | 665.47 | (4 ± 1.2)e-7, marginal |
| P3 | 3.2991 | 0.0422 | 1.49° | 2189.21 | 438.40 | (1.503 ± 0.005)e-4 |
| P4 | 6.2042 | 0.0493 | 1.00° | 5639.8 | 390.55 | (1.0001 ± 0.0005)e-3 |
| P5 | 11.483 | 0.0414 | 2.00° | 14149.8 | 374.94 | (2.500 ± 0.007)e-4 |

Kepler's third law a³/T² = const holds to 1e-5 for P1 and to ~1e-3 for the outer planets (the deviations are
the planetary perturbations, not noise).

**M** circles Tellus (not S) in 27.28966 d (sidereal), 29.4932 d relative to S (phases), on an ellipse with
mean e = 0.056 inclined 5.03° (4.85°–5.17°) to the ecliptic, at mean distance 0.00256 L (1/390 of the
S distance). Its orbit is strongly disturbed by S: osculating e swings between 0.026 and 0.077.

**How things change over the 100 years**

* M: the line of nodes regresses (westward) 1934°/century, a full turn in 6799 d = 18.61 yr; the perigee
  advances (eastward) 4063°/century, a full turn in 3236 d = 8.86 yr. Hence draconic month 27.1806 d,
  anomalistic month 27.5217 d. Both are caused by S's pull.
* Pure fixed ellipses fail for the outer planets: single-ellipse fits leave rms residuals of 1.8′ (P2),
  8.1′ (P3), 2.0′ (P4), 3.1′ (P5) against 1.0′ noise. The planets pull on each other, chiefly P4
  (mass 1/1000 of S). P4 and P5 are close to a 5:2 period commensurability (14150/5640 = 2.509), giving a
  slow large oscillation: P5's osculating period drifts from ≈14270 d (first decade) to ≈14145 d (last decade),
  its a between 11.43 and 11.53 L; P3's a varies 3.297–3.304 L and its e grows 0.0414 → 0.0439.
* Slow turning of the ellipses (perihelion advance per century): Tellus +0.24°, P1 +0.11°, P2 +0.95°,
  P3 ≈ +1.9°, P4 ≈ +0.4°, P5 ≈ +0.7° (the last three are mixed with long-period oscillations).
* The orbit planes tilt slightly: the plane of Tellus's orbit turns by 0.0053° (0.32′) in a century, which is
  why the S great circle is not perfectly fixed; P3's inclination drops 1.499° → 1.485°, P4's rises
  1.000° → 1.008°.
* Mean sizes and periods of Tellus, P1 and M orbits show no secular change (|Δa/a| < 1e-6 per century): no
  sign of drag, tides or any non-gravitational drift.

---

## 2. The laws, as equations

**Law (Newtonian gravitation of 8 point masses).** For each body i ∈ {S, T(ellus), M, P1, P2, P3, P4, P5},
with position **r**_i(t) in an inertial frame (axes = the lon/lat frame of the data, t in days):

    d²r_i/dt² = Σ_{j≠i}  G m_j (r_j − r_i) / |r_j − r_i|³                    (1)

**Observation equation.** The direction of body X on the sky of Tellus at day t is that of **r**_X − **r**_T:

    lon = atan2(d_y, d_x),  lat = asin(d_z/|d|),   d = r_X(t) − r_T(t)       (2)

(no light-travel-time, aberration or parallax correction is needed at the 1′ level; a body is hidden when
within 15° of S).

Consequences contained in (1), each checked on the data: orbits are conics with the attracting body at a focus
and equal areas in equal times; n²a³ = G(m_S + m) (third law); M's nodal regression and apsidal advance;
planet–planet perturbations; Tellus's monthly wobble about the Tellus–M centre of mass.

**Constants.** Only directions are observed, so absolute length and mass are unobservable. I fix the scale by
*defining* the length unit L through

    G m_S = 2.9591e-4 L³/day²  (exact, definition)

Masses relative to S (± 1σ formal):

| quantity | value | ± |
|---|---|---|
| (m_T + m_M)/m_S | 3.039985397e-6 | ~1e-10 formal, tied to M's distance scale; see note |
| μ = m_M/(m_T + m_M) | 0.01180 | 0.0006 (so m_M/m_T = 1/84 ± 4) |
| m_P1/m_S | 0 (fixed) | not detectable; |m| < ~3e-6 |
| m_P2/m_S | 3.97e-7 | 1.2e-7 (marginal) |
| m_P3/m_S | 1.50272e-4 | 4.9e-7 |
| m_P4/m_S | 1.000128e-3 | 5.2e-7 |
| m_P5/m_S | 2.49966e-4 | 6.8e-7 |

Exact values used: m_T = 3.0040984945048e-06, m_M = 3.5886902937516e-08, m_P2 = 3.9734929710287e-07,
m_P3 = 1.5027170880253e-04, m_P4 = 1.0001283624708e-03, m_P5 = 2.4996589240570e-04 (× m_S).

Initial conditions at day 0, heliocentric (S at rest at the origin; the script then shifts to the
centre-of-mass frame), in the data's lon/lat axes, positions in L, velocities in L/day:

| body | x, y, z | vx, vy, vz |
|---|---|---|
| T | −0.394260876772, −0.871803677116, +0.346937334799 | +0.0155124441755, −0.0062150757658, +0.0024728472798 |
| M | −0.393726861089, −0.874292699163, +0.347855124850 | +0.0160590540239, −0.0060879783176, +0.0024729396470 |
| P1 | −0.312574970513, −0.280258502792, +0.111530944174 | +0.0174773583867, −0.0185686179707, +0.0073894647386 |
| P2 | −0.714778209480, +1.517465419472, −0.603879547548 | −0.0114537867214, −0.0048018613780, +0.0019109128914 |
| P3 | −3.056483955071, +1.287910267015, −0.455467332649 | −0.0041175024164, −0.0077065746059, +0.0032820130637 |
| P4 | +6.458854301488, −0.038258785775, −0.019205541576 | +0.0002107786045, +0.0061138588545, −0.0025533072478 |
| P5 | +3.323202144279, +10.245808774887, −4.077414797193 | −0.0049017660339, +0.0012413566517, −0.0003038742665 |

`constants.json` holds exactly these numbers at full double precision and is what `predict.py` reads.

Uncertainties of the initial conditions (1σ, formal, from the least-squares covariance): the fitted quantities
were the Tellus–M centre of mass B, the M-minus-Tellus vector "rel", and the planets
(T = B − μ·rel, M = B + (1−μ)·rel):

| | position σ (L) per component | velocity σ (L/day) per component |
|---|---|---|
| B (Tellus–M centre) | 1.5e-6 – 3.6e-6 | 2.6e-8 – 4.3e-8 |
| rel (M − T) | 0.7e-8 – 1.6e-8 | 1.7e-9 – 3.2e-9 |
| P1 | 3.5e-6 – 6.8e-6 | 2.2e-7 – 3.4e-7 |
| P2 | 3.5e-6 – 6.5e-6 | 2.4e-8 – 3.8e-8 |
| P3 | 0.7e-5 – 1.1e-5 | 1.8e-8 – 3.0e-8 |
| P4 | 1.1e-5 – 2.0e-5 | 1.4e-8 – 2.2e-8 |
| P5 | 2.4e-5 – 4.0e-5 | 0.8e-8 – 1.8e-8 |

(full list per component in `params.txt`, covariance in `pcov.npy`). Note: these are formal errors from a
finite-difference Jacobian; parameters are strongly correlated (e.g. M's distance scale with m_T + m_M, which
are tied by M's period as m ∝ a³), so individual σ's understate what any one number alone is known to.

**Numerical solution used.** Equation (1) is integrated with a fixed-step 6th-order symplectic scheme
(Yoshida composition of 7 leapfrog steps), step 0.25 day, from day 0. The constants above were fitted with this
same scheme and step, so `predict.py` reproduces the fit exactly.

**Fit quality** (rms residual per angle, arcmin; noise is 1.00): S 0.995, M 0.996, P1 1.010, P2 1.001,
P3 0.994, P4 1.003, P5 1.001; all 156 178 sightings 0.9995.

---

## 3. Predictions

`predict.py` (numpy only, one process, about a minute) reads `constants.json`, integrates equation (1) from
day 0 to day 43 815 and writes equation (2) for every row of `predictions.csv` (1701 rows, all filled).
`python predict.py --check` additionally prints the residuals against `observations.csv`.

Expected accuracy. Out-of-sample test actually run (`holdout.py`): the same model fitted to days 0–15 000 only
(41 yr) predicted the unseen days 29 000–36 525 (38–59 yr ahead) with total angular rms 1.40–1.43′ for S and
the planets and 1.51′ for M, against 1.41′ expected from measurement noise alone. So for the 100-yr fit carried
20 yr ahead I expect model errors well below the 1′ noise: ≲0.1′ for S and P1–P5, ≲0.3′ for M (M is the
most sensitive: an error of 1e-8 in its mean motion is 0.25′ after 20 yr).

---

## 4. Self-assessment

| law / claim | confidence | what would prove it wrong |
|---|---|---|
| Tellus and P1–P5 orbit S; M orbits Tellus (heliocentric arrangement) | 99 % | planets' retrograde loops not locked to S's position; P1 reaching large elongations |
| Inverse-square mutual gravitation of 8 point masses, eq. (1), explains everything to 1′ | 97 % | future residuals growing systematically above ~1′ (e.g. an M longitude drift ∝ t², a perihelion turning faster than computed, or a planet drifting as if pulled by an unseen body) |
| Kepler's laws as the first approximation (ellipses, areas, a³ ∝ T²) | 99 % | – (already seen to fail at the few-arcmin level, exactly as eq. (1) predicts) |
| P4, P3, P5 masses (1.000e-3, 1.50e-4, 2.50e-4 of S) | 90 % within 1 % | perturbation amplitudes (e.g. P3's 8′ departure from an ellipse) inconsistent in later decades |
| m_T + m_M = 3.04e-6 m_S and M's distance 0.00256 L | 75 % | they rest on small effects (parallactic inequality of M, S's 0.1′ monthly wobble); a direct distance measure of M (parallax) differing by > few % |
| m_M/m_T ≈ 1/84 | 70 % within 15 % | amplitude of S's 29.49-d wobble different from 0.1′ |
| P1 and P2 masses ≈ 0 | 60 % that P2's is really ~4e-7; P1's is unknown | – these do not matter for the predictions at 1′ |
| No light-time, aberration, relativity, tides or drag visible | 85 % | residuals correlated with distance or velocity of the body; secular change of M's month |
| Predictions for years 101–120 correct to about the 1′ noise | 90 % (planets, S), 80 % (M) | see above |

Caveats, stated plainly:

* The least-squares search converges slowly in its last digits (χ² was still falling by ~10 out of 312 000 per
  iteration when I stopped). The parameters had stopped moving at the level that matters for 1′ predictions,
  and the out-of-sample test above was made with an even less converged fit.
* A free fit preferred a slightly negative P1 mass (−3.2e-6, formally many σ). That is unphysical; refitting
  with m_P1 = 0 raised χ² by only 34 out of 312 000, spread over several bodies with both signs, which is
  inside the convergence slop. I therefore fixed m_P1 = 0 and regard P1's mass as undetermined, and the formal
  σ's of section 2 as optimistic. I could not rule out that this hints at a very small unmodelled effect.
* Absolute distances and masses are not determined by directions alone; only ratios are (hence the defined
  unit L).
* I did not test alternative force laws explicitly (e.g. an exponent different from 2); the evidence for the
  exponent is that one law with fixed constants fits both the heliocentric orbits and the planet–planet and
  S–M perturbations to the noise over 100 years, and extrapolates.
