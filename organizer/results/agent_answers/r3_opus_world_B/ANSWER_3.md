# Stage 3: the new faint body N

Short version: N obeys the stage-2 laws unchanged, with its own coupling **κ_N = 1.3220 ± 0.0022**. It has just passed perihelion on a steep orbit (a ≈ 10.1 AU, e ≈ 0.75, i ≈ 54°). On **day 36967 (early in year 102) it passes 0.032 AU from P4**. That encounter throws it into a smaller, steeper orbit (a ≈ 6.4 AU, i ≈ 62°), which brings it back to perihelion near Tellus's orbit in year ~112. The encounter limits how well N can be predicted after year 102.

Units and frames are as in stages 1 and 2: AU (k = 0.01720209895), days, ecliptic frame, data frame obtained by rotating by ε = 21.7015°.

## 1. Laws and constants used for N

**The laws are the same as in stage 2. N's sightings did not change them.**

- **Law 1 (massive bodies, unchanged).** d²**r**_i/dt² = Σ_{j≠i} GM_j(**r**_j − **r**_i)/|**r**_j − **r**_i|³ for S, P1, P2, T, P3, P4, P5, M. All constants are those of ANSWER_1.md §2 (copied verbatim into `predict_3.py`), for example μ_P4 = GM_P4/GM_S = 9.9992e-4 ± 5e-7.
- **Law 1′ (faint bodies, unchanged).** N is a massless test particle:
  d²**r**_N/dt² = κ_N Σ_j GM_j(**r**_j − **r**_N)/|**r**_j − **r**_N|³, with j over the 8 massive bodies.
- **Law 2 (unchanged).** The seen direction is unit(**r**_N − **r**_T), rotated to the data frame by ε = 21.7015°.

**Why κ_N is needed.** These are N-body fits of all 207 sightings (days 36233–36525), rms per angle:

| κ_N fixed at | 0.70 | 0.90 | 1.00 | 1.10 | 1.20 | 1.25 | 1.30 | **1.322 (free)** | 1.35 | 1.40 |
|---|---|---|---|---|---|---|---|---|---|---|
| rms (′) | 15.1 | 10.2 | 7.76 | 5.37 | 3.06 | 1.99 | 1.16 | **1.033** | 1.22 | 2.09 |

So the stage-1 law (κ = 1) fails for N by 7.8′. With its own κ, N fits to the 1′ noise exactly as the stage-2 bodies did. κ_N = 1.322 lies just above the stage-2 range (0.70–1.30). Nothing in N's sightings calls for a new law.

**Constants for N.** Heliocentric ecliptic state at epoch day 36379.0 (formal 1σ from the fit; residual rms 1.033′):

| | value | ± 1σ |
|---|---|---|
| x (AU) | 0.470141446 | 1.1e-4 |
| y (AU) | 1.611943078 | 1.8e-4 |
| z (AU) | 1.947179881 | 3.1e-4 |
| vx (AU/d) | −0.0156624231 | 9.8e-7 |
| vy (AU/d) | −0.0008985069 | 1.4e-6 |
| vz (AU/d) | 0.0044295782 | 1.5e-6 |
| κ_N | 1.3220 | 0.0022 |

Full-precision values are in `predict_3.py`. The corresponding osculating elements at the epoch:
- a = 10.147 AU, e = 0.7466, i = 53.68°, Ω = 15.26°, ω = 70.63°
- q = 2.571 AU (perihelion is at about the epoch, day ≈ 36379), Q = 17.72 AU
- period = 2π(a³/κ_N k²)^{1/2} = 10,270 d

**Method.**
1. A two-body orbit with κ free gives the starting point.
2. A Levenberg–Marquardt N-body fit follows: test particle in the stage-1 field, DOP853 at rtol = 1e-11, all 7 parameters free.
3. Run on the observed days, `predict_3.py` reproduces the 207 sightings at 1.46′ total angle (the noise alone gives 1.41′).

## 2. Predictions

`predict_3.py` integrates the stage-1 massive bodies (Law 1), then N from its epoch (Law 1′). It applies Law 2 and fills all 1,461 rows of `predictions_3.csv` (days 36530–43830).

## 3. N's path over the next 20 years

1. **Year 100 (day 36525).** N is 3.0 AU from S and 3.3 AU from Tellus, moving outward just after perihelion (q = 2.57 AU on day ≈ 36379). It is on a steeply inclined orbit (54°), seen high in the north of the data frame (lon 117°, lat +22°) and moving east quickly.
2. **Day 36966.75 ± 0.5 (year 101.2).** N passes **0.032 AU from P4** (5–95% range over the uncertainties: 0.026–0.045 AU). P4's pull, scaled by κ_N, transforms the orbit:
   - a drops from 10.15 to ≈ 6.44 AU (5–95%: 6.14–7.24)
   - the period drops from ≈ 28 yr to ≈ 14 yr
   - aphelion drops from 17.7 to ≈ 11 AU
   - the inclination rises from 53.7° to ≈ 62.6° (57.8–66.5°)
   - e stays ≈ 0.74 (0.66–0.75)
   - N stays bound in every Monte Carlo sample.
   Seen from Tellus, it swings within days to the other side of the sky. By year 102 it is at lon ≈ 196°, lat ≈ −6°.
3. **Years 102–108.** N recedes to its greatest distance, r ≈ 10.2 AU around day 38600 (year 105.7), and drifts slowly through lon ≈ 210–236°, moving southward (lat −18° at year 104, −27° at 106, −36° at 108).
4. **Years 109–112.** N falls back toward S, speeding up and diving further south (lat ≈ −49° at year 110). Its perihelion is ≈ 1.7 AU, near day 40850 (year ≈ 111.8). On the way it passes 1.12 AU from P1 (day 40846), 1.53 AU from P3 (day 40773) and 2.0 AU from Tellus (day 40714). Around perihelion it sweeps rapidly across the sky and is back in the north (lon ≈ 85°, lat +22°) by year 112.
5. **Years 112–120.** N goes outward again on the new orbit (r ≈ 7.0 AU at year 114, 9.5 at 116, 10.2 at 118, 9.4 at year 120). Its sky track repeats the year 102–108 pattern: lon ≈ 194° → 235°, lat ≈ −3° → −35°. No further close planetary encounters occur before year 120: the closest is P3 at 1.5 AU.

## 4. Expected errors and what would prove me wrong

I propagated 150 samples drawn from the fit covariance (state and κ_N) through the same equations. The resulting spread of N's predicted direction:

| time after last sighting | median error | 90% bound |
|---|---|---|
| 1 yr (before the P4 encounter) | 0.05° | 0.15° |
| 2 yr (just after the encounter) | 0.9° | 2.1° |
| 5 yr | 2.9° | 8° |
| 10 yr (approaching the second perihelion) | 17° | 43° |
| 20 yr | 13° | 27° |

Rounded expectations:
- **after 1 yr: ≈ 0.1°** (up to 0.2°)
- **after 5 yr: ≈ 3°** (up to ~8°)
- **after 20 yr: ≈ 15°** (up to ~30°)
- Around the second perihelion (years 110–113), N moves fast and close to Tellus, so errors there can reach tens of degrees.

The growth is almost entirely due to the P4 encounter. Errors of a few 1e-4 AU in N's position, and the 0.2% uncertainty in κ_N, become large changes in the deflection. Before the encounter the prediction is good to a fraction of a degree.

The model limits from stage 2 (unexplained arcminute-level drifts of some faint bodies over decades) are negligible compared with this.

**Confidence.**
- N follows Law 1′ with a constant κ_N ≈ 1.32: 85%.
- The P4 encounter happens on day 36967 ± 1 at 0.026–0.045 AU: 85%.
- Errors within the ranges in the table: 75%.

**What would prove me wrong.**
- In year 101 (before day 36960), N strays more than ~0.3° from the predicted track. That would mean the orbit or κ_N is wrong.
- After the encounter, N fails to swing across the sky around day 36967, or keeps moving on its old 28-year orbit, with aphelion ~17.7 AU and no perihelion return near year 112. That would mean the encounter did not happen as computed, i.e. P4's pull on N is not κ_N·GM_P4.
- N is not back to perihelion at ≈ 1.7 AU around years 111–112.
- N shows a sudden change of course with no massive body nearby, which would point to a force outside Laws 1/1′.
