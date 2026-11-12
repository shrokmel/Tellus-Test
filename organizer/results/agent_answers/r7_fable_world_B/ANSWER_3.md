# Stage 3: the new body N

Data: `observations_3.csv` (207 sightings of N, days 36233–36525), plus the earlier files for the constants of the
other bodies. Scripts: `nfit.py` (fit), `nmulti.py` (search for other solutions), `nmc.py`, `nmc2.py` (uncertainty),
`predict_3.py` (prediction).

## 1. Laws used for N

The stage 2 laws, unchanged. N's sightings confirm them and do not call for a revision.

**Law 1 (gravitation with a body-specific factor).** For every body b,

    d²r_b/dt² = (1 + η_b) · Σ_{j≠b} G m_j (r_j − r_b) / |r_j − r_b|³

summed over the eight massive bodies S, P1, Tellus, M, P2, P3, P4, P5; η = 0 for those eight, and N has negligible mass.

**Law 2.** Direction of N seen from Tellus: g = r_N − r_T, lon = atan2(g_y, g_x), lat = arcsin(g_z/|g|).

**Law 3.** M–Tellus relative vector ρ: extra term γ dρ/dt, γ = (4.17 ± 0.58)×10⁻¹³ day⁻¹ (irrelevant for N).

**Constants of the massive bodies**: exactly those of stage 2 (G m_S = 2.959×10⁻⁴ L³/day² defines the length L; mass
ratios T+M (3.0428 ± 0.0030)×10⁻⁶, P1 (3.3 ± 6.0)×10⁻⁷, P2 (6.07 ± 0.10)×10⁻⁷, P3 (1.50113 ± 0.00050)×10⁻⁴,
P4 (1.000021 ± 0.000050)×10⁻³, P5 (2.4961 ± 0.0009)×10⁻⁴; states at day 18262). Full values are in the appendix.

**Constants of N** (new; 1σ from the fit, scaled by the fit's reduced χ²):

    η_N = +0.3220 ± 0.0021

Position (L) and velocity (L/day) relative to S at day 36358:

| | x | y | z | vx | vy | vz |
|---|---|---|---|---|---|---|
| value | 0.796109 | 2.189619 | 1.113723 | −0.01535753 | 0.00186426 | 0.00500808 |
| ± | 0.00011 | 0.00026 | 0.00021 | 9.4e−7 | 2.1e−6 | 1.2e−6 |

(The errors are strongly correlated with η: correlation 0.68–0.82 for the position components.)

**Why I keep the laws.**

| model fitted to N's 207 sightings | rms residual (′) | χ² (414 angles) |
|---|---|---|
| stage 1 law (η = 0), 6 constants | 7.76 | 24,912 |
| stage 2 law, η free, 7 constants | 1.033 | 442 |
| stage 2 law + light delay | 1.032 | 441 (delay 0.018 ± 0.021 day/L: not needed) |

- The stage 1 law fails by a factor 56 in χ² and would describe N as an unbound body (e = 1.26) leaving for good.
  With its own η, N fits to the noise (χ² = 442 for 407 degrees of freedom) and is bound.
- The solution is unique: 36 fits started from distances 0.3–14 L and η from −0.5 to +1.0 all end at the same place.
- η_N = +0.322 is slightly outside the range of the twelve earlier bodies (−0.298 to +0.297) but of the same kind.
- From N alone I cannot tell whether η multiplies all gravity or only that of S (the two fit equally: 1.0333 vs
  1.0332). I use "all gravity" because stage 2 established it to 0.1% on twelve bodies. For N this choice matters
  after year 1.2 (see section 4).

## 2. Predictions

`predict_3.py` integrates the massive bodies from day 18262, inserts N at day 36358 with the constants above, and
integrates Law 1 forward (fixed 0.5-day step, 8th-order Runge–Kutta; halving the step changes N's position in year 3
by 5×10⁻¹² L). It fills `predictions_3.csv` (1461 days, 36530–43830). Run: `python predict_3.py` (about 40 s).
Its ephemeris reproduces the sightings after day 36358 to 1.06′/1.03′ rms and agrees with the fitting code to 10⁻⁵ °.

## 3. N's path

**Now.** N was found just before perihelion: it passed 2.57 L from S around day 36358, moving on an ellipse with
a = 10.15 L, e = 0.747, tilted about 33° to the lat = 0 plane, period 2π[a³/((1+η)G m_S)]^{1/2} = 10,280 d (28 yr),
aphelion 17.7 L. At the end of the data it is 3.0 L from S and 3.3 L from Tellus, receding.

**Year 1.2: a very close passage of P4.** On day 36966 ± 1 (1.21 yr after the last sighting) N passes P4 at
0.033 L (68% range 0.028–0.038 L; all 400 trial orbits between 0.026 and 0.059 L), at a relative speed of
0.0094 L/day. That is 13 times closer than the edge of P4's sphere of dominance (0.43 L). P4 bends the path sharply
and takes away orbital energy.

**After the passage** (nominal orbit; ranges from 400 trial orbits drawn from the fit's uncertainty):

| | before (day 36890) | after (day 37620) | range after |
|---|---|---|---|
| a (L) | 10.35 | 5.84 | 5.45–7.27 |
| e | 0.751 | 0.754 | 0.57–0.77 |
| perihelion (L) | 2.58 | 1.44 | 1.39–2.78 |
| aphelion (L) | 18.1 | 10.2 | — |
| period (d) | 10,570 | 4,480 (12.3 yr) | 4,040–6,220 |

So N does not go out to 18 L and stay away for 28 years. It is captured into a short orbit: it reaches only about
10.2 L (year 5.7, day ≈ 38605), falls back, passes perihelion at about 1.4 L from S near year 12 (within about 2.2 L
of Tellus, crossing the sky fast: 316° → 81° in longitude between years 11 and 12), goes out again to about 10.2 L
near year 18, and at year 20 is 9.4 L from S, coming in again. No trial orbit is ejected, and none hits P4's centre.

Nominal track (distance from S, from Tellus, direction):

| years after day 36525 | day | r from S (L) | from Tellus (L) | lon | lat |
|---|---|---|---|---|---|
| 0.01 | 36530 | 3.01 | 3.32 | 119.0° | +21.1° |
| 0.5 | 36705 | 4.12 | 4.63 | 182.6° | +21.0° |
| 1 | 36890 | 5.40 | 4.89 | 179.5° | +5.8° |
| 2 | 37255 | 7.45 | 6.80 | 196.5° | −5.6° |
| 3 | 37620 | 8.81 | 8.13 | 205.1° | −12.6° |
| 5 | 38350 | 10.14 | 9.45 | 217.4° | −22.7° |
| 7 | 39080 | 9.92 | 9.27 | 229.2° | −31.6° |
| 9 | 39810 | 8.06 | 7.53 | 247.0° | −42.1° |
| 10 | 40175 | 6.35 | 5.95 | 265.5° | −49.0° |
| 11 | 40540 | 3.83 | 3.80 | 316.1° | −52.5° |
| 12 | 40905 | 1.67 | 2.23 | 81.2° | +21.7° |
| 13 | 41270 | 4.79 | 4.27 | 175.0° | +10.1° |
| 15 | 42000 | 8.51 | 7.79 | 203.5° | −11.0° |
| 18 | 43095 | 10.22 | 9.51 | 222.0° | −26.0° |
| 20 | 43830 | 9.38 | 8.77 | 234.5° | −35.2° |

Whether N is bright enough to be seen again I cannot say (no brightness data); by distance, the next good chance is
the perihelion near year 12.

## 4. Self-assessment

Expected error of the predicted direction (angle on the sky between prediction and truth), from 400 trial orbits
consistent with the sightings:

| time after last sighting | typical (median) | 68% below | 95% below |
|---|---|---|---|
| 1 year | 0.05° | 0.07° | 0.15° |
| 2 years | 0.8° | 1.2° | 2.4° |
| 5 years | 2.6° | 4.0° | 10° |
| 10 years | 15° | 23° | 53° |
| 12 years (perihelion, close to Tellus) | 85° | 106° | 159° (i.e. unpredictable) |
| 20 years | 12° | 18° | 32° |

My answer to the question: about **0.05–0.1° after 1 year, about 3° after 5 years (10° would not surprise me), about
15° after 20 years (30° would not surprise me)**. Before the P4 passage the prediction is as good as the fit. After
it, everything depends on the miss distance, known only to ±15%, and the errors are 20–50 times larger. They peak
around the perihelion of year 12, when a timing error of a few weeks swings the direction across the sky.

Additional risks not in that table:
- If N's factor acted on S's pull only and P4 pulled N with the plain stage 1 strength, the path after the passage
  would differ from mine by 1.3° at year 2, 2.7° at year 5 and 14° at year 20. Stage 2 says this is unlikely
  (planet factor equal to S factor to 0.1%), but N alone cannot check it.
- P4 is treated as a point. At 0.033 L any extent, flattening or companions of P4 would matter; I have no data.
- P4's own position and mass are known far better than N's miss distance and add little.

Confidence: N is bound to S and has η ≈ +0.32: 97%. It passes within 0.06 L of P4 near day 36966: 90%. It is
thrown into a shorter orbit with period 11–17 yr and returns to perihelion within the 20 years: 80%.

What would prove me wrong:
- Sightings in the first year more than 0.2° from the predictions: the fit or the law is wrong before any encounter.
- N continuing outward past 11 L after year 6, or not returning near year 12: the P4 passage was not as computed.
  A stage 1 (η = 0) orbit would put N 10° from my prediction after 1 year and 24° after 5.
- A later, longer arc of N giving an η different from +0.322 by more than about 0.006.

## Appendix: all constants at full working precision (identical to `predict_3.py`)

G m_S = 2.959e-4 L^3/day^2, T0 = day 18262. `GM_REL` = G m / G m_S; `GAMMA` in 1/day; `STATE` = x, y, z (L), vx, vy, vz (L/day) at T0 relative to S (M relative to Tellus); `FAINT` = (eta, state of N relative to S at day TN).

```
GM_REL = {   # GM_body / GM_S
    'P1': 3.2741614495e-07,
    'T':  3.0083240800e-06,
    'M':  3.4441577080e-08,
    'P2': 6.0736690970e-07,
    'P3': 1.5011246810e-04,
    'P4': 1.0000211634e-03,
    'P5': 2.4961482571e-04,
}
GAMMA = 4.171622e-13   # 1/day : secular along-track acceleration coefficient of M
STATE = {   # at day T0: x,y,z [L], vx,vy,vz [L/day]; relative to S (M: relative to Tellus)
    'P1': [9.1064625793917428e-02, 4.0711307714642903e-01, -1.6201496553362094e-01, -2.5512649259988664e-02, 3.4643403204493821e-03, -1.3791398583217075e-03],
    'T':  [-4.0776510470894889e-01, -8.6627310270167379e-01, 3.4475232702735803e-01, 1.5420533311223866e-02, -6.4330704887881835e-03, 2.5604599109755742e-03],
    'M':  [2.4967154385874607e-03, -1.7254207205097418e-04, -1.3165789847240892e-04, 4.8481672576106441e-05, 5.6945717821121943e-04, -1.9570359391707365e-04],
    'P2': [1.2113147734079719e+00, -1.0129050205619601e+00, 4.0335357397262367e-01, 9.4886205124626154e-03, 9.2461983448505064e-03, -3.6797702263660116e-03],
    'P3': [6.5768643691123563e-02, -3.1257926417280570e+00, 1.2873631044800671e+00, 9.2255039138841397e-03, 3.5653317767069101e-04, -3.7583546803934739e-04],
    'P4': [1.3175273244394403e+00, 5.7847146868591590e+00, -2.4224862181861733e+00, -6.5860226060181007e-03, 1.0603492951939233e-03, -4.0658173109741046e-04],
    'P5': [-1.0981162632906226e+01, -3.9631404639702456e-01, 5.4566383420789377e-01, 2.2695437333347088e-04, -4.9360099621809659e-03, 1.8977742448029952e-03],
}
TN = 36358      # epoch (day) of the state of N
FAINT = {   # name: (ETA, [x,y,z (L), vx,vy,vz (L/day)] relative to S at day TN)
    'N': (3.2201120915e-01, [7.9610948908865509e-01, 2.1896192239166838e+00, 1.1137234590298801e+00, -1.5357525246208101e-02, 1.8642546631355878e-03, 5.0080826224349973e-03]),
}
```
