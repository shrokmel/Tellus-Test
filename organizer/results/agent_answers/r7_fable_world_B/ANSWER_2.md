# Stage 2: the faint bodies C1–C12

Data used: `observations.csv` (156,178 sightings of S, P1–P5, M) and `observations_2.csv` (39,886 sightings of C1–C12).
Scripts: `nbody2.py`, `cfit*.py`, `cjoint*.py`, `newmass2.py`, `stageD.py` (fitting), `predict_2.py` (prediction),
`check2.py` (verification).

## 1. Test of the stage 1 laws

Stage 1 laws: (1) Newtonian gravitation between all bodies; (2) the observed direction is the geometric direction from
Tellus; (3) a small secular term γ on M.

**What the faint bodies are.** Each goes round S on an elongated, tilted ellipse (e = 0.50–0.95, perihelion 0.7–2.9 L,
aphelion 3–26 L, return period 1200–15,700 d, tilt to Tellus's orbit 0.1°–17°) and is seen only near perihelion.

**Where stage 1 holds.**
- Law 2 holds for all 19 bodies. A light-delay term fitted separately for each faint body gives values of both signs
  (−0.020 to +0.003 day per L, mean −0.003) and does not improve the fits; the direction is geometric.
- Law 1 holds in *form* for every body: inverse-square attraction towards S and towards every planet, summed.
  The planets' pulls on the faint bodies are large (C5 has its orbit changed from q = 2.16 L, tilt 1.3° to q = 2.68 L,
  tilt 18° by a passage near P4) and are reproduced.
- Law 3 is not tested by the faint bodies. The bright bodies fit as before.

**Where stage 1 fails.** Treated as ordinary small bodies under Law 1 (attraction G m_S/r² towards S etc.), the faint
bodies cannot be fitted even within one apparition. Best-fit residuals (rms, arcmin; noise is 1.0), stage 1 law, one
apparition at a time:

| body | stage 1 law, single apparition | with the factor (1+η) below, single apparition |
|---|---|---|
| C1 | 132–168 | 0.95–1.01 |
| C2 | 11–19 | 0.97–1.04 |
| C3 | 21–37 | 0.94–1.00 |
| C4 | 383 (the 480-day apparition; 1.1–1.4 on the two 110–130-day ones) | 0.90–1.02 |
| C5 | 5–24 | 0.97–1.01 |
| C6 | 1.9–4.2 | 0.98–0.99 |
| C7 | 57–79 | 0.97–1.02 |
| C8 | 104–173 | 0.98–1.05 |
| C9 | 20–37 | 0.98–1.01 |
| C10 | 3.9–4.9 (1.0 on a 190-day one) | 0.97–1.01 |
| C11 | 16–30 | 0.99–1.00 |
| C12 | 6–27 | 0.98–1.00 |

The stage 1 law also predicts the wrong return period (C6: 2903 d from its motion near perihelion, 5403 d observed).

**What is wrong.** Every faint body moves exactly as if all gravity acting on it were multiplied by a constant of its
own, 1 + η_b, with η_b between −0.298 and +0.297 (both signs occur). Tests of what that factor multiplies:
- *Only the pull of S* (e.g. a second S-centred 1/r² force): excluded. C6 then fits two apparitions at 1.66′ and
  three at 3.85′ rms; with the factor applied to all attracting bodies it fits all four at 1.000′.
- *A separate factor for the planets' pull*, fitted for each faint body: the extra planetary factor is −0.0012 to
  +0.0002 for all twelve, i.e. the planets' pull carries the same (1+η_b) to within 0.1%.
- *Dependence on the orbit*: none. C4 and C10 have nearly the same orbit (a ≈ 9.1 L, e ≈ 0.91) and η = −0.274 and
  +0.033. η is the same at every return of a body (e.g. C7: −0.2175, −0.2176, −0.2177 in three separate apparitions).
- For the bright bodies η = 0: P1–P5 obey P²/a³ = const to 1–2 parts in 10⁴ once their own masses are allowed for.

**Second, smaller finding.** With η alone the fits over the whole century left 2–13% excess rms for the bodies with
many returns (C8 1.126, C11 1.091, C3 1.085, C2 1.057, C7 1.039). This is not a failure of the law: the faint bodies
pass close to P2, P3 and Tellus and measure their masses better than stage 1 could. Adjusting four masses, all within
1σ of the stage 1 values, removes the excess (all twelve at 0.99–1.013).

## 2. The laws

Bodies: eight massive bodies S, P1, T (Tellus), M, P2, P3, P4, P5, and twelve faint bodies C1–C12 of negligible mass.
Positions **r** in the fixed star frame of the data, t in days.

**Law 1 (generalised).** Every body b is accelerated by the inverse-square attraction of all massive bodies,
multiplied by a constant η-factor belonging to the accelerated body:

    d²r_b/dt² = (1 + η_b) · Σ_{j≠b} G m_j (r_j − r_b) / |r_j − r_b|³

with η_b = 0 for S, P1–P5, T, M (this is then exactly stage 1's Law 1) and η_b ≠ 0 for the faint bodies.
(Equivalent statement: a faint body's "gravitational charge" is (1+η_b) times its inertia; the equality of the two is
not universal.)

**Law 2 (unchanged).** Direction of b seen from Tellus: g = r_b − r_T, lon = atan2(g_y, g_x), lat = arcsin(g_z/|g|).

**Law 3 (unchanged, empirical).** The M–Tellus relative vector ρ has the extra term d²ρ/dt² = … + γ dρ/dt,
γ = (4.17 ± 0.58)×10⁻¹³ day⁻¹ (M's mean motion falls by 3γ per day).

What forced the change: the twelve faint bodies (section 1). The equations for the bright bodies are the same as in
stage 1; only constants were refined.

**Constants.** Length unit L: G m_S = 2.959×10⁻⁴ L³ day⁻² (definition). Epoch t₀ = day 18262.

Masses, G m / G m_S (stage 1 value in brackets; all changes are within 1σ of stage 1):

| body | value | ± | stage 1 |
|---|---|---|---|
| P1 | 3.3×10⁻⁷ | 6.0×10⁻⁷ (not detected) | (2.7 ± 6.1)×10⁻⁷ |
| T | 3.0083×10⁻⁶ | 0.0030×10⁻⁶ | (2.975 ± 0.035)×10⁻⁶ |
| M | 3.444×10⁻⁸ | 0.060×10⁻⁸ | (3.438 ± 0.061)×10⁻⁸ |
| T + M | 3.0428×10⁻⁶ | 0.0030×10⁻⁶ | (3.010 ± 0.035)×10⁻⁶ |
| P2 | 6.07×10⁻⁷ | 0.10×10⁻⁷ | (6.69 ± 1.2)×10⁻⁷ |
| P3 | 1.50113×10⁻⁴ | 0.00050×10⁻⁴ | (1.5026 ± 0.0050)×10⁻⁴ |
| P4 | 1.000021×10⁻³ | 0.000050×10⁻³ | (0.99986 ± 0.00055)×10⁻³ |
| P5 | 2.4961×10⁻⁴ | 0.0009×10⁻⁴ | (2.5002 ± 0.0072)×10⁻⁴ |

(The ± are 2.5 × the formal errors, or the size of the last iteration step if larger, because the individual faint
bodies' mass estimates scatter more than their formal errors: χ² = 349 for 55 degrees of freedom.)

Positions and velocities of the massive bodies at t₀: the stage 1 values re-adjusted for the new masses (shifts below
their stage 1 uncertainties, which still apply); full values in the appendix. The bright bodies' predicted directions
change by at most 0.014′ relative to stage 1.

Faint bodies: η, and position (L) / velocity (L/day) relative to S at t₀ (1σ formal errors; full precision in the
appendix):

| body | η | ± | x | y | z | vx | vy | vz |
|---|---|---|---|---|---|---|---|---|
| C1 | +0.264886 | 5.7e−5 | −9.1189895 | −21.8249640 | 6.5764645 | 0.001299360 | 0.000761017 | −0.000206835 |
| ± | | | 2.9e−4 | 3.8e−4 | 3.6e−4 | 1.9e−8 | 2.0e−8 | 1.6e−8 |
| C2 | −0.099704 | 6.0e−5 | 11.8528580 | 4.9029344 | −1.9736802 | −0.001070153 | 0.001730198 | −0.000687161 |
| ± | | | 2.5e−4 | 1.9e−4 | 1.5e−4 | 3.0e−8 | 4.1e−8 | 2.5e−8 |
| C3 | +0.230589 | 5.5e−5 | −3.1166271 | 0.0923077 | −0.5606415 | 0.006761183 | −0.009323501 | 0.005720685 |
| ± | | | 4.4e−5 | 3.1e−5 | 2.2e−5 | 1.0e−7 | 1.7e−7 | 1.2e−7 |
| C4 | −0.273608 | 2.0e−5 | −22.5997546 | 3.1879265 | −0.1979902 | −0.002256834 | −0.000617743 | 0.000433200 |
| ± | | | 6.1e−3 | 9.7e−4 | 3.0e−3 | 1.7e−6 | 2.5e−7 | 4.6e−7 |
| C5 | +0.171787 | 7.4e−6 | −0.4996732 | 5.1939148 | −1.9212824 | −0.005970569 | 0.002297842 | −0.000905534 |
| ± | | | 2.6e−5 | 1.7e−5 | 3.6e−5 | 1.8e−8 | 2.0e−8 | 1.3e−8 |
| C6 | −0.143973 | 1.3e−4 | −4.4908159 | 7.0123162 | −1.6755726 | −0.002731657 | −0.003047991 | 0.000013174 |
| ± | | | 5.1e−4 | 2.6e−4 | 1.3e−4 | 8.5e−8 | 2.1e−7 | 1.5e−7 |
| C7 | −0.217623 | 3.3e−6 | 2.5395808 | 0.5449561 | 0.2811915 | 0.002557965 | 0.006543047 | −0.003811417 |
| ± | | | 3.7e−6 | 7.0e−6 | 4.6e−6 | 1.9e−8 | 1.7e−8 | 1.4e−8 |
| C8 | −0.297797 | 1.9e−5 | 9.6307288 | −6.5812847 | −0.9675765 | 0.001229938 | 0.001140329 | −0.000517445 |
| ± | | | 9.2e−5 | 1.0e−4 | 1.6e−4 | 1.1e−8 | 1.4e−8 | 1.4e−8 |
| C9 | −0.081762 | 1.2e−5 | 9.0170598 | −1.7633125 | −0.0055133 | 0.000468307 | 0.002101810 | −0.001539323 |
| ± | | | 3.7e−5 | 8.6e−5 | 5.6e−5 | 2.0e−8 | 1.4e−8 | 1.6e−8 |
| C10 | +0.033424 | 1.5e−4 | 6.2530222 | 4.9730520 | −4.1195499 | 0.002043816 | 0.004695711 | −0.002871703 |
| ± | | | 3.5e−4 | 1.7e−4 | 2.8e−4 | 1.7e−7 | 1.9e−7 | 1.9e−7 |
| C11 | +0.158105 | 2.4e−5 | 0.3573218 | −4.6247991 | 1.6614283 | 0.005611509 | 0.000440068 | −0.000340169 |
| ± | | | 2.0e−5 | 3.1e−5 | 2.2e−5 | 4.2e−8 | 1.8e−8 | 2.0e−8 |
| C12 | +0.297253 | 1.3e−4 | 12.3162779 | 4.1618475 | 0.0245814 | 0.002109819 | 0.003273857 | −0.001691519 |
| ± | | | 4.1e−4 | 3.8e−4 | 2.2e−4 | 1.0e−7 | 1.3e−7 | 1.1e−7 |

Orbits near perihelion (they change from return to return under the planets' pulls; periods follow
P = 2π [a³ / ((1+η) G m_S)]^{1/2}, which matches the observed return intervals):

| body | a (L) | e | perihelion (L) | period (d) | tilt | returns seen |
|---|---|---|---|---|---|---|
| C1 | 13.3 | 0.947 | 0.70 | 15,720 | 4.6° | 2 |
| C2 | 7.4 | 0.78 | 1.62 | 7600–7780 | 0.1° | 4 |
| C3 | 5.74 | 0.70 | 1.71 | 4520–4550 | 9.9° | 9 |
| C4 | 9.1 | 0.906 | 0.86 | ≈ 11,800 | 1° | 2 |
| C5 | 4.4 | 0.51 → 0.39 | 2.16 → 2.68 | 3100 | 1.3° → 18° | 12 |
| C6 | 5.67 | 0.50 | 2.85 | 5300–5400 | 16° | 4 |
| C7 | 2.04 | 0.50 | 1.02 | 1200 | 17° | 31 |
| C8 | 6.44 | 0.836 | 1.06 | 7100 | 17° | 6 |
| C9 | 5.19 | 0.764 | 1.22 | 4510 | 15° | 8 |
| C10 | 9.08 | 0.909 | 0.82 | 9820 | 14° | 3 |
| C11 | 3.2 | 0.543 | 1.47 | 1950 | 2.6° | 19 |
| C12 | 9.4 | 0.767 | 2.2 | ≈ 9200 | 16° | 4 |

**Quality of fit.** `check2.py` runs the ephemeris of `predict_2.py` (not the fitting code) against every sighting:
reduced χ² = 0.9998 for all 392,128 measured angles of the 19 bodies (σ = 1′). Per body, rms in arcmin (lon·cos lat,
lat): S 1.000/1.005, P1 0.998/0.997, P2 0.999/1.002, P3 0.996/1.003, P4 0.998/0.997, P5 0.995/0.999, M 0.993/1.002,
C1 0.961/1.023, C2 1.008/1.017, C3 1.000/0.989, C4 1.012/0.985, C5 0.982/0.995, C6 1.028/0.970, C7 1.009/1.003,
C8 1.017/1.007, C9 1.009/0.991, C10 0.987/1.003, C11 1.012/1.005, C12 1.008/1.015.

## 3. Predictions

`predict_2.py` integrates Law 1 (with η) and Law 3 for all 20 bodies from t₀ with the constants above (fixed 0.5-day
step, 8th-order Runge–Kutta) and applies Law 2 on the days in `predictions_2.csv` (243 days × 19 bodies, days
36555–43815). Run: `python predict_2.py` (about a minute).

**Test of predictive power (hold-out).** Each faint body's seven constants were refitted using only days 0–29220 and
the result compared with its sightings in days 29221–36525, which that fit never saw (the planet masses were not
refitted for this test, so it is not fully independent). Residuals on the unseen 20 years, arcmin:

| body | unseen sightings | rms along lon | rms along lat | mean offset (lon, lat) |
|---|---|---|---|---|
| C2 | 280 | 2.46 | 1.25 | −2.17, −0.51 |
| C3 | 382 | 1.02 | 1.03 | −0.03, +0.06 |
| C5 | 754 | 0.97 | 1.03 | +0.06, −0.01 |
| C6 | 162 | 1.05 | 0.97 | −0.27, −0.12 |
| C7 | 3787 | 1.00 | 1.02 | +0.11, +0.03 |
| C8 | 412 | 1.11 | 1.00 | +0.41, −0.05 |
| C9 | 627 | 1.05 | 1.01 | +0.16, −0.02 |
| C11 | 1081 | 1.01 | 1.02 | −0.09, +0.05 |
| C12 | 245 | 1.03 | 0.98 | −0.07, −0.11 |
| C1, C4, C10 | 0 | — | — | not testable |

Nine bodies are predicted 20 years ahead to within 0.1–0.4′ of the truth; C2 (three earlier returns, orbit crossing
P4's) is off by 2.2′. I therefore expect the faint-body predictions for years 101–120 to be good to a few tenths of an
arcminute for C3, C5, C7, C9, C11, C12, about 0.5′ for C6 and C8, 1–3′ for C2, and possibly worse for C1, C4 and C10,
which were last seen on days 24074, 13141 and 27299. The bright bodies' predictions are those of stage 1 to within
0.014′ (stage 1 hold-out: indistinguishable from noise).

## 4. Self-assessment

| law | confidence | what would prove it wrong |
|---|---|---|
| Faint bodies orbit S on perturbed ellipses; observed direction is geometric (Law 2) | 97% | A light-delay or parallax signature above ≈ 0.2′ in a close passage. |
| Inverse-square form of the attraction of S and of each planet, for all bodies | 95% | A faint body whose residuals grow with distance range in one apparition (the 12 bodies cover r = 0.7–3.5 L while observed and 3–26 L unobserved between returns, and all return on time). |
| Each faint body has a constant factor (1+η_b) on *all* gravity (Law 1 generalised) | 85% | η changing from return to return (seen constant to ≤ 2×10⁻⁴ over up to 31 returns); a planet-only factor differing from the S factor by more than 0.2%; a correlation of η with distance from S. I cannot say *why* η differs from body to body, nor exclude a mechanism that happens to scale like total gravity. |
| η = 0 for the eight bright bodies | 90% | P²/a³ of a bright body off by more than 2×10⁻⁴ (not seen). For M and Tellus a small η is not separable from their masses. |
| Refined masses of T+M, P2, P3, P4, P5 | 80% that the true values lie within 2 quoted σ | The faint bodies' individual mass estimates still disagree more than their formal errors (C11 and C8 mainly), so some unmodelled effect at the 0.1–0.2′ level remains for bodies passing close to planets. |
| Law 3 (γ term on M) | 60% form / 85% existence, as in stage 1 | As in stage 1. |
| Predictions for years 101–120 | bright bodies: systematic error < 0.1′ (90%), M < 0.3′ (75%); faint bodies: < 0.3′ for most (see hold-out), but up to several arcmin is possible for C1, C4, C10 | These three were last seen 25–65 years before the end of the data and return (or pass planets) inside the prediction window; a close planetary passage amplifies any error in their η or state. |

## Appendix: all constants at full working precision (identical to `predict_2.py`)

G m_S = 2.959e-4 L^3/day^2, t0 = day 18262. `GM_REL` = G m / G m_S; `GAMMA` in 1/day; `STATE` = x, y, z (L), vx, vy, vz (L/day) at t0 relative to S (M relative to Tellus); `FAINT` = (eta, state relative to S at t0).

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
FAINT = {   # name: (ETA, [x,y,z (L), vx,vy,vz (L/day)] relative to S at day T0)
    'C1':  (2.6488613715e-01, [-9.1189895381480959e+00, -2.1824964008081334e+01, 6.5764644673391404e+00, 1.2993603454929613e-03, 7.6101727343148352e-04, -2.0683458611232272e-04]),
    'C2':  (-9.9704278132e-02, [1.1852857969619459e+01, 4.9029343550001849e+00, -1.9736802361150958e+00, -1.0701525560640036e-03, 1.7301984942718165e-03, -6.8716110673063110e-04]),
    'C3':  (2.3058868636e-01, [-3.1166271320118932e+00, 9.2307651831623039e-02, -5.6064150839779547e-01, 6.7611827020680332e-03, -9.3235005362170900e-03, 5.7206846730579039e-03]),
    'C4':  (-2.7360823636e-01, [-2.2599754638888832e+01, 3.1879265093057629e+00, -1.9799016839606856e-01, -2.2568335523305383e-03, -6.1774297913043245e-04, 4.3320042472451956e-04]),
    'C5':  (1.7178701773e-01, [-4.9967321815203253e-01, 5.1939148312102503e+00, -1.9212824305168206e+00, -5.9705686679097699e-03, 2.2978417189274264e-03, -9.0553387468603013e-04]),
    'C6':  (-1.4397289489e-01, [-4.4908159191976500e+00, 7.0123161731723176e+00, -1.6755725536024453e+00, -2.7316571561804550e-03, -3.0479907310770279e-03, 1.3174370046549965e-05]),
    'C7':  (-2.1762328559e-01, [2.5395807503006296e+00, 5.4495610615521117e-01, 2.8119147286943663e-01, 2.5579651700346436e-03, 6.5430470264162916e-03, -3.8114170937963827e-03]),
    'C8':  (-2.9779706945e-01, [9.6307288410462135e+00, -6.5812846551285054e+00, -9.6757651516755394e-01, 1.2299384741256958e-03, 1.1403285472154823e-03, -5.1744545672735359e-04]),
    'C9':  (-8.1761783009e-02, [9.0170597898079574e+00, -1.7633124544658563e+00, -5.5133496016859969e-03, 4.6830680038379124e-04, 2.1018099049162318e-03, -1.5393230856623433e-03]),
    'C10': (3.3423849529e-02, [6.2530221611912111e+00, 4.9730520226929817e+00, -4.1195498914057023e+00, 2.0438158045139373e-03, 4.6957114516685992e-03, -2.8717027502551009e-03]),
    'C11': (1.5810486804e-01, [3.5732181183038170e-01, -4.6247990579348626e+00, 1.6614283485480703e+00, 5.6115087925891439e-03, 4.4006836088981294e-04, -3.4016949472848206e-04]),
    'C12': (2.9725269557e-01, [1.2316277877206435e+01, 4.1618475198974076e+00, 2.4581374014470282e-02, 2.1098194438991186e-03, 3.2738565334331318e-03, -1.6915189251079917e-03]),
}
```
