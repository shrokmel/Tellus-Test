# The sky of Tellus: motions, laws, predictions

All numbers below come from fitting the 156,178 sightings in `observations.csv` (days 0–36525).
Scripts used: `fit.py`, `nbody.py`, `stageA.py`, `stageB.py`, `stageC.py` (fitting), `predict.py` (prediction).

## 1. Description of the motions

**Frame.** `S` travels once round the sky in 365.2575 days along a great circle tilted 21.70° to the
`lat = 0` plane, crossing it at lon ≈ 0°/180°. All other bodies stay within a few degrees of that circle.

**S.** Its motion along the circle is not uniform: it is the motion of a Kepler ellipse of eccentricity 0.0200,
seen from the focus. Read heliocentrically: Tellus goes round S on an ellipse with a = 0.99998 L (L = length unit,
defined in section 2), e = 0.01997, sidereal period 365.2575 d.

**P1–P5.** Each wanders along the same band, mostly eastward, with periodic retrograde loops.
P1 never gets more than 29° from S (it swings from one side of S to the other every 79 d; synodic period 158 d).
P2–P5 can be anywhere relative to S and retrograde around opposition. All of this is reproduced to the
measurement noise by letting the bodies go round S, not round Tellus:

| body | a (L) | sidereal period (d) | e | tilt of orbit to Tellus's orbit | mass / mass of S |
|---|---|---|---|---|---|
| P1 | 0.44999 | 110.258 | 0.0600 | < 0.01° | (3 ± 6)×10⁻⁷ (not detected) |
| Tellus (+M) | 0.99998 | 365.2575 | 0.0200 | 0 | (3.010 ± 0.035)×10⁻⁶ |
| P2 | 1.7000 | 809.6 | 0.0501 | ≈ 0.01° | (6.7 ± 1.2)×10⁻⁷ |
| P3 | 3.299 | ≈ 2189 | 0.040–0.044 | 1.5° | (1.5026 ± 0.0050)×10⁻⁴ |
| P4 | 6.204 | ≈ 5640 | 0.049 | 1.0° | (9.9986 ± 0.0055)×10⁻⁴ |
| P5 | 11.45–11.52 | ≈ 14170 | 0.041–0.043 | 2.0° | (2.5002 ± 0.0072)×10⁻⁴ |

Period² / a³ is the same for all six (to 0.1%, the level at which the mutual pulls disturb it).

**M.** M goes round *Tellus*: sidereal month 27.2897 d, synodic month 29.49 d, mean distance 0.002552 L
(1/392 of the S–Tellus distance), e ≈ 0.056, orbit tilted 5.0° to Tellus's orbit. Its ellipse is strongly disturbed by S:
the line of apsides advances once in 3237 d (8.86 yr; anomalistic month 27.5217 d), the nodes regress once in 6805 d
(18.63 yr; draconic month 27.1807 d), and the longitude carries large periodic terms (6.31° at 27.52 d, 1.28° at 31.77 d,
0.66° at 14.75 d, 0.22° at one year, 0.035° at 29.49 d …). All of these come out of the law in section 2 without
separate constants.

**Changes over the 100 years.**
- No orbit is a fixed ellipse. Fitting a fixed ellipse to each body leaves systematic residuals of 1.7′ (P2), 8.1′ (P3),
  2.0′ (P4), 3.1′ (P5) rms, against 1.00′ noise; these vanish when the bodies attract one another.
- Tellus's perihelion advances 0.26° per century, P1's 0.11°, P2's 1.09°; Tellus's orbital plane shifts by 0.004° per
  century (tilt to the lat = 0 plane 21.7001° → 21.7036°, node 180.001° → 180.009° between the first and last decade). P3, P4, P5 exchange energy: e.g. P5's osculating
  a drifts 11.515 → 11.447 L and P3's e moves between 0.040 and 0.044 (P4 has 1/1000 of S's mass and dominates this).
- S itself shows a 29.5-day wobble of ≈ 0.1′: Tellus circles the Tellus–M mass centre (M/Tellus mass = 0.0116 ± 0.0002).
- M's month lengthens: its mean motion decreases by a fraction 1.25×10⁻¹² per day (4.6×10⁻⁸ per century), i.e. M falls
  behind a constant-rate motion by ≈ 0.7′ × (centuries)². Small, but detected at 7σ.

## 2. The laws, as equations

Bodies: S, P1, Tellus (T), M, P2, P3, P4, P5, each a point mass. Positions **r**ᵢ(t) in the fixed star frame of the data
(x towards lon = 0°, lat = 0°; z towards lat = +90°), t in days.

**Law 1 (universal gravitation; one law for every body).**

    d²r_i/dt² = Σ_{j≠i}  G m_j (r_j − r_i) / |r_j − r_i|³            for all 8 bodies

**Law 2 (what is observed).** The direction of body b on the sky is the direction of the straight line from Tellus to b at
that instant (no measurable light delay, no parallax from the observer's place on Tellus):

    g = r_b − r_T ,   lon = atan2(g_y, g_x) ,   lat = arcsin(g_z / |g|)

**Law 3 (slow secular term for M only; empirical).** The relative motion **ρ** = r_M − r_T has one extra term,

    d²ρ/dt² = [Law 1] + γ · dρ/dt ,        γ = (4.18 ± 0.58)×10⁻¹³ day⁻¹

(applied so that the T–M mass centre is unaffected: M gets +γ ρ̇ m_T/(m_T+m_M), T gets −γ ρ̇ m_M/(m_T+m_M)).
It makes M's orbit widen by a fraction 2γ per day and its mean motion drop by 3γ per day.

**Consequences (not separate laws).** For one small body round S, Law 1 gives Kepler's laws:
r = a(1−e²)/(1+e cos ν), equal areas in equal times, and n²a³ = G(m_S + m), i.e.
P = 2π [a³ / G(m_S + m)]^{1/2} = 365.2644 d × (a/L)^{3/2} / (1 + m/m_S)^{1/2}. The apsidal/nodal motions of M and the planetary drifts follow from the
same Law 1 with all bodies included.

**Constants.** Unit of length L is defined by G m_S = 2.959×10⁻⁴ L³ day⁻² exactly (angles alone cannot give an absolute
scale; L is within 2×10⁻⁵ of Tellus's mean distance from S). Epoch t₀ = day 18262. Values ± 1σ.

Mass ratios G m / G m_S:

| body | value | ± |
|---|---|---|
| P1 | 2.727879×10⁻⁷ | 6.1×10⁻⁷ |
| T | 2.975294×10⁻⁶ | 3.5×10⁻⁸ |
| M | 3.437505×10⁻⁸ | 6.1×10⁻¹⁰ |
| P2 | 6.688391×10⁻⁷ | 1.2×10⁻⁷ |
| P3 | 1.502585×10⁻⁴ | 5.0×10⁻⁷ |
| P4 | 9.998619×10⁻⁴ | 5.5×10⁻⁷ |
| P5 | 2.500207×10⁻⁴ | 7.2×10⁻⁷ |

Positions (L) and velocities (L/day) at t₀, relative to S (for M: relative to Tellus). The exact 17-digit values used
are in `predict.py`; rounded values and 1σ here:

| body | x | y | z | vx | vy | vz |
|---|---|---|---|---|---|---|
| P1 | 0.09106480 | 0.40711317 | −0.16201501 | −0.025512640 | 0.003464343 | −0.001379139 |
| ± | 6.8e−6 | 2.9e−6 | 3.0e−6 | 1.7e−7 | 2.7e−7 | 1.9e−7 |
| T | −0.40776512 | −0.86627329 | 0.34475248 | 0.015420531 | −0.006433098 | 0.002560469 |
| ± | 2.7e−6 | 1.7e−6 | 1.5e−6 | 2.8e−8 | 1.3e−7 | 4.9e−8 |
| M (rel. T) | 0.0024876286 | −0.0001719213 | −0.0001311762 | 0.00004830638 | 0.00056738506 | −0.00019499166 |
| ± | 9.5e−6 | 6.5e−7 | 5.1e−7 | 1.8e−7 | 2.2e−6 | 7.5e−7 |
| P2 | 1.21131541 | −1.01290400 | 0.40335320 | 0.009488618 | 0.009246203 | −0.003679772 |
| ± | 4.1e−6 | 4.7e−6 | 3.7e−6 | 2.8e−8 | 2.7e−8 | 3.0e−8 |
| P3 | 0.06576793 | −3.12579403 | 1.28736371 | 0.009225502 | 0.000356533 | −0.000375835 |
| ± | 1.1e−5 | 6.3e−6 | 8.3e−6 | 2.1e−8 | 2.1e−8 | 2.3e−8 |
| P4 | 1.31752714 | 5.78471711 | −2.42248731 | −0.006586015 | 0.001060351 | −0.000406582 |
| ± | 2.0e−5 | 1.2e−5 | 1.6e−5 | 2.2e−8 | 1.6e−8 | 1.7e−8 |
| P5 | −10.98116377 | −0.39631146 | 0.54566286 | 0.000226957 | −0.004936010 | 0.001897774 |
| ± | 1.6e−5 | 3.6e−5 | 3.1e−5 | 2.1e−8 | 8.4e−9 | 1.3e−8 |

(The M row's large ± on x and vy is one shared uncertainty: the overall size of M's orbit, known to 0.4%, tied to
G(m_T+m_M) = (3.010 ± 0.035)×10⁻⁶ G m_S. M's *direction* is fixed far more tightly.)

**Quality of fit.** With these 50 constants the model reproduces all 312,356 measured angles with χ² = 311,655
(reduced χ² = 0.998 for σ = 1′). Residual rms per body (arcmin): S 1.002, P1 0.997, P2 1.000, P3 0.999, P4 0.998,
P5 0.997, M 0.998; no trend in any 10-year block exceeds 0.045′.

## 3. Predictions

`predict.py` integrates Laws 1 + 3 from t₀ with the constants above (fixed 0.5-day step, 8th-order Runge–Kutta; M is
carried as the vector **ρ** to keep rounding error below 10⁻⁵ arcmin) and applies Law 2 on the requested days. It fills
`predictions.csv` (243 days × 7 bodies, days 36555–43815). Run: `python predict.py`.

**Test of predictive power (hold-out).** I refitted every constant using only days 0–29220 (the first 80 years) and
predicted the last 20 years of data, which that fit never saw. Residuals of those unseen 20 years (arcmin):

| body | rms along lon, days 29221–32872 | rms along lat | mean offset (lon, lat) | rms lon, days 32873–36525 | rms lat | mean offset (lon, lat) |
|---|---|---|---|---|---|---|
| S | 0.989 | 1.007 | +0.008, +0.015 | 1.000 | 0.996 | +0.001, −0.018 |
| P1 | 0.982 | 1.004 | −0.044, +0.014 | 1.018 | 0.978 | −0.040, +0.018 |
| P2 | 0.990 | 0.999 | +0.021, −0.004 | 1.011 | 1.011 | +0.001, −0.014 |
| P3 | 0.986 | 0.978 | −0.005, −0.011 | 0.998 | 0.982 | +0.029, +0.035 |
| P4 | 0.997 | 0.990 | +0.001, +0.007 | 1.003 | 1.003 | −0.024, +0.039 |
| P5 | 0.976 | 0.990 | −0.077, −0.046 | 1.000 | 0.995 | +0.010, +0.007 |
| M | 1.015 | 0.999 | −0.021, −0.008 | 1.014 | 0.997 | +0.000, −0.001 |

So a 20-year extrapolation is indistinguishable from the 1′ measurement noise for every body, M included; the 80-year
fit gives the same constants within 1σ (e.g. γ = 4.25×10⁻¹³ day⁻¹). The submitted predictions use the fit to all 100
years. Expected systematic error of the predictions for years 101–120: ≲ 0.05′ for S and P1–P5, ≲ 0.2′ for M
(dominated by the uncertainty of γ: ±0.05′ at year 120, and of M's mean motion).

## 4. Self-assessment

| law | confidence | what would prove it wrong |
|---|---|---|
| Bodies P1–P5 and Tellus orbit S; M orbits Tellus | 99% | Any P seen at a place no S-centred orbit allows (e.g. P1 beyond 29° from S). |
| Law 1, inverse-square gravitation between *all* bodies with the masses above | 97% | Residuals of the years 101–120 growing beyond noise in a pattern (e.g. a perihelion advance not produced by the other bodies; a failure of P²∝a³ for a newly found body). P3–P5 are the sharp test: their mutual perturbations reach tens of arcmin. |
| Law 2, geometric direction from Tellus's centre, no light delay | 90% | A residual in P1/P2 with the synodic period and size ∝ distance (light delay), or a parallax term in M depending on where on Tellus the observer stands. The fit leaves no room for any unmodelled term above ≈ 0.1′ rms in any body. |
| Masses: P3, P4, P5, Tellus+M | 95% | — (0.06–1% precision). P1's mass is not detected (consistent with 0); P2's is a 5σ detection. |
| Law 3, the γ term on M | 60% for the form, 85% that some secular slowing of M exists | The term is 7σ (Δχ² = 55 for one constant) but I cannot tell its form: a tidal exchange with Tellus's spin, or any other slow change of M's mean motion, fits equally. If M's longitude in years 101–120 follows the pure Law 1 curve instead (difference grows to ≈ 0.35′ by year 120), Law 3 is wrong. |
| 20-year prediction accuracy | planets and S: systematic error < 0.1′ (90%); M: < 0.3′ (75%) | — |

Known limits: absolute distances are unknown (only ratios to L); the noise (1′) hides any effect below about 0.02′ per
body-decade; bodies that are never observed (if any exist beyond P5) would show first as a slow drift of P5.

## Appendix: the constants at full working precision (identical to `predict.py`)

G m_S = 2.959e-4 L^3/day^2, t0 = day 18262. `GM_REL` = G m / G m_S; `GAMMA` in 1/day; `STATE` = x, y, z (L), vx, vy, vz (L/day) at t0, relative to S (M relative to Tellus).

```
GM_REL = {   # GM_body / GM_S
    'P1': 2.7278790470e-07,
    'T':  2.9752939618e-06,
    'M':  3.4375049071e-08,
    'P2': 6.6883910260e-07,
    'P3': 1.5025853709e-04,
    'P4': 9.9986185473e-04,
    'P5': 2.5002068485e-04,
}
GAMMA = 4.178126e-13   # 1/day : secular along-track acceleration coefficient of M
STATE = {   # at day T0: x,y,z [L], vx,vy,vz [L/day]; relative to S (M: relative to Tellus)
    'P1': [9.1064796662090955e-02, 4.0711317232946975e-01, -1.6201501250526185e-01, -2.5512639971106227e-02, 3.4643427437372246e-03, -1.3791392777226775e-03],
    'T':  [-4.0776511809444843e-01, -8.6627328675652981e-01, 3.4475248009507425e-01, 1.5420530812504218e-02, -6.4330980864530894e-03, 2.5604692349140575e-03],
    'M':  [2.4876285751986205e-03, -1.7192134000170546e-04, -1.3117618115432645e-04, 4.8306378175254615e-05, 5.6738506415462953e-04, -1.9499165628962406e-04],
    'P2': [1.2113154114620732e+00, -1.0129040011801935e+00, 4.0335319792128582e-01, 9.4886175917886645e-03, 9.2462030680732707e-03, -3.6797723924188413e-03],
    'P3': [6.5767934024599795e-02, -3.1257940334883516e+00, 1.2873637087575227e+00, 9.2255020854899832e-03, 3.5653305429762360e-04, -3.7583537456235399e-04],
    'P4': [1.3175271434769387e+00, 5.7847171102032453e+00, -2.4224873113543706e+00, -6.5860153798669718e-03, 1.0603505997303056e-03, -4.0658244403189080e-04],
    'P5': [-1.0981163772060590e+01, -3.9631145766255582e-01, 5.4566285928723568e-01, 2.2695668132832275e-04, -4.9360101399250249e-03, 1.8977743474789051e-03],
}
```
