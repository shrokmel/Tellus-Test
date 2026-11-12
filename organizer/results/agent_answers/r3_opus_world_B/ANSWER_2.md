# Stage 2: the faint bodies C1–C12

Short version:
- **Stage-1 laws, applied exactly, fail badly for every faint body.** No orbit under the Sun's ordinary pull reproduces even a single apparition: the misfit is 1.4′ to 400′.
- **One new constant per faint body removes almost all of this.** I call it κ: each faint body responds to gravity κ times as strongly as a stage-1 body does. The κ values range from 0.70 to 1.30, and κ applies to the pull of the planets as well as of S.
- **With κ, every apparition fits to the 1′ noise.** Every run of 3–4 consecutive apparitions fits to 1.00–1.08′, and 4 of the 12 bodies fit their whole century to ≤1.07′.
- **A smaller problem remains unexplained.** Over the full century some bodies keep 1.2–3′ excess residuals (C7, C8, C9, C11, C2, C3). These are small deflections that I could not attribute to any law I tested.

Units, frames and conventions are as in stage 1: AU defined by k = 0.01720209895, days, ecliptic frame rotated from the data frame by ε = 21.7015°.

---

## 1. Testing the stage-1 laws

**The stage-1 bodies (S, P1–P5, M).** These are unaffected: the faint bodies have no detectable mass. A faint body of even 1e-7 of S's mass passing near Tellus would disturb M and S, and nothing of the kind is seen in the stage-1 fit. Law 1 still fits them at 1.00′ in every decade.

**The faint bodies under the stage-1 law (κ = 1).** For every apparition (a continuous run of sightings) I fitted the best Kepler orbit about S with the stage-1 solar GM. Planetary perturbations are negligible within one apparition, except for C4's close passage of Tellus and P4. The rms misfit per angle (the noise is 1.0′):

| body | rms with κ = 1, by apparition (′) | rms with κ free (′) |
|---|---|---|
| C1 | 134, 168 | 0.97, 1.01 |
| C2 | 4.2, 19.1, 11.2, 18.3 | 0.97–1.03 |
| C3 | 9.0, 22.1, 36.7, 7.4, 10.5, 11.0, 20.5, 6.7, 6.6 | 0.94–1.03 |
| C5 | 23.8, 23.0, 18.3, 20.6, 23.9, 16.2, 17.2, 14.7, 12.1, 5.8, 7.0, 5.1 | 0.94–1.02 |
| C6 | 1.1, 1.9, 4.2, 3.5 | 0.98–1.03 |
| C7 | 7–221 (38 apparitions; median ≈ 57) | 0.97–1.06 |
| C8 | 31, 173, 2.8, 2.2, 108, 104, 32 | 0.96–1.05 |
| C9 | 37, 18, 11, 32, 1.8, 1.1, 9.8, 27, 20 | 0.95–1.02 |
| C10 | 1.1, 3.9, 4.9 | 0.97–1.01 |
| C11 | 1.0–32.7 (19 apparitions; typically 15–31) | 0.95–1.04 |
| C12 | 5.5, 26.6, 9.6, 20.0 | 0.98–1.03 |

C4's apparitions need the planets (see below). When an N-body fit with κ = 1 is forced through all apparitions of a body, it fails completely; for C6, for example, the misfit is 176′.

**The stage-1 law therefore fails for every faint body, by up to hundreds of arcminutes.** The single repair is to let the strength of gravity felt by each faint body be its own constant κ. The fitted κ is the same in every apparition of a given body, to ~1e-3. Examples:
- C7: κ = 0.7825 ± 0.0001 in all 38 apparitions.
- C11: 1.1571–1.1583 over 19 apparitions.
- C5: 1.1709–1.1730 over 12.

κ is also independent of distance from S: C5's perihelion distance moved from 1.9 to 2.7 AU with no change in κ.

**Tests of what κ multiplies.**
- **Sun only, planets at normal strength:** C11's 100-year fit gives 21.9′.
- **All gravity on the body scaled by κ (S and planets):** 1.35′.
- **Separate factors for the Sun and for the planets:** they agree, 1.1593 vs 1.1582.

So κ scales the body's response to all gravity, like a ratio of gravitational to inertial mass. Values above 1 (C1, C3, C5, C10, C11, C12) rule out radiation pressure, which can only weaken the Sun's pull.

**Test of the full law (Newtonian N-body field of the 8 stage-1 bodies, with κ per faint body), whole century, rms per angle:**

| | C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100-yr fit (′) | 1.42 | 1.22 | 1.29 | **1.00** | 1.10 | **1.00** | 2.15 | 2.95 | 1.61 | 1.14 | 1.35 | **1.07** |
| any window of 3 (C11: 4) apparitions (′) | – | 1.03–1.07 | 0.96–1.04 | – | 0.98–1.13 | – | 0.99–1.09 | 0.99–1.57 | 1.00–1.37 | – | 0.98–1.12 | 1.02–1.03 |
| final fit window (′) | 1.42 (all) | 1.06 | 1.03 | 1.00 (all) | 1.00 | 1.00 (all) | 1.08 | 1.02 | 1.02 | 1.14 (all) | 0.99 | 1.07 (all) |

C4's apparitions include passes 0.064 AU from Tellus and 0.059 AU from P4. A pure two-body orbit misses them by 400′, but the N-body law with κ = 0.7264 fits all 501 sightings at 0.998′. This is a strong confirmation of the planets' pull, scaled by κ.

**What is still wrong.** Over a whole century several bodies keep a real excess misfit, beyond what the fit can remove:
- C8's first apparition is predicted 9.3′ off by a fit to its other six.
- C9's κ shifts by 0.2% between the early and late record.
- C11's along-track timing wanders with a period of ~18,000 d, close to its synodic period with P3.
- C2 drifts by 22′ over one revolution.

Where the excess lies, a fit to any 3–4 consecutive apparitions still reaches noise level, so the extra effect accumulates slowly or acts in isolated episodes.

I tested and **rejected** the following explanations; none removes the excess:
1. Outgassing-type forces A₁, A₂, A₃ ∝ 1/r² (fitted values ≲ 1e-8 AU/d², no improvement).
2. A power-law change of the 1/r² exponent for faint bodies (|δ| < 1e-5).
3. Light-time (makes it slightly worse).
4. A wrong mass for P3: the comets disagree on the correction, which ranges from −2.5% to +0.2%.
5. A single unseen planet. A linearized search over circular orbits from 4 to 80 AU at all phases improves χ² by at most 1.5% of the excess, and only for unphysical negative masses.

I therefore report the excess as an unexplained residual.

## 2. The laws (equations) and constants

The laws for the stage-1 bodies are **unchanged** from stage 1. One law is new, forced by the faint bodies.

**Law 1 (massive bodies, unchanged).** For i, j ∈ {S, P1, P2, T, P3, P4, P5, M}:

  d²**r**_i/dt² = Σ_{j≠i} GM_j (**r**_j − **r**_i)/|**r**_j − **r**_i|³

with GM_j = μ_j k². All the constants (the μ_j and the 42 epoch state components) are exactly those of ANSWER_1.md §2, and they are copied verbatim into `predict_2.py`. Main values: μ_P3 = 1.5057e-4, μ_P4 = 9.9992e-4, μ_P5 = 2.5013e-4, μ_T = 2.608337e-6, μ_M = 3.19492e-8.

**Law 1′ (faint bodies, new).** Each faint body C_k is a massless test particle. Its acceleration is the stage-1 gravitational field multiplied by its own constant κ_k:

  d²**r**_k/dt² = κ_k Σ_j GM_j (**r**_j − **r**_k)/|**r**_j − **r**_k|³   (j over the 8 massive bodies)

For every stage-1 body κ = 1. This is checked to ~1e-5 by Kepler's third law in stage 1.

**What forced the change.** With κ = 1, no orbit fits any single apparition of most faint bodies, with misfits up to 400′. With κ free, every apparition fits to noise and κ comes out the same in every apparition of a body.

**Law 2 (what is seen, unchanged).** **u** = unit(**r**_body − **r**_T), rotated to the data frame by ε = 21.7015°.

**Constants for the faint bodies.** For each body I give κ and its heliocentric ecliptic state at the stated epoch (day number). Positions are in AU and velocities in AU/day. The ± values are formal 1σ.

Because of the unexplained residual, each state was fitted to the body's most recent apparitions, so that it predicts the future best:
- the last 3 apparitions: C2, C3, C8, C9
- the last 4: C5, C11
- days ≥ 27000: C7
- all data: C1, C4, C6, C10, C12

| body | epoch | x | y | z | vx | vy | vz | κ (formal) |
|---|---|---|---|---|---|---|---|---|
| C1 | 23898.0 | 0.506404396 ± 7.2e-06 | 0.507513617 ± 7.3e-06 | 0.050463572 ± 1.3e-05 | -0.025816622 ± 1.2e-08 | 0.018607438 ± 5.4e-08 | 0.000870351 ± 3.6e-07 | 1.264859 |
| C2 | 21667.0 | -1.160452616 ± 3.0e-05 | 2.517949935 ± 6.4e-05 | 0.001153898 ± 3.2e-05 | -0.005929916 ± 1.4e-07 | -0.011008997 ± 2.5e-07 | 0.000009893 ± 1.7e-07 | 0.900273 ± 0.000062 |
| C3 | 31750.0 | -2.954550190 ± 8.3e-05 | -0.244958299 ± 7.2e-05 | -0.415541390 ± 4.3e-05 | 0.008345286 ± 3.2e-07 | -0.010307968 ± 3.8e-07 | 0.002106564 ± 2.1e-07 | 1.230460 ± 0.000106 |
| C4 | 12899.0 | 0.843601297 ± 5.5e-06 | -0.163851448 ± 8.3e-06 | 0.001753946 ± 6.8e-06 | 0.003515555 ± 1.7e-07 | 0.021544575 ± 1.6e-07 | -0.000385301 ± 1.5e-07 | 0.726379 ± 0.000011 |
| C5 | 32576.0 | 0.026911592 ± 7.3e-05 | -2.919341485 ± 9.2e-05 | -0.677999780 ± 5.1e-05 | 0.011648860 ± 4.9e-07 | 0.002434986 ± 3.9e-07 | 0.003241944 ± 2.1e-07 | 1.171641 ± 0.000107 |
| C6 | 20446.0 | 1.028443234 ± 7.9e-05 | -2.616210486 ± 1.1e-04 | -0.464018937 ± 3.9e-05 | 0.010707343 ± 5.9e-07 | 0.003545876 ± 4.4e-07 | 0.002532990 ± 3.7e-07 | 0.855992 ± 0.000127 |
| C7 | 31638.0 | 2.547084268 ± 2.0e-05 | 1.571963572 ± 1.9e-05 | 0.193389822 ± 1.5e-05 | -0.002005569 ± 3.0e-08 | 0.005796784 ± 4.0e-08 | -0.001736699 ± 4.0e-08 | 0.782572 ± 0.000017 |
| C8 | 28914.0 | 2.392854199 ± 3.2e-05 | 1.786667281 ± 5.5e-05 | -0.197643775 ± 3.7e-05 | -0.010010764 ± 1.3e-07 | 0.000571485 ± 1.2e-07 | 0.002458391 ± 1.9e-07 | 0.702293 ± 0.000029 |
| C9 | 29539.0 | -0.960356607 ± 2.2e-05 | -1.197709271 ± 5.8e-05 | 0.376715080 ± 2.7e-05 | 0.007395383 ± 3.3e-07 | -0.015200284 ± 3.1e-07 | 0.002373027 ± 2.1e-07 | 0.918179 ± 0.000069 |
| C10 | 17251.0 | -0.132071195 ± 4.8e-05 | -0.813350975 ± 2.2e-05 | 0.100824734 ± 2.0e-05 | 0.024956562 ± 3.0e-07 | -0.007263866 ± 8.5e-07 | -0.005255500 ± 3.3e-07 | 1.033723 ± 0.000004 |
| C11 | 31333.0 | -2.671925774 ± 4.6e-05 | -1.360870359 ± 4.2e-05 | 0.016773817 ± 2.4e-05 | -0.001084205 ± 1.4e-07 | -0.011021185 ± 1.9e-07 | -0.000392072 ± 1.1e-07 | 1.157893 ± 0.000054 |
| C12 | 34671.0 | -1.743761452 ± 3.2e-05 | -1.355056271 ± 3.3e-05 | -0.045492401 ± 2.9e-05 | 0.010435967 ± 2.4e-07 | -0.013247504 ± 2.6e-07 | 0.004751607 ± 2.3e-07 | 1.296555 ± 0.000039 |

`predict_2.py` holds the full-precision values.

**Realistic uncertainty of κ.** Here I take the spread between independent apparition and window fits, which is larger than the formal errors:

| | C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| κ | 1.2649 | 0.9003 | 1.2305 | 0.7264 | 1.1716 | 0.8560 | 0.7826 | 0.7023 | 0.918 | 1.0337 | 1.158 | 1.2966 |
| ± | 0.0005 | 0.0003 | 0.0006 | 0.0002 | 0.0003 | 0.0005 | 0.0002 | 0.0004 | 0.002 | 0.0004 | 0.0008 | 0.0008 |

Other derived quantities:

| body | a (AU) | e | i | period (d) |
|---|---|---|---|---|
| C1 | 13.3 | 0.947 | 3.4° | ≈15,700 |
| C2 | 7.4 | 0.78 | 0.1° | 7,780 |
| C3 | 5.75 | 0.70 | 9.8° | 4,545 |
| C4 | 9.1 | 0.91 | 1.0° | 11,800 |
| C5 | 4.39 | 0.51 | 1.3° | 3,107 |
| C6 | 5.73 | 0.50 | 16° | 5,421 |
| C7 | 2.03 | 0.50 | 17.3° | 1,198 |
| C8 | 6.44 | 0.84 | 16.9° | 7,127 |
| C9 | 5.19 | 0.77 | 15.6° | 4,507 |
| C10 | 9.07 | 0.91 | 14.1° | 9,807 |
| C11 | 3.24 | 0.54 | 2.6° | 1,980 |
| C12 | 9.5 | 0.77 | 15.8° | 9,385 |

Here the period is 2π (a³/κk²)^{1/2}.

**How the constants were obtained.**
1. The 8 massive bodies are integrated from the stage-1 constants and tabulated every 0.5 d (positions are Hermite-interpolated between table points).
2. Each faint body is integrated as a test particle with DOP853 at rtol = 1e-11.
3. Its 6 state components and κ are fitted by Levenberg–Marquardt least squares to its sighting directions.

Applied to synthetic data made from the same model, the fitter recovers 0.994′ rms.

## 3. Predictions

`predict_2.py` is self-contained:
1. It integrates the 8 massive bodies with Law 1, using the stage-1 constants. S, P1–P5 and M come out identical to stage 1.
2. It integrates each faint body with Law 1′ from its epoch.
3. It applies Law 2.

It writes all 4,617 rows of `predictions_2.csv` (days 36555–43815). Run on observed days, the script reproduces each faint body's recent sightings at 1.37–1.79′ total angle (the noise alone gives 1.41′).

**Backtest.** I refitted each body without its last apparition and predicted that apparition 1,960–9,350 d ahead. The mean miss:

| body | mean miss (′) | ahead (d) |
|---|---|---|
| C5 | 1.6 | – |
| C6 | 1.3 | – |
| C7 | 1.25 | – |
| C8 | 1.5 | – |
| C9 | 1.3 | – |
| C11 | 1.25 | – |
| C3 | 3.1 | 4,300 |
| C12 | 5.5 | 9,300 |
| C2 | 22 | 7,800 |

The first six are noise-level, since pure noise alone gives a mean miss of 1.25′. Expected errors for years 101–120:
- S, P1–P5, M: as in stage 1 (≪ 1′).
- Most faint bodies: ~1–3′.
- C2, C3, C12: up to ~5–20′ late in the interval.
- C1 and C4 (no apparition after days 24074 and 13141) and C10 (three apparitions): several arcminutes or more. Their next returns lie far ahead of the data and they pass close to planets.

## 4. Self-assessment

| law | confidence | what would prove it wrong |
|---|---|---|
| Law 1 (stage-1 Newtonian N-body) for S, P1–P5, T, M | 97% | Unchanged from stage 1. Any systematic drift above ~1′ of those bodies in years 101–120. |
| Faint bodies are massless test particles in that field | 95% | A faint body visibly disturbing M, S or a planet during a close pass. |
| Law 1′: each faint body's response to all gravity is scaled by a constant κ_k (0.70–1.30), and κ is the same for the pull of S and of the planets | 85% | κ of one body changing between apparitions by much more than ~1e-3. Or a close planetary pass that requires the planet's pull on that body to be scaled differently from the Sun's. C4's double close pass already supports a common κ. |
| κ is constant in time for each body | 70% | The 100-year excess residuals (C2, C8, C9, C11) may mean κ or the orbit changes slowly or in episodes. A future apparition missing its predicted place by more than ~5′ (for C5–C9, C11) would show the remaining effect is not small. |
| Law 2 unchanged (no light-time or aberration) | 90% | Fits that improve systematically when light-time is included (I found the opposite). |
| The remaining 1–3′ excess is a real, unexplained effect, not fitting error | 80% | A better global optimization fitting all 100 years of C7/C8 to 1.0′ with constant κ would show it was a fitting artifact. |

**What would prove the whole picture wrong:** a faint body whose apparition cannot be fitted to ~1′ by *any* κ with the Sun and planets acting together; or two faint bodies on nearly the same orbit moving with different κ, at a level that couples to some measurable property of them.
