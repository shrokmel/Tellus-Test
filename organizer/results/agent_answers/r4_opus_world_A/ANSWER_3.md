# ANSWER 3 — The new body N

Short version: N obeys the stage-2 laws.
- **Faint-body law:** it is a faint body with its own gravity factor, β = −0.3235 ± 0.0021, and it fits to 0.99′. Pure gravity leaves 7.8′.
- **P4 encounter:** in 1.2 years (day ≈ 36966) it will pass 0.035 L from the giant P4. The encounter throws it onto a larger, steeper orbit that carries it out to about 20 L from S by year 20.
- **Encounter rule:** I now state explicitly the rule from stage 2 that such a close pass resets β to about zero. That rule is the main uncertainty in the prediction after year 1.

Scripts: `nfit.py` (fit), `npath.py`, `nscen.py`, `nB.py` (path, ensembles, encounter scenarios), `nident.py` (is N a known body?), `check_3.py` (checks predict_3.py), `predict_3.py` (predictions).

---

## 1. The laws used for N

They are the stage-2 laws, with the same equations and constants (see ANSWER_2.md, whose constants are in `predict_2.py`, which `predict_3.py` imports):

- **Law 1 (frame):** **e** = R(ε) **u**, ε = 21.7015° ± 0.0010°. Directions are geometric, from Tellus.
- **Law 2 (Newtonian N-body for S, Tellus, P1–P5):**
  d²**r**ᵢ/dt² = −G(M_S+mᵢ)**r**ᵢ/rᵢ³ + Σ_{j≠i} G mⱼ[(**r**ⱼ−**r**ᵢ)/|**r**ⱼ−**r**ᵢ|³ − **r**ⱼ/rⱼ³].
  Constants:
  - G M_S = 2.958466×10⁻⁴ L³/d²
  - m_T/M_S = (3.0 ± 0.5)×10⁻⁶
  - m_P1/M_S = (0.57 ± 0.17)×10⁻⁶
  - m_P2/M_S = (0.94 ± 0.24)×10⁻⁶
  - m_P3/M_S = (1.4974 ± 0.0010)×10⁻⁴
  - m_P4/M_S = (1.00077 ± 0.00011)×10⁻³
  - m_P5/M_S = (2.4884 ± 0.0014)×10⁻⁴

  The day-0 states are tabulated in ANSWER_2.md.
- **Law 4 (faint bodies):**
  d²**r**/dt² = κ[−GM_S **r**/r³ + Σ_{j∈T,P1..P5} G mⱼ(**r**ⱼ−**r**)/|**r**ⱼ−**r**|³] − Σ_{j∈T,P1..P5} G mⱼ **r**ⱼ/rⱼ³, with κ = 1 − β.
- **Law 4c (encounter rule, now stated as an equation).** If a faint body passes within D = 0.1 L of Tellus or a P, then at the moment of closest approach (d|**r** − **r**ⱼ|/dt = 0) its β is replaced by β_after = 0 ± 0.016.

  Evidence (stage 2) — every recorded pass inside 0.07 L reduced |β| to ≤ 5% of its value:
  - C5: 0.054 L from P4, β −0.172 → −0.001
  - C7: 0.07 L from Tellus, 0.105 → 0.0002
  - C7: 0.07 L from Tellus, −0.021 → −0.0001
  - C1: 0.03 L from Tellus, −0.263 → −0.014
  - C7's 0.04 L pass left an already-zero β at zero.
  - A wider pass (0.16 L) only halved β, so the threshold D lies between 0.07 and 0.16 L. I take D = 0.1 L.
  - The ±0.016 is C1's leftover |β|.

**What changed since stage 2, and why.** No equation or constant changed. Law 4c was only described in words in stage 2 ("β jumps at close encounters, size unknown"). I now give it as an explicit rule with numbers, because N's predicted path runs through a 0.035 L pass of P4, so its prediction depends on that rule. The rule comes from the stage-2 data, not from N's sightings, which end before the encounter.

**Constants for N** (fitted to all 207 sightings, days 36233–36525; heliocentric ecliptic state at epoch t₀ = day 36379; ±1σ):

| quantity | value |
|---|---|
| x, y, z [L] | 0.470283 ± 1.1e-4, 1.612052 ± 1.7e-4, 1.947458 ± 2.9e-4 |
| vx, vy, vz [L/day] | −0.01566212 ± 9.3e-7, −0.00089738 ± 1.4e-6, 0.00443080 ± 1.4e-6 |
| β (until the P4 encounter) | −0.3235 ± 0.0021 |
| β after the encounter (Law 4c) | 0 ± 0.016 |

**N's sightings test the laws:**
- Law 4 (one constant β plus the stage-2 field) fits the arc at 0.990′ rms per angle, the noise level.
- Pure gravity (β = 0) gives 7.8′, so N is, like most C's, a body with κ ≠ 1.
- Every one of 6 starting orbits × 3 starting β converged to the same solution.
- N is not a known faint body seen again: every C's stage-2 orbit is 35°–160° from N on N's sighting days.

---

## 2. Predictions

`predict_3.py` reads `predictions_3.csv` and integrates Law 4 for N from t₀ in the Law-2 field of the stage-2 solution (DOP853, relative tolerance 10⁻¹¹). It detects the closest approach to P4 within 0.1 L (found at day 36966.25, 0.035 L) and resets β to 0 there, as Law 4c says. It then writes lon/lat in the observation frame (Law 1) for all 1461 days, in about 9 s. It agrees with the independent fitting code to 0.02′.

---

## 3. N's path

**Now (day 36525).** N is 2.98 L from S and 3.26 L from Tellus, just past perihelion: it was found at r ≈ 2.6–3 L, at the edge of visibility (faint bodies are seen only inside r ≈ 3 L). It is moving outward and southward, at ecliptic lon 123°, lat +41°.

Its orbit now: a = 10.35 L, e = 0.751, q = 2.58 L, Q = 18.1 L, i = 53.9°, period 28.9 yr (for κ = 1.323). It is steep, steeper than any of C1–C12, and runs from just inside 3 L out past P5's distance.

**Year 1.2: close encounter with P4** (day ≈ 36966, at r ≈ 5.9 L). N passes P4 at 0.035 L, about 8% of P4's Hill radius. Across the ensemble of fits the miss distance is 0.026–0.056 L. Under Law 4c its β resets to ≈ 0, and P4's pull reshapes the orbit to:

a = 11.6 L, e = 0.825, q = 2.02 L, Q = 21.1 L, i = 64°, period 39 yr.

**Years 2–20.** N climbs steadily outward:

| time | r | sky position (ecliptic) |
|---|---|---|
| year 2 | 7.7 L | lon 194°, lat −11° |
| year 5 | 12.7 L | lon 207°, lat −30° |
| year 10 | 17.7 L | lon 218°, lat −43° |
| year 15 | 20.1 L | lon 227°, lat −50° |
| year 20 | 20.5 L (near aphelion, around day 43260+) | lon 237°, lat −56° |

On the sky it drifts slowly east and further south. Its motion is dominated by Tellus's yearly parallax loop, whose size shrinks from tens of degrees to about 3° as N recedes.
- Solar elongation alternates between about 40° and 130° each year, so N passes near conjunction every year. The predictions are given anyway, but N will in any case be too far from S (> 3 L) to be seen after about day 36600.
- After the encounter N passes no other body closer than 1.0 L (P4 again, on its way out). It does not return to the inner system within the 20 years: perihelion is in about year 40.

Without the β reset (if κ stayed 1.323), the same encounter would instead send N back to perihelion at about 1.5 L in year 12.5. That path differs from the prediction by 17–25° after year 10.

---

## 4. Self-assessment

The expected direction error (1σ) combines four parts: the formal fit errors (40-member ensemble drawn from the covariance), the switching time (closest approach vs. entering 0.1 L), β_after = ±0.016, and integration (< 0.001°).

| after | formal | switching time | β_after | expected error | if Law 4c is wrong (β stays −0.32) |
|---|---|---|---|---|---|
| 1 yr (before the encounter) | 0.08° | 0 | 0 | **≈ 0.1°** | same |
| 5 yr | 2.8° | 1.7° | 0.1° | **≈ 3–4°** | ≈ 3–4° (2.5° from the nominal) |
| 20 yr | 8.7° | 2.5° | 2° | **≈ 10°** | ≈ 17–25° off |

The large growth after year 1 comes from the P4 encounter: a 0.035 L pass amplifies the 2×10⁻⁴ L uncertainty of today's position into degrees.

Confidence:
- N obeys Law 4 up to the encounter: 85%.
- The P4 encounter happens closer than 0.06 L: 90%.
- β resets to ≈ 0 at it (Law 4c): 60%.

What would prove me wrong:
1. N's position in the next months departing from `predictions_3.csv` by more than about 0.3° before day 36966. That would mean Law 4 or the fitted β is wrong.
2. After the encounter, N following the no-reset path (back toward the inner system at about year 12) rather than the outbound one. That would falsify Law 4c.
3. N found far (> 10°) from the predicted track in years 5–20 in the direction of the A/B difference, i.e. a different β after the encounter.
4. Any brightening or return of N to r < 3 L before year 20. The predicted orbit does not allow it.
