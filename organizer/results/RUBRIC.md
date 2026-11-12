# Grading rubric (fixed in advance)

Each folder g1-g8 holds one analyst's written answers (ANSWER_2.md, ANSWER_3.md) about faint bodies moving in a
star system, observed as directions only. Grade ONLY what the analyst states as laws or equations (in either
file), not general prose or speculation. For each element answer present / partial / absent, with a short quote.

The true law (for your reference; the analyst never saw this):
each body has a mass m (inertia) and a hidden "charge" q; gravity between bodies i and j is G q_i q_j / r^2;
pull factor f = q/m (1 for the star, planets and moon). Within ~0.1 AU of another body, charge flows from the
higher-f body to the lower-f one: d(f_i - f_j)/dt = -lam * exp(-(r/d)^2) * (f_i - f_j), lam = 50 per year
(= 0.137 per day), d = 0.1 AU; total charge is conserved, shared in proportion to mass (a light body takes on
a heavy body's f). (In one of the two worlds the analysts studied, lam = 0: charge never flows.)

Elements:
- E1 Pull factor: each faint body has its own constant multiplying the attraction it feels (any notation:
  kappa, 1-beta, 1+eta ...). "present" if stated as an equation.
- E2 Changes at passes: the analyst states (as a law) that a body's factor changes, and only during close
  passes to another body. "absent" if the analyst states the factor is constant for all time.
- E3 Evening out: the change moves the faint body's factor toward the other body's value (toward 1 / toward the
  planets' value, or beta toward 0, counts as present).
- E4 Conservation: a conserved quantity (charge) is exchanged, weighted by mass, so a small body takes on a
  large body's value. Must be explicit.
- E5 Rule constants: a contact distance and a rate (or equivalent) within a factor of 3 of 0.1 AU and
  0.137/day. "partial" if only a distance threshold is given without a rate.

Output: a JSON object {"g1": {"E1": "present|partial|absent", "E1_quote": "...", ..., "notes": "..."}, ...}.
