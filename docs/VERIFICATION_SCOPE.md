# Verification scope

The main result is an exact computer-assisted theorem. That phrase has a
specific boundary here: the paper contains the infinite reductions and
coverage proof, while programs perform finite exact arithmetic and verify
large symbolic certificates. There is not yet a proof-assistant
formalization of the whole theorem.

## The theorem

For every pair of integers `r >= t >= 2`,

```text
Omega_t(r) >= Omega_2(1350),
```

with equality only at `(r,t)=(1350,2)`. The sharp raw bridge factor is

```text
2 Omega_2(1350) = 2.655018313913361199726406175676...
```

Equivalently, the raw coefficient range is `0<=q<=r-2`. The excluded `t=1`
edge is not covered and the corresponding unrestricted sentence would be
false; for example, `Omega_1(6)=0.9824439139...<1`.

The target is a reduced rational. Its canonical `numerator/denominator`
fingerprint is

```text
2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691
```

## What is proved symbolically

- The wall recurrence and reverse-row identities.
- Positivity and orientation of the square-tail decrement, inherited from
  Paper 3.
- The lossless transfer
  `Omega_t(r) >= min_(t<=j<=r) rho(j,t)`.
- The two terminal formulas for `rho`.
- The top-region estimate and its exact re-anchoring at row 504.
- The top/joint dichotomy and the polynomial identity
  `(t-2)^2-5t=(t-9)^2+9(t-9)+4`.
- The logical assembly and uniqueness of equality.

These are half-line or all-parameter arguments. They do not rest on a finite
sample.

## What is proved by finite exact computation

- The `t=2` bridge prefix on `2 <= r <= 3001`, with equality only at 1350.
- The `t=3` bridge prefix on `3 <= r <= 3001`.
- The `t=4` square-tail prefix on `4 <= r <= 1000`.
- The 123,753-cell core
  `7 <= r <= 503`, `5 <= t <= r-2`.
- The fixed-column prefixes for `t=5,6,7,8`.
- Exact Wallis-level crossing and anchor gates for `t=9,10`.
- The D-LAW finite ladder on `60000 <= m <= 262144`, all seven indices,
  with a complete outward interval anchor at the final row.

Each claim is limited to its displayed range. The unbounded complement is
handled separately.

## What is proved by exact half-line certificates

- Positive-coefficient tail polynomials for the `t=2` and `t=3` bridge
  stripes.
- Positive-coefficient square-tail tails for columns `t=4,...,10`.
- Exact Sturm certificates for the joint core at `t=11,...,16`.
- One 105-monomial nonnegative bivariate certificate for every `t>=17`.
- One single-anchor D-LAW envelope certificate over 57 contiguous windows,
  with lower-envelope positivity checked on every window and a rejected
  `4/5` negative control.

All polynomial coefficients and Sturm signs are rational or integer. A
decimal capacity is a rendering of an exact rational after the sign test; it
is never the gate.

## What is not claimed

- No higher fixed stripe is claimed to have sharp minimum `4/3`.
- The reported numerical locations of possible higher-stripe minima are not
  theorems.
- Pointwise monotonicity in `t` is false, with an exact witness at row 45082.
- D-LAW is proved only on its stated index and level window. The limiting
  profile, saturation constant, and uniform-in-`t` transfer remain open.
- The theorem is not described as Lean-verified or kernel-checked.
- The two proof chains for each of `t=9` and `t=10` share definitions and
  upstream facts; they are not claimed to be fully independent
  implementations.

## Master join

The final join has 22 named gates. It rebuilds the extremal value, verifies
the exact `4/3` comparison, checks the closed `t=2,3,4` records, checks every
regional artifact and splice, validates both D-LAW legs and their negative
control, and checks symbolic coverage arithmetic.
Current status:

```text
22/22 PASS
conditional legs: none
failures: none
```

## Interval-constructor disclosure

An inherited fixed-point constructor used a floor where an upper endpoint
required a ceiling. The replay detects and repairs the one-unit defect before
building certificates. Every present consumer immediately divides the
interval by two; exact outward rounding makes the corrected and inherited
inputs identical at `pi/2`, so every downstream coefficient delta is exactly
zero. This is a source-contract defect with zero theorem blast radius, not a
numerical approximation argument.
