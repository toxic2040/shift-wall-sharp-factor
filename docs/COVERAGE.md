# Final coverage ledger for `rho(j,t)>4/3`, `j>=t>=5`

Status: **complete and unconditional**.

Every integer cell `(j,t)` with `j>=t>=5` lies in exactly one of the cases
below.

## 1. Terminal pair

- `j=t`: `rho(t,t)=t^2/2 >= 25/2 > 4/3`.
- `j=t+1`: `rho(t+1,t)>2(t+1)^2/9 >= 8 > 4/3`.

These are all-parameter symbolic theorems.

## 2. Finite core

For

```text
7 <= j <= 503,
5 <= t <= j-2,
```

`core/A_CORE_SWEEP_SUMMARY.json` records 123,753 exact cells, zero
failures, and zero degenerate cells. Its minimum is at `(j,t)=(503,5)`:

```text
rho(503,5) = 3.577109931649850292311076... > 4/3.
```

The core starts immediately after the terminal pair. For `5<=t<=501` it is
the interval `[t+2,503]`; for `t>=502` it is empty.

## 3. Fixed columns `t=5,...,8`

| `t` | exact prefix | tail | join record |
| ---: | --- | --- | --- |
| 5 | `[5,1000]` | `j>=1001` | `columns/H2_COL_T5.json` |
| 6 | `[6,4006]` | `j>=4007` | `columns/H2_COL_T6.json` |
| 7 | `[7,12398]` | `j>=12399` | `columns/H2_COL_T7.json` |
| 8 | `[8,19999]` | `j>=20000` | `columns/H2_COL_T8.json` |

Every prefix is an exact integer stream. Every tail is a strict
positive-coefficient rational polynomial at target `4/3`. The prefix and tail
overlap at the first tail cell, so each whole column is closed.

## 4. Columns `t=9,10`

The rows below 504 are terminal cells or finite core cells. For `j>=504`, use
the top region through the anchor and the exact tail above it.

### Column 9

- The level ladder proves `a_68133 < 49 < a_68134`.
- Strict increase gives `a_j<=49=(9-2)^2` for every `504<=j<=68133`, so the
  top theorem applies.
- `columns/H2_TAIL_T9_M68133.json` proves `rho(j,9)>4/3` for every
  `j>=68134`.
- `columns/H2_COL_T9.json` validates the join.

The sharpened-corridor record plus the proved D-LAW gives a second exact
certificate chain across the same column. It is corroborating evidence and
is not needed for the coverage join above.

### Column 10

- The level ladder proves `a_200000<64=(10-2)^2`.
- Strict increase gives the top region for every `504<=j<=200000`.
- `columns/H2_TAIL_T10_M200000.json` proves `rho(j,10)>4/3` for every
  `j>=200001`.
- `columns/H2_COL_T10.json` validates the join.

The sharpened-corridor Sturm record gives a second exact certificate chain
for the `t=10` joint range. It is corroborating evidence and not needed for
coverage. Both second chains share definitions and upstream Paper 3 facts
with the primary proof; they are not claimed as fully independent
implementations.

## 5. Columns `t>=11`

Assume `j>=504` and `t>=11`. Let `a_j` be the increasing Wallis level.

- If `a_j <= (t-2)^2`, the top theorem gives

  ```text
  rho(j,t) > 4100936855929/2058631521408 > 4/3.
  ```

- If `a_j > (t-2)^2`, then

  ```text
  (t-2)^2 - 5t = (t-9)^2 + 9(t-9) + 4 > 0,
  ```

  so the joint hypotheses hold. The `27/20` core and the row-504 attached
  factor give

  ```text
  rho(j,t) >= 254783667121/190614029760 > 4/3.
  ```

  Exact Sturm sequences cover `t=11,...,16`; the 105-monomial bivariate
  certificate covers every `t>=17`.

The top condition includes equality and the joint condition is its strict
complement, so the split is exhaustive and disjoint.

## Adjacency at the only changing boundary

For `t<=501`, the ordered integer intervals are

```text
[t,t+1] | [t+2,503] | [504,infinity).
```

At `t=501`, `t+2=503`. At `t=502`, `t+2=504`, the finite core becomes empty,
and the interior starts immediately after the terminal pair. Thus no integer
cell is duplicated or omitted when the core disappears.

The unclaimed set is empty.
