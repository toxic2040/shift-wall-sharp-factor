# A15 — top end and gap-free integer coverage

**Date:** 2026-08-05  
**Status:** theorem. This closes the top-end and integer-coverage obligations
and, together with A1 and A10--A14, proves `(SQ)` on its full domain.

## 1. Top-end theorem

A10-TOP proves

```text
r psi >= r eps+(1-eps)[2q-q^2/r],
q=(t-2)sqrt(a+57/20)/a.                                  (TOP)
```

> **A15-TOP [theorem].** For every canonical cell with
>
> ```text
> r>=4,       t>=sqrt(a_r)+2,
> ```
>
> `(SQ)` holds.

The threshold gives `q>=1`. Also `L1_lo<L1<=r/(r-1)` and
`t<=r-2`, so `q<r`. The concave quadratic `2q-q^2/r` takes its minimum on
`[1,r]` at an endpoint and is therefore at least `2-1/r`. From A14 `(EPS)`,

```text
r psi > [1-2/(4r-3)](2-1/r).
```

Both factors increase with `r`. At `r=4`, the right side is `77/52`, and

```text
(77/52)^2-2=521/2704>0.
```

Thus `r psi>sqrt(2)`, which is A10's top gate. The same coarsened scalar gate
fails at `r=3`, outside the `(SQ)` domain.

## 2. Exact splice anchors

The normalized ladder has

```text
a_r=6E_1(r-1),
E_1(m)=E_1(m-1)+D_1(m-1)/R_m.
```

All terms are positive, so `a_r` is strictly increasing. A 192-bit
outward-rounded replay of A3's cancellation-free ladder proves the exact
integer comparisons

```text
a_2000000 < 121,       a_2000001 > 49.                   (ANCHOR)
```

The canonical SHA-256 of the two dyadic interval rows is

```text
f28e67c4395fc0e926b5734162a0ba099e15da7b705fc6091deafe8c709b6356.
```

These are point anchors with analytic monotone propagation, not sampled
uniformity claims.

## 3. Frozen finite coverage

`A15_top_coverage.py` validates both upstream manifests before reading their
verdict records:

- 17-entry `sech_instrument_2026-08-04/SHA256SUMS`;
- 21-entry `A5_SHA256SUMS.txt`.

It then checks every JSONL row sequentially: contiguous `r`, the exact `t_hi`,
zero failures and indeterminates, no error record, and the exact pass count.
The finite frontier needed here is

| integer rows | certified columns |
|---|---:|
| `4..600` | full row, 178,503 exact cells |
| `601..1450` | full row, 869,125 exact cells |
| `1451..8000` | `t<=55` |
| `8001..20000` | `t<=30` |
| `20001..60000` | `t<=20` |
| `60001..2000000` | `t<=12` |

A1 independently closes `t=2,3,4,5` on their entire half-lines.

## 4. Integer partition

There are two splice regions.

### `4<=r<=2000000`

The two full-row legs cover `r<=1450`. Above that, any cell not in the finite
packet has `t>=13`. By monotonicity and the left anchor,

```text
a_r<121,       hence sqrt(a_r)+2<13<=t.
```

Every such residual cell is covered by A15-TOP.

### `r>=2000001`

A1 covers `t<=5`. For `t>=6`, split at the exact real threshold
`sqrt(a_r)+2`.

- If `t>=sqrt(a_r)+2`, A15-TOP applies.
- If `t<sqrt(a_r)+2`, the right anchor and monotonicity give `a_r>49`.
  Put `y=sqrt(a_r)>7`. Then

  ```text
  y^2-5y-10 > 7^2-5*7-10 = 4,
  ```

  so `a_r>5(sqrt(a_r)+2)>5t`. Also `r>504`. All hypotheses of
  A14-UNIFORM hold.

The cases are complementary, the finite row intervals are adjacent, and the
splices occur between the consecutive integers `2000000` and `2000001`.
There is no uncovered boundary cell.

> **A15-SQ [theorem].** For every pair of integers
>
> ```text
> r>=4,       2<=t<=r-2,
> ```
>
> the square-tail inequality `(SQ)` holds.

## 5. Replay

```bash
python3 A15_top_coverage.py
```

The anchor computation and all new analytic gates are exact. The large finite
boxes are not recomputed; their immutable files and source are hash-validated,
then every stored row verdict is audited.
