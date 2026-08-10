# Open questions

The square-tail inequality is closed on the natural range

```text
r >= 2,    2 <= t <= r.
```

The adjacent stripe `t=1` is not open: it is false, and it is characterized
exactly on `4<=r<=3000` — failure set exactly `[6,41]`, minimum at `r=13`,
strictly increasing recovery from `r=42` (`A45_t1_stripe.py`). The upper cells
in the natural range are not open either: they have explicit formulas and
ratios at least 2.

## The uniform constant is now determined

Earlier releases listed the exact uniform constant as the main open question.
It is closed. With

```text
rho(r,t) = r^2 d_(r,t)^2 / (2 S_(r,t)V_(r,t)),
```

the minimum of `rho` over the whole range is attained, at one cell only:

```text
min rho = rho(129,2) = 1.014576164564870742...,   attained only at (129,2).
```

Equivalently the sharp uniform constant replacing the factor 2 is

```text
C* = 2 rho(129,2) = 2.029152329129741484617172393669...,
```

an explicit rational whose reduced numerator and denominator have 1,233 digits
each. The runner-up is the neighbouring cell `(128,2)`, separated from the
minimum by `1.5878e-6`. The replayed chain is `global_infimum_tails.json`,
`global_infimum_reduction.json`, `global_infimum_windows.json`,
`core_uniqueness.json`, `dyadic_screen.json` and `global_infimum_join.json`.

The factor 2 is retained as the structural theorem. Its proof is the one that
carries the ideas, and the slack it leaves is uniform and explicit: `rho >= 1`
has strict margin `1.4576%` everywhere except at `(129,2)`.

## What is still open

**The shift-wall constant.** `C*` propagates through the shift-wall bridge:
the published `(2r)^2 L(G)_q >= 2 L(C)_q` holds with `C*` in place of 2 on
`0 <= q <= r-2`, by the same induction with the rescaled radical, and no new
certificate is needed. `C*` is *not* sharp there. The bridge's own extremal
cell is not the image of `(129,2)`, and its exact optimal constant is unknown.
Determining it is a separate extremal problem.

**Shape of the column minima away from the core.** Each fixed column `t` has an
interior minimum in `r`, and on the certified core `4 <= r <= 503` the 500
column minima increase strictly in `t`. Whether that monotonicity persists for
every `t` is not proved. Ordinary floating-point exploration suggests the
argmins grow geometrically, but nothing here bounds its error, and no such scan
is offered as evidence.

**Formalization.** Formalizing the manuscript in Lean is unfinished. That is an
implementation frontier rather than an uncertainty in the theorem proved here.
