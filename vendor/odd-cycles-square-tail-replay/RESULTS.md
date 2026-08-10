# Certificate index

The compact files under `results/` are deterministic summaries. Resumable
exact work records are regenerated in a temporary directory and are not
shipped.

| Certificate | Exact role |
|---|---|
| `fixed_tail_t2.json` through `fixed_tail_t8.json` | Seven fixed-column half-lines, anchored at `M_0=312,115,249,699,2003,6199,19999`; polynomial degrees `22,30,38,46,54,62,70`. |
| `finite_core.json` | Every interior cell on `4<=r<=503`: 125,250 exact comparisons in 20 blocks. |
| `fixed_prefixes.json` | The four residual column gaps: 26,888 exact comparisons in one staircase. |
| `terminal_cells.json` | Symbolic proofs of `rho(r,r)=r^2/2` and `rho(r,r-1)>2r^2/9`, plus a 21-cell transcription replay. |
| `t2_orientation.json` | The 309 finite signs `d_(r,2)>0` on `4<=r<=312`; the fixed-tail certificate supplies `r>=313`. |
| `symbolic_join.json` | Binds the seven tails, finite core, fixed prefixes, terminal checks, and the `t>=9` top/joint dichotomy into the complete square-tail coverage statement. |
| `global_infimum_tails.json` | The four fixed columns `t=2,3,4,5` retargeted from `rho>=1` to `rho>=203/200`; anchors `M_0=503,503,503,1000`, degrees `22,30,38,46`. Records the rejection of the superseded `t=2` anchor `M_0=312`. |
| `global_infimum_reduction.json` | The row-free core constant `513/500`, the row-corrected joint floor at the worst row `r=504`, the top floor, and the five exact `a_r` level crossings that bound the `t=6,7,8` windows. |
| `global_infimum_windows.json` | The four residual bounded windows: 19,146 exact comparisons against `203/200` on 16,483 rows. |
| `core_uniqueness.json` | All 125,250 core cells compared exactly with `rho(129,2)`. Exactly one lies at or below it. Authoritative. |
| `dyadic_screen.json` | The same 125,250 verdicts via a directed-rounding dyadic interval screen with exact fallback. Accelerator only; its sign stream is gated cell-by-cell against `core_uniqueness.json`. |
| `global_infimum_join.json` | Binds the six certificates above into `min rho = rho(129,2)`, attained only at `(129,2)`. |

The gaps below the tail anchors in columns `t=2,3,4` are absorbed by the finite
core. Only columns `t=5,6,7,8` contribute separate prefix cells. The sum

```text
125,250 + 26,888 = 152,138
```

is the complete finite SQ residue consumed by the factor-2 proof. The
sharp-constant proof reuses the 125,250-cell core and adds the 19,146 window
comparisons, for

```text
125,250 + 19,146 = 144,396
```

exact comparisons in its own chain.

Each certificate records `passed=true`. `verify_all.py` validates its schema
and load-bearing counts, regenerates it from the corresponding driver, and
compares the bytes with the banked file. `SHA256SUMS` binds the certificate to
the rest of the release tree.
