# Odd cycles, sech spectra, and a square-tail inequality

This is the exact-arithmetic replay accompanying J. Councilman's manuscript
*Odd cycles, sech spectra, and a square-tail inequality*.
Archived at doi:10.5281/zenodo.21853982.

For the odd-cycle probability polynomial

```text
Pi_r(u) = E[u^((o(pi)-1)/2)],    pi uniform in S_(2r-1),
```

write `L(f)_t=f_t^2-f_(t-1)f_(t+1)` and

```text
S_(r,t) = L((1+u)Pi_r)_t,
V_(r,t) = L(Pi_r)_(t-1),
d_(r,t) = S_(r,t) - ((2r-1)/(2r))^2 S_(r-1,t).
```

The paper proves

```text
d_(r,t)^2 >= (2/r^2) S_(r,t)V_(r,t),
r >= 2,  2 <= t <= r.
```

and, as a second theorem, that the factor 2 is not extremal: the sharp uniform
constant is

```text
C* = 2 rho(129,2) = 2.029152329129741484...,
d_(r,t)^2 >= (C*/r^2) S_(r,t)V_(r,t),
```

with equality only at `(r,t)=(129,2)`. The second statement implies the first;
both are proved by separate certificate chains and replayed independently.

It also proves `d_(r,t)>0` on this range. The lower boundary is exact:
`(r,t)=(2,1)` fails. The two upper cells are elementary:
`rho(r,r)=r^2/2` for `r>=2`, and `rho(r,r-1)>2r^2/9` for `r>=3`.

## Verification

Run the complete release replay from this directory:

```bash
python3.14t -I -B verify_all.py
```

The runner is fail-closed. It checks the sorted release manifest, rejects
optimized Python, runs every exact gate, regenerates each banked certificate in
a temporary directory, and requires byte-for-byte agreement with the release
record. Missing or extra files, dependency drift, exceptions, timeouts, failed
gates, surviving named mutations, or certificate differences produce a
nonzero exit.

The replay reconstructs:

- the coefficient-profile curvature and connection-ratio signs;
- the first-lowering and uniform-decrement bounds;
- the linear-extremal and joint-window certificates;
- seven fixed-column tail polynomials, with 329 rational coefficients total;
- 125,250 exact cells in the full-row core and 26,888 exact cells in the four
  fixed-column gaps;
- the two terminal cells and the remaining decrement orientation;
- a finite definition-level transcription check of the shift-wall bridge;
- the four fixed columns retargeted to the rational barrier `203/200`, the
  exact `a_r` level crossings, and 19,146 exact comparisons in the four
  bounded windows; and
- all 125,250 core cells compared exactly against `rho(129,2)`, establishing
  that exactly one of them attains it.

The seven separate tail certificates run concurrently, one single-worker
job per column. The finite core then uses one outer `ProcessPoolExecutor` with
the runtime CPU count. The fixed-column staircase is sequential because its row
recurrence is sequential. Both finite generators write and flush resumable
JSONL records as work is completed. The release runner creates those work files
under its own temporary directory.

Every sign decision is made with exact integer, rational, interval, symbolic
polynomial, or Sturm arithmetic. Reporting decimals do not determine a gate.

## What this packet is

The release is classified as a `REPRODUCIBLE_MANUSCRIPT_REPLAY`. The paper is
the proof. The code checks its long algebraic identities, coefficient signs,
finite residues, and named transcription boundaries. It is not a proof
assistant development, and a clean run does not convert an ordinary analytic
argument in the manuscript into kernel-checked code.

The separate Lean development is intentionally absent. It formalizes related
structure and remains useful work, but it is not a dependency of this theorem
or this release. The exact boundary is recorded in
[VERIFICATION_SCOPE.md](VERIFICATION_SCOPE.md); the remaining open questions
are stated in [OPEN.md](OPEN.md).

## Contents

- `paper/`: manuscript source and PDF;
- `verification/`: exact replay drivers;
- `results/`: compact banked certificates;
- `verify_all.py`: authoritative release runner; and
- `SHA256SUMS`: sorted manifest for the complete release tree.

The older two-million-row interval scan is not included. It remains useful as
a regression archive. The factor-2 proof consumes 152,138 finite cells; the
sharp-constant proof reuses the 125,250-cell core and adds the 19,146 window
comparisons.

`core_uniqueness.json` is the authoritative uniqueness certificate: it compares
every core cell with `rho(129,2)` by exact cross-multiplication of the cancelled
integer form. `dyadic_screen.json` reaches the same 125,250 verdicts through a
directed-rounding dyadic interval screen with exact fallback; it is an
accelerator, its entire sign stream is gated cell-by-cell against the exact
certificate, and its only permitted fallback is the equality at `(129,2)`.

## Citation

J. Councilman, *Odd cycles, sech spectra, and a square-tail inequality*, draft,
2026. Machine-readable metadata is in `CITATION.cff`. A DOI will be added only
after the release record exists.

## License

Verification code is MIT licensed. The manuscript, documentation, and
certificate records are CC BY 4.0.
