# Verification scope

The manuscript proves the square-tail inequality on its natural coefficient
range. This package is its deterministic exact-arithmetic replay, classified
as `REPRODUCIBLE_MANUSCRIPT_REPLAY`.

That label is deliberate. A clean run proves that the released source tree,
implemented identities, exact finite certificates, and named mutation controls
agree. It does not mean that every proof in the paper has been translated into
a proof assistant.

## Proof branches

| Branch | Executable coverage | Boundary |
|---|---|---|
| Odd-cycle, Jacobi, determinant, and DPP descriptions | Definition-level recurrence, determinant, normalization, and finite transcription identities are replayed where they feed later gates. | The all-row spectral and combinatorial arguments are proofs in the paper. |
| Coefficient-ratio signs | Connection formulas, determinant orientations, positive connection-ratio polynomials, and exact row checks are reconstructed. | Total positivity of the coefficient triangle and its planar-network proof remain manuscript arguments. |
| Adjacent-row profile curvature | The three- and four-column determinant factors, 36 positive connection coefficients, and exact profile rows are checked. | The passage from total-positive minors to all-row curvature is the displayed manuscript proof. |
| First lowering and decrement bound | Finite exact seeds, recurrence identities, rational interval envelopes, and nonnegative tail coefficients are checked. | Elementary induction and monotonicity surrounding the displayed identities remain paper reasoning. |
| Linear extremality | All 192 Bernstein coefficients, denominator signs, row-ratio monotonicity, and top-comparison coefficients are exact symbolic gates. | The pure-defect log-concavity argument and the reduction to the five-variable polytope are written in the paper. |
| Joint window | Positive polynomial coefficients, the two exceptional Sturm comparisons, and the final boundary inequality are exact. | The coverage implication using those gates is the manuscript argument. |
| Fixed columns | Seven exact rational tail certificates are rebuilt. Their degrees are `22,30,38,46,54,62,70`, with 329 coefficients total. | Each certificate proves only its stated half-line; the finite gaps are separate. |
| Finite residues | The factor-2 chain uses the full core `4<=r<=503`, `2<=t<=r-2` (125,250 exact comparisons) and the fixed-column staircase (26,888). The sharp-constant chain reuses the same core and adds 19,146 window comparisons. | No larger scan is used as theorem evidence. |
| Terminal cells and decrement orientation | The two terminal columns are symbolic. Column `t=2` has 309 exact finite signs followed by an oriented tail polynomial; `t>=3` follows from the paper's scalar comparison. | The profile comparison proving the `t>=3` orientation is an analytic manuscript step. |
| Symbolic join | The runner binds all seven tails, both finite certificates, terminal gates, and `(t-2)^2-5t>0` for every integer `t>=9`. | This is an exact coverage join, not a replacement for the regional proofs it cites. |
| Sharp constant (Theorem, `thm:sharp`) | The retargeted tails for `t=2..5` at target `203/200` (anchors `503,503,503,1000`), the exact `a_r` level crossings at `25,30,35,36,40`, the four bounded windows totalling 19,146 exact comparisons on 16,483 rows, the 125,250 exact core comparisons against `rho(129,2)`, and the join that binds them are all rebuilt. | The chain proves `min rho = rho(129,2)` attained only at `(129,2)`. It implies the factor-2 theorem but does not replace its proof, which is retained and replayed independently. |
| Core uniqueness authority | `A43_core_uniqueness.py` is the authoritative exact certificate: every core cell is compared with `rho(129,2)` by exact cross-multiplication of the cancelled integer form, with no `Fraction` reduction in the comparison. | `A44_dyadic_screen.py` is an accelerator only. Its entire 125,250-cell sign stream is gated against A43's, and its sole permitted exact fallback is the equality at `(129,2)`. A44 is never the authority for a verdict. |
| Shift-wall consequence | Recurrences, normalizers, the radical inequality, induction indices, reversal identities, and negative controls are replayed on a finite transcription range. The quoted lower-stripe failure window is checked separately on `4<=r<=69`, and the stripe is characterized exactly on `4<=r<=3000` (failure set exactly `[6,41]`, strictly increasing recovery). | The proper-position argument and the all-row induction are proofs in the paper; no all-row claim is made for the lower stripe. |

## Fail-closed boundary

`verify_all.py` fails on a missing or extra release file, an unsafe path, a
manifest mismatch, optimized Python, a dependency mismatch, a component
exception or timeout, a malformed or failed certificate, a rebuilt-certificate
difference, or a surviving named mutation.

Mutation controls test their stated failure modes. They do not prove that every
possible source mutation would be caught independently of the manifest, code
review, and rebuilt-certificate comparison.

Two files are shipped as hashed provenance for the reduction certificate and are
not executed by the runner: `verification/A15_top_coverage.py` and
`verification/A15_TOP_COVERAGE.md`. `A40_global_infimum_reduction.py` records
their SHA-256 among its dependencies; `A15_top_coverage.py` refers to a sibling
capsule directory that is not part of this release and will not run here. They
are present so that the reduction's dependency fingerprints can be checked, not
so that they can be re-run.

The Lean development is a separate formalization project and is not consulted
by this runner. Conversely, this replay makes no claim about Lean's kernel or
axiom boundary.
