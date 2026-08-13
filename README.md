# The sharp constant in the shift-wall bridge

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21866366.svg)](https://doi.org/10.5281/zenodo.21866366)

The paper *The sharp constant in the shift-wall bridge for odd-cycle
polynomials* determines the best possible constant in the bridge between
consecutive shift-wall rows. This release carries that manuscript together
with the exact-arithmetic replay that reproduces every decision-bearing number
in it.

For the odd-cycle recurrence

```text
M_0(z)=1,
M_1(z)=1+z,
M_n(z)=M_(n-1)(z)+n^2 z M_(n-2)(z),
```

write `L(f)_q=f_q^2-f_(q-1)f_(q+1)` and

```text
Omega_t(r) =
    (2r)^2 L((1+z)M_(2r-2))_(r-t)
    --------------------------------
          2 L(M_(2r-1))_(r-t)
```

The exact computer-assisted theorem is

```text
Omega_t(r) >= Omega_2(1350)  for all integers r>=t>=2,
```

with equality only at `(r,t)=(1350,2)`. Hence the sharp raw bridge factor is

```text
2 Omega_2(1350) = 2.655018313913361199726406175676...
```

The value is stored and compared as a reduced rational. The SHA-256 digest of
its canonical `numerator/denominator` encoding is

```text
2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691
```

The domain matters. In the raw coefficient indexing the bridge used here is
proved for `0<=q<=r-2`, equivalently `2<=t<=r`. The omitted `t=1` edge is not
part of the theorem; for example, `Omega_1(6)=0.9824439139...<1`.

## Read the result

- [Paper PDF](paper/shift_wall_sharp_factor.pdf)
- [Verification scope](docs/VERIFICATION_SCOPE.md)
- [Complete coverage ledger](docs/COVERAGE.md)
- [Certificate source map](docs/SOURCE_MAP.md)
- [Relation to the three preceding papers](docs/PROGRAM_MAP.md)
- [Open problems](docs/OPEN_PROBLEMS.md)
- [Fixed-column tail certificate, t=2](docs/FIXED_COLUMN_TAIL_T2.md)

The main join has 22 exact gates, no conditional legs, and no failures. The
finite core contains 123,753 exact cells. The D-LAW companion certificate has
57 contiguous passing windows and a rejected `4/5` negative control. D-LAW is
a second chain for the ninth stripe; the primary global coverage does not
depend on it.

## Reproduce the release

Every number the paper quotes is rebuilt from source by the replay below; the
banked records are compared, not trusted.

Use CPython 3.14.5t with free threading enabled and SymPy 1.14.0. The work
directory must be outside this checkout because the runner verifies the
release inventory before and after copying it.

```bash
python3.14t -m pip install -r requirements.txt
python3.14t -I -B verify_all.py \
  --work-dir /absolute/path/to/shift-wall-replay-work
```

The Paper 4 finite-core and prefix sweeps use one six-process outer pool. The
eleven fixed-column tail attempts use six concurrent single-process jobs; four
are required to fail at their deliberately undersized anchors. Each completed
tail attempt is appended immediately to a resumable JSONL record. The vendored
Paper 3 replay runs first under its own pinned contract. No replay step uses
the network.

Success requires all of the following:

- the complete vendored Paper 3 replay passes;
- all four lower-stripe theorem producers and the exact row-45082
  monotonicity-witness producer pass;
- all 123,753 finite-core cells and four fixed-column prefixes replay;
- the seven passing and four failing tail attempts match their contracts;
- the two D-LAW legs, 57 windows, and negative control match;
- the master join reports 22/22 PASS;
- 35 deterministic records reproduce byte-for-byte, the 497 core rows
  reproduce modulo reporting-only elapsed time, and all 43 release-anchor
  gates pass.

The run is resumable from the same `--work-dir`. See
[REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the full execution and comparison
contract.

## Contents

```text
paper/          manuscript source and reader PDF
docs/           claim boundaries, coverage, program map, and source map
verification/   Paper 4 producers and banked exact records
vendor/         the exact Paper 3 v0.1.0 replay consumed by the transfer step
verify_all.py   isolated full producer replay
SHA256SUMS      sorted release inventory and digests
```

One warning about `vendor/`. `vendor/odd-cycles-square-tail-replay/` is a
complete second deposit — the Paper 3 v0.1.0 release, carried verbatim so the
transfer step replays with nothing fetched. It brings its own `README.md`,
`OPEN.md`, `CITATION.cff`, `SHA256SUMS`, and license files, and every one of
them describes Paper 3, not this work. Read the top-level files here for
anything about the sharp constant.

## Claim boundary

This release supports an exact computer-assisted theorem. It is not described
as proof-assistant verified or kernel checked. Decimals are reader renderings;
all decisive comparisons use integers, rational arithmetic, coefficient
signs, or exact Sturm sequences. The speculative limiting decrement profile
in `docs/OPEN_PROBLEMS.md` is explicitly not part of the theorem.

## Citation and licenses

Cite the paper. Its DOI is
[doi:10.5281/zenodo.21866366](https://doi.org/10.5281/zenodo.21866366), and
machine-readable citation metadata are in [CITATION.cff](CITATION.cff).

The manuscript, documentation, and certificate records are CC BY 4.0; see
[LICENSE-CC-BY-4.0](LICENSE-CC-BY-4.0). Verification code is MIT licensed; see
[LICENSE](LICENSE). The vendored Paper 3 release retains the license notices
in its own directory.
