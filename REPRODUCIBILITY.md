# Reproducibility contract

## Frozen inputs

The release is self-contained. It vendors the exact public Paper 3 replay at
tag `v0.1.0` under `vendor/odd-cycles-square-tail-replay/`, including that
release's own sorted manifest and full runner. Every Paper 3 source consumed
by a Paper 4 producer is read from that vendored tree; nothing is read from
outside the release.

During packaging, absolute development paths were replaced by release-relative
paths. Reporting-only elapsed-time fields were removed from deterministic
certificate records. Mathematical inputs, exact rational fields, coefficient
streams, gates, and verdicts were not changed.

## Command

From the release root:

```bash
python3.14t -I -B verify_all.py \
  --work-dir /absolute/path/to/shift-wall-replay-work
```

The work directory is persistent and must lie outside the checkout. On its
first run, the runner validates `SHA256SUMS`, scans the text inventory for
nonportable paths, copies the release to `WORK/tree`, and removes every Paper
4 banked output from that copy. All producers then run against the copy. On a
later invocation, completed stages and tail attempts are validated by digest
and skipped; partial prefix and core JSONL files are resumed by their own
producers.

## Concurrency

- The finite core uses one `ProcessPoolExecutor(max_workers=6)`.
- Each fixed-prefix phase uses one outer process pool with at most six workers.
- Eleven fixed-column tail attempts are submitted to one six-entry scheduler.
  Each attempt is a single-process producer, so no process pool is nested.
- The vendored Paper 3 replay finishes before Paper 4 begins and follows its
  own released CPU-count contract.

The tail batch submits the largest anchors first. Every completed attempt,
including expected failures, is appended immediately to
`WORK/tail_progress.jsonl` with a digest and exact exit/status contract. One
failed attempt does not terminate the remaining jobs.

## Expected tail contracts

| Column | Anchor `m0` | Expected status |
| ---: | ---: | --- |
| 4 | 1000 | PASS, independent crosscheck |
| 5 | 1000 | PASS |
| 6 | 2003 | FAIL, undersized control |
| 6 | 4006 | PASS |
| 7 | 6199 | FAIL, undersized control |
| 7 | 12398 | PASS |
| 8 | 19999 | PASS |
| 9 | 55000 | FAIL, undersized control |
| 9 | 68133 | PASS |
| 10 | 100000 | FAIL, undersized control |
| 10 | 200000 | PASS |

An expected failure is evidence only when the producer exits nonzero and the
record says `status: FAIL`. A missing record, an unexpected pass, or an
unexpected failure stops the release replay.

## Acceptance criteria

The runner accepts only if:

1. the release manifest and inventory are exact and portable;
2. the complete vendored Paper 3 replay passes;
3. every Paper 4 producer completes under its stated contract;
4. all 34 deterministic Paper 4 outputs reproduce byte-for-byte;
5. each of the 497 finite-core row records matches after removing only its
   reporting-only `secs` field;
6. the master theorem join reports 22/22 PASS with no conditional leg;
7. the release-anchor verifier reports 40/40 PASS.

The final work record is `WORK/REPLAY_REPORT.json`. A failed or interrupted run
does not modify the release checkout and does not authorize publication.

## PDF build

The reader PDF can be rebuilt separately:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  shift_wall_sharp_factor.tex
```

PDF byte identity is not an arithmetic gate because TeX installations can
embed different metadata. The archived PDF is visually inspected page by
page and is covered by the release manifest.
