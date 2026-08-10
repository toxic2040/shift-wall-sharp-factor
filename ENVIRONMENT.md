# Reference environment

The release replay was frozen on:

```text
Linux x86_64
CPython 3.14.5 free-threading build
SymPy 1.14.0
16 logical CPUs
```

Install the sole third-party Paper 4 dependency with

```bash
python3.14t -m pip install -r requirements.txt
```

Then run

```bash
python3.14t -I -B verify_all.py \
  --work-dir /absolute/path/to/shift-wall-replay-work
```

The runner refuses optimized mode, a GIL-enabled interpreter, the wrong
Python minor version, or a different SymPy version. Several symbolic
components use Python assertions as exact gates. The runner does not use the
network and writes only to the explicit external work directory.

Paper 4's finite-core and fixed-prefix producers use one outer pool of six
processes. The tail batch uses one six-entry thread scheduler to supervise six
single-process arithmetic jobs; the jobs do not contain inner pools. Results
are appended to `tail_progress.jsonl` with `flush()` and `fsync()` as each job
finishes. The vendored Paper 3 replay runs sequentially before Paper 4 and may
use the runtime CPU count under its own no-nesting contract.

The included PDF was built with `latexmk` and pdfTeX. LaTeX is not required
for the arithmetic replay, and the runner does not demand a byte-identical PDF
rebuild across TeX installations.
