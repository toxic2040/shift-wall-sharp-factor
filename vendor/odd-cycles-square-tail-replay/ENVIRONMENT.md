# Reference environment

The release replay was frozen on:

```text
Linux x86_64
CPython 3.14.5 free-threading build
SymPy 1.14.0
16 logical CPUs
```

Install the sole third-party replay dependency with

```bash
python3.14t -m pip install -r requirements.txt
```

Then run

```bash
python3.14t -I -B verify_all.py
```

The runner refuses optimized mode because several symbolic components use
Python assertions as exact gates. It does not use the network. Each component
receives a fixed timeout, and all generated certificates and resumable JSONL
work files are placed in a temporary directory.

The seven separate fixed-column tails run concurrently with at most seven
single-worker jobs. The finite core then uses one outer process pool with
`os.cpu_count()` workers. No component starts a second process pool inside
those workers. The fixed-column prefix recurrence runs sequentially;
parallelizing its dependent rows would change the algorithm rather than speed
it up.

The included PDF was built with `latexmk` and pdfTeX from TeX Live 2023. LaTeX
is not required for the arithmetic replay, and the runner does not demand a
byte-identical PDF rebuild across TeX installations.
