#!/usr/bin/env python3
"""Finite-core sweep of the 4/3 barrier on the t>=5 bulk.

Claim proved by this run, exactly:

    rho(r,t) > 4/3  for every  7 <= r <= 503  and  5 <= t <= r-2.

Definition-level exact-integer route: rho is rebuilt here from the wall
recurrence and the L-operator alone, with no reuse of the certificate
machinery that proves the half-line legs.  The (2r-1)! normalisations
cancel identically, so with

    Lg_r = L((1+u)F_r)_t,   Lv_r = L(F_r)_{t-1},
    M    = r^2 Lg_r - (2r-1)^4 (r-1)^2 Lg_{r-1},

    rho(r,t) = M^2 / (2 r^2 Lg_r Lv_r)

and the 4/3 test is the exact integer sign of  3*M^2 - 4*(2 r^2 Lg Lv).
No floats anywhere; the only rationals are the final reported minima.

Plan: iterate rows r = 7..503; per row compute the wall row once and
evaluate ALL columns t = 5..r-2 from it.  Rows are grouped into contiguous
chunks; each chunk streams the Wall recurrence from n=0 with a rolling
window (rows never travel as task arguments).  One ProcessPoolExecutor,
6 workers, no nesting.  Every finished row is written immediately to
A_core_sweep.jsonl with flush(); a rerun skips completed rows.  Every row
body is wrapped in try/except so one failure cannot kill the batch.

A cell with 3*M^2 - 4*den <= 0 is a KILL finding and is recorded with its
exact value; the summary gate then reads FAILED.
"""
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as Frac

HERE = os.path.dirname(os.path.abspath(__file__))
JSONL = os.path.join(HERE, "A_core_sweep.jsonl")
SUMMARY = os.path.join(HERE, "A_CORE_SWEEP_SUMMARY.json")

R_LO, R_HI = 7, 503
T_LO = 5
CHUNK = 8
WORKERS = 6
STR_DIGITS = 1_000_000


# ------------------------------------------------- definitions (from A24)
def wall_stream(n_max):
    """Yield (n, M_n coefficient list) with a rolling window."""
    prev2, prev = None, None
    for n in range(0, n_max + 1):
        if n == 0:
            row = [1]
        elif n == 1:
            row = [1, 1]
        else:
            row = list(prev)
            for i, c in enumerate(prev2):
                while len(row) <= i + 1:
                    row.append(0)
                row[i + 1] += n * n * c
        yield n, row
        prev2, prev = prev, row


def F_from_row(row, r):
    """[u^k] F_r = m_{2r-2, r-1-k}."""
    return [row[r - 1 - k] if 0 <= r - 1 - k < len(row) else 0 for k in range(r)]


def shift_add(f):
    """(1+u) * f"""
    out = [0] * (len(f) + 1)
    for i, c in enumerate(f):
        out[i] += c
        out[i + 1] += c
    return out


def Lop(f, t):
    g = lambda i: f[i] if 0 <= i < len(f) else 0
    return g(t) * g(t) - g(t - 1) * g(t + 1)


# ------------------------------------------------------------- workers
def _init_worker():
    sys.set_int_max_str_digits(STR_DIGITS)


def one_row(Fcur, Fprev, r):
    """All columns t = 5..r-2 for one row.  Exact integers throughout."""
    t0 = time.time()
    Gcur = shift_add(Fcur)
    Gprev = shift_add(Fprev)
    shift_const = (2 * r - 1) ** 4 * (r - 1) ** 2
    cmps = []
    fails = []
    degenerate = []
    min_t = None
    min_num = min_den = None
    for t in range(T_LO, r - 1):
        Lg = Lop(Gcur, t)
        Lv = Lop(Fcur, t - 1)
        Lg1 = Lop(Gprev, t)
        if Lg <= 0 or Lv <= 0:
            degenerate.append({"r": r, "t": t, "Lg_sign": (Lg > 0) - (Lg < 0),
                               "Lv_sign": (Lv > 0) - (Lv < 0)})
            cmps.append("!")
            continue
        M = r * r * Lg - shift_const * Lg1
        num = M * M
        den = 2 * r * r * Lg * Lv
        diff = 3 * num - 4 * den
        cmp_char = ">" if diff > 0 else "=" if diff == 0 else "<"
        cmps.append(cmp_char)
        if diff <= 0:
            fails.append({"r": r, "t": t, "cmp": cmp_char,
                          "num": str(num), "den": str(den)})
        if min_num is None or num * min_den < min_num * den:
            min_t, min_num, min_den = t, num, den
    return {"r": r, "t_lo": T_LO, "t_hi": r - 2, "cells": r - 2 - T_LO + 1,
            "cmps": "".join(cmps), "fails": fails, "degenerate": degenerate,
            "min_t": min_t, "min_num": str(min_num), "min_den": str(min_den),
            "secs": round(time.time() - t0, 3)}


def run_chunk(bounds):
    """Contiguous r-range [r_lo, r_hi]; streams the Wall from n=0.

    Returns a list of per-row records.  Each row body is individually
    wrapped so one row's failure cannot kill the chunk.
    """
    r_lo, r_hi = bounds
    out = []
    Fprev = None
    for n, row in wall_stream(2 * r_hi - 2):
        if n % 2 or n < 2 * (r_lo - 1) - 2:
            continue
        r_here = (n + 2) // 2
        Fcur = F_from_row(row, r_here)
        if r_lo <= r_here <= r_hi and Fprev is not None:
            try:
                out.append(one_row(Fcur, Fprev, r_here))
            except Exception as exc:            # noqa: BLE001
                out.append({"r": r_here, "error": repr(exc)})
        Fprev = Fcur
    return out


# --------------------------------------------------------------- driver
def load_done():
    done = {}
    if os.path.exists(JSONL):
        with open(JSONL) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                if "error" not in rec:
                    done[rec["r"]] = rec       # last good record wins
    return done


def contiguous_chunks(rs, size):
    chunks = []
    run = []
    for r in rs:
        if run and r != run[-1] + 1:
            chunks.append(run)
            run = []
        run.append(r)
        if len(run) == size:
            chunks.append(run)
            run = []
    if run:
        chunks.append(run)
    return [(c[0], c[-1]) for c in chunks]


def main():
    sys.set_int_max_str_digits(STR_DIGITS)
    done = load_done()
    todo = [r for r in range(R_LO, R_HI + 1) if r not in done]
    chunks = contiguous_chunks(todo, CHUNK)
    print("A core sweep: rows %d..%d, %d already done, %d to do in %d chunks "
          "on %d workers (one outer pool, no nesting)"
          % (R_LO, R_HI, len(done), len(todo), len(chunks), WORKERS), flush=True)

    t0 = time.time()
    if chunks:
        sink = open(JSONL, "a", buffering=1)
        with ProcessPoolExecutor(max_workers=WORKERS,
                                 initializer=_init_worker) as ex:
            futs = {ex.submit(run_chunk, c): c for c in chunks}
            k = 0
            for fut in as_completed(futs):
                bounds = futs[fut]
                try:
                    records = fut.result()
                except Exception as exc:        # one chunk cannot kill the batch
                    records = [
                        {"r": r, "error": repr(exc)}
                        for r in range(bounds[0], bounds[1] + 1)
                    ]
                for rec in records:
                    sink.write(json.dumps(rec) + "\n")
                    sink.flush()
                    if "error" not in rec:
                        done[rec["r"]] = rec
                k += 1
                if k % 8 == 0 or k == len(chunks):
                    print("  ... %d/%d chunks, %d rows recorded, %.0fs"
                          % (k, len(chunks), len(done), time.time() - t0),
                          flush=True)
        sink.close()

    # ---------------------------------------------------- consolidation
    missing = [r for r in range(R_LO, R_HI + 1) if r not in done]
    if missing:
        print("INCOMPLETE: %d rows missing (first: %s) - rerun to resume"
              % (len(missing), missing[:10]))
        return 1

    total_cells = 0
    fails = []
    degenerate = []
    stream = hashlib.sha256()
    gmin = None                        # (num, den, r, t) exact
    for r in range(R_LO, R_HI + 1):
        rec = done[r]
        span = rec["t_hi"] - rec["t_lo"] + 1
        assert rec["t_lo"] == T_LO and rec["t_hi"] == r - 2
        assert len(rec["cmps"]) == span == rec["cells"]
        total_cells += span
        fails.extend(rec["fails"])
        degenerate.extend(rec["degenerate"])
        for t, c in zip(range(T_LO, r - 1), rec["cmps"]):
            stream.update(f"{r}:{t}:{c}\n".encode("ascii"))
        num, den = int(rec["min_num"]), int(rec["min_den"])
        if gmin is None or num * gmin[1] < gmin[0] * den:
            gmin = (num, den, r, rec["min_t"])

    num, den, gr, gt = gmin
    value = Frac(num, den)
    literal = f"{value.numerator}/{value.denominator}"
    margin = value - Frac(4, 3)

    def truncate(v, places=24):
        sign = "-" if v < 0 else ""
        whole, rem = divmod(abs(v.numerator), v.denominator)
        digits = []
        for _ in range(places):
            rem *= 10
            d, rem = divmod(rem, v.denominator)
            digits.append(str(d))
        return f"{sign}{whole}." + "".join(digits)

    ok = not fails and not degenerate and total_cells == sum(
        r - 6 for r in range(R_LO, R_HI + 1))
    summary = {
        "schema": "h3-core-sweep-four-thirds-v1",
        "claim": "rho(r,t) > 4/3 for every 7 <= r <= 503 and 5 <= t <= r-2",
        "route": "definition-level exact integers, (2r-1)! cancelled; "
                 "sign test = 3*M^2 - 4*(2 r^2 Lg Lv)",
        "status": "PASS" if ok else "FAILED",
        "rows": R_HI - R_LO + 1,
        "total_cells": total_cells,
        "expected_cells": sum(r - 6 for r in range(R_LO, R_HI + 1)),
        "failures": fails,
        "degenerate": degenerate,
        "comparison_stream_encoding": "ASCII r:t:cmp newline, r=7..503 "
                                      "ascending, t=5..r-2 ascending",
        "comparison_stream_sha256": stream.hexdigest(),
        "global_minimum": {
            "r": gr, "t": gt,
            "sha256_of_num_slash_den": hashlib.sha256(
                literal.encode("ascii")).hexdigest(),
            "numerator_digits": len(str(num)),
            "denominator_digits": len(str(den)),
            "decimal_truncated_24": truncate(value),
            "margin_over_four_thirds_decimal_truncated_24": truncate(margin),
        },
        "execution": {
            "interpreter": sys.version,
            "workers": WORKERS,
            "chunk_rows": CHUNK,
        },
    }
    with open(SUMMARY, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    print(json.dumps({k: summary[k] for k in
                      ("status", "total_cells", "comparison_stream_sha256")},
                     indent=2))
    print("global minimum: r=%d t=%d  rho=%s  (margin over 4/3: %s)"
          % (gr, gt, summary["global_minimum"]["decimal_truncated_24"],
             summary["global_minimum"]
             ["margin_over_four_thirds_decimal_truncated_24"]))
    if not ok:
        print("RESULT: SWEEP FAILED - see failures/degenerate in summary")
        return 1
    print("RESULT: all %d cells clear 4/3 strictly" % total_cells)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
