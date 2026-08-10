#!/usr/bin/env python3
"""Exact fixed-column prefix certification of rho(r,t) > 4/3.

Recomputes the exact factorial-cancelled rho stream (verbatim recurrence from
verify_sq_t4_transfer.py, t-general) in r-range chunks under ONE outer
ProcessPoolExecutor.  The A3_cert_col*.jsonl records carry only strip minima
over t in [2, t_hi] (argmin at t=2..3, below 4/3), so they cannot certify a
fixed column at 4/3 and a recompute is required.

run mode:      --plan "6:8013:4,8:20000:12" [--r-start R] [--workers W]
               each item is t:r_max:chunks; rows r in [r_start(default t), r_max]
assemble mode: --assemble "6:8012" (t:m0) -> H2_PREFIX_T{t}.json

Chunks tile with one shared boundary row; the shared row's raw fingerprint
must agree across chunks.  Chunk bodies are [r_lo, r_hi-1]; the top row of the
final chunk (r = m0+1) is the overlap cell with the tail certificate.
Progress is appended to h2_prefix_progress.jsonl (flush+fsync per chunk,
resume by skip).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parents[2]
JSONL = HERE / "h2_prefix_progress.jsonl"
ENGINE_PATH = (RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay"
               / "verification" / "engine.py")
SCHEMA = "h2-fixed-column-four-thirds-prefix-v1"
WORK_EXPONENT = 2.585
ENCLOSURE_BITS = 192


def turan(row: list[int], index: int) -> int:
    def coefficient(j: int) -> int:
        return row[j] if 0 <= j < len(row) else 0

    return coefficient(index) ** 2 - coefficient(index - 1) * coefficient(index + 1)


def shifted(row: list[int]) -> list[int]:
    return [row[k] + (row[k - 1] if k else 0) for k in range(len(row))]


def rho_stream_range(t: int, r_lo: int, r_hi: int):
    """Yield (r, numerator, denominator, decrement) for r in [r_lo, r_hi].

    The recurrence is run from n=2 exactly as in verify_sq_t4_transfer.py's
    rho_stream; only the quadratic-cost window products are gated on r >= r_lo.
    """
    k_max = t + 1
    previous_two = [1] + [0] * k_max
    previous = [1, 1] + [0] * (k_max - 1)
    previous_even = [1] + [0] * k_max
    for n in range(2, 2 * r_hi - 1):
        if n % 2 == 0:
            current = [
                previous[k] + n * n * previous_two[k]
                for k in range(k_max + 1)
            ]
        else:
            current = [
                (previous[k - 1] if k else 0) + n * n * previous_two[k]
                for k in range(k_max + 1)
            ]
        if n % 2 == 0:
            r = n // 2 + 1
            if r >= t and r >= r_lo:
                current_s = turan(shifted(current), t)
                predecessor_s = turan(shifted(previous_even), t)
                current_v = turan(current, t - 1)
                decrement = (
                    r * r * current_s
                    - (2 * r - 1) ** 4 * (r - 1) ** 2 * predecessor_s
                )
                yield (
                    r,
                    decrement * decrement,
                    2 * r * r * current_s * current_v,
                    decrement,
                )
            previous_even = current
        previous_two, previous = previous, current


def enclosure(numerator: int, denominator: int) -> tuple[F, F]:
    """Cheap outward bracket of numerator/denominator using top bits only."""
    n_shift = max(numerator.bit_length() - ENCLOSURE_BITS, 0)
    d_shift = max(denominator.bit_length() - ENCLOSURE_BITS, 0)
    a = numerator >> n_shift
    b = denominator >> d_shift
    scale = F(2) ** (n_shift - d_shift)
    return F(a, b + 1) * scale, F(a + 1, b) * scale


def hex_fingerprint(numerator: int, denominator: int) -> str:
    return hashlib.sha256(
        f"{numerator:x}/{denominator:x}".encode("ascii")
    ).hexdigest()


def decimal_truncation(value: F, places: int) -> str:
    sign = "-" if value < 0 else ""
    whole, remainder = divmod(abs(value.numerator), value.denominator)
    digits: list[str] = []
    for _ in range(places):
        remainder *= 10
        digit, remainder = divmod(remainder, value.denominator)
        digits.append(str(digit))
    return f"{sign}{whole}." + "".join(digits)


def reduced_fingerprint(numerator: int, denominator: int) -> dict:
    value = F(numerator, denominator)
    literal = f"{value.numerator}/{value.denominator}"
    return {
        "sha256": hashlib.sha256(literal.encode("ascii")).hexdigest(),
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
        "decimal_truncated_50": decimal_truncation(value, 50),
        "margin_over_four_thirds_24": decimal_truncation(value - F(4, 3), 24),
    }


def worker_init() -> None:
    sys.set_int_max_str_digits(1_000_000)


def load_engine():
    specification = importlib.util.spec_from_file_location(
        "h2_engine", ENGINE_PATH
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def certify_chunk(task: tuple[int, int, int]) -> dict:
    t, r_lo, r_hi = task
    started = time.time()
    engine = None
    wall = None
    if r_lo <= 60:
        engine = load_engine()
        wall = engine.wall_M(2 * 60)

    comparisons: list[str] = []
    strict_failures: list[int] = []
    nonpositive_denominators: list[int] = []
    decrement_positive = decrement_negative = decrement_zero = 0
    engine_cells = 0
    engine_failures: list[int] = []
    terminal_rows: list[dict] = []
    exact_min_comparisons = 0
    best: dict | None = None
    first_row: dict | None = None
    top_row: dict | None = None

    for r, numerator, denominator, decrement in rho_stream_range(
        t, r_lo, r_hi
    ):
        if denominator <= 0:
            nonpositive_denominators.append(r)
            continue
        if decrement > 0:
            decrement_positive += 1
        elif decrement < 0:
            decrement_negative += 1
        else:
            decrement_zero += 1

        if engine is not None and r <= 60:
            if F(numerator, denominator) == engine.exact_rho(r, t, wall):
                engine_cells += 1
            else:
                engine_failures.append(r)
        if r <= t + 1:
            record = {
                "r": r,
                "value_decimal_6": decimal_truncation(
                    F(numerator, denominator), 6
                ),
            }
            if r == t:
                record["equals_half_r_squared"] = (
                    F(numerator, denominator) == F(t * t, 2)
                )
            terminal_rows.append(record)

        if first_row is None:
            first_row = {"r": r, "raw_sha256": hex_fingerprint(
                numerator, denominator)}
        if r == r_hi:
            value = F(numerator, denominator)
            top_row = {
                "r": r,
                "raw_sha256": hex_fingerprint(numerator, denominator),
                "comparison": ">" if 3 * numerator - 4 * denominator > 0
                else "=" if 3 * numerator - 4 * denominator == 0 else "<",
                **reduced_fingerprint(numerator, denominator),
            }
            break

        difference = 3 * numerator - 4 * denominator
        comparison = ">" if difference > 0 else "=" if difference == 0 else "<"
        comparisons.append(comparison)
        if difference <= 0:
            strict_failures.append(r)

        low, high = enclosure(numerator, denominator)
        if best is None:
            best = {"r": r, "num": numerator, "den": denominator,
                    "lo": low, "hi": high}
        elif low > best["hi"]:
            pass
        elif high < best["lo"]:
            best = {"r": r, "num": numerator, "den": denominator,
                    "lo": low, "hi": high}
        else:
            exact_min_comparisons += 1
            if numerator * best["den"] < best["num"] * denominator:
                best = {"r": r, "num": numerator, "den": denominator,
                        "lo": low, "hi": high}

    if top_row is None or first_row is None or best is None:
        raise RuntimeError(f"chunk ({t},{r_lo},{r_hi}) produced no rows")
    minimum = {"r": best["r"],
               **reduced_fingerprint(best["num"], best["den"])}
    return {
        "t": t,
        "r_lo": r_lo,
        "r_hi": r_hi,
        "rows_body": len(comparisons),
        "comparisons": "".join(comparisons),
        "strict_failures": strict_failures,
        "nonpositive_denominators": nonpositive_denominators,
        "decrement_positive": decrement_positive,
        "decrement_negative": decrement_negative,
        "decrement_zero": decrement_zero,
        "engine_crosscheck_cells": engine_cells,
        "engine_crosscheck_failures": engine_failures,
        "terminal_rows": terminal_rows,
        "exact_min_comparisons": exact_min_comparisons,
        "minimum": minimum,
        "first_row": first_row,
        "top_row": top_row,
        "elapsed_seconds": round(time.time() - started, 1),
    }


def chunk_boundaries(r_start: int, r_max: int, chunks: int) -> list[int]:
    low = float(max(r_start, 2)) ** WORK_EXPONENT
    high = float(r_max) ** WORK_EXPONENT
    cuts = [r_start]
    for index in range(1, chunks):
        cut = round((low + (high - low) * index / chunks)
                    ** (1.0 / WORK_EXPONENT))
        cut = min(max(cut, cuts[-1] + 2), r_max - 2)
        cuts.append(cut)
    cuts.append(r_max)
    return cuts


def read_progress() -> dict[tuple[int, int, int], dict]:
    records: dict[tuple[int, int, int], dict] = {}
    if not JSONL.exists():
        return records
    with JSONL.open("r", encoding="utf-8") as source:
        for line in source:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            records[(record["t"], record["r_lo"], record["r_hi"])] = record
    return records


def run_mode(plan: str, r_start_override: int | None, workers: int) -> int:
    done = read_progress()
    tasks: list[tuple[int, int, int]] = []
    for item in plan.split(","):
        t_text, r_max_text, chunk_text = item.split(":")
        t, r_max, chunks = int(t_text), int(r_max_text), int(chunk_text)
        r_start = r_start_override if r_start_override is not None else t
        cuts = chunk_boundaries(r_start, r_max, chunks)
        for lower, upper in zip(cuts, cuts[1:]):
            if (t, lower, upper) not in done:
                tasks.append((t, lower, upper))
        print(f"t={t}: cuts {cuts}")
    if not tasks:
        print("all requested chunks already recorded")
        return 0
    print(f"submitting {len(tasks)} chunks with {workers} workers")
    failures = 0
    with JSONL.open("a", encoding="utf-8") as sink:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_init
        ) as pool:
            futures = {pool.submit(certify_chunk, task): task
                       for task in tasks}
            for future in as_completed(futures):
                task = futures[future]
                try:
                    record = future.result()
                except Exception as error:  # noqa: BLE001 one failure never kills the batch
                    failures += 1
                    print(f"[ERROR] chunk {task}: {error!r}")
                    continue
                sink.write(json.dumps(record) + "\n")
                sink.flush()
                os.fsync(sink.fileno())
                print(
                    f"[done] t={record['t']} [{record['r_lo']},"
                    f"{record['r_hi']}] rows={record['rows_body']} "
                    f"min@r={record['minimum']['r']} "
                    f"min={record['minimum']['decimal_truncated_50'][:14]} "
                    f"fails={len(record['strict_failures'])} "
                    f"{record['elapsed_seconds']}s"
                )
    print(f"run complete, failures={failures}")
    return 1 if failures else 0


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def assemble_mode(spec: str) -> int:
    t_text, m0_text = spec.split(":")
    t, m0 = int(t_text), int(m0_text)
    records = [record for record in read_progress().values()
               if record["t"] == t]
    records.sort(key=lambda record: record["r_lo"])
    if not records:
        print(f"no chunks recorded for t={t}")
        return 1

    problems: list[str] = []
    if records[0]["r_lo"] != t:
        problems.append(f"coverage starts at {records[0]['r_lo']} not {t}")
    for before, after in zip(records, records[1:]):
        if before["r_hi"] != after["r_lo"]:
            problems.append(
                f"gap or overlap between chunk ending {before['r_hi']} "
                f"and chunk starting {after['r_lo']}"
            )
        elif before["top_row"]["raw_sha256"] != after["first_row"]["raw_sha256"]:
            problems.append(
                f"boundary row {before['r_hi']} fingerprint mismatch"
            )
    if records[-1]["r_hi"] != m0 + 1:
        problems.append(
            f"coverage ends at {records[-1]['r_hi']} not m0+1={m0 + 1}"
        )

    digest = hashlib.sha256()
    cells = 0
    strict_failures: list[int] = []
    denominator_bad: list[int] = []
    decrement = {"positive": 0, "negative": 0, "zero": 0}
    engine_cells = 0
    engine_failures: list[int] = []
    terminal_rows: list[dict] = []
    for record in records:
        expected_rows = record["r_hi"] - record["r_lo"]
        if record["rows_body"] != expected_rows:
            problems.append(
                f"chunk [{record['r_lo']},{record['r_hi']}] has "
                f"{record['rows_body']} body rows, expected {expected_rows}"
            )
        for offset, comparison in enumerate(record["comparisons"]):
            digest.update(
                f"{record['r_lo'] + offset}:{comparison}\n".encode("ascii")
            )
        cells += record["rows_body"]
        strict_failures.extend(record["strict_failures"])
        denominator_bad.extend(record["nonpositive_denominators"])
        decrement["positive"] += record["decrement_positive"]
        decrement["negative"] += record["decrement_negative"]
        decrement["zero"] += record["decrement_zero"]
        engine_cells += record["engine_crosscheck_cells"]
        engine_failures.extend(record["engine_crosscheck_failures"])
        terminal_rows.extend(record["terminal_rows"])

    minima = [record["minimum"] for record in records]
    minima.sort(key=lambda item: item["decimal_truncated_50"])
    if (len(minima) > 1 and minima[0]["decimal_truncated_50"]
            == minima[1]["decimal_truncated_50"]):
        problems.append("ambiguous global minimum at 50 decimal digits")
    minimum = minima[0]
    overlap = records[-1]["top_row"]

    if strict_failures:
        problems.append(f"strict failures at rows {strict_failures[:8]}")
    if denominator_bad:
        problems.append(f"nonpositive denominators at {denominator_bad[:8]}")
    if engine_failures:
        problems.append(f"engine crosscheck failed at {engine_failures[:8]}")
    if cells != m0 - t + 1:
        problems.append(f"cell count {cells} != {m0 - t + 1}")

    passed = not problems
    payload = {
        "schema": SCHEMA,
        "t": t,
        "m0": m0,
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "problems": problems,
        "claim": f"rho(r,{t})>4/3 for every integer r in [{t},{m0}]",
        "range": [t, m0],
        "cell_count": cells,
        "comparison_stream_encoding": "ASCII r:<|=|> followed by newline",
        "comparison_sha256": digest.hexdigest(),
        "strict_failures": strict_failures,
        "decrement_signs": decrement,
        "definition_crosscheck": {
            "engine": ENGINE_PATH.relative_to(RELEASE_ROOT).as_posix(),
            "range": [t, 60],
            "cells": engine_cells,
            "failures": engine_failures,
        },
        "terminal_rows": terminal_rows,
        "finite_minimum_r": minimum["r"],
        "finite_minimum": {k: v for k, v in minimum.items() if k != "r"},
        "minimum_at_right_endpoint": minimum["r"] == m0,
        "overlap_row": overlap,
        "raw_row_fingerprint_encoding":
            "ASCII lowercase-hex numerator/denominator",
        "chunks": [
            {
                "r_lo": record["r_lo"],
                "r_hi": record["r_hi"],
                "exact_min_comparisons": record["exact_min_comparisons"],
            }
            for record in records
        ],
        "workers_note":
            "one outer ProcessPoolExecutor over r-range chunks; "
            "no nested pools",
    }
    output = HERE / f"H2_PREFIX_T{t}.json"
    write_json(output, payload)
    print(f"[{'PASS' if passed else 'FAIL'}] t={t} prefix [{t},{m0}] "
          f"cells={cells} min@r={minimum['r']} "
          f"min={minimum['decimal_truncated_50'][:16]} -> {output.name}")
    if problems:
        for problem in problems:
            print(f"  problem: {problem}")
    return 0 if passed else 1


def main() -> int:
    sys.set_int_max_str_digits(1_000_000)
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=str)
    parser.add_argument("--r-start", type=int, default=None)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--assemble", type=str)
    args = parser.parse_args()
    if bool(args.plan) == bool(args.assemble):
        parser.error("exactly one of --plan or --assemble is required")
    if args.plan:
        return run_mode(args.plan, args.r_start, min(args.workers, 6))
    return assemble_mode(args.assemble)


if __name__ == "__main__":
    raise SystemExit(main())
