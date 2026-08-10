#!/usr/bin/env python3
"""Strict uniqueness of the global minimiser (129,2) on the proof-consumed core.

Everything outside the finite core is already bounded below by 203/200 by the
retargeted tails (A39), the bounded windows (A41), the top/joint dichotomy and
the terminal columns, and rho(129,2) < 203/200.  Uniqueness of the global
minimiser therefore reduces entirely to the 125,250 interior cells

    4 <= r <= 503,  2 <= t <= r-2,

which is what this file settles.  For every such cell the exact cancelled
integer pair (numerator, denominator) of rho(r,t) is rebuilt from the Wall
recurrence and compared with rho(129,2) by cross-multiplication only; no
Fraction reduction enters the comparison.

The certificate records the number of cells at or below the witness (which must
be exactly one), the runner-up cell and its exact separation from the witness,
and the exact per-column minima over the core.

Blocks are independent, run in one outer process pool, and flush immediately to
a resumable JSONL work file.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFY = HERE
RESULTS = HERE.parent / "results"
DEFAULT_JSONL = HERE / "A43_CORE_UNIQUENESS.work.jsonl"
DEFAULT_OUTPUT = RESULTS / "core_uniqueness.json"

R_LO, R_HI, BLOCK_SIZE = 4, 503, 25
EXPECTED_CELLS = 125_250
WITNESS = (129, 2)
ALGORITHM = "square-tail-core-uniqueness-block-v1"

# Bound externally: both published fingerprints of rho(129,2).
A40_DECIMAL_SHA256 = "9c09b2864de0e742c291698d4c63ad631d2d872c122f62b6c0fdc8c10b47e373"
A38_BINARY_SHA256 = "63929d5845e7dfff87de264d7945ddb50c2e58f0d9ba73cb3484b46fd1e7b744"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_a38():
    import importlib.util

    path = VERIFY / "A38_targeted_core.py"
    require(path.is_file(), f"missing producer {path}")
    spec = importlib.util.spec_from_file_location("_a38_core", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A38 = load_a38()
reverse_rows = A38.reverse_rows
rho_parts = A38.rho_parts
definitional_rho = A38.definitional_rho
integer_encoding = A38.integer_encoding
decimal_truncation = A38.decimal_truncation

_STAR_NUM = _STAR_DEN = None


def witness_fraction() -> tuple[int, int]:
    """rho(129,2) in lowest terms, from the same cancelled integer machinery."""
    r, t = WITNESS
    rows = reverse_rows(2 * r - 2, r)
    numerator, denominator, _ = rho_parts(rows, r, t)
    require(
        F(numerator, denominator) == definitional_rho(rows, r, t),
        "witness cancelled identity mismatch",
    )
    value = F(numerator, denominator)
    return value.numerator, value.denominator


def initializer(star_num: int, star_den: int) -> None:
    global _STAR_NUM, _STAR_DEN
    _STAR_NUM, _STAR_DEN = star_num, star_den


def evaluate_block(block: tuple[int, int]) -> dict:
    start, end = block
    star_num, star_den = _STAR_NUM, _STAR_DEN
    rows = reverse_rows(2 * end - 2, end)
    cells = 0
    at_or_below: list[list[int]] = []
    equal_cells: list[list[int]] = []
    runner_cell = None
    runner_num = runner_den = None      # smallest rho strictly above the witness
    column_min: dict[int, tuple[int, int, int]] = {}   # t -> (num, den, r)
    digest = hashlib.sha256()
    sign_digest = hashlib.sha256()

    for r in range(start, end + 1):
        for t in range(2, r - 1):
            numerator, denominator, _ = rho_parts(rows, r, t)
            require(denominator > 0, f"non-positive denominator at ({r},{t})")
            cells += 1
            # rho(r,t) ? rho*   <=>   numerator*star_den ? denominator*star_num
            left = numerator * star_den
            right = denominator * star_num
            for value in (r, t, left - right):
                digest.update(integer_encoding(value))
            sign_digest.update(integer_encoding(r))
            sign_digest.update(integer_encoding(t))
            sign_digest.update(b">" if left > right else (b"=" if left == right else b"<"))
            if left < right:
                at_or_below.append([r, t])
            elif left == right:
                at_or_below.append([r, t])
                equal_cells.append([r, t])
            else:
                if runner_num is None or numerator * runner_den < runner_num * denominator:
                    runner_num, runner_den, runner_cell = numerator, denominator, [r, t]
            previous = column_min.get(t)
            if previous is None or numerator * previous[1] < previous[0] * denominator:
                column_min[t] = (numerator, denominator, r)

    return {
        "block": [start, end],
        "cells": cells,
        "at_or_below_witness": at_or_below,
        "equal_to_witness": equal_cells,
        "runner_up": {
            "cell": runner_cell,
            "numerator_hex": format(runner_num, "x") if runner_num is not None else None,
            "denominator_hex": format(runner_den, "x") if runner_den is not None else None,
        },
        "column_minima": {
            str(t): {"r": r, "numerator_hex": format(n, "x"), "denominator_hex": format(d, "x")}
            for t, (n, d, r) in sorted(column_min.items())
        },
        "block_sha256": digest.hexdigest(),
        "sign_sha256": sign_digest.hexdigest(),
    }


def identity_and_mutation_gates(star_num: int, star_den: int) -> dict:
    """Independent re-derivation of the witness plus two negative controls."""
    rows = reverse_rows(2 * WITNESS[0] - 2, WITNESS[0])

    # (1) both published fingerprints of the same rational
    decimal_digest = hashlib.sha256(f"{star_num}/{star_den}".encode()).hexdigest()
    binary_digest = hashlib.sha256(
        integer_encoding(star_num) + integer_encoding(star_den)
    ).hexdigest()

    # (2) a mutated reverse diagonal must break the comparison
    mutated = dict(rows)
    mutated[100] = list(rows[100])
    mutated[100][3] += 1
    bad_num, bad_den, _ = rho_parts(mutated, 51, 2)
    good_num, good_den, _ = rho_parts(rows, 51, 2)
    mutation_detected = bad_num * good_den != good_num * bad_den

    # (3) a false witness must be rejected: rho(128,2) and rho(130,2) both exceed rho*
    wider = reverse_rows(2 * 130 - 2, 130)
    neighbours = {}
    for r in (128, 130):
        numerator, denominator, _ = rho_parts(wider, r, 2)
        neighbours[str(r)] = numerator * star_den > denominator * star_num
    false_witness_rejected = all(neighbours.values())

    return {
        "witness_decimal_sha256": decimal_digest,
        "witness_decimal_sha256_matches_A40": decimal_digest == A40_DECIMAL_SHA256,
        "witness_binary_sha256": binary_digest,
        "witness_binary_sha256_matches_A38": binary_digest == A38_BINARY_SHA256,
        "mutated_reverse_diagonal_detected": mutation_detected,
        "false_witness_neighbours_rejected": false_witness_rejected,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--require-fresh", action="store_true")
    arguments = parser.parse_args()

    star_num, star_den = witness_fraction()
    gates = identity_and_mutation_gates(star_num, star_den)
    for name in (
        "witness_decimal_sha256_matches_A40",
        "witness_binary_sha256_matches_A38",
        "mutated_reverse_diagonal_detected",
        "false_witness_neighbours_rejected",
    ):
        require(gates[name], f"gate failed: {name}")

    blocks = [
        (start, min(start + BLOCK_SIZE - 1, R_HI))
        for start in range(R_LO, R_HI + 1, BLOCK_SIZE)
    ]

    done: dict[str, dict] = {}
    if arguments.jsonl.is_file() and not arguments.require_fresh:
        for line in arguments.jsonl.read_text().splitlines():
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if record.get("algorithm") == ALGORITHM and record.get("witness_den") == star_den:
                done[json.dumps(record["block"])] = record

    pending = [b for b in blocks if json.dumps(list(b)) not in done]
    print(
        f"core uniqueness: {len(blocks)} blocks, {len(pending)} to run, "
        f"{len(done)} resumed",
        flush=True,
    )

    if pending:
        with arguments.jsonl.open("a", encoding="utf-8") as handle:
            with ProcessPoolExecutor(
                max_workers=max(1, min(arguments.workers, len(pending))),
                initializer=initializer,
                initargs=(star_num, star_den),
            ) as pool:
                futures = {pool.submit(evaluate_block, b): b for b in pending}
                for future in as_completed(futures):
                    block = futures[future]
                    try:
                        record = future.result()
                    except Exception as error:            # one block never kills the batch
                        record = {"block": list(block), "error": repr(error)}
                    record["algorithm"] = ALGORITHM
                    record["witness_den"] = star_den
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    handle.flush()
                    done[json.dumps(list(block))] = record
                    print(f"  block {block} done", flush=True)

    records = [done[json.dumps(list(b))] for b in blocks]
    errors = [r for r in records if "error" in r]
    require(not errors, f"block errors: {errors}")

    cells = sum(r["cells"] for r in records)
    at_or_below = [c for r in records for c in r["at_or_below_witness"]]
    equal = [c for r in records for c in r["equal_to_witness"]]

    runner_num = runner_den = None
    runner_cell = None
    for record in records:
        candidate = record["runner_up"]
        if candidate["cell"] is None:
            continue
        n = int(candidate["numerator_hex"], 16)
        d = int(candidate["denominator_hex"], 16)
        if runner_num is None or n * runner_den < runner_num * d:
            runner_num, runner_den, runner_cell = n, d, candidate["cell"]

    columns: dict[int, tuple[int, int, int]] = {}
    for record in records:
        for key, entry in record["column_minima"].items():
            t = int(key)
            n = int(entry["numerator_hex"], 16)
            d = int(entry["denominator_hex"], 16)
            r = entry["r"]
            previous = columns.get(t)
            if previous is None or n * previous[1] < previous[0] * d:
                columns[t] = (n, d, r)

    ordered = sorted(columns.items())
    column_increasing = all(
        a[1][0] * b[1][1] < b[1][0] * a[1][1] for a, b in zip(ordered, ordered[1:])
    )

    separation = F(runner_num, runner_den) - F(star_num, star_den)
    # Core-only. Global uniqueness is the join's assertion, not this file's.
    unique = cells == EXPECTED_CELLS and at_or_below == [list(WITNESS)] and equal == [list(WITNESS)]

    aggregate = hashlib.sha256()
    sign_aggregate = hashlib.sha256()
    for record in records:
        aggregate.update(record["block_sha256"].encode())
        sign_aggregate.update(record["sign_sha256"].encode())

    payload = {
        "algorithm": ALGORITHM,
        "schema": "square-tail-core-uniqueness-v1",
        "claim": (
            "rho(r,t) > rho(129,2) for every core cell 4<=r<=503, 2<=t<=r-2 "
            "other than (129,2) itself"
        ),
        "scope": (
            "finite core only; every cell outside it is already above 203/200 "
            "> rho(129,2) by the retargeted tails, bounded windows, top/joint "
            "dichotomy and terminal columns"
        ),
        "cells": cells,
        "expected_cells": EXPECTED_CELLS,
        "blocks": len(blocks),
        "witness": {
            "cell": list(WITNESS),
            "numerator_digits": len(str(star_num)),
            "denominator_digits": len(str(star_den)),
            "decimal_truncated": decimal_truncation(star_num, star_den, 30),
        },
        "cells_at_or_below_witness": at_or_below,
        "cells_equal_to_witness": equal,
        "runner_up": {
            "cell": runner_cell,
            "decimal_truncated": decimal_truncation(runner_num, runner_den, 30),
            "separation_from_witness_decimal": decimal_truncation(
                separation.numerator, separation.denominator, 24
            )
            if separation.numerator >= 0
            else None,
            "separation_numerator_digits": len(str(separation.numerator)),
            "separation_denominator_digits": len(str(separation.denominator)),
            "separation_positive": separation > 0,
        },
        "column_minima": {
            str(t): {
                "argmin_r": r,
                "decimal_truncated": decimal_truncation(n, d, 24),
            }
            for t, (n, d, r) in ordered
        },
        "column_minima_strictly_increasing_in_t": column_increasing,
        "witness_fingerprints": {
            "decimal_sha256": gates.pop("witness_decimal_sha256"),
            "binary_sha256": gates.pop("witness_binary_sha256"),
            "decimal_encoding": "sha256 of the ASCII string \"numerator/denominator\"",
            "binary_encoding": "sha256 of integer_encoding(numerator)+integer_encoding(denominator)",
        },
        "gates": {
            **gates,
            "cell count matches core": cells == EXPECTED_CELLS,
            "exactly one cell at or below witness": at_or_below == [list(WITNESS)],
            "the witness attains equality": equal == [list(WITNESS)],
            "runner-up strictly above witness": separation > 0,
            "column minima strictly increase in t": column_increasing,
        },
        "unique_minimiser_on_core": unique,
        "not_certified_here": [
            "global uniqueness; that needs the outside-core floor, supplied by the "
            "retargeted tails, the bounded windows, the top/joint dichotomy and the "
            "terminal columns, and asserted only by the join"
        ],
        "canonical_blocks_sha256": aggregate.hexdigest(),
        "canonical_sign_sha256": sign_aggregate.hexdigest(),
        "configuration": {
            "worker_policy": "one outer ProcessPoolExecutor; no nested pool",
            "block_size": BLOCK_SIZE,
        },
        "passed": unique and all(
            v for v in (
                gates["witness_decimal_sha256_matches_A40"],
                gates["witness_binary_sha256_matches_A38"],
                gates["mutated_reverse_diagonal_detected"],
                gates["false_witness_neighbours_rejected"],
                separation > 0,
            )
        ),
    }
    payload["status"] = "PASS" if payload["passed"] else "FAIL"
    arguments.output.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: payload[k] for k in ("status", "cells", "unique_minimiser_on_core")}))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
