#!/usr/bin/env python3
"""Exact definition-level certificate for the proof-consumed 500-row core.

Every interior cell on 4<=r<=503 is rebuilt from the Wall recurrence.  Blocks
are independent, run in one outer process pool, and are flushed immediately to
a resumable JSONL work file.  The compact JSON certificate is assembled from
that work file after all blocks finish.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_JSONL = HERE / "A38_TARGETED_CORE.work.jsonl"
DEFAULT_OUTPUT = HERE / "A38_TARGETED_CORE.json"
R_LO = 4
R_HI = 503
BLOCK_SIZE = 25
EXPECTED_CELLS = 125_250
HEX_DIGITS = frozenset("0123456789abcdef")
ALGORITHM = "square-tail-finite-core-block-v2"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reverse_rows(n_max: int, kmax: int) -> dict[int, list[int]]:
    previous_two = [1] + [0] * kmax
    previous = [1, 1] + [0] * (kmax - 1)
    rows = {0: previous_two, 1: previous}
    for n in range(2, n_max + 1):
        if n % 2 == 0:
            current = [
                previous[k] + n * n * previous_two[k]
                for k in range(kmax + 1)
            ]
        else:
            current = [
                (previous[k - 1] if k else 0) + n * n * previous_two[k]
                for k in range(kmax + 1)
            ]
        rows[n] = current
        previous_two, previous = previous, current
    return rows


def shifted(row):
    return [row[k] + (row[k - 1] if k else 0) for k in range(len(row))]


def turan(row, k):
    zero = F(0) if row and isinstance(row[0], F) else 0

    def coefficient(j):
        return row[j] if 0 <= j < len(row) else zero

    return coefficient(k) ** 2 - coefficient(k - 1) * coefficient(k + 1)


def rho_parts(rows, r: int, t: int):
    current = rows[2 * r - 2]
    predecessor = rows[2 * r - 4]
    current_s = turan(shifted(current), t)
    predecessor_s = turan(shifted(predecessor), t)
    current_v = turan(current, t - 1)
    decrement = (
        r * r * current_s
        - (2 * r - 1) ** 4 * (r - 1) ** 2 * predecessor_s
    )
    return (
        decrement * decrement,
        2 * r * r * current_s * current_v,
        decrement,
    )


def integer_encoding(value: int) -> bytes:
    """Canonical binary encoding, independent of decimal-string limits."""
    magnitude = abs(value).to_bytes(max(1, (abs(value).bit_length() + 7) // 8), "big")
    return (
        (b"-" if value < 0 else b"+")
        + len(magnitude).to_bytes(8, "big")
        + magnitude
    )


def update_integer(digest, value: int) -> None:
    digest.update(integer_encoding(value))


def decimal_truncation(numerator: int, denominator: int, places: int = 18) -> str:
    """Return a deterministic, explicitly truncated decimal expansion."""
    require(denominator > 0, "decimal denominator must be positive")
    sign = "-" if numerator < 0 else ""
    whole, remainder = divmod(abs(numerator), denominator)
    digits = []
    for _ in range(places):
        remainder *= 10
        digit, remainder = divmod(remainder, denominator)
        digits.append(str(digit))
    return f"{sign}{whole}." + "".join(digits)


def definitional_rho(rows, r: int, t: int) -> F:
    def pi(rr):
        return [
            F(value, math.factorial(2 * rr - 1))
            for value in rows[2 * rr - 2]
        ]

    current = pi(r)
    predecessor = pi(r - 1)
    current_s = turan(shifted(current), t)
    predecessor_s = turan(shifted(predecessor), t)
    current_v = turan(current, t - 1)
    decrement = (
        current_s - F((2 * r - 1) ** 2, (2 * r) ** 2) * predecessor_s
    )
    return F(r * r) * decrement * decrement / (2 * current_s * current_v)


def identity_gate() -> int:
    rows = reverse_rows(34, 18)
    checked = 0
    for r in range(4, 19):
        for t in range(2, r - 1):
            numerator, denominator, _ = rho_parts(rows, r, t)
            require(
                F(numerator, denominator) == definitional_rho(rows, r, t),
                f"cancelled identity mismatch at (r,t)=({r},{t})",
            )
            checked += 1

    mutated = dict(rows)
    mutated[14] = list(rows[14])
    mutated[14][2] += 1
    numerator, denominator, _ = rho_parts(mutated, 8, 2)
    require(
        F(numerator, denominator) != definitional_rho(rows, 8, 2),
        "mutated reverse diagonal was not detected",
    )
    return checked


def evaluate_block(block: tuple[int, int]) -> dict:
    start, end = block
    rows = reverse_rows(2 * end - 2, end)
    cells = failures = positive = 0
    digest = hashlib.sha256()
    worst_num = worst_den = None
    worst_cell = None
    for r in range(start, end + 1):
        try:
            for t in range(2, r - 1):
                numerator, denominator, decrement = rho_parts(rows, r, t)
                for value in (r, t, numerator, denominator, decrement):
                    update_integer(digest, value)
                cells += 1
                if denominator <= 0 or numerator <= denominator:
                    failures += 1
                if decrement > 0:
                    positive += 1
                if (
                    worst_num is None
                    or numerator * worst_den < worst_num * denominator
                ):
                    worst_num, worst_den = numerator, denominator
                    worst_cell = [r, t]
        except Exception as error:
            return {
                "block": [start, end],
                "error": f"row {r}: {error!r}",
            }
    return {
        "algorithm": ALGORITHM,
        "block": [start, end],
        "cells": cells,
        "passes": cells - failures,
        "failures": failures,
        "positive_decrements": positive,
        "worst_cell": worst_cell,
        "worst_numerator_hex": format(worst_num, "x"),
        "worst_denominator_hex": format(worst_den, "x"),
        "exact_sha256": digest.hexdigest(),
    }


def blocks() -> list[tuple[int, int]]:
    return [
        (start, min(start + BLOCK_SIZE - 1, R_HI))
        for start in range(R_LO, R_HI + 1, BLOCK_SIZE)
    ]


def load_records(path: Path) -> dict[tuple[int, int], dict]:
    records: dict[tuple[int, int], dict] = {}
    if not path.exists():
        return records
    with path.open(encoding="utf-8") as source:
        for line_number, literal in enumerate(source, 1):
            try:
                record = json.loads(literal)
            except json.JSONDecodeError as error:
                raise RuntimeError(
                    f"malformed JSONL line {line_number}: {error}"
                ) from error
            raw_block = record.get("block")
            require(
                isinstance(raw_block, list) and len(raw_block) == 2,
                f"missing block at JSONL line {line_number}",
            )
            records[tuple(raw_block)] = record
    return records


def repair_trailing_jsonl(path: Path) -> None:
    """Keep complete flushed lines and repair only a torn final write."""
    if not path.exists():
        return
    data = path.read_bytes()
    if not data or data.endswith(b"\n"):
        return
    start = data.rfind(b"\n") + 1
    prefix, tail = data[:start], data[start:]
    try:
        json.loads(tail.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        repaired = prefix
    else:
        repaired = prefix + tail + b"\n"
    temporary = path.with_suffix(path.suffix + ".repair.tmp")
    with temporary.open("wb") as sink:
        sink.write(repaired)
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def record_is_complete(record: dict, block: tuple[int, int]) -> bool:
    start, end = block
    cells = sum(r - 3 for r in range(start, end + 1))
    worst_cell = record.get("worst_cell")
    numerator_hex = record.get("worst_numerator_hex")
    denominator_hex = record.get("worst_denominator_hex")
    digest = record.get("exact_sha256")
    return (
        "error" not in record
        and record.get("algorithm") == ALGORITHM
        and record.get("block") == [start, end]
        and record.get("cells") == cells
        and record.get("passes") == cells
        and record.get("failures") == 0
        and record.get("positive_decrements") == cells
        and isinstance(worst_cell, list)
        and len(worst_cell) == 2
        and all(isinstance(value, int) for value in worst_cell)
        and start <= worst_cell[0] <= end
        and 2 <= worst_cell[1] <= worst_cell[0] - 2
        and isinstance(numerator_hex, str)
        and bool(numerator_hex)
        and set(numerator_hex) <= HEX_DIGITS
        and int(numerator_hex, 16) > 0
        and format(int(numerator_hex, 16), "x") == numerator_hex
        and isinstance(denominator_hex, str)
        and bool(denominator_hex)
        and set(denominator_hex) <= HEX_DIGITS
        and int(denominator_hex, 16) > 0
        and format(int(denominator_hex, 16), "x") == denominator_hex
        and isinstance(digest, str)
        and len(digest) == 64
        and set(digest) <= HEX_DIGITS
    )


def write_summary(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def canonical_digest(records: dict[tuple[int, int], dict]) -> str:
    digest = hashlib.sha256()
    for block in blocks():
        literal = json.dumps(records[block], sort_keys=True, separators=(",", ":"))
        digest.update(literal.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def run(
    jsonl_path: Path, output_path: Path, workers: int, require_fresh: bool = False
) -> int:
    require(workers >= 1, "worker count must be positive")
    require(
        jsonl_path.resolve() != output_path.resolve(),
        "work JSONL and compact certificate must use different paths",
    )
    identity_cells = identity_gate()
    expected_blocks = blocks()
    if require_fresh:
        require(
            not jsonl_path.exists() or jsonl_path.stat().st_size == 0,
            "--require-fresh refuses a nonempty work JSONL",
        )
    repair_trailing_jsonl(jsonl_path)
    existing = load_records(jsonl_path)
    todo = [
        block
        for block in expected_blocks
        if not record_is_complete(existing.get(block, {}), block)
    ]
    reused_blocks = len(expected_blocks) - len(todo)

    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    with jsonl_path.open("a", encoding="utf-8", buffering=1) as sink:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(evaluate_block, block): block for block in todo}
            for future in as_completed(futures):
                block = futures[future]
                try:
                    record = future.result()
                except Exception as error:
                    record = {"block": list(block), "error": repr(error)}
                sink.write(json.dumps(record, sort_keys=True) + "\n")
                sink.flush()

    records = load_records(jsonl_path)
    malformed = [
        list(block)
        for block in expected_blocks
        if not record_is_complete(records.get(block, {}), block)
    ]
    extras = [
        list(block) for block in records if block not in set(expected_blocks)
    ]
    cells = sum(
        records.get(block, {}).get("cells", 0) for block in expected_blocks
    )
    passed = not malformed and not extras and cells == EXPECTED_CELLS
    worst = None
    if passed:
        candidates = [
            (
                F(
                    int(record["worst_numerator_hex"], 16),
                    int(record["worst_denominator_hex"], 16),
                ),
                record["worst_cell"],
            )
            for block, record in records.items()
            if block in set(expected_blocks) and record_is_complete(record, block)
        ]
        worst = min(candidates)
        witness_literal = (
            integer_encoding(worst[0].numerator)
            + integer_encoding(worst[0].denominator)
        )
    payload = {
        "schema": "square-tail-finite-core-v1",
        "algorithm": ALGORITHM,
        "passed": passed,
        "claim": "rho(r,t)>1 for 4<=r<=503 and 2<=t<=r-2",
        "configuration": {
            "worker_policy": (
                "one outer ProcessPoolExecutor; default workers=os.cpu_count()"
            ),
            "row_range": [R_LO, R_HI],
            "block_size": BLOCK_SIZE,
        },
        "identity_replay_cells": identity_cells,
        "blocks": len(expected_blocks),
        "cells": cells,
        "malformed_blocks": malformed,
        "unexpected_blocks": extras,
        "finite_box_witness": {
            "cell": worst[1] if worst else None,
            "numerator_hex": format(worst[0].numerator, "x") if worst else None,
            "denominator_hex": format(worst[0].denominator, "x") if worst else None,
            "reporting_decimal_truncated": (
                decimal_truncation(worst[0].numerator, worst[0].denominator)
                if worst else None
            ),
            "exact_fraction_sha256": (
                hashlib.sha256(witness_literal).hexdigest() if worst else None
            ),
            "role": "minimum among these replayed cells, not a global infimum",
        },
        "canonical_blocks_sha256": canonical_digest(records) if passed else None,
        "resume_policy": (
            "structurally checked local recovery cache; authoritative release replay "
            "uses --require-fresh"
        ),
        "mutation_controls": {"reverse_diagonal_corruption_rejected": True},
    }
    write_summary(output_path, payload)
    print(f"[PASS] cancelled identity replay: {identity_cells} cells")
    print(f"[INFO] workers: {workers}")
    print(f"[INFO] resume blocks accepted: {reused_blocks}")
    print(
        f"[{'PASS' if passed else 'FAIL'}] finite core: "
        f"{cells}/{EXPECTED_CELLS} exact comparisons in {len(expected_blocks)} blocks"
    )
    if worst:
        print(f"[INFO] finite-core witness: {worst[0]} at {tuple(worst[1])}")
    print(f"work JSONL: {jsonl_path}")
    print(f"certificate: {output_path}")
    return 0 if passed else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--require-fresh", action="store_true")
    args = parser.parse_args()
    return run(args.jsonl, args.output, args.workers, args.require_fresh)


if __name__ == "__main__":
    raise SystemExit(main())
