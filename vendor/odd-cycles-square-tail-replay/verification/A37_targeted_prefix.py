#!/usr/bin/env python3
"""Exact finite certificate for the four proof-consumed fixed-column gaps.

The square-tail proof needs only the following staircase beyond its full-row
core:

    t=5,  504<=r<=699       196 cells
    t=6,  504<=r<=2003     1500 cells
    t=7,  504<=r<=6199     5696 cells
    t=8,  504<=r<=19999   19496 cells

The recurrence is sequential in the row index.  Each row is written and
flushed to a resumable JSONL work file as soon as it is complete; the compact
JSON certificate is assembled from that work file at the end.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
DEFAULT_JSONL = HERE / "A37_TARGETED_PREFIX.work.jsonl"
DEFAULT_OUTPUT = HERE / "A37_TARGETED_PREFIX.json"
R_LO = 504
R_HI = 19_999
EXPECTED_CELLS = 26_888
HEX_DIGITS = frozenset("0123456789abcdef")
ALGORITHM = "square-tail-fixed-prefix-row-v2"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reverse_step(n: int, previous: list[int], previous_two: list[int]) -> list[int]:
    if n % 2 == 0:
        return [
            previous[k] + n * n * previous_two[k]
            for k in range(len(previous))
        ]
    return [
        (previous[k - 1] if k else 0) + n * n * previous_two[k]
        for k in range(len(previous))
    ]


def shifted(row: list[int]) -> list[int]:
    return [row[k] + (row[k - 1] if k else 0) for k in range(len(row))]


def turan(row: list[int], k: int) -> int:
    def coefficient(j: int) -> int:
        return row[j] if 0 <= j < len(row) else 0

    return coefficient(k) ** 2 - coefficient(k - 1) * coefficient(k + 1)


def rho_parts(r: int, current: list[int], predecessor: list[int], t: int):
    current_s = turan(shifted(current), t)
    predecessor_s = turan(shifted(predecessor), t)
    current_v = turan(current, t - 1)
    decrement = (
        r * r * current_s
        - (2 * r - 1) ** 4 * (r - 1) ** 2 * predecessor_s
    )
    numerator = decrement * decrement
    denominator = 2 * r * r * current_s * current_v
    return numerator, denominator, decrement


def update_integer(digest, value: int) -> None:
    """Hash an integer without Python's decimal-string size limit."""
    magnitude = abs(value).to_bytes(max(1, (abs(value).bit_length() + 7) // 8), "big")
    digest.update(b"-" if value < 0 else b"+")
    digest.update(len(magnitude).to_bytes(8, "big"))
    digest.update(magnitude)


def evaluate_row(
    r: int, current: list[int], predecessor: list[int]
) -> dict:
    t_lo = required_t_lo(r)
    digest = hashlib.sha256()
    passes = failures = positive = 0
    for t in range(t_lo, 9):
        numerator, denominator, decrement = rho_parts(
            r, current, predecessor, t
        )
        for value in (t, numerator, denominator, decrement):
            update_integer(digest, value)
        if denominator > 0 and numerator > denominator:
            passes += 1
        else:
            failures += 1
        if decrement > 0:
            positive += 1
    return {
        "algorithm": ALGORITHM,
        "r": r,
        "t_lo": t_lo,
        "t_hi": 8,
        "cells": 9 - t_lo,
        "passes": passes,
        "failures": failures,
        "positive_decrements": positive,
        "exact_sha256": digest.hexdigest(),
    }


def required_t_lo(r: int) -> int:
    if r <= 699:
        return 5
    if r <= 2_003:
        return 6
    if r <= 6_199:
        return 7
    return 8


def exact_identity_gate() -> int:
    sys.path.insert(0, str(HERE))
    from engine import exact_rho, wall_M

    kmax = 9
    previous_two = [1] + [0] * kmax
    previous = [1, 1] + [0] * (kmax - 1)
    even_rows = {1: previous_two}
    checked = 0
    walls = wall_M(40)
    for n in range(2, 35):
        current = reverse_step(n, previous, previous_two)
        previous_two, previous = previous, current
        if n % 2:
            continue
        r = n // 2 + 1
        even_rows[r] = current
        if r < 4:
            continue
        for t in range(2, min(r - 2, 8) + 1):
            numerator, denominator, _ = rho_parts(
                r, current, even_rows[r - 1], t
            )
            require(
                F(numerator, denominator) == exact_rho(r, t, walls),
                f"cancelled identity mismatch at (r,t)=({r},{t})",
            )
            checked += 1

    current = list(even_rows[8])
    current[2] += 1
    numerator, denominator, _ = rho_parts(8, current, even_rows[7], 2)
    require(
        F(numerator, denominator) != exact_rho(8, 2, walls),
        "mutated reverse diagonal was not detected",
    )
    return checked


def load_records(path: Path) -> dict[int, dict]:
    records: dict[int, dict] = {}
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
            r = record.get("r")
            require(isinstance(r, int), f"row missing at JSONL line {line_number}")
            records[r] = record
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


def record_is_complete(record: dict) -> bool:
    r = record.get("r")
    if not isinstance(r, int) or not R_LO <= r <= R_HI:
        return False
    t_lo = required_t_lo(r)
    return (
        "error" not in record
        and record.get("algorithm") == ALGORITHM
        and record.get("t_lo") == t_lo
        and record.get("t_hi") == 8
        and record.get("cells") == 9 - t_lo
        and record.get("passes") == 9 - t_lo
        and record.get("failures") == 0
        and record.get("positive_decrements") == 9 - t_lo
        and isinstance(record.get("exact_sha256"), str)
        and len(record["exact_sha256"]) == 64
        and set(record["exact_sha256"]) <= HEX_DIGITS
    )


def canonical_digest(records: dict[int, dict]) -> str:
    digest = hashlib.sha256()
    for r in range(R_LO, R_HI + 1):
        record = records[r]
        literal = json.dumps(record, sort_keys=True, separators=(",", ":"))
        digest.update(literal.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def write_summary(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def run(jsonl_path: Path, output_path: Path) -> int:
    require(
        jsonl_path.resolve() != output_path.resolve(),
        "work JSONL and compact certificate must use different paths",
    )
    identity_cells = exact_identity_gate()
    repair_trailing_jsonl(jsonl_path)
    existing = load_records(jsonl_path)
    done = {r for r, record in existing.items() if record_is_complete(record)}

    kmax = 9
    previous_two = [1] + [0] * kmax
    previous = [1, 1] + [0] * (kmax - 1)
    predecessor_even = previous_two
    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    with jsonl_path.open("a", encoding="utf-8", buffering=1) as sink:
        for n in range(2, 2 * R_HI - 1):
            current = reverse_step(n, previous, previous_two)
            previous_two, previous = previous, current
            if n % 2:
                continue
            r = n // 2 + 1
            if r < R_LO:
                predecessor_even = current
                continue
            if r > R_HI:
                break

            try:
                record = evaluate_row(r, current, predecessor_even)
            except Exception as error:
                record = {"r": r, "error": repr(error)}
            if r in done:
                require(
                    record == existing[r],
                    f"resume record disagrees with exact replay at r={r}",
                )
            else:
                sink.write(json.dumps(record, sort_keys=True) + "\n")
                sink.flush()
            predecessor_even = current

    records = load_records(jsonl_path)
    expected_rows = set(range(R_LO, R_HI + 1))
    missing = sorted(expected_rows - set(records))
    extras = sorted(set(records) - expected_rows)
    malformed = sorted(
        r for r in expected_rows & set(records) if not record_is_complete(records[r])
    )
    cells = sum(records[r].get("cells", 0) for r in expected_rows & set(records))
    ranges = [
        {"t": 5, "r_lo": 504, "r_hi": 699, "cells": 196},
        {"t": 6, "r_lo": 504, "r_hi": 2_003, "cells": 1_500},
        {"t": 7, "r_lo": 504, "r_hi": 6_199, "cells": 5_696},
        {"t": 8, "r_lo": 504, "r_hi": 19_999, "cells": 19_496},
    ]
    passed = not missing and not extras and not malformed and cells == EXPECTED_CELLS
    payload = {
        "schema": "square-tail-fixed-prefix-v1",
        "algorithm": ALGORITHM,
        "passed": passed,
        "claim": "rho(r,t)>1 on the four proof-consumed fixed-column gaps",
        "configuration": {
            "workers": 1,
            "row_range": [R_LO, R_HI],
        },
        "identity_replay_cells": identity_cells,
        "ranges": ranges,
        "rows": len(expected_rows),
        "cells": cells,
        "missing_rows": missing,
        "unexpected_rows": extras,
        "malformed_rows": malformed,
        "canonical_rows_sha256": canonical_digest(records) if passed else None,
        "resume_policy": "every reused row is recomputed before acceptance",
        "mutation_controls": {"reverse_diagonal_corruption_rejected": True},
    }
    write_summary(output_path, payload)
    print(f"[PASS] cancelled identity replay: {identity_cells} cells")
    print(
        f"[{'PASS' if passed else 'FAIL'}] targeted staircase: "
        f"{cells}/{EXPECTED_CELLS} exact comparisons on {len(expected_rows)} rows"
    )
    print(f"[INFO] resume rows revalidated: {len(done)}")
    print(f"work JSONL: {jsonl_path}")
    print(f"certificate: {output_path}")
    return 0 if passed else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    return run(args.jsonl, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
