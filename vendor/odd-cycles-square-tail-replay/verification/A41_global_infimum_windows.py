#!/usr/bin/env python3
"""Exact finite comparisons for the residual global-infimum windows.

The analytic reduction leaves one finite prefix before the t=5 tail and three
bounded windows before the t=6,7,8 cells enter the existing top/joint regions.
Every comparison is the exact cancelled integer inequality rho>203/200.

One compact row record is appended and flushed immediately.  A matching,
structurally complete row from the same source digest is skipped on resume.
The consolidated JSON certificate is rebuilt from the JSONL at the end.
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
RESULTS = HERE.parent / "results"
DEFAULT_JSONL = HERE / "A41_GLOBAL_INFIMUM_WINDOWS.work.jsonl"
DEFAULT_OUTPUT = RESULTS / "global_infimum_windows.json"
SCHEMA = "square-tail-global-infimum-windows-v1"
ALGORITHM = "square-tail-global-infimum-window-row-v1"
TARGET = F(203, 200)
RANGES = {
    5: (504, 1_000),
    6: (504, 3_633),
    7: (1_468, 8_413),
    8: (9_884, 18_456),
}
EXPECTED_CELLS = sum(end - start + 1 for start, end in RANGES.values())
EXPECTED_ROWS = frozenset(
    r
    for start, end in RANGES.values()
    for r in range(start, end + 1)
)
HEX_DIGITS = frozenset("0123456789abcdef")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def integer_encoding(value: int) -> bytes:
    magnitude = abs(value).to_bytes(max(1, (abs(value).bit_length() + 7) // 8), "big")
    return (
        (b"-" if value < 0 else b"+")
        + len(magnitude).to_bytes(8, "big")
        + magnitude
    )


def update_integer(digest, value: int) -> None:
    digest.update(integer_encoding(value))


def required_ts(r: int) -> list[int]:
    return [
        t for t, (start, end) in RANGES.items()
        if start <= r <= end
    ]


def expected_rows() -> set[int]:
    return set(EXPECTED_ROWS)


def evaluate_row(
    r: int,
    current: list[int],
    predecessor: list[int],
    source_digest: str,
) -> dict:
    columns = required_ts(r)
    require(columns, f"row {r} is outside every requested range")
    digest = hashlib.sha256()
    passes = failures = positive = 0
    for t in columns:
        numerator, denominator, decrement = rho_parts(r, current, predecessor, t)
        margin = TARGET.denominator * numerator - TARGET.numerator * denominator
        for value in (r, t, numerator, denominator, decrement, margin):
            update_integer(digest, value)
        if denominator > 0 and margin > 0:
            passes += 1
        else:
            failures += 1
        if decrement > 0:
            positive += 1
    return {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "source_sha256": source_digest,
        "target": f"{TARGET.numerator}/{TARGET.denominator}",
        "r": r,
        "columns": columns,
        "cells": len(columns),
        "passes": passes,
        "failures": failures,
        "positive_decrements": positive,
        "exact_sha256": digest.hexdigest(),
    }


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
            numerator, denominator, _ = rho_parts(r, current, even_rows[r - 1], t)
            require(
                F(numerator, denominator) == exact_rho(r, t, walls),
                f"cancelled identity mismatch at (r,t)=({r},{t})",
            )
            checked += 1

    mutated = list(even_rows[8])
    mutated[2] += 1
    numerator, denominator, _ = rho_parts(8, mutated, even_rows[7], 2)
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
            require(isinstance(r, int), f"missing row at JSONL line {line_number}")
            records[r] = record
    return records


def repair_trailing_jsonl(path: Path) -> None:
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


def record_complete(record: dict, source_digest: str) -> bool:
    r = record.get("r")
    if not isinstance(r, int) or r not in EXPECTED_ROWS:
        return False
    columns = required_ts(r)
    digest = record.get("exact_sha256")
    return (
        record.get("schema") == SCHEMA
        and record.get("algorithm") == ALGORITHM
        and record.get("source_sha256") == source_digest
        and record.get("target") == f"{TARGET.numerator}/{TARGET.denominator}"
        and record.get("columns") == columns
        and record.get("cells") == len(columns)
        and record.get("passes") == len(columns)
        and record.get("failures") == 0
        and record.get("positive_decrements") == len(columns)
        and isinstance(digest, str)
        and len(digest) == 64
        and set(digest) <= HEX_DIGITS
    )


def canonical_digest(records: dict[int, dict], rows: set[int]) -> str:
    digest = hashlib.sha256()
    for r in sorted(rows):
        literal = json.dumps(records[r], sort_keys=True, separators=(",", ":"))
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
    require(jsonl_path.resolve() != output_path.resolve(), "JSONL and summary paths must differ")
    source_digest = sha256(Path(__file__).resolve())
    identity_cells = exact_identity_gate()
    repair_trailing_jsonl(jsonl_path)
    existing = load_records(jsonl_path)
    rows = expected_rows()
    done = {
        r for r in rows
        if record_complete(existing.get(r, {}), source_digest)
    }

    kmax = 9
    previous_two = [1] + [0] * kmax
    previous = [1, 1] + [0] * (kmax - 1)
    predecessor_even = previous_two
    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    with jsonl_path.open("a", encoding="utf-8", buffering=1) as sink:
        for n in range(2, 2 * max(rows) - 1):
            current = reverse_step(n, previous, previous_two)
            previous_two, previous = previous, current
            if n % 2:
                continue
            r = n // 2 + 1
            if r > max(rows):
                break
            if r in rows and r not in done:
                try:
                    record = evaluate_row(r, current, predecessor_even, source_digest)
                except Exception as error:
                    record = {
                        "schema": SCHEMA,
                        "algorithm": ALGORITHM,
                        "source_sha256": source_digest,
                        "target": f"{TARGET.numerator}/{TARGET.denominator}",
                        "r": r,
                        "error": repr(error),
                    }
                sink.write(json.dumps(record, sort_keys=True) + "\n")
                sink.flush()
            predecessor_even = current

    records = load_records(jsonl_path)
    missing = sorted(rows - set(records))
    unexpected = sorted(set(records) - rows)
    malformed = sorted(
        r for r in rows & set(records)
        if not record_complete(records[r], source_digest)
    )
    cells = sum(records[r].get("cells", 0) for r in rows & set(records))
    passed = not missing and not unexpected and not malformed and cells == EXPECTED_CELLS
    ranges = [
        {"t": t, "r_lo": start, "r_hi": end, "cells": end - start + 1}
        for t, (start, end) in sorted(RANGES.items())
    ]
    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "claim": "rho(r,t)>203/200 on the four residual finite column ranges",
        "target": f"{TARGET.numerator}/{TARGET.denominator}",
        "source_sha256": source_digest,
        "configuration": {
            "workers": 1,
            "worker_reason": "serial row recurrence; each row comparison is sub-second",
        },
        "identity_replay_cells": identity_cells,
        "ranges": ranges,
        "rows": len(rows),
        "cells": cells,
        "expected_cells": EXPECTED_CELLS,
        "resume_rows_reused": len(done),
        "missing_rows": missing,
        "unexpected_rows": unexpected,
        "malformed_rows": malformed,
        "canonical_rows_sha256": canonical_digest(records, rows) if passed else None,
        "mutation_controls": {"reverse_diagonal_corruption_rejected": True},
    }
    write_summary(output_path, payload)
    print(f"[PASS] cancelled identity replay: {identity_cells} cells")
    print(f"[INFO] interpreter: python3.14t; workers: 1; resumed rows: {len(done)}")
    print(
        f"[{'PASS' if passed else 'FAIL'}] residual windows: "
        f"{cells}/{EXPECTED_CELLS} exact comparisons on {len(rows)} rows"
    )
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
