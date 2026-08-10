#!/usr/bin/env python3
"""Exact sign gate for the square-tail decrement in column t=2.

The analytic profile argument proves d_(r,t)>0 for every interior t>=3,
and the two terminal columns are explicit.  This file binds the oriented
fixed-column tail from r=313 onward to the only finite gap:

    d_(r,2) > 0,  4 <= r <= 312.

The gate uses the cancelled integer numerator M and checks its normalization
against the literal Fraction definition on every row through r=18.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from fractions import Fraction as F
from pathlib import Path


R_LO = 4
R_HI = 312
T = 2
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "A35_T2_DECREMENT_SIGN.json"
sys.path.insert(0, str(HERE))
import A1_verify_tail as a1


def wall_reverse_rows(n_max: int, k_max: int) -> dict[int, list[int]]:
    """Return the top reverse diagonals used by the cancelled SQ formula."""
    previous_two = [1] + [0] * k_max
    previous_one = [1, 1] + [0] * (k_max - 1)
    rows = {0: previous_two, 1: previous_one}
    for n in range(2, n_max + 1):
        if n % 2 == 0:
            current = [
                previous_one[k] + n * n * previous_two[k]
                for k in range(k_max + 1)
            ]
        else:
            current = [
                (previous_one[k - 1] if k else 0) + n * n * previous_two[k]
                for k in range(k_max + 1)
            ]
        rows[n] = current
        previous_two, previous_one = previous_one, current
    return rows


def shifted(row: list[int]) -> list[int]:
    return [row[k] + (row[k - 1] if k else 0) for k in range(len(row))]


def turan(row: list[int] | list[F], k: int):
    zero = F(0) if row and isinstance(row[0], F) else 0

    def coefficient(index: int):
        return row[index] if 0 <= index < len(row) else zero

    return coefficient(k) ** 2 - coefficient(k - 1) * coefficient(k + 1)


def cancelled_numerator(rows: dict[int, list[int]], r: int) -> int:
    current = rows[2 * r - 2]
    predecessor = rows[2 * r - 4]
    current_s = turan(shifted(current), T)
    predecessor_s = turan(shifted(predecessor), T)
    return (
        r * r * current_s
        - (2 * r - 1) ** 4 * (r - 1) ** 2 * predecessor_s
    )


def definitional_decrement(rows: dict[int, list[int]], r: int) -> F:
    current_factorial = math.factorial(2 * r - 1)
    predecessor_factorial = math.factorial(2 * r - 3)
    current = [F(value, current_factorial) for value in rows[2 * r - 2]]
    predecessor = [
        F(value, predecessor_factorial) for value in rows[2 * r - 4]
    ]
    current_s = turan(shifted(current), T)
    predecessor_s = turan(shifted(predecessor), T)
    return current_s - F((2 * r - 1) ** 2, (2 * r) ** 2) * predecessor_s


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(payload: dict, output: Path) -> None:
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, output)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Certify the finite t=2 decrement-sign gap exactly."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT,
        help=f"certificate path (default: {OUTPUT.name})",
    )
    args = parser.parse_args()

    rows = wall_reverse_rows(2 * R_HI - 2, 3)

    normalization_failures = []
    for r in range(R_LO, 19):
        numerator = cancelled_numerator(rows, r)
        expected = F(
            numerator,
            r * r * math.factorial(2 * r - 1) ** 2,
        )
        actual = definitional_decrement(rows, r)
        if actual != expected:
            normalization_failures.append(r)

    values = [(cancelled_numerator(rows, r), r) for r in range(R_LO, R_HI + 1)]
    sign_failures = [r for value, r in values if value <= 0]
    minimum_value, minimum_row = min(values)
    tail_orientation, _, _ = a1.section3_tail(
        312, t=T, verbose=False, record=False
    )
    passed = not normalization_failures and not sign_failures and tail_orientation

    payload = {
        "schema": "square-tail-t2-orientation-v1",
        "version": 1,
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "claim": "d_(r,2)>0 for every integer r>=4",
        "cells": len(values),
        "normalization_replay": {
            "range": "4<=r<=18",
            "failures": normalization_failures,
        },
        "sign_failures": sign_failures,
        "tail_orientation": {
            "anchor_m0": 312,
            "range": "integer r>=313",
            "passed": tail_orientation,
            "gate": "V*G-corr has nonnegative coefficients and positive constant term",
        },
        "smallest_integer_numerator": {
            "r": minimum_row,
            "M": str(minimum_value),
        },
        "configuration": {
            "workers": 1,
        },
    }
    write_json(payload, args.output)

    print("A35 -- finite t=2 decrement sign")
    print(
        f"  [{'PASS' if not normalization_failures else 'FAIL'}] "
        "cancelled numerator equals the definitional decrement, r=4..18"
    )
    print(
        f"  [{'PASS' if not sign_failures else 'FAIL'}] "
        f"M>0 on {len(values)} rows, r={R_LO}..{R_HI}"
    )
    print(
        f"  [{'PASS' if tail_orientation else 'FAIL'}] "
        "oriented analytic tail, r>=313"
    )
    print(f"  smallest raw M={minimum_value} at r={minimum_row}")
    print(f"  output: {args.output}  sha256={sha256(args.output)}")
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
