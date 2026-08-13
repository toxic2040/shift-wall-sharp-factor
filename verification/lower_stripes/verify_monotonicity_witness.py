#!/usr/bin/env python3
"""Rebuild the exact row-45082 counterexample to pointwise monotonicity.

For ``L(f)_q = f_q^2 - f_(q-1) f_(q+1)``, put

    C = M_(2r-1),  G = (1+z) M_(2r-2),
    D = L(G)_(r-3) L(C)_(r-2) - L(G)_(r-2) L(C)_(r-3).

The factors ``L(C)_(r-2)`` and ``L(C)_(r-3)`` are positive, so the sign of
``D`` is the sign of ``Omega_3(r) - Omega_2(r)``.  Only a fixed-width window
at the top of each wall polynomial is needed to compute these four minors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "MONOTONICITY_WITNESS.json"
SCHEMA = "shift-wall-monotonicity-witness-v1"
R_WITNESS = 45_082
WINDOW_WIDTH = 8
SMALL_R_MIN = 5
SMALL_R_MAX = 40
EXPECTED_BIT_LENGTH = 5_416_154
EXPECTED_DECIMAL_DIGITS = 1_630_425
EXPECTED_BYTE_LENGTH = 677_020
EXPECTED_SHA256 = (
    "0aed686f60e5eaf6056db6cb8e1f59aa99610dc26293b7f9f57fb53d2f8ed47e"
)


class VerificationError(RuntimeError):
    """A fail-closed witness gate did not pass."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def degree(n: int) -> int:
    return (n + 1) // 2


def top_rows(
    n_max: int,
    *,
    width: int = WINDOW_WIDTH,
    mutation: tuple[int, int, int] | None = None,
) -> tuple[list[int], list[int]]:
    """Return top-coefficient windows for ``M_(n_max-1), M_n_max``.

    Window position ``i`` stores the coefficient at ``degree(M_n)-i``.
    ``mutation`` is used only by the deterministic negative control.
    """

    require(n_max >= 1, "top-window recurrence needs n_max >= 1")
    require(width >= 5, "top-window width is too short for the witness")
    previous_two = [1] + [0] * (width - 1)
    previous = [1, 1] + [0] * (width - 2)
    for n in range(2, n_max + 1):
        offset = degree(n) - degree(n - 1)
        square = n * n
        current = [0] * width
        for index in range(width):
            value = square * previous_two[index]
            previous_index = index - offset
            if 0 <= previous_index < width:
                value += previous[previous_index]
            current[index] = value
        if mutation is not None and n == mutation[0]:
            _, index, delta = mutation
            require(0 <= index < width, "mutation index is outside the window")
            current[index] += delta
        previous_two, previous = previous, current
    return previous_two, previous


def full_rows(n_max: int) -> list[list[int]]:
    """Build ascending-coefficient wall rows for small definition checks."""

    rows = [[1], [1, 1]]
    for n in range(2, n_max + 1):
        prior = rows[n - 1]
        prior_two = rows[n - 2]
        current = [0] * max(len(prior), len(prior_two) + 1)
        for index, value in enumerate(prior):
            current[index] += value
        for index, value in enumerate(prior_two):
            current[index + 1] += n * n * value
        rows.append(current)
    return rows


def coefficient(values: list[int], index: int) -> int:
    return values[index] if 0 <= index < len(values) else 0


def turan(values: list[int], q: int) -> int:
    return (
        coefficient(values, q) ** 2
        - coefficient(values, q - 1) * coefficient(values, q + 1)
    )


def full_witness(rows: list[list[int]], r: int) -> tuple[int, tuple[int, ...]]:
    even = rows[2 * r - 2]
    shifted = [0] * (len(even) + 1)
    for index, value in enumerate(even):
        shifted[index] += value
        shifted[index + 1] += value
    odd = rows[2 * r - 1]
    lg2 = turan(shifted, r - 2)
    lg3 = turan(shifted, r - 3)
    lc2 = turan(odd, r - 2)
    lc3 = turan(odd, r - 3)
    return lg3 * lc2 - lg2 * lc3, (lg2, lg3, lc2, lc3)


def window_witness(
    r: int,
    *,
    mutation: tuple[int, int, int] | None = None,
) -> tuple[int, tuple[int, ...]]:
    even, odd = top_rows(2 * r - 1, mutation=mutation)

    def shifted(index: int) -> int:
        high = even[index - 1] if 0 <= index - 1 < len(even) else 0
        low = even[index] if 0 <= index < len(even) else 0
        return high + low

    def odd_coefficient(index: int) -> int:
        return odd[index] if 0 <= index < len(odd) else 0

    def top_turan(accessor, offset: int) -> int:
        return (
            accessor(offset) ** 2
            - accessor(offset + 1) * accessor(offset - 1)
        )

    lg2 = top_turan(shifted, 2)
    lg3 = top_turan(shifted, 3)
    lc2 = top_turan(odd_coefficient, 2)
    lc3 = top_turan(odd_coefficient, 3)
    return lg3 * lc2 - lg2 * lc3, (lg2, lg3, lc2, lc3)


def definition_checks() -> dict:
    rows = full_rows(2 * SMALL_R_MAX - 1)
    for r in range(SMALL_R_MIN, SMALL_R_MAX + 1):
        full = full_witness(rows, r)
        windowed = window_witness(r)
        require(windowed == full, f"top-window/full-definition mismatch at r={r}")

    control_r = 8
    original = window_witness(control_r)
    mutated = window_witness(control_r, mutation=(14, 2, 1))
    require(mutated != original, "top-window recurrence mutation was not rejected")
    return {
        "full_definition_range": [SMALL_R_MIN, SMALL_R_MAX],
        "full_definition_cells": SMALL_R_MAX - SMALL_R_MIN + 1,
        "top_window_matches_full_definition": True,
        "recurrence_mutation": "add 1 at n=14, top-window position 2",
        "recurrence_mutation_rejected": True,
    }


def magnitude_digest(value: int) -> tuple[int, str]:
    magnitude = abs(value)
    byte_length = (magnitude.bit_length() + 7) // 8
    payload = magnitude.to_bytes(byte_length, "big")
    return byte_length, hashlib.sha256(payload).hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def build_record() -> dict:
    definition = definition_checks()
    witness, (lg2, lg3, lc2, lc3) = window_witness(R_WITNESS)
    sign = (witness > 0) - (witness < 0)
    bit_length = abs(witness).bit_length()
    byte_length, witness_sha256 = magnitude_digest(witness)
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(2_000_000)
    decimal_digits = len(str(abs(witness)))

    require(lc2 > 0 and lc3 > 0, "the two denominator minors are not positive")
    require(sign == -1, "the claimed strict inequality did not reproduce")
    require(bit_length == EXPECTED_BIT_LENGTH, "witness bit length changed")
    require(decimal_digits == EXPECTED_DECIMAL_DIGITS,
            "witness decimal length changed")
    require(byte_length == EXPECTED_BYTE_LENGTH, "witness byte length changed")
    require(witness_sha256 == EXPECTED_SHA256, "witness digest changed")

    return {
        "schema": SCHEMA,
        "status": "PASS",
        "passed": True,
        "claim": "Omega_3(45082) < Omega_2(45082)",
        "claim_boundary": (
            "This is one exact counterexample to pointwise monotonicity in t; "
            "it is not a theorem about any other row or stripe."
        ),
        "r": R_WITNESS,
        "witness_definition": (
            "D=L(G)_(r-3)L(C)_(r-2)-L(G)_(r-2)L(C)_(r-3), "
            "G=(1+z)M_(2r-2), C=M_(2r-1), "
            "L(f)_q=f_q^2-f_(q-1)f_(q+1)"
        ),
        "sign_D": sign,
        "strict_inequality_holds": True,
        "denominator_minors_positive": {"L_C_r_minus_2": lc2 > 0,
                                         "L_C_r_minus_3": lc3 > 0},
        "component_signs": {
            "L_G_r_minus_2": (lg2 > 0) - (lg2 < 0),
            "L_G_r_minus_3": (lg3 > 0) - (lg3 < 0),
            "L_C_r_minus_2": 1,
            "L_C_r_minus_3": 1,
        },
        "window_width": WINDOW_WIDTH,
        "definition_checks": definition,
        "witness_bit_length": bit_length,
        "witness_decimal_digits": decimal_digits,
        "witness_bytes_length": byte_length,
        "digest_convention": (
            "SHA-256 of the minimal-length big-endian byte encoding of |D|"
        ),
        "witness_sha256": witness_sha256,
        "expected_values_match": {
            "bit_length": True,
            "decimal_digits": True,
            "byte_length": True,
            "sha256": True,
        },
        "execution": {
            "arithmetic": "exact integers",
            "parallelism": "one process, no pool",
            "resume": "not applicable; one fixed-width recurrence",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    try:
        record = build_record()
        atomic_json(arguments.output, record)
    except (OSError, VerificationError, ValueError) as error:
        print(f"[FAIL] {error}", file=sys.stderr)
        return 1
    print(
        f"[PASS] exact monotonicity witness at r={record['r']}; "
        f"bits={record['witness_bit_length']}; "
        f"sha256={record['witness_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
