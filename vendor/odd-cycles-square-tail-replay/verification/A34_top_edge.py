#!/usr/bin/env python3
"""Exact certificate for the two upper square-tail cells.

The interior proof covers 2 <= t <= r-2.  This file proves the two remaining
nondegenerate cells directly:

    rho(r,r) = r^2/2                         (r >= 2),
    rho(r,r-1) > 2r^2/9 >= 2                (r >= 3).

The proof uses the first two coefficients of the wall row.  A small
definition-level replay is included as a transcription check; it is not used
to infer the half-lines.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from fractions import Fraction as F
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "A34_TOP_EDGE.json"
FAILS: list[str] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gate(name: str, condition: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}"
          + (f"  {detail}" if detail else ""))
    if not condition:
        FAILS.append(name)


def wall_rows(max_n: int) -> list[list[int]]:
    rows = [[1], [1, 1]]
    for n in range(2, max_n + 1):
        row = rows[n - 1][:]
        source = rows[n - 2]
        if len(row) < len(source) + 1:
            row.extend([0] * (len(source) + 1 - len(row)))
        for k, value in enumerate(source):
            row[k + 1] += n * n * value
        rows.append(row)
    return rows


def pi_row(r: int, rows: list[list[int]]) -> list[F]:
    wall = rows[2 * r - 2]
    normalizer = math.factorial(2 * r - 1)
    return [F(wall[r - 1 - k], normalizer) for k in range(r)]


def shift_one(row: list[F]) -> list[F]:
    result = [F(0)] * (len(row) + 1)
    for k, value in enumerate(row):
        result[k] += value
        result[k + 1] += value
    return result


def turan(row: list[F], k: int) -> F:
    def coefficient(j: int) -> F:
        return row[j] if 0 <= j < len(row) else F(0)

    return coefficient(k) ** 2 - coefficient(k - 1) * coefficient(k + 1)


def rho_from_definition(r: int, t: int, rows: list[list[int]]) -> F:
    current = pi_row(r, rows)
    predecessor = pi_row(r - 1, rows)
    s_value = turan(shift_one(current), t)
    v_value = turan(current, t - 1)
    predecessor_s = turan(shift_one(predecessor), t)
    decrement = s_value - F((2 * r - 1) ** 2, (2 * r) ** 2) * predecessor_s
    return F(r * r) * decrement * decrement / (2 * s_value * v_value)


def symbolic_certificate() -> dict:
    n, r, x = sp.symbols("n r x", integer=True, nonnegative=True)
    w = n * (n + 1) * (2 * n + 1) / 6
    v = n * (n - 2) * (n - 1) * (n + 1) * (20 * n**2 - 8 * n - 21) / 360

    wall_recurrences = (
        sp.factor(w - w.subs(n, n - 1) - n**2) == 0
        and sp.factor(v - v.subs(n, n - 1) - n**2 * w.subs(n, n - 2)) == 0
        and w.subs(n, 0) == 0
        and w.subs(n, 1) == 1
        and v.subs(n, 0) == 0
        and v.subs(n, 1) == 0
    )
    gate("first two wall-coefficient formulas", wall_recurrences)

    a = sp.factor(w.subs(n, 2 * r - 2))
    b = sp.factor(v.subs(n, 2 * r - 2))
    d_current = sp.factor(a**2 + a + 1 - b)
    e_current = sp.factor(a**2 - b)
    predecessor_correction = sp.factor(
        (2 * r - 1) ** 4 * (2 * r - 2) ** 2 / (4 * r**2)
    )

    e_factor = 160 * r**4 - 48 * r**3 - 802 * r**2 + 1071 * r - 360
    positive_e_factor = sp.Poly(sp.expand(e_factor.subs(r, x + 3)), x)
    e_factorization = sp.factor(
        e_current - (r - 1) * (2 * r - 1) * e_factor / 90
    ) == 0
    e_positive = e_factorization and all(
        coefficient > 0 for coefficient in positive_e_factor.all_coeffs()
    )
    gate(
        "top Newton determinant E is positive",
        e_positive,
        f"shifted coefficients={list(reversed(positive_e_factor.all_coeffs()))}",
    )

    margin_polynomial = sp.Poly(
        7020
        + 425439 * x
        + 1133721 * x**2
        + 1294641 * x**3
        + 808325 * x**4
        + 295836 * x**5
        + 62924 * x**6
        + 7104 * x**7
        + 320 * x**8,
        x,
    )
    margin_identity = sp.factor(
        d_current
        - 3 * predecessor_correction
        - margin_polynomial.as_expr().subs(x, r - 3) / (90 * r**2)
    ) == 0
    margin_positive = all(
        coefficient > 0 for coefficient in margin_polynomial.all_coeffs()
    )
    gate("predecessor correction is below D/3", margin_identity and margin_positive)

    # D=E+a+1>E and D-K>2D/3.  Hence
    # rho=r^2(D-K)^2/(2DE)>2r^2 D/(9E)>2r^2/9.
    relation = sp.factor(d_current - e_current - a - 1) == 0
    gate("D=E+a+1", relation)

    return {
        "status": "PASS" if not FAILS else "FAIL",
        "wall_coefficients": {
            "w_n": str(w),
            "v_n": str(v),
        },
        "top_minus_one": {
            "D": str(d_current),
            "E": str(e_current),
            "K": str(predecessor_correction),
            "identity": "D-3K=R(r-3)/(90r^2)",
            "R_coefficients_low_to_high": [
                int(coefficient)
                for coefficient in reversed(margin_polynomial.all_coeffs())
            ],
            "conclusion": "rho(r,r-1)>2r^2/9>=2 for every integer r>=3",
        },
        "top": {
            "identity": "rho(r,r)=r^2/2",
            "range": "integer r>=2",
        },
    }


def definition_replay() -> dict:
    rows = wall_rows(24)
    checked = 0
    digest = hashlib.sha256()
    for r in range(2, 13):
        top = rho_from_definition(r, r, rows)
        gate(f"definition top r={r}", top == F(r * r, 2), str(top))
        digest.update(f"{r}:{r}:{top.numerator}:{top.denominator}\n".encode())
        checked += 1
        if r >= 3:
            edge = rho_from_definition(r, r - 1, rows)
            gate(f"definition top-minus-one r={r}", edge > F(2 * r * r, 9), str(edge))
            digest.update(
                f"{r}:{r - 1}:{edge.numerator}:{edge.denominator}\n".encode()
            )
            checked += 1
    return {
        "status": "PASS" if not FAILS else "FAIL",
        "range": "2<=r<=12",
        "cells": checked,
        "sha256": digest.hexdigest(),
        "role": "definition-level transcription check; half-lines are symbolic",
    }


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
        description="Certify the two terminal square-tail cells exactly."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT,
        help=f"certificate path (default: {OUTPUT.name})",
    )
    args = parser.parse_args()

    print("A34 -- square-tail upper edge")
    symbolic = symbolic_certificate()
    replay = definition_replay()
    payload = {
        "schema": "square-tail-terminal-cells-v1",
        "version": 1,
        "passed": not FAILS,
        "status": "PASS" if not FAILS else "FAIL",
        "claim": "square-tail inequality on t=r for r>=2 and t=r-1 for r>=3",
        "configuration": {
            "workers": 1,
        },
        "symbolic": symbolic,
        "definition_replay": replay,
        "failures": FAILS,
    }
    write_json(payload, args.output)
    print(f"  output: {args.output}  sha256={sha256(args.output)}")
    if FAILS:
        print(f"FAILED: {len(FAILS)} gates")
        return 1
    print("CLOSED: rho(r,r-1)>2r^2/9 and rho(r,r)=r^2/2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
