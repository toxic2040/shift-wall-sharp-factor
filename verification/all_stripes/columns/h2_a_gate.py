#!/usr/bin/env python3
"""Certified outward-rounded ladder for the Wallis variable a_r = 6 E_1(r-1).

Recurrence (A2_stable.py, rewrite (G), all terms positive, no cancellation):
    R_1 = 5,  R_m = 1 + (2m/(2m-1))^2 R_{m-1}
    E_1(m) = E_1(m-1) + D_1(m-1)/R_m
    D_1(m) = D_1(m-1)(1 - 1/R_m) + R_m/(2m+1)^2        [E_0(m) = 1 for all m]
    a_r = 6 E_1(r-1)

Interval arithmetic is the deposited A1 fixed-point Iv class (256 bits,
outward rounding), so every recorded bound is certified.  D_1(m) > 0 at every
step certifies that a_r is strictly increasing.

Gates (prefix-free t=9 and t=10 closures; TOP region for column t is
a_r <= (t-2)^2):
    a_68133  <= 49 = (9-2)^2
    a_200000 <= 64 = (10-2)^2

Crosschecks: exact_a(r) from A2_stable.py must land inside the interval at
r in {60, 200, 500}; the recorded exact crossings a=25@1468, 30@3634,
35@8414, 36@9884, 40@18457 must be bracketed.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parents[2]
A1_PATH = (RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay"
           / "verification" / "A1_verify_tail_independent.py")
A2_PATH = (RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay"
           / "verification" / "A2_stable.py")
OUTPUT = HERE / "H2_A_GATE.json"
SCHEMA = "h2-wallis-a-gate-v1"
R_MAX = 200_001
GATES = {55_000: 49, 68_133: 49, 100_000: 64, 200_000: 64}
CROSSINGS = {25: 1_468, 30: 3_634, 35: 8_414, 36: 9_884, 40: 18_457}
LOCATE_THRESHOLDS = [45, 49, 50, 64]
EXTRA_CHECKPOINTS = [1_001, 2_004, 4_007, 6_200, 12_399, 20_000, 27_000]
EXACT_CROSSCHECK_ROWS = [60, 200, 500]


def load_module(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decimal_truncation(value: F, places: int) -> str:
    sign = "-" if value < 0 else ""
    whole, remainder = divmod(abs(value.numerator), value.denominator)
    digits: list[str] = []
    for _ in range(places):
        remainder *= 10
        digit, remainder = divmod(remainder, value.denominator)
        digits.append(str(digit))
    return f"{sign}{whole}." + "".join(digits)


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def main() -> int:
    sys.set_int_max_str_digits(1_000_000)
    started = time.time()
    a1 = load_module(A1_PATH, "h2_agate_iv_source")
    a2 = load_module(A2_PATH, "h2_agate_a2_stable")
    Iv, ONE = a1.Iv, a1.ONE

    wanted_rows: set[int] = set(GATES) | set(EXTRA_CHECKPOINTS)
    wanted_rows |= set(EXACT_CROSSCHECK_ROWS)
    for crossing_r in CROSSINGS.values():
        wanted_rows |= {crossing_r - 1, crossing_r}
    wanted_m = {row - 1 for row in wanted_rows}

    R = Iv.exact(5)
    E1 = Iv.exact(0)
    D1 = Iv.exact(1)
    # m = 0 state: E_1(0) = 0, D_1(0) = 1 (delta_{k1}).
    monotone_failures: list[int] = []
    bounds: dict[int, tuple[F, F]] = {}
    located: dict[int, dict] = {}
    previous_upper: F | None = None
    for m in range(1, R_MAX):
        if m == 1:
            Rm = Iv.exact(5)
        else:
            Rm = ONE + R.scal(4 * m * m, (2 * m - 1) ** 2)
        R = Rm
        inv_R = Rm.inv()
        xm = Rm.scal(1, (2 * m + 1) ** 2)
        if D1.lo <= 0:
            monotone_failures.append(m)
        E1 = E1 + D1 * inv_R
        D1 = D1 * (ONE - inv_R) + xm
        a = E1.scal(6)
        low, high = a.flo(), a.fhi()
        if m in wanted_m:
            bounds[m + 1] = (low, high)
        for threshold in LOCATE_THRESHOLDS:
            if threshold not in located and low >= threshold:
                located[threshold] = {
                    "threshold": threshold,
                    "first_r_at_or_above": m + 1,
                    "previous_upper_below": previous_upper is not None
                    and previous_upper < threshold,
                    "a_previous_upper": decimal_truncation(
                        previous_upper, 12
                    ) if previous_upper is not None else None,
                    "a_at_lower": decimal_truncation(low, 12),
                    "max_gate_anchor": m,
                }
        previous_upper = high

    failures: list[str] = []
    if monotone_failures:
        failures.append(
            f"D_1 lower bound nonpositive at m={monotone_failures[:5]}"
        )

    exact_checks = []
    for row in EXACT_CROSSCHECK_ROWS:
        exact_value = a2.exact_a(row)
        low, high = bounds[row]
        ok = low <= exact_value <= high
        if not ok:
            failures.append(f"exact_a({row}) escaped the interval")
        exact_checks.append({
            "r": row,
            "exact_decimal_20": decimal_truncation(exact_value, 20),
            "inside_interval": ok,
        })

    crossing_records = []
    for threshold, row in CROSSINGS.items():
        below_high = bounds[row - 1][1]
        at_low = bounds[row][0]
        ok = below_high < threshold <= at_low
        if not ok:
            failures.append(f"crossing a={threshold}@r={row} not bracketed")
        crossing_records.append({
            "threshold": threshold,
            "first_r_at_or_above": row,
            "a_below_upper": decimal_truncation(below_high, 12),
            "a_at_lower": decimal_truncation(at_low, 12),
            "bracketed": ok,
        })

    gate_records = {}
    for row, limit in GATES.items():
        low, high = bounds[row]
        ok = high <= limit
        if not ok:
            failures.append(f"gate a_{row} <= {limit} failed")
        gate_records[str(row)] = {
            "limit": limit,
            "a_lower": decimal_truncation(low, 18),
            "a_upper": decimal_truncation(high, 18),
            "certified": ok,
            "geometry": f"a_r <= {limit} = (t-2)^2 puts (r,t) in the TOP "
                        "region for the matching column",
        }

    for threshold in LOCATE_THRESHOLDS:
        record = located.get(threshold)
        if record is None:
            located[threshold] = {
                "threshold": threshold,
                "first_r_at_or_above": None,
                "note": f"not reached: a_r stays below {threshold} "
                        f"through r={R_MAX}",
            }
        elif record["previous_upper_below"] is not True:
            failures.append(
                f"crossing location for a={threshold} not certified"
            )

    checkpoints = {
        str(row): {
            "a_lower": decimal_truncation(bounds[row][0], 15),
            "a_upper": decimal_truncation(bounds[row][1], 15),
        }
        for row in sorted(wanted_rows)
    }
    top_row = max(bounds)
    width_at_top = bounds[top_row][1] - bounds[top_row][0]

    passed = not failures
    payload = {
        "schema": SCHEMA,
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "claim": (
            "a_r = 6 E_1(r-1) is strictly increasing; "
            "a_68133 <= 49 and a_200000 <= 64, certified outward"
        ),
        "method": (
            "256-bit fixed-point outward-rounded interval ladder over the "
            "cancellation-free (G) recurrence of A2_stable.py, using the "
            "deposited A1 Iv class; strict increase certified by D_1 > 0 "
            "at every step"
        ),
        "configuration": {
            "r_max": R_MAX,
            "interval_bits": a1.PREC,
            "workers": 1,
        },
        "sources": {
            A1_PATH.relative_to(RELEASE_ROOT).as_posix(): sha256_file(A1_PATH),
            A2_PATH.relative_to(RELEASE_ROOT).as_posix(): sha256_file(A2_PATH),
            "producer": sha256_file(Path(__file__).resolve()),
        },
        "monotonicity": {
            "certified_strictly_increasing": not monotone_failures,
            "d1_positive_steps": R_MAX - 1 - len(monotone_failures),
        },
        "gates": gate_records,
        "crossings": crossing_records,
        "threshold_locations": {
            str(threshold): located[threshold]
            for threshold in LOCATE_THRESHOLDS
        },
        "exact_crosschecks": exact_checks,
        "checkpoints": checkpoints,
        "interval_width_at_r_max": decimal_truncation(width_at_top, 60),
        "failures": failures,
    }
    write_json(OUTPUT, payload)
    print(f"[{'PASS' if passed else 'FAIL'}] a-gate ladder to r={R_MAX}: "
          f"a_68133<={GATES[68_133]}: {gate_records['68133']['certified']}, "
          f"a_200000<={GATES[200_000]}: "
          f"{gate_records['200000']['certified']}, "
          f"runtime={round(time.time() - started, 1)}s")
    for failure in failures:
        print(f"  failure: {failure}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
