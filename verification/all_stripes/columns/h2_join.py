#!/usr/bin/env python3
"""Fail-closed join of one fixed column: exact prefix + retargeted tail.

Full mode (t with a passing tail anchor m0):
  prefix H2_PREFIX_T{t}.json     proves rho(r,t)>4/3 for t <= r <= m0,
  tail   H2_TAIL_T{t}_M{m0}.json proves rho(r,t)>4/3 for r >= m0+1,
  the overlap cell r=m0+1 is verified exactly on both sides.

Prefix-only mode (--prefix-only): closes just the prefix leg [t, m0]; rows
above m0 are then covered by the joint-region certificate rather than by a
fixed-column tail.  Recorded tail attempts are kept as capacity-slope data.

Writes H2_COL_T{t}.json (atomic).

Usage: h2_join.py t m0 [--attempts m0a,m0b,...] [--prefix-only]
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
SCHEMA = "h2-fixed-column-four-thirds-join-v1"
FAILURES: list[str] = []
GATES: dict[str, bool] = {}


def gate(name: str, condition: bool) -> bool:
    passed = condition is True
    GATES[name] = passed
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
    if not passed:
        FAILURES.append(name)
    return passed


def positive_margin(text: str) -> bool:
    return (not text.startswith("-")
            and any(digit != "0" for digit in text if digit.isdigit()))


def decimal_to_fraction(text: str) -> F:
    whole, _, fractional = text.partition(".")
    return F(int(whole + fractional), 10 ** len(fractional))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    parser = argparse.ArgumentParser()
    parser.add_argument("t", type=int)
    parser.add_argument("m0", type=int)
    parser.add_argument("--attempts", type=str, default="")
    parser.add_argument("--prefix-only", action="store_true")
    args = parser.parse_args()
    t, m0 = args.t, args.m0

    print(f"H2 join, column t={t}, "
          f"{'prefix leg only, R0 anchor' if args.prefix_only else 'final anchor'} "
          f"m0={m0}")
    prefix_path = HERE / f"H2_PREFIX_T{t}.json"
    prefix = json.loads(prefix_path.read_text(encoding="utf-8"))

    gate("prefix passed", prefix.get("passed") is True
         and prefix.get("status") == "PASS")
    gate("prefix schema", prefix.get("schema")
         == "h2-fixed-column-four-thirds-prefix-v1")
    gate("prefix column and range", prefix.get("t") == t
         and prefix.get("m0") == m0 and prefix.get("range") == [t, m0])
    gate("prefix cell count", prefix.get("cell_count") == m0 - t + 1)
    gate("prefix strict", prefix.get("strict_failures") == []
         and prefix.get("problems") == [])
    gate("prefix definition crosscheck",
         prefix.get("definition_crosscheck", {}).get("failures") == []
         and prefix.get("definition_crosscheck", {}).get("cells") == 61 - t)
    minimum = prefix.get("finite_minimum", {})
    gate("prefix minimum above target",
         positive_margin(minimum.get("margin_over_four_thirds_24", "-")))

    overlap = prefix.get("overlap_row", {})
    gate("overlap cell is m0 plus one", overlap.get("r") == m0 + 1)
    gate("overlap cell exact comparison", overlap.get("comparison") == ">")
    gate("overlap margin positive",
         positive_margin(overlap.get("margin_over_four_thirds_24", "-")))

    tail_path = HERE / f"H2_TAIL_T{t}_M{m0}.json"
    tail_summary = None
    if not args.prefix_only:
        tail = json.loads(tail_path.read_text(encoding="utf-8"))
        gate("tail passed", tail.get("passed") is True
             and tail.get("status") == "PASS")
        gate("tail schema", tail.get("schema")
             == "h2-fixed-column-four-thirds-tail-v1")
        configuration = tail.get("configuration", {})
        gate("tail column and anchor", configuration.get("t") == t
             and configuration.get("anchor_m0") == m0
             and configuration.get("claim_start_r") == m0 + 1
             and configuration.get("target") == "4/3")
        gate("tail base gates", tail.get("base_certificate_passed") is True
             and all(tail.get("base_certificate_gates", {}).values())
             and tail.get("source_interval_ladder_crosscheck") is True)
        gate("tail reconstruction and denominator",
             tail.get("target_one_reconstruction_exact") is True
             and tail.get("comparison_denominator_nonnegative") is True)
        gate("tail strictly positive",
             tail.get("retargeted_negative_indices") == []
             and tail.get("retargeted_zero_indices") == []
             and tail.get("retargeted_positive_coefficients")
             == tail.get("retargeted_degree") + 1)
        capacity_text = tail.get("coefficient_capacity", {}).get(
            "decimal_truncated", "0.0")
        gate("tail capacity above four thirds",
             decimal_to_fraction(capacity_text) > F(4, 3))
        gate("tail pi repair", tail.get("a1_pi_interval_repair", {}).get(
            "corrected_reference_tripwire_passed") is True)
        gate("tail mutation control", bool(
            tail.get("mutation_control", {})
            and tail["mutation_control"]["first_failing_negative_indices"]))
        three_halves = tail.get("three_halves_control", {})
        gate("three halves control recorded",
             three_halves.get("consistent_with_capacity") is True)
        tail_summary = {
            "file": tail_path.name,
            "sha256": sha256_file(tail_path),
            "anchor_m0": m0,
            "retargeted_degree": tail.get("retargeted_degree"),
            "retargeted_sha256": tail.get("retargeted_sha256"),
            "coefficient_capacity": tail.get("coefficient_capacity"),
            "capacity_limiting_indices":
                tail.get("capacity_limiting_indices"),
            "mutation_control": tail.get("mutation_control"),
            "three_halves_control_strictly_positive":
                three_halves.get("strictly_positive"),
        }

    attempts = []
    for item in args.attempts.split(","):
        if not item:
            continue
        attempt_m0 = int(item)
        path = HERE / f"H2_TAIL_T{t}_M{attempt_m0}.json"
        record = json.loads(path.read_text(encoding="utf-8"))
        attempts.append({
            "anchor_m0": attempt_m0,
            "status": record.get("status"),
            "capacity_decimal": record.get(
                "coefficient_capacity", {}).get("decimal_truncated"),
            "anchor_ratio_decimal": record.get(
                "anchor_ratio", {}).get("decimal_truncated"),
            "retargeted_negative_indices": record.get(
                "retargeted_negative_indices"),
            "file": path.name,
            "sha256": sha256_file(path),
        })
    if attempts:
        if args.prefix_only:
            gate("all recorded tail attempts failed",
                 all(item["status"] == "FAIL" for item in attempts))
        else:
            gate("final attempt is the passing one",
                 attempts[-1]["anchor_m0"] == m0
                 and attempts[-1]["status"] == "PASS"
                 and all(item["status"] == "FAIL"
                         for item in attempts[:-1]))

    passed = not FAILURES
    if args.prefix_only:
        claim = (f"rho(j,{t})>4/3 for every integer j in [{t},{m0}]; "
                 f"j>={m0 + 1} covered by the joint-region certificate")
        status = "PREFIX_LEG_PROVED" if passed else "FAIL"
        coverage = {
            "prefix": [t, m0],
            "first_row_above_prefix": m0 + 1,
            "tail": "covered by the joint-region certificate",
        }
    else:
        claim = f"rho(j,{t})>4/3 for every integer j>={t}"
        status = "PASS" if passed else "FAIL"
        coverage = {
            "prefix": [t, m0],
            "tail": f"r>={m0 + 1}",
            "overlap_cell": m0 + 1,
        }
    payload = {
        "schema": SCHEMA,
        "passed": passed,
        "status": status,
        "claim": claim,
        "column": t,
        "coverage": coverage,
        "prefix": {
            "file": prefix_path.name,
            "sha256": sha256_file(prefix_path),
            "cell_count": prefix.get("cell_count"),
            "comparison_sha256": prefix.get("comparison_sha256"),
            "finite_minimum_r": prefix.get("finite_minimum_r"),
            "finite_minimum": prefix.get("finite_minimum"),
            "minimum_at_right_endpoint":
                prefix.get("minimum_at_right_endpoint"),
        },
        "tail": tail_summary,
        "overlap_row": overlap,
        "tail_attempts": attempts,
        "gates": GATES,
        "failures": FAILURES,
    }
    output = HERE / f"H2_COL_T{t}.json"
    write_json(output, payload)
    print(f"[{status}] column t={t}: prefix [{t},{m0}]"
          + ("" if args.prefix_only else f" + tail r>={m0 + 1}")
          + f" -> {output.name}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
