#!/usr/bin/env python3
"""Prefix-free column join: retargeted tail + certified Wallis a-gate.

For column t with tail anchor m0 and TOP-region bound L = (t-2)^2:
  tail   H2_TAIL_T{t}_M{m0}.json proves rho(r,t)>4/3 for r >= m0+1;
  a-gate H2_A_GATE.json certifies a_{m0} <= L and a_r strictly increasing,
         so every 504 <= r <= m0 has a_r <= L, i.e. (r,t) lies in the TOP
         region, which is closed at 4/3 elsewhere (not in this lane);
  cells with r < 504 lie in the finite core, closed elsewhere.
Writes H2_COL_T{t}.json (atomic).

Usage: h2_join_prefix_free.py t m0
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
SCHEMA = "h2-fixed-column-four-thirds-prefix-free-join-v1"
FAILURES: list[str] = []
GATES: dict[str, bool] = {}


def gate(name: str, condition: bool) -> bool:
    passed = condition is True
    GATES[name] = passed
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
    if not passed:
        FAILURES.append(name)
    return passed


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
    args = parser.parse_args()
    t, m0 = args.t, args.m0
    limit = (t - 2) ** 2

    print(f"H2 prefix-free join, column t={t}, anchor m0={m0}, "
          f"TOP bound a<={limit}")
    tail_path = HERE / f"H2_TAIL_T{t}_M{m0}.json"
    gate_path = HERE / "H2_A_GATE.json"
    tail = json.loads(tail_path.read_text(encoding="utf-8"))
    a_gate = json.loads(gate_path.read_text(encoding="utf-8"))

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

    gate("a-gate passed", a_gate.get("passed") is True)
    record = a_gate.get("gates", {}).get(str(m0), {})
    gate("a-gate row matches anchor", bool(record))
    gate("a-gate limit matches (t-2)^2", record.get("limit") == limit)
    gate("a-gate certified", record.get("certified") is True)
    gate("a-gate monotonicity certified", a_gate.get(
        "monotonicity", {}).get("certified_strictly_increasing") is True)
    gate("a-gate crossings bracketed", all(
        item.get("bracketed") is True
        for item in a_gate.get("crossings", [])))
    gate("a-gate exact crosschecks", all(
        item.get("inside_interval") is True
        for item in a_gate.get("exact_crosschecks", [])))

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
            "base_gate_failures": [
                name for name, ok in record.get(
                    "base_certificate_gates", {}).items() if not ok],
            "file": path.name,
            "sha256": sha256_file(path),
        })
    if attempts:
        gate("final attempt is the passing one",
             attempts[-1]["anchor_m0"] == m0
             and attempts[-1]["status"] == "PASS"
             and all(item["status"] == "FAIL" for item in attempts[:-1]))

    passed = not FAILURES
    payload = {
        "schema": SCHEMA,
        "passed": passed,
        "status": "TAIL_PLUS_TOP_GATE_PROVED" if passed else "FAIL",
        "claim": (
            f"rho(j,{t})>4/3 for every integer j>={m0 + 1} (tail leg, "
            f"proved here); for 504<=j<={m0}, a_j<={limit}=(t-2)^2 is "
            "certified here, placing every such cell in the TOP region "
            "whose 4/3 closure is outside this lane; cells j<504 belong "
            "to the finite core, also outside this lane"
        ),
        "column": t,
        "coverage": {
            "tail": f"r>={m0 + 1}",
            "top_gate": f"a_r<={limit} certified for all r<={m0} "
                        "(monotone + endpoint)",
            "external_dependencies": [
                "TOP-region 4/3 closure for a_r<=(t-2)^2",
                "finite core r<504",
            ],
        },
        "tail": {
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
        },
        "a_gate": {
            "file": gate_path.name,
            "sha256": sha256_file(gate_path),
            "row": m0,
            "limit": limit,
            "bounds": record,
            "interval_width_at_r_max":
                a_gate.get("interval_width_at_r_max"),
        },
        "tail_attempts": attempts,
        "gates": GATES,
        "failures": FAILURES,
    }
    output = HERE / f"H2_COL_T{t}.json"
    write_json(output, payload)
    print(f"[{payload['status']}] column t={t}: tail r>={m0 + 1} + "
          f"a_{m0}<={limit} -> {output.name}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
