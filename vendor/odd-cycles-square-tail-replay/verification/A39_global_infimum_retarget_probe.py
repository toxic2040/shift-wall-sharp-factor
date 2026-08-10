#!/usr/bin/env python3
"""Scratch exact probe for strengthening fixed-column tail bounds.

This leaves the frozen replay packet untouched.  It reuses the live
``A1_verify_tail_independent`` construction, verifies that the reconstructed
target-1 polynomial is exactly the polynomial returned by that construction,
and then changes only the final comparison target from 1 to a supplied exact
rational number.

The source verifier does not expose its numerator and denominator envelope
polynomials.  This probe captures those local exact objects at the function's
return boundary.  That makes this a diagnostic/reconstruction artifact, not a
replacement for a clean generalized verifier in a release packet.

Each completed column is flushed to JSONL immediately.  Existing matching rows
are reused only after structural validation; the compact JSON summary is
always rebuilt from the JSONL.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "A1_verify_tail_independent.py"
RESULTS = HERE.parent / "results"
DEFAULT_JSONL = HERE / "A39_GLOBAL_INFIMUM_RETARGET.work.jsonl"
DEFAULT_OUTPUT = RESULTS / "global_infimum_tails.json"
SCHEMA = "square-tail-global-infimum-retarget-probe-v3"
ALGORITHM = "square-tail-global-infimum-retarget-components-v3"
DEFAULT_TARGET = F(203, 200)
DEFAULT_ANCHORS = {
    2: 503,
    3: 503,
    4: 503,
    5: 1_000,
}
ORIGINAL_T2_CONTROL_ANCHOR = 312


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_fraction(literal: str) -> F:
    try:
        value = F(literal)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    if value <= 0:
        raise argparse.ArgumentTypeError("target must be positive")
    return value


def parse_anchor(literal: str) -> tuple[int, int]:
    try:
        raw_t, raw_m0 = literal.split("=", 1)
        t, m0 = int(raw_t), int(raw_m0)
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError("anchor must have form t=m0") from error
    if not 2 <= t <= 8 or m0 < 2:
        raise argparse.ArgumentTypeError("require 2<=t<=8 and m0>=2")
    return t, m0


def fraction_fingerprint(value: F) -> dict[str, str | int]:
    literal = f"{value.numerator}/{value.denominator}"
    return {
        "sha256": hashlib.sha256(literal.encode("ascii")).hexdigest(),
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
    }


def polynomial_fingerprint(values: list[F]) -> str:
    digest = hashlib.sha256()
    for index, value in enumerate(values):
        if index:
            digest.update(b"\n")
        digest.update(str(value.numerator).encode("ascii"))
        digest.update(b"/")
        digest.update(str(value.denominator).encode("ascii"))
    return digest.hexdigest()


def decimal_truncation(value: F, places: int = 15) -> str:
    sign = "-" if value < 0 else ""
    numerator = abs(value.numerator)
    whole, remainder = divmod(numerator, value.denominator)
    digits: list[str] = []
    for _ in range(places):
        remainder *= 10
        digit, remainder = divmod(remainder, value.denominator)
        digits.append(str(digit))
    return f"{sign}{whole}." + "".join(digits)


def load_source_module():
    sys.path.insert(0, str(HERE))
    import A1_verify_tail_independent as source

    return source


def build_with_components(t: int, m0: int):
    source = load_source_module()
    captured: dict[str, object] = {}
    target_code = source.build_certificate.__code__

    def trace(frame, event, argument):
        if frame.f_code is target_code:
            if event == "return":
                captured.update(frame.f_locals)
            return trace
        return None

    prior_trace = sys.gettrace()
    try:
        sys.settrace(trace)
        result = source.build_certificate(m0=m0, t=t)
    finally:
        sys.settrace(prior_trace)

    required = {"u_lower", "p_upper", "q_upper", "v_index"}
    require(required <= captured.keys(), "source locals were not captured")
    return source, result, captured


def evaluate(task: tuple[int, int, int, int, str, str]) -> dict:
    sys.set_int_max_str_digits(1_000_000)
    t, m0, target_num, target_den, source_digest, producer_digest = task
    target = F(target_num, target_den)
    try:
        require(sha256(SOURCE) == source_digest, "source changed after launch")
        require(sha256(Path(__file__).resolve()) == producer_digest,
                "retarget producer changed after launch")
        source, result, captured = build_with_components(t, m0)
        u_squared = source.pmul(captured["u_lower"], captured["u_lower"])
        denominator = source.pmul(
            [F(8)],
            source.pmul(
                captured["p_upper"],
                source.pmul(
                    captured["q_upper"],
                    source.pmul(captured["v_index"], captured["v_index"]),
                ),
            ),
        )
        polynomial_length = max(len(u_squared), len(denominator))
        u_squared += [F(0)] * (polynomial_length - len(u_squared))
        denominator += [F(0)] * (polynomial_length - len(denominator))
        base = source.padd(u_squared, denominator, -1)
        retargeted = source.padd(
            u_squared,
            [target * coefficient for coefficient in denominator],
            -1,
        )
        require(len(u_squared) == len(denominator) == len(retargeted),
                "padded polynomial lengths disagree")

        denominator_nonnegative = all(value >= 0 for value in denominator)
        ratios = [
            u_squared[index] / value
            for index, value in enumerate(denominator)
            if value > 0
        ]
        require(ratios, "comparison denominator polynomial is zero")
        capacity = min(ratios)
        limiting_indices = [
            index
            for index, value in enumerate(denominator)
            if value > 0 and u_squared[index] / value == capacity
        ]
        negative_indices = [
            index for index, value in enumerate(retargeted) if value < 0
        ]
        zero_indices = [
            index for index, value in enumerate(retargeted) if value == 0
        ]
        mutation_target = capacity + F(1, capacity.denominator)
        mutated = source.padd(
            u_squared,
            [mutation_target * coefficient for coefficient in denominator],
            -1,
        )
        gates = {
            "source_digest_stable": sha256(SOURCE) == source_digest,
            "producer_digest_stable": sha256(Path(__file__).resolve())
            == producer_digest,
            "source_certificate_passed": result.get("ok") is True,
            "source_interval_ladder_replay": source.check_interval_ladder(t),
            "target_one_reconstruction_exact": base
            == result.get("certificate_poly"),
            "comparison_denominator_nonnegative": denominator_nonnegative,
            "retargeted_coefficients_strictly_positive": all(
                value > 0 for value in retargeted
            ),
            "capacity_exceeds_target": capacity > target,
            "above_capacity_mutation_rejected": any(
                coefficient < 0 for coefficient in mutated
            ),
        }
        passed = all(gates.values())
        return {
            "schema": SCHEMA,
            "passed": passed,
            "status": "PASS" if passed else "FAIL",
            "t": t,
            "anchor_m0": m0,
            "claim_start_r": m0 + 1,
            "target": f"{target.numerator}/{target.denominator}",
            "source_sha256": source_digest,
            "producer_sha256": producer_digest,
            "algorithm": ALGORITHM,
            "gates": gates,
            "degree": len(retargeted) - 1,
            "coefficient_count": len(retargeted),
            "retargeted_polynomial_sha256": polynomial_fingerprint(retargeted),
            "negative_coefficient_indices": negative_indices,
            "zero_coefficient_indices": zero_indices,
            "coefficient_capacity": fraction_fingerprint(capacity),
            "coefficient_capacity_decimal_truncated": decimal_truncation(capacity),
            "coefficient_capacity_vs_target": (
                "above" if capacity > target else "equal" if capacity == target else "below"
            ),
            "capacity_limiting_indices": limiting_indices,
            "above_capacity_mutation": fraction_fingerprint(mutation_target),
            "anchor_ratio": fraction_fingerprint(result["anchor_ratio"]),
            "anchor_ratio_decimal_truncated": decimal_truncation(
                result["anchor_ratio"]
            ),
        }
    except Exception as error:
        return {
            "schema": SCHEMA,
            "passed": False,
            "status": "ERROR",
            "t": t,
            "anchor_m0": m0,
            "claim_start_r": m0 + 1,
            "target": f"{target.numerator}/{target.denominator}",
            "source_sha256": source_digest,
            "producer_sha256": producer_digest,
            "algorithm": ALGORITHM,
            "error": repr(error),
        }


def record_key(record: dict) -> tuple[int, int, str, str, str] | None:
    values = (
        record.get("t"),
        record.get("anchor_m0"),
        record.get("target"),
        record.get("source_sha256"),
        record.get("producer_sha256"),
    )
    if (
        isinstance(values[0], int)
        and isinstance(values[1], int)
        and isinstance(values[2], str)
        and isinstance(values[3], str)
        and isinstance(values[4], str)
    ):
        return values
    return None


def record_complete(record: dict) -> bool:
    return (
        record.get("schema") == SCHEMA
        and record.get("algorithm") == ALGORITHM
        and record.get("status") in {"PASS", "FAIL"}
        and isinstance(record.get("passed"), bool)
        and record_key(record) is not None
        and isinstance(record.get("gates"), dict)
        and isinstance(record.get("degree"), int)
        and isinstance(record.get("coefficient_count"), int)
        and isinstance(record.get("retargeted_polynomial_sha256"), str)
        and len(record["retargeted_polynomial_sha256"]) == 64
    )


def load_records(path: Path) -> dict[tuple[int, int, str, str, str], dict]:
    records: dict[tuple[int, int, str, str, str], dict] = {}
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
            if record.get("schema") != SCHEMA:
                continue
            key = record_key(record)
            require(key is not None, f"invalid key at JSONL line {line_number}")
            records[key] = record
    return records


def write_summary(path: Path, payload: dict) -> None:
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
    parser.add_argument("--target", type=parse_fraction, default=DEFAULT_TARGET)
    parser.add_argument(
        "--anchor",
        action="append",
        type=parse_anchor,
        help="override one or more defaults with t=m0",
    )
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(args.jsonl.resolve() != args.output.resolve(),
            "JSONL and summary paths must differ")

    anchors = dict(DEFAULT_ANCHORS)
    if args.anchor:
        anchors.update(args.anchor)
    source_digest = sha256(SOURCE)
    producer_digest = sha256(Path(__file__).resolve())
    target_literal = f"{args.target.numerator}/{args.target.denominator}"
    main_tasks = [
        (
            t,
            anchors[t],
            args.target.numerator,
            args.target.denominator,
            source_digest,
            producer_digest,
        )
        for t in sorted(anchors)
    ]
    control_task = (
        2,
        ORIGINAL_T2_CONTROL_ANCHOR,
        args.target.numerator,
        args.target.denominator,
        source_digest,
        producer_digest,
    )
    tasks = main_tasks + [control_task]
    main_keys = {
        (t, m0, target_literal, source_digest, producer_digest)
        for t, m0, _, _, _, _ in main_tasks
    }
    control_key = (
        2,
        ORIGINAL_T2_CONTROL_ANCHOR,
        target_literal,
        source_digest,
        producer_digest,
    )
    expected_keys = main_keys | {control_key}

    existing = load_records(args.jsonl)
    todo = [
        task
        for task in tasks
        if not record_complete(
            existing.get(
                (task[0], task[1], target_literal, source_digest, producer_digest),
                {},
            )
        )
    ]
    reused = len(tasks) - len(todo)
    args.jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.jsonl.open("a", encoding="utf-8", buffering=1) as sink:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(evaluate, task): task for task in todo}
            for future in as_completed(futures):
                task = futures[future]
                try:
                    record = future.result()
                except Exception as error:
                    record = {
                        "schema": SCHEMA,
                        "passed": False,
                        "status": "ERROR",
                        "t": task[0],
                        "anchor_m0": task[1],
                        "claim_start_r": task[1] + 1,
                        "target": target_literal,
                        "source_sha256": source_digest,
                        "producer_sha256": producer_digest,
                        "algorithm": ALGORITHM,
                        "error": repr(error),
                    }
                sink.write(json.dumps(record, sort_keys=True) + "\n")
                sink.flush()

    records = load_records(args.jsonl)
    selected = {key: records.get(key) for key in expected_keys}
    selected_main = {key: selected[key] for key in main_keys}
    control = selected[control_key]
    missing = [list(key[:3]) for key, value in selected.items() if value is None]
    errors = [
        [*key[:3], value.get("error")]
        for key, value in selected.items()
        if value is not None and value.get("status") == "ERROR"
    ]
    columns = {
        str(key[0]): value
        for key, value in sorted(selected_main.items())
        if value is not None
    }
    control_gates = control.get("gates", {}) if isinstance(control, dict) else {}
    control_passed = (
        isinstance(control, dict)
        and control.get("status") == "FAIL"
        and control.get("passed") is False
        and control.get("coefficient_capacity_vs_target") == "below"
        and control.get("negative_coefficient_indices") == [0, 1]
        and control_gates.get("source_digest_stable") is True
        and control_gates.get("producer_digest_stable") is True
        and control_gates.get("source_certificate_passed") is True
        and control_gates.get("source_interval_ladder_replay") is True
        and control_gates.get("target_one_reconstruction_exact") is True
        and control_gates.get("comparison_denominator_nonnegative") is True
        and control_gates.get("retargeted_coefficients_strictly_positive") is False
        and control_gates.get("capacity_exceeds_target") is False
        and control_gates.get("above_capacity_mutation_rejected") is True
    )
    passed = (
        not missing
        and not errors
        and len(columns) == len(anchors)
        and all(value.get("passed") is True for value in columns.values())
        and control_passed
    )
    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "scope": "four unbounded fixed-column analytic tails t=2..5 only; finite windows and joint/top regions are separate",
        "target": target_literal,
        "source_file": SOURCE.name,
        "source_sha256": source_digest,
        "producer_sha256": producer_digest,
        "missing": missing,
        "errors": errors,
        "columns": columns,
        "controls": {
            "original_t2_anchor_rejected": control_passed,
            "record": control,
        },
    }
    write_summary(args.output, payload)

    for t in sorted(columns, key=int):
        record = columns[t]
        print(
            f"[{'PASS' if record['passed'] else 'FAIL'}] t={t} "
            f"m0={record['anchor_m0']} capacity="
            f"{record.get('coefficient_capacity_decimal_truncated', 'error')} "
            f"limiting={record.get('capacity_limiting_indices', [])}"
        )
    print(
        f"[{'PASS' if control_passed else 'FAIL'}] control t=2 "
        f"m0={ORIGINAL_T2_CONTROL_ANCHOR} rejected at 203/200"
    )
    print(f"[INFO] workers: {args.workers}; resumed rows: {reused}")
    print(f"work JSONL: {args.jsonl}")
    print(f"summary: {args.output}")
    print("ALL RETARGETED TAILS PASS" if passed else "RETARGET PROBE INCOMPLETE OR FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
