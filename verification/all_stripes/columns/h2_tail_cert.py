#!/usr/bin/env python3
"""Fixed-column tail certificate at target 4/3 for one (t, m0) attempt.

Same construction as the t=4 transfer certificate in
verification/lower_stripes/verify_sq_t4_transfer.py: load the deposited
corrected-independent builder, apply the mandatory one-ulp outward PI repair,
capture u_lower/p_upper/q_upper/v_index from build_certificate(m0, t), and
test strict coefficientwise positivity of

    H = u_lower^2 - target * 8 * p_upper * q_upper * v_index^2.

All-positive H proves rho(r,t) > target for every integer r >= m0+1.
Writes one JSON per attempt (atomic), pass or fail.
"""

from __future__ import annotations

import argparse
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
REPO = RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay"
INDEPENDENT_PATH = REPO / "verification" / "A1_verify_tail_independent.py"
COMMON_PATH = (RELEASE_ROOT / "verification"
               / "lower_stripes"
               / "verify_bridge_t2_minimum.py")
SCHEMA = "h2-fixed-column-four-thirds-tail-v1"
TARGET = F(4, 3)


def load_module(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def indexed_polynomial_fingerprint(common, values: list[F]) -> str:
    digest = hashlib.sha256()
    for index, value in enumerate(values):
        digest.update(
            f"{index}:{common.fraction_literal(value)}\n".encode("ascii")
        )
    return digest.hexdigest()


def repair_independent_pi(common, source) -> dict:
    """Verbatim from verify_sq_t4_transfer.py (mandatory one-ulp repair)."""
    scale = source.SCALE
    lo = source.PI_LO.numerator * scale // source.PI_LO.denominator
    upper_numerator = source.PI_HI.numerator * scale
    bad_hi = upper_numerator // source.PI_HI.denominator
    good_hi = -((-upper_numerator) // source.PI_HI.denominator)
    common.require(lo == bad_hi and good_hi == bad_hi + 1,
                   "independent PI defect shape changed")
    bad_half = (lo // 2, -((-bad_hi) // 2))
    good_half = (lo // 2, -((-good_hi) // 2))
    common.require(bad_half == good_half,
                   "PI repair survives the half scaling")
    reference_denominator = 10**100
    reference_numerator = int(common.PI_100_DECIMAL_DIGITS)
    reference_lo = F(reference_numerator, reference_denominator)
    reference_hi = F(reference_numerator + 1, reference_denominator)
    bad = source.Iv(lo, bad_hi)
    good = source.Iv(lo, good_hi)
    common.require(not (bad.flo() < reference_lo and bad.fhi() > reference_hi),
                   "bad independent PI passed the reference tripwire")
    common.require(good.flo() < reference_lo and good.fhi() > reference_hi,
                   "corrected independent PI failed the reference tripwire")
    source.PI = good
    return {
        "defective_width_ulps": 0,
        "corrected_width_ulps": 1,
        "defective_reference_tripwire_rejected": True,
        "corrected_reference_tripwire_passed": True,
        "half_pi_exactly_invariant": True,
    }


def build_with_components(source, m0: int, t: int):
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
    if not required <= captured.keys():
        raise RuntimeError("independent tail components were not captured")
    return result, captured


def coefficientwise_verdict(u2: list[F], den: list[F], target: F) -> dict:
    values = [u2[i] - target * den[i] for i in range(len(u2))]
    negative = [i for i, v in enumerate(values) if v < 0]
    zero = [i for i, v in enumerate(values) if v == 0]
    return {
        "strictly_positive": not negative and not zero,
        "negative_indices": negative,
        "zero_indices": zero,
    }


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
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    t, m0 = args.t, args.m0
    output = args.output or HERE / f"H2_TAIL_T{t}_M{m0}.json"

    started = time.time()
    common = load_module(COMMON_PATH, f"h2_common_t{t}_{m0}")
    source = load_module(INDEPENDENT_PATH, f"h2_a1_independent_t{t}_{m0}")
    independent_sha = common.sha256_file(INDEPENDENT_PATH)

    pi_repair = repair_independent_pi(common, source)
    ladder_ok = source.check_interval_ladder(t)

    result, captured = build_with_components(source, m0, t)
    base_ok = result["ok"] is True

    payload: dict = {
        "schema": SCHEMA,
        "claim": f"rho(r,{t})>4/3 for every integer r>={m0 + 1}",
        "configuration": {
            "t": t,
            "anchor_m0": m0,
            "claim_start_r": m0 + 1,
            "interval_bits": source.PREC,
            "kmax": t + 1,
            "workers": 1,
            "target": "4/3",
        },
        "source": {
            "independent_path": INDEPENDENT_PATH.relative_to(
                RELEASE_ROOT).as_posix(),
            "independent_sha256": independent_sha,
            "common_verifier_sha256": common.sha256_file(COMMON_PATH),
            "producer_sha256": common.sha256_file(Path(__file__).resolve()),
        },
        "a1_pi_interval_repair": pi_repair,
        "source_interval_ladder_crosscheck": ladder_ok,
        "base_certificate_gates": {k: bool(v) for k, v in result["gates"].items()},
        "base_certificate_passed": base_ok,
        "tail_log_moment_max_degree": result["tail_degree"],
        "anchor_ratio": {
            **common.fraction_fingerprint(result["anchor_ratio"]),
            "decimal_truncated": common.decimal_truncation(
                result["anchor_ratio"], 24
            ),
        },
        "anchor_width": common.fraction_fingerprint(result["anchor_width"]),
    }

    if not base_ok:
        payload.update({"passed": False, "status": "FAIL",
                        "failure": "base target-one certificate failed"})
        write_json(output, payload)
        print(f"[FAIL] t={t} m0={m0}: base certificate failed")
        return 1

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
    length = max(len(u_squared), len(denominator))
    u_squared += [F(0)] * (length - len(u_squared))
    denominator += [F(0)] * (length - len(denominator))
    reconstructed_one = source.padd(u_squared, denominator, -1)
    common.require(reconstructed_one == result["certificate_poly"],
                   "target-one tail reconstruction mismatch")
    common.require(all(value >= 0 for value in denominator),
                   "tail comparison denominator polynomial is negative")

    h_poly = source.padd(
        u_squared, [TARGET * value for value in denominator], -1
    )
    verdict = coefficientwise_verdict(u_squared, denominator, TARGET)
    passed = verdict["strictly_positive"]

    capacities = [
        (u_squared[index] / value, index)
        for index, value in enumerate(denominator)
        if value > 0
    ]
    capacity = min(value for value, _ in capacities)
    limiting_indices = [
        index for value, index in capacities if value == capacity
    ]
    zero_denominator_indices = [
        index for index, value in enumerate(denominator) if value == 0
    ]

    # Negative control: target 3/2, recorded with capacity consistency.
    three_halves = coefficientwise_verdict(u_squared, denominator, F(3, 2))
    three_halves["consistent_with_capacity"] = (
        three_halves["strictly_positive"] == (capacity > F(3, 2))
    )

    # Mutation control: first failing target on the 1/100 grid above 4/3,
    # demonstrated by direct coefficientwise evaluation on both sides.
    mutation = None
    if passed:
        steps = (capacity - TARGET) * 100
        k_fail = int(steps) + 1 if steps == int(steps) else int(steps) + 1
        target_fail = TARGET + F(k_fail, 100)
        target_prev = TARGET + F(k_fail - 1, 100)
        fail_verdict = coefficientwise_verdict(
            u_squared, denominator, target_fail
        )
        prev_verdict = coefficientwise_verdict(
            u_squared, denominator, target_prev
        )
        common.require(bool(fail_verdict["negative_indices"]),
                       "mutation grid target did not go negative")
        common.require(not prev_verdict["negative_indices"],
                       "pre-mutation grid target already negative")
        mutation = {
            "grid_step": "1/100",
            "first_failing_target": common.fraction_literal(target_fail),
            "first_failing_target_decimal": common.decimal_truncation(
                target_fail, 6
            ),
            "first_failing_negative_indices":
                fail_verdict["negative_indices"],
            "previous_grid_target": common.fraction_literal(target_prev),
            "previous_grid_target_no_negative": True,
        }

    h_sha = indexed_polynomial_fingerprint(common, h_poly)
    payload.update({
        "passed": bool(passed),
        "status": "PASS" if passed else "FAIL",
        "target_one_reconstruction_exact": True,
        "comparison_denominator_nonnegative": True,
        "retargeted_degree": len(h_poly) - 1,
        "retargeted_coefficient_count": len(h_poly),
        "retargeted_positive_coefficients": sum(v > 0 for v in h_poly),
        "retargeted_negative_indices": verdict["negative_indices"],
        "retargeted_zero_indices": verdict["zero_indices"],
        "retargeted_sha256": h_sha,
        "coefficient_fingerprint_encoding":
            "ASCII index:numerator/denominator followed by newline",
        "coefficient_capacity": {
            **common.fraction_fingerprint(capacity),
            "decimal_truncated": common.decimal_truncation(capacity, 24),
        },
        "capacity_limiting_indices": limiting_indices,
        "zero_denominator_indices": zero_denominator_indices,
        "three_halves_control": three_halves,
        "mutation_control": mutation,
    })
    write_json(output, payload)
    capacity_decimal = common.decimal_truncation(capacity, 12)
    print(
        f"[{'PASS' if passed else 'FAIL'}] t={t} m0={m0}: "
        f"degree={len(h_poly) - 1}, capacity={capacity_decimal}, "
        f"negatives={verdict['negative_indices'][:8]}, "
        f"runtime={round(time.time() - started, 1)}s"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
