#!/usr/bin/env python3
"""Exact certificate for the normalized shift-wall t=2 minimum.

The normalized bridge ratio is

    Omega_t(r) = (2r)^2 L((1+z)M_{2r-2})_{r-t}
                 / (2 L(M_{2r-1})_{r-t}).

This verifier proves that Omega_2 has its unique global minimum at r=1350.
It combines a definition-level exact prefix with the two-sided A1 continuum
envelopes for the half-line.  All verdicts use integers or Fractions.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parents[1]
DEFAULT_REPO = RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay"
DEFAULT_OUTPUT = HERE / "BRIDGE_T2_MINIMUM.json"

SCHEMA = "shift-wall-omega-t2-minimum-v1"
ALGORITHM = "exact-prefix-plus-a1-envelope-v1"
T = 2
R_STAR = 1_350
M0 = 3_000
PREFIX_R_MAX = M0 + 1

EXPECTED_SOURCE_SHA256 = {
    "OPEN.md": "e4537f1fa12de71931e66929ad0eef8e138e42941c83ccf6e0b4496e9030b6da",
    "paper/PANEL3_odd_cycles_sech_square_tail.tex":
        "e60836828ba688bc8601837d829c1e7becd0454209f2200524d12622f25c49bc",
    "verification/A1_verify_tail.py":
        "7cc90454d04a71b2eecc6123eee98c82896f2809ff1514d0262125732bb50958",
    "verification/A36_bw_bridge.py":
        "2bf61ca5f20f47b91d9a0ab94804ac64cc81ea158447a08faade940ae7de2581",
    "verification/A38_targeted_core.py":
        "4655e738438b9c9671eeefb8f0d4275e628fb3d86de615e405ec4d8b1acd6cef",
    "verification/engine.py":
        "7f90e47d835d0737819a8631781a502dc8e70d9af95e12aab5774cef45783825",
}

EXPECTED_TARGET_SHA256 = (
    "2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691"
)
EXPECTED_RAW_CONSTANT_SHA256 = (
    "e5e91c6fedb65c563e9e9680bba5b5a54752d5f5a46be155712b5fb6d1a4d8c3"
)
EXPECTED_PREFIX_COMPARISON_SHA256 = (
    "943bebfb449b9ec178562598433b1f76dd9f3cb0cab0e8203bfce784882245ae"
)
EXPECTED_CLOSEST_GAP_SHA256 = (
    "6fe47fdd07e7bff4f0e8e0aec3fb9867a720a77e942f156911c37184de621c8d"
)
EXPECTED_TAIL_POLYNOMIAL_SHA256 = (
    "d8ba995d1f52b231ff5e198bd881ec9d1bd0def96c8b9c1c6efbf1f95130482a"
)
PI_100_DECIMAL_DIGITS = (
    "314159265358979323846264338327950288419716939937510"
    "58209749445923078164062862089986280348253421170679"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fraction_literal(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def fraction_fingerprint(value: F) -> dict[str, str | int]:
    literal = fraction_literal(value)
    return {
        "sha256": hashlib.sha256(literal.encode("ascii")).hexdigest(),
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
    }


def polynomial_fingerprint(values: list[F]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(fraction_literal(value).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def decimal_truncation(value: F, places: int = 18) -> str:
    sign = "-" if value < 0 else ""
    whole, remainder = divmod(abs(value.numerator), value.denominator)
    digits: list[str] = []
    for _ in range(places):
        remainder *= 10
        digit, remainder = divmod(remainder, value.denominator)
        digits.append(str(digit))
    return f"{sign}{whole}." + "".join(digits)


def turan(row: list[int] | list[F], index: int):
    zero = 0 * row[0]

    def coefficient(j: int):
        return row[j] if 0 <= j < len(row) else zero

    return coefficient(index) ** 2 - coefficient(index - 1) * coefficient(index + 1)


def shifted_reverse(row: list[int]) -> list[int]:
    return [row[k] + (row[k - 1] if k else 0) for k in range(len(row))]


def reverse_bridge_stream(r_max: int, t: int = T, mutation=None):
    """Yield the raw numerator and denominator of Omega_t(r).

    Rows are the top coefficients of M_n in reverse order.  The recurrence is
    the Wall recurrence itself, with no normalization or floating point.
    """

    k_max = t + 1
    previous_two = [1] + [0] * k_max
    previous = [1, 1] + [0] * (k_max - 1)
    for n in range(2, 2 * r_max):
        if n % 2 == 0:
            current = [
                previous[k] + n * n * previous_two[k]
                for k in range(k_max + 1)
            ]
        else:
            current = [
                (previous[k - 1] if k else 0) + n * n * previous_two[k]
                for k in range(k_max + 1)
            ]
        if mutation is not None and n == mutation[0]:
            _, index, delta = mutation
            current[index] += delta
        if n % 2:
            r = (n + 1) // 2
            p_value = turan(shifted_reverse(previous), t)
            q_value = turan(current, t)
            yield r, (2 * r) ** 2 * p_value, 2 * q_value
        previous_two, previous = previous, current


def one_plus_z(row: list[int]) -> list[int]:
    result = [0] * (len(row) + 1)
    for index, value in enumerate(row):
        result[index] += value
        result[index + 1] += value
    return result


def forward_wall_rows(n_max: int) -> list[list[int]]:
    rows = [[1], [1, 1]]
    for n in range(2, n_max + 1):
        result = rows[-1] + [0] * max(0, len(rows[-2]) + 1 - len(rows[-1]))
        for index, value in enumerate(rows[-2]):
            result[index + 1] += n * n * value
        rows.append(result)
    return rows


def forward_bridge_value(rows: list[list[int]], r: int, t: int = T) -> F:
    q = r - t
    numerator = (2 * r) ** 2 * turan(one_plus_z(rows[2 * r - 2]), q)
    denominator = 2 * turan(rows[2 * r - 1], q)
    require(denominator > 0, f"forward denominator failed at r={r}")
    return F(numerator, denominator)


def ladder_values(r_max: int, wanted: set[int], t: int = T) -> dict[int, F]:
    """Independent exact Wallis-ladder reconstruction of selected cells."""

    h = F(1)
    total = F(1)
    b_row = [F(1)] + [F(0)] * (t + 1)
    c_row = [F(1), F(1)] + [F(0)] * t
    values: dict[int, F] = {}
    for m in range(1, r_max):
        h *= F((2 * m) ** 2, (2 * m - 1) ** 2)
        q_step = 1 / h
        v_step = h / F((2 * m + 1) ** 2)
        old_b, old_c = b_row, c_row
        total += q_step
        b_row = [total] + [
            old_b[k] + q_step * old_c[k] for k in range(1, t + 2)
        ]
        c_row = [F(1)] + [
            old_c[k] + v_step * b_row[k - 1] for k in range(1, t + 2)
        ]
        r = m + 1
        if r in wanted:
            p_value = turan(
                [b_row[k] + (b_row[k - 1] if k else F(0))
                 for k in range(len(b_row))],
                t,
            )
            q_value = turan(c_row, t)
            require(q_value > 0, f"ladder denominator failed at r={r}")
            values[r] = F(2 * r * r) * v_step * v_step * p_value / q_value
    require(values.keys() == wanted, "ladder did not produce every requested row")
    return values


def definition_checks() -> tuple[F, dict[str, int | bool]]:
    wanted = set(range(2, 61)) | {R_STAR}
    direct: dict[int, F] = {}
    for r, numerator, denominator in reverse_bridge_stream(R_STAR):
        if r in wanted:
            require(denominator > 0, f"reverse denominator failed at r={r}")
            direct[r] = F(numerator, denominator)
    require(direct.keys() == wanted, "reverse stream missed requested rows")

    forward_rows = forward_wall_rows(2 * 30 - 1)
    for r in range(2, 31):
        require(
            direct[r] == forward_bridge_value(forward_rows, r),
            f"forward/reverse definition mismatch at r={r}",
        )

    ladder = ladder_values(R_STAR, wanted)
    for r in wanted:
        require(direct[r] == ladder[r], f"definition/ladder mismatch at r={r}")

    unmutated_r8 = direct[8]
    mutated_r8 = None
    for r, numerator, denominator in reverse_bridge_stream(
        8, mutation=(14, 2, 1)
    ):
        if r == 8:
            require(denominator > 0, "mutated denominator unexpectedly failed")
            mutated_r8 = F(numerator, denominator)
    require(mutated_r8 is not None and mutated_r8 != unmutated_r8,
            "reverse-row mutation was not detected")

    target = direct[R_STAR]
    return target, {
        "forward_reverse_cells": 29,
        "ladder_reverse_cells": len(wanted),
        "target_definition_ladder_exact": target == ladder[R_STAR],
        "reverse_row_mutation_rejected": mutated_r8 != unmutated_r8,
    }


def prefix_certificate(target: F) -> dict:
    comparison_digest = hashlib.sha256()
    equal_rows: list[int] = []
    less_rows: list[int] = []
    greater_count = 0
    denominator_positive_count = 0
    numerator_positive_count = 0
    closest: tuple[int, int, int] | None = None
    selected: dict[int, F] = {}

    for r, numerator, denominator in reverse_bridge_stream(PREFIX_R_MAX):
        require(denominator > 0, f"prefix denominator failed at r={r}")
        require(numerator > 0, f"prefix numerator failed at r={r}")
        denominator_positive_count += 1
        numerator_positive_count += 1
        difference = numerator * target.denominator - target.numerator * denominator
        if difference < 0:
            comparison = "<"
            less_rows.append(r)
        elif difference > 0:
            comparison = ">"
            greater_count += 1
            gap_denominator = denominator * target.denominator
            if (
                closest is None
                or difference * closest[2] < closest[1] * gap_denominator
            ):
                closest = (r, difference, gap_denominator)
        else:
            comparison = "="
            equal_rows.append(r)
        comparison_digest.update(f"{r}:{comparison}\n".encode("ascii"))
        if r in {R_STAR - 1, R_STAR, R_STAR + 1}:
            selected[r] = F(numerator, denominator)

    require(equal_rows == [R_STAR], f"prefix equality set changed: {equal_rows}")
    require(not less_rows, f"prefix has values below target: {less_rows[:8]}")
    require(closest is not None and closest[0] == R_STAR - 1,
            "closest prefix competitor changed")
    digest = comparison_digest.hexdigest()
    require(digest == EXPECTED_PREFIX_COMPARISON_SHA256,
            f"prefix digest changed: {digest}")
    gap = F(closest[1], closest[2])
    gap_fp = fraction_fingerprint(gap)
    require(gap_fp["sha256"] == EXPECTED_CLOSEST_GAP_SHA256,
            "closest-gap fingerprint changed")

    return {
        "range": [2, PREFIX_R_MAX],
        "cell_count": PREFIX_R_MAX - 1,
        "denominator_positive_count": denominator_positive_count,
        "numerator_positive_count": numerator_positive_count,
        "less_rows": less_rows,
        "equal_rows": equal_rows,
        "greater_count": greater_count,
        "comparison_stream_encoding": "ASCII r:<|=|> followed by newline",
        "comparison_sha256": digest,
        "closest_competitor_r": closest[0],
        "closest_gap": gap_fp,
        "closest_gap_decimal_truncated": decimal_truncation(gap, 24),
        "neighbors": {
            str(r): {
                **fraction_fingerprint(value),
                "decimal_truncated": decimal_truncation(value, 24),
            }
            for r, value in sorted(selected.items())
        },
    }


def load_a1(repo: Path):
    verification = repo / "verification"
    sys.path.insert(0, str(verification))
    try:
        return importlib.import_module("A1_verify_tail")
    finally:
        sys.path.pop(0)


def repair_a1_pi_interval(a1) -> dict[str, int | bool]:
    """Replace A1's truncating PI upper endpoint by an outward ceiling.

    The pinned source uses ``-int(-PI_HI*D)``.  For a positive Fraction,
    ``int`` truncates toward zero, so that expression is another floor, not a
    ceiling.  The corrected one-ulp-wider interval is installed before any A1
    envelope is constructed.
    """

    original_lo = a1.PI.lo
    original_hi = a1.PI.hi
    require(len(PI_100_DECIMAL_DIGITS) == 101,
            "hard-coded pi reference must have 100 decimal places")
    reference_denominator = 10**100
    reference_lower = F(int(PI_100_DECIMAL_DIGITS), reference_denominator)
    reference_upper = F(
        int(PI_100_DECIMAL_DIGITS) + 1, reference_denominator
    )

    def reference_tripwire(interval) -> bool:
        return (
            interval.flo() < reference_lower
            and interval.fhi() > reference_upper
        )

    inherited_reference_tripwire = reference_tripwire(a1.PI)
    lower_numerator = a1.PI_LO.numerator * a1.D
    upper_numerator = a1.PI_HI.numerator * a1.D
    corrected_lo = lower_numerator // a1.PI_LO.denominator
    corrected_hi = -((-upper_numerator) // a1.PI_HI.denominator)
    inherited_upper_misses_lower_bound = a1.PI.fhi() < a1.PI_LO
    require(inherited_upper_misses_lower_bound,
            "the pinned A1 PI-rounding defect was not reproduced")
    require(F(corrected_lo, 1) <= F(lower_numerator, a1.PI_LO.denominator),
            "corrected PI lower endpoint is not outward")
    require(F(corrected_hi, a1.D) >= a1.PI_HI,
            "corrected PI upper endpoint is not outward")
    a1.PI = a1.Iv(corrected_lo, corrected_hi)
    require(a1.PI.flo() <= a1.PI_LO and a1.PI.fhi() >= a1.PI_HI,
            "corrected PI interval does not contain the rational bracket")
    corrected_reference_tripwire = reference_tripwire(a1.PI)
    require(not inherited_reference_tripwire,
            "defective PI singleton unexpectedly passed the reference tripwire")
    require(corrected_reference_tripwire,
            "corrected PI interval failed the reference tripwire")
    return {
        "inherited_upper_misses_certified_pi_lower_bound":
            inherited_upper_misses_lower_bound,
        "original_width_ulps": original_hi - original_lo,
        "corrected_width_ulps": corrected_hi - corrected_lo,
        "upper_endpoint_delta_ulps": corrected_hi - original_hi,
        "corrected_interval_contains_pi_bracket": True,
        "reference_decimal_places": 100,
        "defective_reference_tripwire_rejected":
            not inherited_reference_tripwire,
        "corrected_reference_tripwire_passed": corrected_reference_tripwire,
    }


def capture_a1_envelopes(a1):
    captured: dict[str, object] = {}
    target_code = a1.section3_tail.__code__

    def trace(frame, event, argument):
        if frame.f_code is target_code:
            if event == "return":
                captured.update(frame.f_locals)
            return trace
        return None

    prior_trace = sys.gettrace()
    try:
        sys.settrace(trace)
        result = a1.section3_tail(M0, t=T, verbose=False, record=False)
    finally:
        sys.settrace(prior_trace)
    require(result[0] is True, "A1 source envelope certificate failed")
    required = {"blo", "bhi", "clo", "chi", "tau0"}
    require(required <= captured.keys(), "A1 envelope locals were not captured")
    return result, captured


def tail_certificate(a1, captured: dict[str, object], target: F) -> dict:
    blo = captured["blo"]
    bhi = captured["bhi"]
    clo = captured["clo"]
    chi = captured["chi"]
    zero = [F(0)]

    def get(rows, index):
        return rows[index] if 0 <= index < len(rows) else zero

    x_lower = a1.padd(get(blo, T - 1), get(blo, T))
    y_upper = a1.padd(get(bhi, T - 2), get(bhi, T - 1))
    z_upper = a1.padd(get(bhi, T), get(bhi, T + 1))
    p_lower = a1.padd(
        a1.pmul(x_lower, x_lower),
        a1.pmul(y_upper, z_upper),
        -1,
    )
    q_upper = a1.padd(
        a1.pmul(get(chi, T), get(chi, T)),
        a1.pmul(get(clo, T - 1), get(clo, T + 1)),
        -1,
    )
    q_lower = a1.padd(
        a1.pmul(get(clo, T), get(clo, T)),
        a1.pmul(get(chi, T - 1), get(chi, T + 1)),
        -1,
    )
    q_width = a1.padd(q_upper, q_lower, -1)
    h_poly = a1.padd(
        p_lower,
        [2 * target * coefficient for coefficient in q_upper],
        -1,
    )

    require(all(coefficient >= 0 for row in blo for coefficient in row),
            "a lower b envelope has a negative coefficient")
    require(all(coefficient >= 0 for row in bhi for coefficient in row),
            "an upper b envelope has a negative coefficient")
    require(all(coefficient >= 0 for row in clo for coefficient in row),
            "a lower C envelope has a negative coefficient")
    require(all(coefficient >= 0 for row in chi for coefficient in row),
            "an upper C envelope has a negative coefficient")
    require(
        all(
            coefficient >= 0
            for lower, upper in zip(blo, bhi, strict=True)
            for coefficient in a1.padd(upper, lower, -1)
        ),
        "a b envelope is inverted",
    )
    require(
        all(
            coefficient >= 0
            for lower, upper in zip(clo, chi, strict=True)
            for coefficient in a1.padd(upper, lower, -1)
        ),
        "a C envelope is inverted",
    )
    require(all(coefficient >= 0 for coefficient in x_lower),
            "X lower polynomial is not nonnegative")
    require(all(coefficient >= 0 for coefficient in y_upper + z_upper),
            "Y/Z upper polynomial is not nonnegative")
    require(all(coefficient > 0 for coefficient in p_lower),
            "P lower polynomial is not strictly coefficient-positive")
    require(all(coefficient > 0 for coefficient in q_upper),
            "Q upper polynomial is not strictly coefficient-positive")
    require(all(coefficient > 0 for coefficient in q_lower),
            "Q lower polynomial is not strictly coefficient-positive")
    require(all(coefficient >= 0 for coefficient in q_width),
            "Q lower/upper polynomial order is inverted")
    require(all(coefficient > 0 for coefficient in h_poly),
            "retargeted bridge polynomial is not strictly coefficient-positive")

    h_digest = polynomial_fingerprint(h_poly)
    require(h_digest == EXPECTED_TAIL_POLYNOMIAL_SHA256,
            f"tail polynomial digest changed: {h_digest}")
    p_anchor = a1.pev(p_lower, F(0))
    q_anchor = a1.pev(q_upper, F(0))
    q_lower_anchor = a1.pev(q_lower, F(0))
    require(q_anchor > 0 and q_lower_anchor > 0,
            "tail denominator anchor is not positive")
    anchor_lower_bound = p_anchor / (2 * q_anchor)
    require(anchor_lower_bound > target, "tail anchor does not beat target")

    capacities = [
        p_lower[index] / (2 * coefficient)
        for index, coefficient in enumerate(q_upper)
        if coefficient > 0
    ]
    capacity = min(capacities)
    require(capacity > target, "coefficient capacity does not beat target")

    mutated_target = target * F(101, 100)
    mutated_poly = a1.padd(
        p_lower,
        [2 * mutated_target * coefficient for coefficient in q_upper],
        -1,
    )
    negative_mutation_indices = [
        index for index, coefficient in enumerate(mutated_poly) if coefficient < 0
    ]
    require(negative_mutation_indices,
            "one-percent target mutation was not rejected")

    tau0 = captured["tau0"]
    return {
        "anchor_m0": M0,
        "claim_start_r": M0 + 1,
        "variable": "x=tau_m-tau_M0_lower >= 0",
        "tau0_lower": fraction_fingerprint(tau0.flo()),
        "tau0_width": fraction_fingerprint(tau0.wid()),
        "p_lower_degree": len(p_lower) - 1,
        "p_lower_positive_coefficients": sum(value > 0 for value in p_lower),
        "p_lower_sha256": polynomial_fingerprint(p_lower),
        "q_upper_degree": len(q_upper) - 1,
        "q_upper_positive_coefficients": sum(value > 0 for value in q_upper),
        "q_upper_sha256": polynomial_fingerprint(q_upper),
        "q_lower_degree": len(q_lower) - 1,
        "q_lower_positive_coefficients": sum(value > 0 for value in q_lower),
        "q_lower_sha256": polynomial_fingerprint(q_lower),
        "q_lower_anchor_positive": q_lower_anchor > 0,
        "q_upper_minus_lower_nonnegative_coefficients": sum(
            value >= 0 for value in q_width
        ),
        "q_upper_minus_lower_strict_coefficients": sum(
            value > 0 for value in q_width
        ),
        "retargeted_degree": len(h_poly) - 1,
        "retargeted_positive_coefficients": sum(value > 0 for value in h_poly),
        "retargeted_sha256": h_digest,
        "anchor_lower_bound": fraction_fingerprint(anchor_lower_bound),
        "anchor_lower_bound_decimal_truncated": decimal_truncation(
            anchor_lower_bound, 24
        ),
        "coefficient_capacity": fraction_fingerprint(capacity),
        "coefficient_capacity_decimal_truncated": decimal_truncation(capacity, 24),
        "one_percent_target_mutation_negative_indices": negative_mutation_indices,
        "wallis_factor_gate": {
            "statement": "k_m=8(m+1)^2 v_m^2/pi^2 > 1/2 for m>=1",
            "algebraic_difference":
                "(m+1)(4m+1)-(2m+1)^2=m>0",
            "passed": M0 >= 1,
        },
    }


def source_hashes(repo: Path) -> dict[str, str]:
    actual: dict[str, str] = {}
    for relative, expected in EXPECTED_SOURCE_SHA256.items():
        path = repo / relative
        require(path.is_file(), f"missing source: {path}")
        digest = sha256_file(path)
        require(digest == expected, f"source hash mismatch: {relative}")
        actual[relative] = digest
    return actual


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
    parser.add_argument("--source-repo", type=Path, default=DEFAULT_REPO)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    repo = args.source_repo.resolve()
    require(repo.is_relative_to(RELEASE_ROOT),
            "source repository must be inside the release tree")

    pinned_sources = source_hashes(repo)
    target, definition = definition_checks()
    target_fp = fraction_fingerprint(target)
    require(target_fp["sha256"] == EXPECTED_TARGET_SHA256,
            "target fingerprint changed")
    require(target_fp["numerator_digits"] == 6_733, "target numerator size changed")
    require(target_fp["denominator_digits"] == 6_732,
            "target denominator size changed")
    raw_constant_fp = fraction_fingerprint(2 * target)
    require(raw_constant_fp["sha256"] == EXPECTED_RAW_CONSTANT_SHA256,
            "raw bridge constant fingerprint changed")

    prefix = prefix_certificate(target)
    a1 = load_a1(repo)
    pi_interval_repair = repair_a1_pi_interval(a1)
    source_result, captured = capture_a1_envelopes(a1)
    tail = tail_certificate(a1, captured, target)

    gates = {
        "source_hashes_pinned": pinned_sources == EXPECTED_SOURCE_SHA256,
        "definition_crosschecks_passed": all(
            bool(value) for key, value in definition.items() if key.endswith("exact")
            or key.endswith("rejected")
        ),
        "target_fingerprint_matches": target_fp["sha256"]
        == EXPECTED_TARGET_SHA256,
        "prefix_unique_minimum": prefix["equal_rows"] == [R_STAR]
        and not prefix["less_rows"],
        "tail_source_certificate_passed": source_result[0] is True,
        "tail_strict_coefficient_positivity":
            tail["retargeted_positive_coefficients"]
            == tail["retargeted_degree"] + 1,
        "tail_denominator_positive": tail["q_lower_anchor_positive"] is True,
        "wallis_factor_strict": tail["wallis_factor_gate"]["passed"] is True,
        "negative_controls_passed": definition["reverse_row_mutation_rejected"]
        and bool(tail["one_percent_target_mutation_negative_indices"])
        and pi_interval_repair[
            "inherited_upper_misses_certified_pi_lower_bound"
        ]
        and pi_interval_repair["corrected_interval_contains_pi_bracket"]
        and pi_interval_repair["defective_reference_tripwire_rejected"]
        and pi_interval_repair["corrected_reference_tripwire_passed"],
    }
    require(all(gates.values()), "one or more final gates failed")

    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "status": "PASS",
        "passed": True,
        "claim": "Omega_2(r) has its unique global minimum at r=1350",
        "claim_domain": "integers r>=2 at reverse index t=2",
        "normalization":
            "Omega_t(r)=(2r)^2 L((1+z)M_(2r-2))_(r-t)/(2 L(M_(2r-1))_(r-t))",
        "source_repo": repo.relative_to(RELEASE_ROOT).as_posix(),
        "source_hashes": pinned_sources,
        "producer_sha256": sha256_file(Path(__file__).resolve()),
        "execution": {
            "arithmetic": "exact integers and fractions",
            "parallelism": "one process, no pool",
            "prefix_resume": "not applicable; 3000 cells",
        },
        "gates": gates,
        "definition_checks": definition,
        "a1_pi_interval_repair": pi_interval_repair,
        "minimizer_r": R_STAR,
        "target": {
            **target_fp,
            "decimal_truncated": decimal_truncation(target, 30),
        },
        "raw_optimal_bridge_constant": {
            **raw_constant_fp,
            "relation": "2*Omega_2(1350)",
            "decimal_truncated": decimal_truncation(2 * target, 30),
        },
        "prefix": prefix,
        "tail": tail,
        "proof_join": {
            "prefix": f"exact on 2<=r<={PREFIX_R_MAX}",
            "tail": f"strict for r>={M0 + 1}",
            "overlap_r": M0 + 1,
            "conclusion": "unique global minimum on every integer r>=2",
        },
        "claim_boundary": (
            "This certificate closes only the t=2 bridge stripe. "
            "It does not compare t>=3 stripes or formalize the proof in Lean."
        ),
    }
    write_json(args.output.resolve(), payload)
    print(
        "[PASS] unique global t=2 minimum at r=1350; "
        f"target_sha256={target_fp['sha256']}"
    )
    print(
        f"[PASS] exact prefix r=2..{PREFIX_R_MAX}; "
        f"closest r={prefix['closest_competitor_r']}; "
        f"comparison_sha256={prefix['comparison_sha256']}"
    )
    print(
        f"[PASS] tail r>={M0 + 1}; degree={tail['retargeted_degree']}; "
        f"positive_coefficients={tail['retargeted_positive_coefficients']}; "
        f"sha256={tail['retargeted_sha256']}"
    )
    print(f"[PASS] wrote {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
