#!/usr/bin/env python3
"""Exact certificate that Omega_3(r) > 4/3 for every integer r >= 3."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
COMMON_PATH = HERE / "verify_bridge_t2_minimum.py"
DEFAULT_OUTPUT = HERE / "BRIDGE_T3_BARRIER.json"
SCHEMA = "shift-wall-omega-t3-four-thirds-v1"
ALGORITHM = "exact-prefix-plus-corrected-a1-envelope-v1"
T = 3
M0 = 3_000
PREFIX_R_MAX = M0 + 1
BARRIER = F(4, 3)

EXPECTED_CSTAR_SHA256 = (
    "2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691"
)
EXPECTED_PREFIX_SHA256 = (
    "3d38ad91d746b02ff761d823ed72a1456bcb86843ea9d2010bb774aabf51f5af"
)
EXPECTED_PREFIX_MINIMUM_SHA256 = (
    "fece2e66ac7cdf95c7f3e19471cc6f09524b2479a6f736e6b5a45d95e4d3def8"
)
EXPECTED_TAIL_SHA256 = (
    "5832976deb2acc02cc1d2b43aabb09a5d97e9d672b1a35d59c4d3535fd934bd9"
)
EXPECTED_Q_LOWER_SHA256 = (
    "a1fa8edef32e9dc33b1bcdfd06b17d435b72c9d7bcb8fab35d54f1b4f66a0eb7"
)


def indexed_polynomial_fingerprint(common, values: list[F]) -> str:
    digest = hashlib.sha256()
    for index, value in enumerate(values):
        digest.update(
            f"{index}:{common.fraction_literal(value)}\n".encode("ascii")
        )
    return digest.hexdigest()


def load_common():
    specification = importlib.util.spec_from_file_location("bridge_t2_common", COMMON_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load the t=2 verifier module")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def capture_envelopes(a1, *, dscale: int = 1):
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
        result = a1.section3_tail(
            M0, t=T, verbose=False, record=False, dscale=dscale
        )
    finally:
        sys.settrace(prior_trace)
    required = {"blo", "bhi", "clo", "chi", "tau0"}
    if not required <= captured.keys():
        raise RuntimeError("A1 envelope locals were not captured")
    return result, captured


def envelope_polynomials(common, a1, captured: dict[str, object], target: F):
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
    return {
        "blo": blo,
        "bhi": bhi,
        "clo": clo,
        "chi": chi,
        "x_lower": x_lower,
        "y_upper": y_upper,
        "z_upper": z_upper,
        "p_lower": p_lower,
        "q_upper": q_upper,
        "q_lower": q_lower,
        "q_width": q_width,
        "h_poly": h_poly,
    }


def definition_checks(common) -> tuple[F, dict]:
    wanted = set(range(T, 61))
    direct: dict[int, F] = {}
    for r, numerator, denominator in common.reverse_bridge_stream(60, t=T):
        if r >= T:
            common.require(denominator > 0, f"reverse denominator failed at r={r}")
            direct[r] = F(numerator, denominator)
    common.require(direct.keys() == wanted, "reverse stream missed a t=3 cell")

    rows = common.forward_wall_rows(2 * 60 - 1)
    for r in range(T, 61):
        common.require(
            direct[r] == common.forward_bridge_value(rows, r, t=T),
            f"forward/reverse mismatch at r={r}",
        )
    ladder = common.ladder_values(60, wanted, t=T)
    for r in wanted:
        common.require(direct[r] == ladder[r], f"ladder mismatch at r={r}")

    c_star = None
    for r, numerator, denominator in common.reverse_bridge_stream(
        common.R_STAR, t=2
    ):
        if r == common.R_STAR:
            c_star = F(numerator, denominator)
    common.require(c_star is not None, "could not rebuild c_star")
    c_star_fp = common.fraction_fingerprint(c_star)
    common.require(c_star_fp["sha256"] == EXPECTED_CSTAR_SHA256,
                   "c_star fingerprint changed")
    common.require(3 * c_star.numerator < 4 * c_star.denominator,
                   "c_star is not below 4/3")

    # A top window one coefficient too short silently substitutes B_4=0 in P_3.
    # It must disagree with the full definition.
    short_value = None
    k_max = T
    previous_two = [1] + [0] * k_max
    previous = [1, 1] + [0] * (k_max - 1)
    for n in range(2, 2 * 10):
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
        if n == 19:
            p_value = common.turan(common.shifted_reverse(previous), T)
            q_value = common.turan(current, T)
            short_value = F(20 * 20 * p_value, 2 * q_value)
        previous_two, previous = previous, current
    common.require(short_value is not None and short_value != direct[10],
                   "undersized top-window mutation was not rejected")

    return c_star, {
        "forward_reverse_cells": len(wanted),
        "ladder_reverse_cells": len(wanted),
        "c_star_rebuilt": c_star_fp,
        "c_star_below_four_thirds_exact": True,
        "undersized_top_window_rejected": short_value != direct[10],
    }


def prefix_certificate(common) -> dict:
    digest = hashlib.sha256()
    cell_count = 0
    denominator_positive_count = 0
    minimum: tuple[int, int, int] | None = None
    failures: list[int] = []
    for r, numerator, denominator in common.reverse_bridge_stream(
        PREFIX_R_MAX, t=T
    ):
        if r < T:
            continue
        common.require(numerator > 0, f"prefix numerator failed at r={r}")
        common.require(denominator > 0, f"prefix denominator failed at r={r}")
        cell_count += 1
        denominator_positive_count += 1
        difference = 3 * numerator - 4 * denominator
        comparison = ">" if difference > 0 else "=" if difference == 0 else "<"
        if difference <= 0:
            failures.append(r)
        digest.update(f"{r}:{comparison}\n".encode("ascii"))
        if (
            minimum is None
            or numerator * minimum[2] < minimum[1] * denominator
        ):
            minimum = (r, numerator, denominator)

    common.require(cell_count == 2_999, "prefix cell count changed")
    common.require(not failures, f"t=3 prefix failed at {failures[:8]}")
    common.require(minimum is not None and minimum[0] == PREFIX_R_MAX,
                   "finite t=3 minimum row changed")
    comparison_sha = digest.hexdigest()
    common.require(comparison_sha == EXPECTED_PREFIX_SHA256,
                   f"t=3 prefix digest changed: {comparison_sha}")
    minimum_value = F(minimum[1], minimum[2])
    minimum_fp = common.fraction_fingerprint(minimum_value)
    common.require(minimum_fp["sha256"] == EXPECTED_PREFIX_MINIMUM_SHA256,
                   "finite t=3 minimum fingerprint changed")
    return {
        "range": [T, PREFIX_R_MAX],
        "cell_count": cell_count,
        "denominator_positive_count": denominator_positive_count,
        "non_strict_rows": failures,
        "comparison_stream_encoding": "ASCII r:<|=|> followed by newline",
        "comparison_sha256": comparison_sha,
        "finite_minimum_r": minimum[0],
        "finite_minimum": {
            **minimum_fp,
            "decimal_truncated": common.decimal_truncation(minimum_value, 24),
        },
        "minimum_margin_over_four_thirds_decimal_truncated":
            common.decimal_truncation(minimum_value - BARRIER, 24),
    }


def tail_certificate(common, a1) -> tuple[dict, dict]:
    source_result, captured = capture_envelopes(a1)
    common.require(source_result[0] is True, "corrected A1 t=3 source gate failed")
    data = envelope_polynomials(common, a1, captured, BARRIER)

    for lower_name, upper_name in (("blo", "bhi"), ("clo", "chi")):
        lower_rows, upper_rows = data[lower_name], data[upper_name]
        common.require(
            all(value >= 0 for row in lower_rows for value in row),
            f"{lower_name} has a negative coefficient",
        )
        common.require(
            all(value >= 0 for row in upper_rows for value in row),
            f"{upper_name} has a negative coefficient",
        )
        common.require(
            all(
                value >= 0
                for lower, upper in zip(lower_rows, upper_rows, strict=True)
                for value in a1.padd(upper, lower, -1)
            ),
            f"{lower_name}/{upper_name} order is inverted",
        )

    for name in ("x_lower", "y_upper", "z_upper", "p_lower", "q_upper"):
        common.require(all(value >= 0 for value in data[name]),
                       f"{name} has a negative coefficient")
    common.require(all(value >= 0 for value in data["q_lower"]),
                   "q_lower has a negative coefficient")
    common.require(data["q_lower"][0] > 0,
                   "q_lower does not have a positive constant")
    common.require(all(value >= 0 for value in data["q_width"]),
                   "q_lower/q_upper order is inverted")
    common.require(all(value > 0 for value in data["h_poly"]),
                   "four-thirds tail polynomial is not strictly positive")

    h_sha = indexed_polynomial_fingerprint(common, data["h_poly"])
    q_lower_sha = indexed_polynomial_fingerprint(common, data["q_lower"])
    common.require(h_sha == EXPECTED_TAIL_SHA256,
                   f"t=3 tail digest changed: {h_sha}")
    common.require(q_lower_sha == EXPECTED_Q_LOWER_SHA256,
                   f"t=3 q_lower digest changed: {q_lower_sha}")

    capacities = [
        data["p_lower"][index] / (2 * coefficient)
        for index, coefficient in enumerate(data["q_upper"])
        if coefficient > 0
    ]
    capacity = min(capacities)
    limiting_indices = [
        index
        for index, coefficient in enumerate(data["q_upper"])
        if coefficient > 0
        and data["p_lower"][index] / (2 * coefficient) == capacity
    ]
    common.require(capacity > BARRIER, "tail coefficient capacity <=4/3")

    target_mutation = F(3, 2)
    mutated_h = a1.padd(
        data["p_lower"],
        [2 * target_mutation * value for value in data["q_upper"]],
        -1,
    )
    mutation_negative_indices = [
        index for index, value in enumerate(mutated_h) if value < 0
    ]
    common.require(mutation_negative_indices == [2, 3, 4, 5],
                   "3/2 mutation witness changed")

    bad_result, bad_captured = capture_envelopes(a1, dscale=2_000)
    bad_data = envelope_polynomials(common, a1, bad_captured, BARRIER)
    bad_lower_indices = [
        (family, level, index)
        for family in ("blo", "clo")
        for level, row in enumerate(bad_data[family])
        for index, value in enumerate(row)
        if value < 0
    ]
    bad_h_indices = [
        index for index, value in enumerate(bad_data["h_poly"]) if value < 0
    ]
    common.require(bad_result[0] is False, "INJ x2000 source mutation passed")
    common.require(bad_lower_indices and bad_h_indices,
                   "INJ x2000 mutation was not rejected")

    return {
        "anchor_m0": M0,
        "claim_start_r": M0 + 1,
        "variable": "x=tau_m-tau_M0_lower >= 0",
        "p_lower_degree": len(data["p_lower"]) - 1,
        "p_lower_nonnegative_coefficients": sum(
            value >= 0 for value in data["p_lower"]
        ),
        "q_upper_degree": len(data["q_upper"]) - 1,
        "q_upper_nonnegative_coefficients": sum(
            value >= 0 for value in data["q_upper"]
        ),
        "q_lower_degree": len(data["q_lower"]) - 1,
        "q_lower_nonnegative_coefficients": sum(
            value >= 0 for value in data["q_lower"]
        ),
        "q_lower_positive_constant": data["q_lower"][0] > 0,
        "q_lower_sha256": q_lower_sha,
        "coefficient_fingerprint_encoding":
            "ASCII index:numerator/denominator followed by newline",
        "retargeted_degree": len(data["h_poly"]) - 1,
        "retargeted_positive_coefficients": sum(
            value > 0 for value in data["h_poly"]
        ),
        "retargeted_sha256": h_sha,
        "coefficient_capacity": common.fraction_fingerprint(capacity),
        "coefficient_capacity_decimal_truncated":
            common.decimal_truncation(capacity, 24),
        "capacity_limiting_indices": limiting_indices,
        "three_halves_mutation_negative_indices": mutation_negative_indices,
        "inj_x2000_source_rejected": bad_result[0] is False,
        "inj_x2000_negative_lower_count": len(bad_lower_indices),
        "inj_x2000_negative_h_indices": bad_h_indices,
        "wallis_factor_gate": {
            "statement": "k_m=8(m+1)^2 v_m^2/pi^2 > 1/2 for m>=1",
            "algebraic_difference":
                "(m+1)(4m+1)-(2m+1)^2=m>0",
            "passed": M0 >= 1,
        },
    }, {
        "source_certificate_passed": source_result[0] is True,
        "strict_h": all(value > 0 for value in data["h_poly"]),
        "positive_denominator_lower": data["q_lower"][0] > 0
        and all(value >= 0 for value in data["q_lower"]),
        "mutations_rejected": mutation_negative_indices == [2, 3, 4, 5]
        and bad_result[0] is False
        and bool(bad_lower_indices)
        and bool(bad_h_indices),
    }


def main() -> int:
    sys.set_int_max_str_digits(1_000_000)
    common = load_common()
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-repo", type=Path, default=common.DEFAULT_REPO)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    repo = args.source_repo.resolve()
    common.require(repo.is_relative_to(common.RELEASE_ROOT),
                   "source repository must be inside the release tree")

    sources = common.source_hashes(repo)
    c_star, definition = definition_checks(common)
    prefix = prefix_certificate(common)
    a1 = common.load_a1(repo)
    pi_repair = common.repair_a1_pi_interval(a1)
    tail, tail_gates = tail_certificate(common, a1)

    gates = {
        "source_hashes_pinned": sources == common.EXPECTED_SOURCE_SHA256,
        "definition_checks_passed": definition["c_star_below_four_thirds_exact"]
        and definition["undersized_top_window_rejected"],
        "prefix_strictly_above_four_thirds": not prefix["non_strict_rows"],
        "corrected_pi_interval":
            pi_repair["corrected_interval_contains_pi_bracket"]
            and pi_repair["defective_reference_tripwire_rejected"]
            and pi_repair["corrected_reference_tripwire_passed"],
        "tail_source_certificate_passed":
            tail_gates["source_certificate_passed"],
        "tail_strict_polynomial": tail_gates["strict_h"],
        "tail_denominator_positive": tail_gates["positive_denominator_lower"],
        "wallis_factor_strict": tail["wallis_factor_gate"]["passed"],
        "negative_controls_passed": tail_gates["mutations_rejected"],
    }
    common.require(all(gates.values()), "one or more final t=3 gates failed")

    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "status": "PASS",
        "passed": True,
        "claim": "Omega_3(r)>4/3>Omega_2(1350) for every integer r>=3",
        "claim_domain": "integers r>=3 at reverse index t=3",
        "normalization":
            "Omega_t(r)=(2r)^2 L((1+z)M_(2r-2))_(r-t)/(2 L(M_(2r-1))_(r-t))",
        "source_repo": repo.relative_to(common.RELEASE_ROOT).as_posix(),
        "source_hashes": sources,
        "common_verifier_sha256": common.sha256_file(COMMON_PATH),
        "producer_sha256": common.sha256_file(Path(__file__).resolve()),
        "execution": {
            "arithmetic": "exact integers and fractions",
            "parallelism": "one process, no pool",
            "prefix_resume": "not applicable; 2999 cells",
        },
        "gates": gates,
        "definition_checks": definition,
        "barrier": "4/3",
        "c_star": {
            **common.fraction_fingerprint(c_star),
            "decimal_truncated": common.decimal_truncation(c_star, 30),
            "exactly_below_barrier": True,
            "barrier_gap_decimal_truncated":
                common.decimal_truncation(BARRIER - c_star, 24),
        },
        "a1_pi_interval_repair": pi_repair,
        "prefix": prefix,
        "tail": tail,
        "proof_join": {
            "prefix": f"strict on {T}<=r<={PREFIX_R_MAX}",
            "tail": f"strict on r>={M0 + 1}",
            "overlap_r": M0 + 1,
            "conclusion": "Omega_3(r)>4/3 for every integer r>=3",
        },
        "claim_boundary": (
            "This certificate closes the t=3 bridge stripe above 4/3. "
            "It does not find the sharp t=3 minimum, compare t>=4, "
            "or assert pointwise monotonicity in t."
        ),
    }
    common.write_json(args.output.resolve(), payload)
    print(
        f"[PASS] Omega_3(r)>4/3 for every r>=3; "
        f"prefix_sha256={prefix['comparison_sha256']}"
    )
    print(
        f"[PASS] tail r>={M0 + 1}; degree={tail['retargeted_degree']}; "
        f"positive_coefficients={tail['retargeted_positive_coefficients']}; "
        f"sha256={tail['retargeted_sha256']}"
    )
    print(
        f"[PASS] exact 4/3-c_star="
        f"{common.decimal_truncation(BARRIER-c_star, 24)}...; "
        f"wrote {args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
