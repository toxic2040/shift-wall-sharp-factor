#!/usr/bin/env python3
"""Exact SQ t=4 barrier and its lossless transfer to the shift-wall bridge."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
COMMON_PATH = HERE / "verify_bridge_t2_minimum.py"
DEFAULT_OUTPUT = HERE / "SQ_T4_BRIDGE_TRANSFER.json"
SCHEMA = "sq-t4-four-thirds-bridge-transfer-v1"
ALGORITHM = "exact-rho-prefix-plus-corrected-independent-tail-v1"
T = 4
M0 = 1_000
PREFIX_R_MAX = M0
TARGET = F(4, 3)

EXPECTED_PREFIX_SHA256 = (
    "6b7fae02acecf8e08aa78d7aa333f902aaacd2def400fe1681efee88ad3f7532"
)
EXPECTED_PREFIX_MINIMUM_SHA256 = (
    "7e04dcafe265b3e6f748c056e4fd60b9a97156d6be50f74e95f2ec28625600a0"
)
EXPECTED_TAIL_SHA256 = (
    "a11810ccd429452dd13e90a969dc76fee7be5060c7e006a65fce525ecc2ed34a"
)
EXPECTED_INDEPENDENT_SHA256 = (
    "ad7afc6ac798c7c2db565b2358325fcf9a68de3926b1596d24a16f358a1a6ec6"
)


def load_module(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def load_common():
    return load_module(COMMON_PATH, "bridge_t4_common")


def indexed_polynomial_fingerprint(common, values: list[F]) -> str:
    digest = hashlib.sha256()
    for index, value in enumerate(values):
        digest.update(
            f"{index}:{common.fraction_literal(value)}\n".encode("ascii")
        )
    return digest.hexdigest()


def turan(row: list[int], index: int) -> int:
    def coefficient(j: int) -> int:
        return row[j] if 0 <= j < len(row) else 0

    return coefficient(index) ** 2 - coefficient(index - 1) * coefficient(index + 1)


def shifted(row: list[int]) -> list[int]:
    return [row[k] + (row[k - 1] if k else 0) for k in range(len(row))]


def rho_stream(r_max: int, t: int = T):
    """Exact factorial-cancelled top-window stream for rho(r,t)."""

    k_max = t + 1
    previous_two = [1] + [0] * k_max
    previous = [1, 1] + [0] * (k_max - 1)
    previous_even = [1] + [0] * k_max
    for n in range(2, 2 * r_max - 1):
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
        if n % 2 == 0:
            r = n // 2 + 1
            if r >= t:
                current_s = turan(shifted(current), t)
                predecessor_s = turan(shifted(previous_even), t)
                current_v = turan(current, t - 1)
                decrement = (
                    r * r * current_s
                    - (2 * r - 1) ** 4 * (r - 1) ** 2 * predecessor_s
                )
                yield (
                    r,
                    decrement * decrement,
                    2 * r * r * current_s * current_v,
                    decrement,
                )
            previous_even = current
        previous_two, previous = previous, current


def definition_and_prefix(common, repo: Path) -> dict:
    verification = repo / "verification"
    sys.path.insert(0, str(verification))
    try:
        import engine
    finally:
        sys.path.pop(0)
    wall_rows = engine.wall_M(2 * 60)

    digest = hashlib.sha256()
    failures: list[int] = []
    denominator_positive = 0
    decrement_positive = 0
    minimum: tuple[int, int, int] | None = None
    crosschecked = 0
    for r, numerator, denominator, decrement in rho_stream(PREFIX_R_MAX):
        common.require(numerator > 0, f"rho numerator failed at r={r}")
        common.require(denominator > 0, f"rho denominator failed at r={r}")
        common.require(decrement > 0, f"rho decrement failed at r={r}")
        denominator_positive += 1
        decrement_positive += 1
        if r <= 60:
            common.require(F(numerator, denominator) == engine.exact_rho(r, T, wall_rows),
                           f"cancelled/full rho mismatch at r={r}")
            crosschecked += 1
        difference = 3 * numerator - 4 * denominator
        comparison = ">" if difference > 0 else "=" if difference == 0 else "<"
        if difference <= 0:
            failures.append(r)
        digest.update(f"{r}:{comparison}\n".encode("ascii"))
        if minimum is None or numerator * minimum[2] < minimum[1] * denominator:
            minimum = (r, numerator, denominator)

    common.require(denominator_positive == 997, "t=4 prefix cell count changed")
    common.require(not failures, f"rho t=4 prefix failed: {failures[:8]}")
    common.require(minimum is not None and minimum[0] == PREFIX_R_MAX,
                   "rho t=4 finite minimum row changed")
    comparison_sha = digest.hexdigest()
    common.require(comparison_sha == EXPECTED_PREFIX_SHA256,
                   f"rho t=4 prefix digest changed: {comparison_sha}")
    minimum_value = F(minimum[1], minimum[2])
    minimum_fp = common.fraction_fingerprint(minimum_value)
    common.require(minimum_fp["sha256"] == EXPECTED_PREFIX_MINIMUM_SHA256,
                   "rho t=4 finite minimum fingerprint changed")
    return {
        "range": [T, PREFIX_R_MAX],
        "cell_count": denominator_positive,
        "definition_crosscheck_range": [T, 60],
        "definition_crosscheck_cells": crosschecked,
        "positive_denominator_count": denominator_positive,
        "positive_decrement_count": decrement_positive,
        "non_strict_rows": failures,
        "comparison_stream_encoding": "ASCII r:<|=|> followed by newline",
        "comparison_sha256": comparison_sha,
        "finite_minimum_r": minimum[0],
        "finite_minimum": {
            **minimum_fp,
            "decimal_truncated": common.decimal_truncation(minimum_value, 24),
        },
        "minimum_margin_over_four_thirds_decimal_truncated":
            common.decimal_truncation(minimum_value - TARGET, 24),
    }


def repair_independent_pi(common, source) -> dict:
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
        "downstream_actual_delta": "0/1",
    }


def build_with_components(source):
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
        result = source.build_certificate(m0=M0, t=T)
    finally:
        sys.settrace(prior_trace)
    required = {"u_lower", "p_upper", "q_upper", "v_index"}
    if not required <= captured.keys():
        raise RuntimeError("independent tail components were not captured")
    return result, captured


def tail_certificate(common, source) -> dict:
    result, captured = build_with_components(source)
    common.require(result["ok"] is True, "source t=4 tail certificate failed")
    common.require(source.check_interval_ladder(T),
                   "independent interval ladder crosscheck failed")
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
    common.require(all(value > 0 for value in h_poly),
                   "rho t=4 four-thirds tail polynomial is not strictly positive")
    h_sha = indexed_polynomial_fingerprint(common, h_poly)
    common.require(h_sha == EXPECTED_TAIL_SHA256,
                   f"rho t=4 tail digest changed: {h_sha}")

    capacities = [
        u_squared[index] / value
        for index, value in enumerate(denominator)
        if value > 0
    ]
    capacity = min(capacities)
    limiting_indices = [
        index
        for index, value in enumerate(denominator)
        if value > 0 and u_squared[index] / value == capacity
    ]
    common.require(capacity > TARGET, "rho t=4 capacity does not beat 4/3")
    mutated = source.padd(
        u_squared, [F(3, 2) * value for value in denominator], -1
    )
    negative_mutation_indices = [
        index for index, value in enumerate(mutated) if value < 0
    ]
    common.require(negative_mutation_indices == list(range(9, 22)),
                   "rho t=4 target mutation witness changed")

    return {
        "anchor_m0": M0,
        "claim_start_r": M0 + 1,
        "source_certificate_passed": True,
        "source_interval_ladder_crosscheck": True,
        "target_one_reconstruction_exact": True,
        "comparison_denominator_nonnegative": True,
        "retargeted_degree": len(h_poly) - 1,
        "retargeted_positive_coefficients": sum(value > 0 for value in h_poly),
        "retargeted_sha256": h_sha,
        "coefficient_fingerprint_encoding":
            "ASCII index:numerator/denominator followed by newline",
        "coefficient_capacity": common.fraction_fingerprint(capacity),
        "coefficient_capacity_decimal_truncated":
            common.decimal_truncation(capacity, 24),
        "capacity_limiting_indices": limiting_indices,
        "three_halves_mutation_negative_indices": negative_mutation_indices,
    }


def rebuild_c_star(common) -> F:
    target = None
    for r, numerator, denominator in common.reverse_bridge_stream(
        common.R_STAR, t=2
    ):
        if r == common.R_STAR:
            target = F(numerator, denominator)
    common.require(target is not None, "could not rebuild c_star")
    common.require(
        common.fraction_fingerprint(target)["sha256"]
        == common.EXPECTED_TARGET_SHA256,
        "c_star fingerprint changed",
    )
    common.require(target < TARGET, "c_star is not below 4/3")
    return target


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

    source_hashes = common.source_hashes(repo)
    independent_path = repo / "verification" / "A1_verify_tail_independent.py"
    common.require(common.sha256_file(independent_path) == EXPECTED_INDEPENDENT_SHA256,
                   "independent A1 source hash changed")
    prefix = definition_and_prefix(common, repo)
    source = load_module(independent_path, "a1_t4_independent")
    pi_repair = repair_independent_pi(common, source)
    tail = tail_certificate(common, source)
    c_star = rebuild_c_star(common)

    gates = {
        "source_hashes_pinned": source_hashes == common.EXPECTED_SOURCE_SHA256,
        "independent_source_pinned":
            common.sha256_file(independent_path) == EXPECTED_INDEPENDENT_SHA256,
        "definition_and_prefix_passed": not prefix["non_strict_rows"],
        "corrected_pi_tripwire_passed":
            pi_repair["defective_reference_tripwire_rejected"]
            and pi_repair["corrected_reference_tripwire_passed"],
        "tail_strictly_above_four_thirds":
            tail["retargeted_positive_coefficients"]
            == tail["retargeted_degree"] + 1,
        "negative_control_passed": bool(
            tail["three_halves_mutation_negative_indices"]
        ),
        "c_star_strictly_below_four_thirds": c_star < TARGET,
    }
    common.require(all(gates.values()), "one or more t=4 transfer gates failed")

    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "status": "PASS",
        "passed": True,
        "claim": (
            "rho(r,4)>4/3 for every integer r>=4; consequently "
            "Omega_4(r)>4/3>Omega_2(1350) for every integer r>=4"
        ),
        "source_repo": repo.relative_to(common.RELEASE_ROOT).as_posix(),
        "source_hashes": {
            **source_hashes,
            "verification/A1_verify_tail_independent.py":
                EXPECTED_INDEPENDENT_SHA256,
        },
        "common_verifier_sha256": common.sha256_file(COMMON_PATH),
        "producer_sha256": common.sha256_file(Path(__file__).resolve()),
        "execution": {
            "arithmetic": "exact integers, rational intervals, and fractions",
            "parallelism": "one process, no pool",
            "prefix_resume": "not applicable; 997 cells",
        },
        "gates": gates,
        "barrier": "4/3",
        "c_star": {
            **common.fraction_fingerprint(c_star),
            "decimal_truncated": common.decimal_truncation(c_star, 30),
        },
        "a1_pi_interval_repair": pi_repair,
        "sq_prefix": prefix,
        "sq_tail": tail,
        "transfer": {
            "theorem": (
                "Omega_t(r) >= (sum_{j=t}^r lambda_(j,t)/sqrt(rho(j,t)))^-2 "
                ">= min_{t<=j<=r} rho(j,t)"
            ),
            "weights":
                "lambda_(j,t)=(a_j-a_(j-1))/a_r >0 and sum lambda=1",
            "input": "rho(j,4)>4/3 for every integer j>=4",
            "conclusion": "Omega_4(r)>4/3 for every integer r>=4",
        },
        "claim_boundary": (
            "This certificate closes only the fixed t=4 SQ and bridge stripes. "
            "It does not prove a fixed-column profile for t>=5 or an all-t theorem."
        ),
    }
    common.write_json(args.output.resolve(), payload)
    print(
        f"[PASS] rho(r,4)>4/3 for every r>=4; "
        f"prefix_sha256={prefix['comparison_sha256']}"
    )
    print(
        f"[PASS] tail r>={M0 + 1}; degree={tail['retargeted_degree']}; "
        f"positive_coefficients={tail['retargeted_positive_coefficients']}; "
        f"sha256={tail['retargeted_sha256']}"
    )
    print(
        f"[PASS] lossless SQ-to-bridge transfer gives Omega_4(r)>4/3>c_star; "
        f"wrote {args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
