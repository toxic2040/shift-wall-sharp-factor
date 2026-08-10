#!/usr/bin/env python3
"""Exact reduction of the global-infimum problem to four tails and three windows.

This is a scratch theorem-boundary certificate.  It verifies the strengthened
row-free joint scalar, retains the finite-row correction factors, computes the
top-region floor, and locates the five exact level-parameter crossings that
replace the eventual t=6,7,8 tails by the existing top and joint theorems.

It does not certify the four fixed-column tails or the three bounded finite
windows.  Those are deliberately separate evidence streams.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
OUTPUT = RESULTS / "global_infimum_reduction.json"
CORE_RECORD = RESULTS / "finite_core.json"
CORE_RECORD_LABEL = "results/finite_core.json"
SCHEMA = "square-tail-global-infimum-reduction-v1"
TARGET = F(203, 200)
JOINT_CORE = F(513, 500)
THRESHOLD_ROWS = {25: 1_468, 30: 3_634, 35: 8_414, 36: 9_884, 40: 18_457}
FAILURES: list[str] = []
GATES: dict[str, bool] = {}


def gate(name: str, condition: bool, detail: str = "") -> bool:
    passed = bool(condition)
    GATES[name] = passed
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)
    return passed


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fraction_fingerprint(value: F) -> dict[str, str | int]:
    literal = f"{value.numerator}/{value.denominator}"
    return {
        "sha256": hashlib.sha256(literal.encode("ascii")).hexdigest(),
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
    }


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


def integer_encoding(value: int) -> bytes:
    magnitude = abs(value).to_bytes(max(1, (abs(value).bit_length() + 7) // 8), "big")
    return (
        (b"-" if value < 0 else b"+")
        + len(magnitude).to_bytes(8, "big")
        + magnitude
    )


def sign_variations(values: list[sp.Expr]) -> int:
    signs = [sp.sign(value) for value in values if value != 0]
    if not all(value in (-1, 1) for value in signs):
        raise RuntimeError("Sturm endpoint has an indeterminate sign")
    return sum(left != right for left, right in zip(signs, signs[1:]))


def primitive_sturm_digest(sequence: list[sp.Expr], variable: sp.Symbol) -> str:
    rows: list[str] = []
    for expression in sequence:
        polynomial = sp.Poly(expression, variable, domain=sp.QQ)
        _, primitive = polynomial.primitive()
        if primitive.LC() < 0:
            primitive = -primitive
        rows.append(",".join(str(value) for value in primitive.all_coeffs()))
    return hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()


def joint_core_certificate() -> dict:
    Q, T, S = sp.symbols("Q T S", nonnegative=True)
    tt = T + 6
    y = Q + 10
    a_upper = y + 5 * (tt - 2)
    alpha = (tt - 1) * (2 * tt - 1) / (tt * (2 * tt + 1))
    beta = tt * (2 * tt + 1) / ((tt + 1) * (2 * tt + 3))
    x0 = (y - 5) / (2 * tt * (2 * tt + 1))
    u0 = alpha * (1 - 5 / y)
    r0 = sp.factor(1 + x0 + (x0 * u0 + x0**2) * (1 - beta) / (1 - u0))
    weight = 2 * tt - 1
    expression = sp.together(
        (weight * r0 - 1) ** 2
        - 2 * sp.Rational(JOINT_CORE.numerator, JOINT_CORE.denominator) * a_upper * r0
    )
    numerator, denominator = expression.as_numer_denom()
    polynomial = sp.Poly(numerator, Q, T).primitive()[1]

    denominator_shifted = sp.Poly(sp.expand(denominator.subs(T, S)), Q, S)
    gate(
        "joint denominator has positive coefficients",
        all(value > 0 for value in denominator_shifted.coeffs()),
        f"terms={len(denominator_shifted.terms())}",
    )

    shifted = sp.Poly(sp.expand(polynomial.as_expr().subs(T, S + 2)), Q, S)
    payload = "\n".join(
        f"{monomial}:{coefficient}" for monomial, coefficient in shifted.terms()
    )
    shifted_digest = hashlib.sha256(payload.encode("ascii")).hexdigest()
    gate(
        "row-free 513/500 core for t>=8",
        len(shifted.terms()) == 74 and all(value > 0 for value in shifted.coeffs()),
        f"bidegree={shifted.degree_list()}, positive_terms={len(shifted.terms())}",
    )

    sturm_records: dict[str, dict] = {}
    for tv in (0, 1):
        univariate = sp.Poly(polynomial.as_expr().subs(T, tv), Q, domain=sp.QQ)
        sequence = sp.sturm(univariate.as_expr(), Q)
        at_zero = [sp.Poly(row, Q).eval(0) for row in sequence]
        at_infinity = [sp.Poly(row, Q).LC() for row in sequence]
        v_zero = sign_variations(at_zero)
        v_infinity = sign_variations(at_infinity)
        digest = primitive_sturm_digest(sequence, Q)
        passed = univariate.eval(0) > 0 and v_zero == v_infinity
        gate(
            f"row-free 513/500 core at t={tv + 6}",
            passed,
            f"Sturm variations={v_zero}->{v_infinity}, length={len(sequence)}",
        )
        sturm_records[str(tv + 6)] = {
            "value_at_zero_positive": bool(univariate.eval(0) > 0),
            "variations_zero": v_zero,
            "variations_infinity": v_infinity,
            "length": len(sequence),
            "primitive_sequence_sha256": digest,
        }

    mutation = sp.together(
        (weight * r0 - 1) ** 2 - 2 * sp.Rational(103, 100) * a_upper * r0
    )
    mutation_value = mutation.as_numer_denom()[0].subs({T: 0, Q: 95})
    gate("103/100 mutation rejected", mutation_value < 0, f"numerator={mutation_value}")
    gate("513/500 strengthens 101/100", JOINT_CORE > F(101, 100))
    return {
        "constant": f"{JOINT_CORE.numerator}/{JOINT_CORE.denominator}",
        "role": "row-free core only",
        "t_ge_8": {
            "bidegree": list(shifted.degree_list()),
            "positive_monomial_terms": len(shifted.terms()),
            "coefficient_sha256": shifted_digest,
        },
        "sturm": sturm_records,
        "mutation_103_over_100_numerator": str(mutation_value),
    }


def corrected_joint_floor() -> dict:
    r0 = 504
    normalizer_factor = 1 - F(2, 4 * r0 - 3)
    quadratic_factor = 1 - F(2, r0)
    factor = normalizer_factor * quadratic_factor
    floor = JOINT_CORE * factor**2
    margin = floor - TARGET

    r = sp.symbols("r", integer=True, positive=True)
    first = 1 - sp.Rational(2, 1) / (4 * r - 3)
    second = 1 - sp.Rational(2, 1) / r
    first_step_num = sp.factor(sp.together(first.subs(r, r + 1) - first).as_numer_denom()[0])
    second_step_num = sp.factor(sp.together(second.subs(r, r + 1) - second).as_numer_denom()[0])
    gate("normalizer correction increases in r", first_step_num == 8)
    gate("quadratic correction increases in r", second_step_num == 2)
    gate("joint floor retains row corrections", factor == F(504_761, 507_276), str(factor))
    gate("corrected joint floor exceeds 203/200", margin > 0, str(margin))
    return {
        "row": r0,
        "row_free_core": f"{JOINT_CORE.numerator}/{JOINT_CORE.denominator}",
        "normalizer_factor": f"{normalizer_factor.numerator}/{normalizer_factor.denominator}",
        "quadratic_factor": f"{quadratic_factor.numerator}/{quadratic_factor.denominator}",
        "attached_factor": f"{factor.numerator}/{factor.denominator}",
        "application": "row_free_core * attached_factor^2",
        "floor": f"{floor.numerator}/{floor.denominator}",
        "floor_decimal_truncated": decimal_truncation(floor),
        "margin_over_target": f"{margin.numerator}/{margin.denominator}",
    }


def top_floor_certificate() -> dict:
    # A15-TOP gives r*psi above this product throughout the top region.
    # Both factors increase, so the integer half-line r>=4 is anchored at 4.
    r = sp.symbols("r", integer=True, positive=True)
    normalizer = 1 - sp.Rational(2, 1) / (4 * r - 3)
    quadratic = 2 - 1 / r
    normalizer_step = sp.factor(
        sp.together(normalizer.subs(r, r + 1) - normalizer).as_numer_denom()[0]
    )
    quadratic_step = sp.factor(
        sp.together(quadratic.subs(r, r + 1) - quadratic).as_numer_denom()[0]
    )
    bridge = F(normalizer.subs(r, 4) * quadratic.subs(r, 4))
    gate("top normalizer factor increases in r", normalizer_step == 8)
    gate("top quadratic factor increases in r", quadratic_step == 1)
    gate("top bridge anchored at r=4", bridge == F(77, 52), str(bridge))
    floor = bridge**2 / 2
    gate("top floor exact value", floor == F(5_929, 5_408), str(floor))
    gate("top floor exceeds 203/200", floor > TARGET, str(floor - TARGET))
    return {
        "dependency": "A15-TOP exact scalar lower bound",
        "anchor_row": 4,
        "r_psi_floor": f"{bridge.numerator}/{bridge.denominator}",
        "rho_floor": f"{floor.numerator}/{floor.denominator}",
        "rho_floor_decimal_truncated": decimal_truncation(floor),
    }


def level_crossing_certificate() -> dict:
    # P_0=1, P_1=u+5 and
    # P_(m+1)=(u+8m^2+12m+5)P_m-[2m(2m+1)]^2P_(m-1).
    # For row r=m+1, a_r=6[u]P_m/[1]P_m.
    old = [1, 0]
    row = [5, 1]
    m = 1
    found: dict[int, dict] = {}
    positive = True
    recurrence_identity = True
    recurrence_identity_rows = 0
    a504_gt_16 = None
    sys.path.insert(0, str(HERE))
    from A2_stable import exact_a

    while len(found) < len(THRESHOLD_ROWS):
        r = m + 1
        positive = positive and row[0] > 0 and row[1] > 0
        if r <= 12:
            recurrence_identity = recurrence_identity and F(6 * row[1], row[0]) == exact_a(r)
            recurrence_identity_rows += 1
        if r == 504:
            a504_gt_16 = 6 * row[1] - 16 * row[0] > 0
        for threshold in THRESHOLD_ROWS:
            difference = 6 * row[1] - threshold * row[0]
            if threshold not in found and difference >= 0:
                previous_difference = 6 * old[1] - threshold * old[0]
                digest = hashlib.sha256()
                digest.update(threshold.to_bytes(2, "big"))
                digest.update(r.to_bytes(8, "big"))
                digest.update(integer_encoding(previous_difference))
                digest.update(integer_encoding(difference))
                found[threshold] = {
                    "first_r_with_a_ge_threshold": r,
                    "previous_difference_negative": previous_difference < 0,
                    "current_difference_positive": difference > 0,
                    "previous_difference_digits": len(str(abs(previous_difference))),
                    "current_difference_digits": len(str(abs(difference))),
                    "witness_sha256": digest.hexdigest(),
                    "a_at_crossing_decimal_truncated": decimal_truncation(
                        F(6 * row[1], row[0]), 12
                    ),
                }
        if len(found) == len(THRESHOLD_ROWS):
            break
        diagonal = 8 * m * m + 12 * m + 5
        back = (2 * m * (2 * m + 1)) ** 2
        new = [
            diagonal * row[0] - back * old[0],
            row[0] + diagonal * row[1] - back * old[1],
        ]
        old, row = row, new
        m += 1

    gate("level recurrence stays positive through final crossing", positive)
    gate(
        "level recurrence matches normalized-row engine",
        recurrence_identity and recurrence_identity_rows == 11,
        f"rows={recurrence_identity_rows}",
    )
    gate("a_504 exceeds 16", a504_gt_16 is True)
    for threshold, expected_r in THRESHOLD_ROWS.items():
        record = found[threshold]
        gate(
            f"exact a_r crossing at {threshold}",
            record["first_r_with_a_ge_threshold"] == expected_r
            and record["previous_difference_negative"] is True
            and record["current_difference_positive"] is True,
            f"r={record['first_r_with_a_ge_threshold']}",
        )

    windows = {
        "6": {"r_lo": 504, "r_hi": found[30]["first_r_with_a_ge_threshold"] - 1},
        "7": {
            "r_lo": found[25]["first_r_with_a_ge_threshold"],
            "r_hi": found[35]["first_r_with_a_ge_threshold"] - 1,
        },
        "8": {
            "r_lo": found[36]["first_r_with_a_ge_threshold"],
            "r_hi": found[40]["first_r_with_a_ge_threshold"] - 1,
        },
    }
    gate(
        "three bounded windows",
        windows
        == {
            "6": {"r_lo": 504, "r_hi": 3_633},
            "7": {"r_lo": 1_468, "r_hi": 8_413},
            "8": {"r_lo": 9_884, "r_hi": 18_456},
        },
    )
    return {
        "dependency": "all-row strict increase of a_r from the proved first-lowering theorem",
        "crossings": {str(key): value for key, value in found.items()},
        "bounded_windows": windows,
    }


def finite_core_witness() -> dict:
    record = json.loads(CORE_RECORD.read_text(encoding="utf-8"))
    witness = record.get("finite_box_witness", {})
    valid = (
        record.get("passed") is True
        and record.get("cells") == 125_250
        and witness.get("cell") == [129, 2]
        and isinstance(witness.get("numerator_hex"), str)
        and isinstance(witness.get("denominator_hex"), str)
    )
    gate("existing exact core witness record", valid, CORE_RECORD.name)
    if not valid:
        return {}
    value = F(int(witness["numerator_hex"], 16), int(witness["denominator_hex"], 16))
    gate("core witness lies below 203/200", value < TARGET, decimal_truncation(value, 18))
    return {
        "cell": [129, 2],
        "value": fraction_fingerprint(value),
        "decimal_truncated": decimal_truncation(value, 18),
        "source_file": CORE_RECORD_LABEL,
        "source_sha256": sha256(CORE_RECORD),
        "role": "exact minimum of the replayed r<=503 core",
    }


def cross_row_quotient_lemma() -> dict:
    # mu_(k-1)=M, successive normalized L-increments p>0 and p+q with q>=0.
    # Then mu_k=M-p and mu_(k+1)=M-2p-q.
    M, p = sp.symbols("M p", positive=True)
    q = sp.symbols("q", nonnegative=True)
    difference = sp.expand((M - p) ** 2 - M * (M - 2 * p - q))
    gate("cross-row quotient strictness identity", difference == p**2 + M * q)
    return {
        "identity": "mu_k^2-mu_(k-1)mu_(k+1)=p^2+mu_(k-1)q>0",
        "ratio_relation": "X_(r,k)/X_(r-1,k)=mu_(k-1)/mu_k",
        "conclusion": "X_(r,k)/X_(r-1,k) is strictly increasing in k on its common domain",
        "dependencies": [
            "positive first lowering and convex L give p>0 and q>=0",
            "positive predecessor coefficients give mu_k>0 on the common domain",
        ],
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
    global OUTPUT, CORE_RECORD
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", type=Path, default=CORE_RECORD)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    CORE_RECORD = arguments.core
    OUTPUT = arguments.output
    sys.set_int_max_str_digits(1_000_000)
    print("A40 -- exact global-infimum reduction")
    joint_core = joint_core_certificate()
    joint_floor = corrected_joint_floor()
    top_floor = top_floor_certificate()
    crossings = level_crossing_certificate()
    core_witness = finite_core_witness()
    quotient = cross_row_quotient_lemma()
    passed = not FAILURES
    payload = {
        "schema": SCHEMA,
        "producer_file": Path(__file__).name,
        "producer_sha256": sha256(Path(__file__).resolve()),
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "target": f"{TARGET.numerator}/{TARGET.denominator}",
        "conclusion": (
            "global infimum reduced to four unbounded tails t=2..5 and three bounded exact windows"
            if passed
            else "reduction gate failure"
        ),
        "not_certified_here": [
            "rho>=203/200 on the four unbounded tails t=2..5",
            "rho>=203/200 on the three bounded t=6,7,8 windows",
            "final global-infimum join",
        ],
        "joint_row_free_core": joint_core,
        "joint_floor_with_row_corrections": joint_floor,
        "top_floor": top_floor,
        "level_crossings": crossings,
        "finite_core_witness": core_witness,
        "cross_row_quotient_lemma": quotient,
        "dependencies": {
            "A2_stable.py": sha256(HERE / "A2_stable.py"),
            "A13_delta5.py": sha256(HERE / "A13_delta5.py"),
            "A14_uniform_sc.py": sha256(HERE / "A14_uniform_sc.py"),
            "A15_TOP_COVERAGE.md": sha256(HERE / "A15_TOP_COVERAGE.md"),
            "A15_top_coverage.py": sha256(HERE / "A15_top_coverage.py"),
            "A12_ratio_sign.py": sha256(HERE / "A12_ratio_sign.py"),
            "A8_l_curvature.py": sha256(HERE / "A8_l_curvature.py"),
        },
        "gates": GATES,
        "failures": FAILURES,
    }
    write_json(OUTPUT, payload)
    print(f"output: {OUTPUT}  sha256={sha256(OUTPUT)}")
    if not passed:
        print(f"FAILED: {len(FAILURES)} gates")
        return 1
    print("REDUCED: four unbounded tails plus three bounded exact windows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
