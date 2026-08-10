#!/usr/bin/env python3
"""Fail-closed exact join for the square-tail global infimum.

The join consumes the strengthened four-tail certificate, the exact
top/joint reduction, the residual finite-window stream, and the frozen finite
core and terminal-cell records.  It does not modify the frozen replay package.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
OUTPUT = RESULTS / "global_infimum_join.json"
INPUTS = {
    "retargeted_tails": RESULTS / "global_infimum_tails.json",
    "top_joint_reduction": RESULTS / "global_infimum_reduction.json",
    "finite_windows": RESULTS / "global_infimum_windows.json",
    "finite_core": RESULTS / "finite_core.json",
    "terminal_cells": RESULTS / "terminal_cells.json",
    "core_uniqueness": RESULTS / "core_uniqueness.json",
}
SCHEMA = "square-tail-global-infimum-join-v1"
TARGET = F(203, 200)
EXPECTED_TAILS = {
    2: (503, 504),
    3: (503, 504),
    4: (503, 504),
    5: (1_000, 1_001),
}
EXPECTED_WINDOWS = {
    5: (504, 1_000),
    6: (504, 3_633),
    7: (1_468, 8_413),
    8: (9_884, 18_456),
}
EXPECTED_REDUCTION_WINDOWS = {
    6: (504, 3_633),
    7: (1_468, 8_413),
    8: (9_884, 18_456),
}
EXPECTED_CROSSINGS = {
    25: 1_468,
    30: 3_634,
    35: 8_414,
    36: 9_884,
    40: 18_457,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def integer_encoding(value: int) -> bytes:
    magnitude = abs(value).to_bytes(
        max(1, (abs(value).bit_length() + 7) // 8), "big"
    )
    return (
        (b"-" if value < 0 else b"+")
        + len(magnitude).to_bytes(8, "big")
        + magnitude
    )


def load(path: Path) -> dict:
    require(path.is_file(), f"missing input: {path}")
    with path.open(encoding="utf-8") as source:
        payload = json.load(source)
    require(isinstance(payload, dict), f"non-object input: {path}")
    return payload


def fraction(literal: str) -> F:
    require(isinstance(literal, str), f"expected rational string, got {literal!r}")
    return F(literal)


def require_all_true(gates: dict, label: str) -> None:
    require(isinstance(gates, dict) and gates, f"{label}: missing gates")
    false_gates = sorted(name for name, value in gates.items() if value is not True)
    require(not false_gates, f"{label}: failed gates: {false_gates}")


def validate_tails(payload: dict) -> None:
    require(payload.get("schema") == "square-tail-global-infimum-retarget-probe-v3", "tail schema")
    require(payload.get("algorithm") == "square-tail-global-infimum-retarget-components-v3", "tail algorithm")
    require(payload.get("producer_sha256") == sha256(HERE / "A39_global_infimum_retarget_probe.py"), "tail producer hash")
    require(payload.get("target") == "203/200", "tail target")
    require(payload.get("passed") is True and payload.get("status") == "PASS", "tail status")
    require(payload.get("errors") == [] and payload.get("missing") == [], "tail completeness")
    columns = payload.get("columns")
    require(isinstance(columns, dict), "tail columns")
    require(set(columns) == {str(t) for t in EXPECTED_TAILS}, "tail column set")
    for t, (anchor, claim_start) in EXPECTED_TAILS.items():
        record = columns[str(t)]
        require(record.get("t") == t, f"tail t={t}: column identity")
        require(record.get("target") == "203/200", f"tail t={t}: target")
        require(record.get("algorithm") == payload.get("algorithm"), f"tail t={t}: algorithm")
        require(record.get("producer_sha256") == payload.get("producer_sha256"), f"tail t={t}: producer")
        require(record.get("anchor_m0") == anchor, f"tail t={t}: anchor")
        require(record.get("claim_start_r") == claim_start, f"tail t={t}: start")
        require(record.get("passed") is True and record.get("status") == "PASS", f"tail t={t}: status")
        require(record.get("negative_coefficient_indices") == [], f"tail t={t}: negative coefficients")
        require(record.get("zero_coefficient_indices") == [], f"tail t={t}: zero coefficients")
        require(record.get("coefficient_capacity_vs_target") == "above", f"tail t={t}: capacity comparison")
        require(record.get("coefficient_count", 0) > 0, f"tail t={t}: coefficient count")
        require_all_true(record.get("gates"), f"tail t={t}")
    controls = payload.get("controls", {})
    require(controls.get("original_t2_anchor_rejected") is True, "original t=2 anchor control")
    control = controls.get("record", {})
    require(control.get("t") == 2 and control.get("anchor_m0") == 312, "original t=2 control identity")
    require(control.get("status") == "FAIL" and control.get("passed") is False, "original t=2 control status")
    require(control.get("coefficient_capacity_vs_target") == "below", "original t=2 control capacity")
    require(control.get("negative_coefficient_indices") == [0, 1], "original t=2 control coefficients")


def validate_reduction(payload: dict) -> None:
    require(payload.get("schema") == "square-tail-global-infimum-reduction-v1", "reduction schema")
    require(payload.get("producer_file") == "A40_global_infimum_reduction.py", "reduction producer file")
    require(payload.get("producer_sha256") == sha256(HERE / "A40_global_infimum_reduction.py"), "reduction producer hash")
    require(payload.get("target") == "203/200", "reduction target")
    require(payload.get("passed") is True and payload.get("status") == "PASS", "reduction status")
    require(payload.get("failures") == [], "reduction failures")
    require_all_true(payload.get("gates"), "reduction")
    dependencies = payload.get("dependencies")
    expected_dependencies = {
        name: sha256(HERE / name)
        for name in (
            "A2_stable.py",
            "A8_l_curvature.py",
            "A12_ratio_sign.py",
            "A13_delta5.py",
            "A14_uniform_sc.py",
            "A15_top_coverage.py",
            "A15_TOP_COVERAGE.md",
        )
    }
    require(dependencies == expected_dependencies, "reduction dependency hashes")

    joint = payload.get("joint_floor_with_row_corrections", {})
    normalizer_factor = 1 - F(2, 4 * 504 - 3)
    quadratic_factor = 1 - F(2, 504)
    expected_factor = normalizer_factor * quadratic_factor
    require(fraction(joint.get("row_free_core")) == F(513, 500), "row-free core")
    require(fraction(joint.get("normalizer_factor")) == normalizer_factor, "normalizer row correction")
    require(fraction(joint.get("quadratic_factor")) == quadratic_factor, "quadratic row correction")
    require(fraction(joint.get("attached_factor")) == expected_factor, "attached row correction")
    require(joint.get("application") == "row_free_core * attached_factor^2", "row-correction exponent")
    joint_floor = fraction(joint.get("floor"))
    require(joint_floor == F(513, 500) * expected_factor**2, "corrected joint floor identity")
    require(joint_floor > TARGET, "corrected joint floor target")
    require(fraction(joint.get("margin_over_target")) == joint_floor - TARGET, "joint margin")

    top = payload.get("top_floor", {})
    require(fraction(top.get("r_psi_floor")) == F(77, 52), "top rPsi floor")
    require(fraction(top.get("rho_floor")) == F(5929, 5408), "top rho floor")
    require(fraction(top.get("rho_floor")) > TARGET, "top target")

    level = payload.get("level_crossings", {})
    crossings = level.get("crossings")
    require(isinstance(crossings, dict), "crossing records")
    require(set(crossings) == {str(x) for x in EXPECTED_CROSSINGS}, "crossing set")
    for threshold, first_r in EXPECTED_CROSSINGS.items():
        record = crossings[str(threshold)]
        require(record.get("first_r_with_a_ge_threshold") == first_r, f"crossing {threshold}: row")
        require(record.get("previous_difference_negative") is True, f"crossing {threshold}: predecessor")
        require(record.get("current_difference_positive") is True, f"crossing {threshold}: current")

    windows = level.get("bounded_windows")
    require(isinstance(windows, dict), "reduction windows")
    require(set(windows) == {str(t) for t in EXPECTED_REDUCTION_WINDOWS}, "reduction window set")
    for t, (lo, hi) in EXPECTED_REDUCTION_WINDOWS.items():
        require(windows[str(t)] == {"r_lo": lo, "r_hi": hi}, f"reduction window t={t}")


def validate_windows(payload: dict) -> None:
    require(payload.get("schema") == "square-tail-global-infimum-windows-v1", "window schema")
    require(payload.get("algorithm") == "square-tail-global-infimum-window-row-v1", "window algorithm")
    require(payload.get("source_sha256") == sha256(HERE / "A41_global_infimum_windows.py"), "window producer hash")
    require(payload.get("target") == "203/200", "window target")
    require(payload.get("passed") is True and payload.get("status") == "PASS", "window status")
    require(payload.get("missing_rows") == [], "window missing rows")
    require(payload.get("unexpected_rows") == [], "window unexpected rows")
    require(payload.get("malformed_rows") == [], "window malformed rows")
    require(payload.get("mutation_controls", {}).get("reverse_diagonal_corruption_rejected") is True, "window mutation control")
    require(payload.get("identity_replay_cells") == 84, "window identity replay count")
    require(
        payload.get("configuration")
        == {
            "workers": 1,
            "worker_reason": "serial row recurrence; each row comparison is sub-second",
        },
        "window execution configuration",
    )

    ranges = payload.get("ranges")
    require(isinstance(ranges, list), "window ranges")
    observed = {
        record.get("t"): (record.get("r_lo"), record.get("r_hi"), record.get("cells"))
        for record in ranges
    }
    expected = {
        t: (lo, hi, hi - lo + 1)
        for t, (lo, hi) in EXPECTED_WINDOWS.items()
    }
    require(observed == expected, "window interfaces")
    expected_cells = sum(cells for _, _, cells in expected.values())
    expected_rows = len({r for lo, hi in EXPECTED_WINDOWS.values() for r in range(lo, hi + 1)})
    require(payload.get("cells") == expected_cells, "window cell total")
    require(payload.get("expected_cells") == expected_cells, "window expected-cell total")
    require(payload.get("rows") == expected_rows, "window row total")
    digest = payload.get("canonical_rows_sha256")
    require(isinstance(digest, str) and len(digest) == 64, "window canonical digest")


def validate_uniqueness(payload: dict, core_value: F) -> None:
    """Bind A43: the core minimum is attained at (129,2) and nowhere else."""
    require(payload.get("schema") == "square-tail-core-uniqueness-v1", "uniqueness schema")
    require(payload.get("algorithm") == "square-tail-core-uniqueness-block-v1", "uniqueness algorithm")
    require(payload.get("passed") is True, "uniqueness status")
    require(payload.get("cells") == payload.get("expected_cells") == 125_250, "uniqueness cell count")
    require(payload.get("cells_at_or_below_witness") == [[129, 2]], "uniqueness at-or-below set")
    require(payload.get("cells_equal_to_witness") == [[129, 2]], "uniqueness equality set")
    require(payload.get("unique_minimiser_on_core") is True, "core uniqueness verdict")
    require_all_true(payload.get("gates"), "uniqueness")
    runner = payload.get("runner_up", {})
    require(runner.get("separation_positive") is True, "uniqueness runner-up separation")
    require(runner.get("cell") == [128, 2], "uniqueness runner-up cell")
    require(
        payload.get("witness_fingerprints", {}).get("decimal_sha256")
        == hashlib.sha256(
            f"{core_value.numerator}/{core_value.denominator}".encode()
        ).hexdigest(),
        "uniqueness witness is the same exact rational as the core witness",
    )


def validate_core(payload: dict) -> F:
    require(payload.get("schema") == "square-tail-finite-core-v1", "core schema")
    require(payload.get("algorithm") == "square-tail-finite-core-block-v2", "core algorithm")
    require(payload.get("passed") is True, "core status")
    require(payload.get("configuration", {}).get("row_range") == [4, 503], "core row range")
    require(payload.get("cells") == sum(r - 3 for r in range(4, 504)), "core cell count")
    require(payload.get("malformed_blocks") == [] and payload.get("unexpected_blocks") == [], "core completeness")
    require(payload.get("mutation_controls", {}).get("reverse_diagonal_corruption_rejected") is True, "core mutation control")
    witness = payload.get("finite_box_witness", {})
    require(witness.get("cell") == [129, 2], "core witness cell")
    require("minimum among these replayed cells" in witness.get("role", ""), "core witness role")
    numerator = int(witness.get("numerator_hex"), 16)
    denominator = int(witness.get("denominator_hex"), 16)
    require(denominator > 0, "core witness denominator")
    value = F(numerator, denominator)
    require(
        hashlib.sha256(
            integer_encoding(value.numerator)
            + integer_encoding(value.denominator)
        ).hexdigest()
        == witness.get("exact_fraction_sha256"),
        "core witness fraction fingerprint",
    )
    require(F(1) < value < TARGET, "core witness separation")
    return value


def validate_terminal(payload: dict) -> None:
    require(payload.get("schema") == "square-tail-terminal-cells-v1", "terminal schema")
    require(payload.get("passed") is True and payload.get("status") == "PASS", "terminal status")
    require(payload.get("failures") == [], "terminal failures")
    symbolic = payload.get("symbolic", {})
    require(symbolic.get("status") == "PASS", "terminal symbolic status")
    require(symbolic.get("top", {}).get("range") == "integer r>=2", "top terminal range")
    require(symbolic.get("top", {}).get("identity") == "rho(r,r)=r^2/2", "top terminal identity")
    require(symbolic.get("top_minus_one", {}).get("conclusion") == "rho(r,r-1)>2r^2/9>=2 for every integer r>=3", "top-minus-one terminal bound")


def write_summary(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def run() -> int:
    payloads = {name: load(path) for name, path in INPUTS.items()}
    validate_tails(payloads["retargeted_tails"])
    validate_reduction(payloads["top_joint_reduction"])
    validate_windows(payloads["finite_windows"])
    core_value = validate_core(payloads["finite_core"])
    validate_terminal(payloads["terminal_cells"])
    validate_uniqueness(payloads["core_uniqueness"], core_value)

    mutation = deepcopy(payloads["finite_windows"])
    for record in mutation["ranges"]:
        if record.get("t") == 8:
            record["r_lo"] += 1
            record["cells"] -= 1
    try:
        validate_windows(mutation)
    except RuntimeError:
        gap_mutation_rejected = True
    else:
        gap_mutation_rejected = False
    require(gap_mutation_rejected, "coverage-gap mutation was accepted")

    result = {
        "schema": SCHEMA,
        "producer_file": Path(__file__).name,
        "producer_sha256": sha256(Path(__file__).resolve()),
        "status": "PASS",
        "passed": True,
        "domain": "integers r>=2 and 2<=t<=r",
        "target": "203/200",
        "theorem": (
            "min rho(r,t) = rho(129,2), attained at (129,2) and nowhere else; "
            "equivalently the sharp constant in d^2 >= (c/r^2) S V is "
            "c = 2 rho(129,2), attained at a unique cell"
        ),
        "uniqueness": (
            "global, and asserted only here: core_uniqueness.json supplies "
            "uniqueness on 4<=r<=503 and explicitly refuses the global claim; "
            "this join supplies the outside-core floor rho>203/200>rho(129,2) "
            "and composes the two. The dependency is one-way -- the core "
            "certificate does not consume this join"
        ),
        "minimum": {
            "cell": [129, 2],
            "exact_fraction_sha256": payloads["finite_core"]["finite_box_witness"]["exact_fraction_sha256"],
            "reporting_decimal_truncated": payloads["finite_core"]["finite_box_witness"]["reporting_decimal_truncated"],
            "strictly_below_target": core_value < TARGET,
        },
        "coverage": {
            "terminal_cells": "t=r for r>=2; t=r-1 for r>=3",
            "finite_core": "4<=r<=503 and 2<=t<=r-2",
            "fixed_columns": {
                "2": "analytic tail r>=504",
                "3": "analytic tail r>=504",
                "4": "analytic tail r>=504",
                "5": "exact window 504<=r<=1000; analytic tail r>=1001",
                "6": "exact window 504<=r<=3633; joint region r>=3634",
                "7": "top region 504<=r<=1467; exact window 1468<=r<=8413; joint region r>=8414",
                "8": "top region 504<=r<=9883; exact window 9884<=r<=18456; joint region r>=18457",
            },
            "remaining_interior": "for t>=9, the exact top/joint partition is exhaustive",
            "outside_core_bound": "rho(r,t)>203/200>rho(129,2)",
        },
        "exact_floors": {
            "joint_with_row_corrections": payloads["top_joint_reduction"]["joint_floor_with_row_corrections"]["floor"],
            "top": payloads["top_joint_reduction"]["top_floor"]["rho_floor"],
        },
        "counts": {
            "finite_core_cells": payloads["finite_core"]["cells"],
            "residual_window_cells": payloads["finite_windows"]["cells"],
            "residual_window_rows": payloads["finite_windows"]["rows"],
            "unbounded_fixed_column_tails": 4,
            "bounded_reduction_windows": 3,
        },
        "inputs": {
            name: {
                "schema": payloads[name].get("schema"),
                "sha256": sha256(path),
                "passed": payloads[name].get("passed"),
            }
            for name, path in INPUTS.items()
        },
        "mutation_controls": {
            "t8_window_gap_rejected": gap_mutation_rejected,
        },
        "not_certified_here": [
            "the separate trigamma-tail route",
        ],
    }
    write_summary(OUTPUT, result)
    print("[PASS] global infimum join")
    print("minimum: rho(129,2) = 1.014576164564870742... (truncated)")
    print(f"certificate: {OUTPUT}")
    return 0


def main() -> int:
    global OUTPUT
    parser = argparse.ArgumentParser()
    for name, default in INPUTS.items():
        parser.add_argument(f"--{name.replace('_', '-')}", type=Path, default=default)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    for name in list(INPUTS):
        INPUTS[name] = getattr(arguments, name)
    OUTPUT = arguments.output
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
