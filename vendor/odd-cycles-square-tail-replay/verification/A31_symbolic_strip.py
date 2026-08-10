#!/usr/bin/env python3
"""Fail-closed exact join for the complete square-tail coverage proof.

The join binds seven fixed-column tails, the targeted finite staircase, the
full-row core, the two terminal cells, and the t=2 decrement orientation.  For
r>=504 and t>=9, the remaining top/joint dichotomy is symbolic:

    t < sqrt(a_r)+2  =>  a_r > (t-2)^2 > 5t.

No legacy interval strip or unvalidated prior-column string is consumed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re

import sympy as sp


TAIL_SPECS = {
    2: {"anchor": 312, "degree": 22, "controls": True},
    3: {"anchor": 115, "degree": 30, "controls": False},
    4: {"anchor": 249, "degree": 38, "controls": False},
    5: {"anchor": 699, "degree": 46, "controls": False},
    6: {"anchor": 2_003, "degree": 54, "controls": False},
    7: {"anchor": 6_199, "degree": 62, "controls": False},
    8: {"anchor": 19_999, "degree": 70, "controls": False},
}
TAIL_GATES = {
    "P and Q upper",
    "anchor factor",
    "anchor width",
    "certificate",
    "certificate strict",
    "continuum lower",
    "index polynomial",
    "interval_ladder",
    "numerator lower",
    "optional_controls",
    "pi bracket",
    "predecessor lower",
}
PREFIX_RANGES = [
    {"t": 5, "r_lo": 504, "r_hi": 699, "cells": 196},
    {"t": 6, "r_lo": 504, "r_hi": 2_003, "cells": 1_500},
    {"t": 7, "r_lo": 504, "r_hi": 6_199, "cells": 5_696},
    {"t": 8, "r_lo": 504, "r_hi": 19_999, "cells": 19_496},
]
HEX64 = re.compile(r"[0-9a-f]{64}")
FAILURES: list[str] = []
GATES: dict[str, bool] = {}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def valid_digest(value) -> bool:
    return isinstance(value, str) and HEX64.fullmatch(value) is not None


def valid_positive_hex_integer(value) -> bool:
    if not isinstance(value, str) or not value or value.startswith("0"):
        return False
    try:
        decoded = int(value, 16)
    except ValueError:
        return False
    return decoded > 0 and format(decoded, "x") == value


def gate(name: str, condition: bool, detail: str = "") -> bool:
    passed = condition is True
    GATES[name] = passed
    print(
        f"  [{'PASS' if passed else 'FAIL'}] {name}"
        + (f"  {detail}" if detail else "")
    )
    if not passed:
        FAILURES.append(name)
    return passed


def load_record(path: Path, label: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        gate(f"read {label}", False, str(error))
        return {}
    gate(f"read {label}", isinstance(value, dict), path.name)
    return value if isinstance(value, dict) else {}


def descriptor(path: Path, record: dict) -> dict:
    return {
        "file": path.name,
        "sha256": sha256(path) if path.is_file() else None,
        "schema": record.get("schema"),
        "passed": record.get("passed"),
    }


def validate_tail(path: Path, t: int, spec: dict) -> tuple[dict, int]:
    record = load_record(path, f"fixed tail t={t}")
    config = record.get("configuration")
    polynomial = record.get("certificate_polynomial")
    tail_gates = record.get("gates")
    config = config if isinstance(config, dict) else {}
    polynomial = polynomial if isinstance(polynomial, dict) else {}
    tail_gates = tail_gates if isinstance(tail_gates, dict) else {}
    ok = (
        record.get("schema") == "square-tail-fixed-column-v1"
        and record.get("version") == 1
        and record.get("passed") is True
        and record.get("status") == "PASS"
        and config.get("t") == t
        and config.get("anchor_m0") == spec["anchor"]
        and config.get("interval_bits") == 256
        and config.get("workers") == 1
        and config.get("controls") is spec["controls"]
        and set(tail_gates) == TAIL_GATES
        and all(value is True for value in tail_gates.values())
        and polynomial.get("degree") == spec["degree"]
        and polynomial.get("coefficient_count") == spec["degree"] + 1
        and valid_digest(polynomial.get("sha256"))
        and record.get("tail_log_moment_max_degree") == 2 * t + 1
    )
    gate(
        f"fixed tail t={t}",
        ok,
        f"r>={spec['anchor'] + 1}, degree={polynomial.get('degree')}",
    )
    return descriptor(path, record), polynomial.get("coefficient_count", 0)


def validate_prefix(path: Path) -> tuple[dict, int]:
    record = load_record(path, "fixed-column prefixes")
    config = record.get("configuration")
    config = config if isinstance(config, dict) else {}
    ok = (
        record.get("schema") == "square-tail-fixed-prefix-v1"
        and record.get("algorithm") == "square-tail-fixed-prefix-row-v2"
        and record.get("passed") is True
        and config.get("workers") == 1
        and config.get("row_range") == [504, 19_999]
        and record.get("identity_replay_cells") == 84
        and record.get("ranges") == PREFIX_RANGES
        and record.get("rows") == 19_496
        and record.get("cells") == 26_888
        and record.get("missing_rows") == []
        and record.get("unexpected_rows") == []
        and record.get("malformed_rows") == []
        and valid_digest(record.get("canonical_rows_sha256"))
        and record.get("mutation_controls")
        == {"reverse_diagonal_corruption_rejected": True}
    )
    gate("fixed-column prefix staircase", ok, "26,888 cells")
    return descriptor(path, record), record.get("cells", 0)


def validate_core(path: Path) -> tuple[dict, int]:
    record = load_record(path, "finite core")
    config = record.get("configuration")
    witness = record.get("finite_box_witness")
    config = config if isinstance(config, dict) else {}
    witness = witness if isinstance(witness, dict) else {}
    ok = (
        record.get("schema") == "square-tail-finite-core-v1"
        and record.get("algorithm") == "square-tail-finite-core-block-v2"
        and record.get("passed") is True
        and config.get("row_range") == [4, 503]
        and config.get("block_size") == 25
        and record.get("identity_replay_cells") == 120
        and record.get("blocks") == 20
        and record.get("cells") == 125_250
        and record.get("malformed_blocks") == []
        and record.get("unexpected_blocks") == []
        and witness.get("cell") == [129, 2]
        and valid_positive_hex_integer(witness.get("numerator_hex"))
        and valid_positive_hex_integer(witness.get("denominator_hex"))
        and valid_digest(witness.get("exact_fraction_sha256"))
        and valid_digest(record.get("canonical_blocks_sha256"))
        and record.get("mutation_controls")
        == {"reverse_diagonal_corruption_rejected": True}
    )
    gate("finite full-row core", ok, "125,250 cells")
    return descriptor(path, record), record.get("cells", 0)


def validate_terminal(path: Path) -> dict:
    record = load_record(path, "terminal cells")
    symbolic = record.get("symbolic")
    replay = record.get("definition_replay")
    symbolic = symbolic if isinstance(symbolic, dict) else {}
    replay = replay if isinstance(replay, dict) else {}
    top = symbolic.get("top") if isinstance(symbolic.get("top"), dict) else {}
    edge = (
        symbolic.get("top_minus_one")
        if isinstance(symbolic.get("top_minus_one"), dict)
        else {}
    )
    ok = (
        record.get("schema") == "square-tail-terminal-cells-v1"
        and record.get("version") == 1
        and record.get("passed") is True
        and record.get("status") == "PASS"
        and record.get("failures") == []
        and symbolic.get("status") == "PASS"
        and top.get("identity") == "rho(r,r)=r^2/2"
        and top.get("range") == "integer r>=2"
        and edge.get("identity") == "D-3K=R(r-3)/(90r^2)"
        and edge.get("conclusion")
        == "rho(r,r-1)>2r^2/9>=2 for every integer r>=3"
        and replay.get("status") == "PASS"
        and replay.get("cells") == 21
    )
    gate("two terminal cells", ok, "symbolic half-lines; 21 replay cells")
    return descriptor(path, record)


def validate_orientation(path: Path) -> dict:
    record = load_record(path, "t=2 decrement orientation")
    normalization = record.get("normalization_replay")
    tail = record.get("tail_orientation")
    minimum = record.get("smallest_integer_numerator")
    normalization = normalization if isinstance(normalization, dict) else {}
    tail = tail if isinstance(tail, dict) else {}
    minimum = minimum if isinstance(minimum, dict) else {}
    ok = (
        record.get("schema") == "square-tail-t2-orientation-v1"
        and record.get("version") == 1
        and record.get("passed") is True
        and record.get("status") == "PASS"
        and record.get("cells") == 309
        and normalization.get("range") == "4<=r<=18"
        and normalization.get("failures") == []
        and record.get("sign_failures") == []
        and minimum == {"r": 4, "M": "15995366"}
        and tail.get("anchor_m0") == 312
        and tail.get("range") == "integer r>=313"
        and tail.get("passed") is True
        and tail.get("gate")
        == "V*G-corr has nonnegative coefficients and positive constant term"
    )
    gate("oriented t=2 decrement", ok, "309 finite signs plus r>=313 tail")
    return descriptor(path, record)


def symbolic_join() -> tuple[dict, bool]:
    x, t = sp.symbols("x t", integer=True, nonnegative=True)
    difference = sp.expand((t - 2) ** 2 - 5 * t)
    shifted = sp.Poly(sp.expand(difference.subs(t, x + 9)), x)
    identity = shifted.as_expr() == x**2 + 9 * x + 4
    positivity = identity and all(value > 0 for value in shifted.all_coeffs())
    gate("t>=9 top/joint partition", positivity, str(shifted.as_expr()))
    boundary = difference.subs(t, 8) == -4
    gate("boundary mutation t=8 rejected", boundary, str(difference.subs(t, 8)))
    return {
        "identity": "(t-2)^2-5t=(t-9)^2+9(t-9)+4",
        "range": "integer t>=9",
        "boundary_t8": "-4",
    }, boundary


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Bind every exact square-tail coverage certificate."
    )
    parser.add_argument("--tails-dir", type=Path, required=True)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--terminal", type=Path, required=True)
    parser.add_argument("--orientation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    tail_paths = {
        t: args.tails_dir / f"fixed_tail_t{t}.json" for t in TAIL_SPECS
    }
    input_paths = [
        *tail_paths.values(),
        args.prefix,
        args.core,
        args.terminal,
        args.orientation,
    ]
    if args.output.resolve() in {path.resolve() for path in input_paths}:
        print("FAILED: output path aliases an input certificate")
        return 1

    print("A31 -- complete symbolic coverage join")
    inputs: dict[str, object] = {"fixed_tails": {}}
    tail_coefficients = 0
    for t, spec in TAIL_SPECS.items():
        item, count = validate_tail(tail_paths[t], t, spec)
        inputs["fixed_tails"][str(t)] = item
        tail_coefficients += count

    inputs["fixed_prefixes"], prefix_cells = validate_prefix(args.prefix)
    inputs["finite_core"], core_cells = validate_core(args.core)
    inputs["terminal_cells"] = validate_terminal(args.terminal)
    inputs["t2_orientation"] = validate_orientation(args.orientation)

    finite_cells = prefix_cells + core_cells
    gate("tail coefficient total", tail_coefficients == 329, str(tail_coefficients))
    gate("finite SQ comparison total", finite_cells == 152_138, str(finite_cells))
    columns_complete = (
        all(TAIL_SPECS[t]["anchor"] + 1 <= 504 for t in (2, 3, 4))
        and all(
            item["r_lo"] == 504
            and item["r_hi"] == TAIL_SPECS[item["t"]]["anchor"]
            for item in PREFIX_RANGES
        )
    )
    gate("seven fixed columns have no finite gap", columns_complete)
    symbolic, boundary = symbolic_join()

    passed = not FAILURES
    payload = {
        "schema": "square-tail-symbolic-join-v2",
        "version": 2,
        "passed": passed,
        "status": "PASS" if passed else "FAIL",
        "claim": "complete square-tail coverage for integers r>=2 and 2<=t<=r",
        "columns": list(TAIL_SPECS),
        "tail_coefficients": tail_coefficients,
        "finite_sq_cells": finite_cells,
        "coverage": {
            "finite_core": "4<=r<=503, 2<=t<=r-2",
            "fixed_prefixes": PREFIX_RANGES,
            "fixed_tails": {
                str(t): f"r>={spec['anchor'] + 1}"
                for t, spec in TAIL_SPECS.items()
            },
            "terminal_cells": "t=r for r>=2; t=r-1 for r>=3",
            "top_or_joint": (
                "r>=504, t>=9: top if t>=sqrt(a_r)+2; otherwise "
                "a_r>(t-2)^2>5t and the joint theorem applies"
            ),
            "decrement_orientation": (
                "terminal cells explicit; t>=3 by profile comparison; "
                "t=2 by 309 finite signs plus its oriented tail"
            ),
        },
        "symbolic_join": symbolic,
        "inputs": inputs,
        "gates": GATES,
        "mutation_controls": {"boundary_t8_rejected": boundary},
        "failures": FAILURES,
    }
    write_json(args.output, payload)
    print(f"  output: {args.output}  sha256={sha256(args.output)}")
    if not passed:
        print(f"FAILED: {len(FAILURES)} checks")
        return 1
    print("CLOSED: core + seven columns + terminals + top/joint")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
