#!/usr/bin/env python3
"""Direction and blast-radius audit for the inherited A1 PI endpoint defect.

The pinned A1 constructors round both fixed-point endpoints down.  This audit
proves that the missing upper ulp is absorbed by the first operation in every
current arithmetic consumer, records coefficient-level endpoint-direction
tags for the new t=2 and t=3 bridge certificates, and installs an external
100-decimal-digit regression tripwire.
"""

from __future__ import annotations

import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
COMMON_PATH = HERE / "verify_bridge_t2_minimum.py"
DEFAULT_OUTPUT = HERE / "A1_PI_ENDPOINT_AUDIT.json"
SCHEMA = "a1-pi-endpoint-direction-audit-v1"
ALGORITHM = "ast-consumer-inventory-plus-exact-endpoint-transfer-v1"
M0 = 3_000

EXPECTED_SOURCE_SHA256 = {
    "verify_all.py":
        "fd568bada17e9dce4b9248135aa6c41b10a780993d11f1e3e2c52c8aa7b67535",
    "verification/A1_verify_tail.py":
        "7cc90454d04a71b2eecc6123eee98c82896f2809ff1514d0262125732bb50958",
    "verification/A1_verify_tail_independent.py":
        "ad7afc6ac798c7c2db565b2358325fcf9a68de3926b1596d24a16f358a1a6ec6",
    "verification/A13_delta5.py":
        "b15c9ac43b051f8eca00ca37fef13be52e920deed8dd710b0ac2672a9d52f98c",
    "verification/A35_t2_decrement_sign.py":
        "34649e7ef042176e9bc56255ab9469958d6758bd9734ff5cbb4371c9dbc9104b",
    "verification/A39_global_infimum_retarget_probe.py":
        "475a3d9a3dc2e569ef255ff122dda7ed0b9874ec0eb47fac2c8197ba50d6e425",
}

BOUND_ROLES = {
    "blo": "lower envelope; an upward error can be unsafe",
    "bhi": "upper envelope; a downward error can be unsafe",
    "clo": "lower envelope; an upward error can be unsafe",
    "chi": "upper envelope; a downward error can be unsafe",
    "Plo": "lower numerator; an upward error can be unsafe",
    "Qlo": "positive lower denominator; invalid upward movement can be unsafe",
    "Qhi": "upper denominator; a downward error can be unsafe",
    "H": "lower surplus Plo-2*c*Qhi; an upward Plo or downward Qhi is unsafe",
}


def load_module(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def load_common():
    return load_module(COMMON_PATH, "bridge_endpoint_common")


def ceil_fraction_scaled(value: F, scale: int) -> int:
    numerator = value.numerator * scale
    return -((-numerator) // value.denominator)


def floor_fraction_scaled(value: F, scale: int) -> int:
    return value.numerator * scale // value.denominator


def half_scaled_bounds(lo: int, hi: int) -> tuple[int, int]:
    return lo // 2, -((-hi) // 2)


def pi_reference_cell(common) -> tuple[F, F]:
    common.require(len(common.PI_100_DECIMAL_DIGITS) == 101,
                   "pi reference does not have 100 decimal places")
    denominator = 10**100
    numerator = int(common.PI_100_DECIMAL_DIGITS)
    return F(numerator, denominator), F(numerator + 1, denominator)


def constructor_audit(common, module, scale_name: str) -> dict:
    scale = getattr(module, scale_name)
    defective_lo = floor_fraction_scaled(module.PI_LO, scale)
    defective_hi = floor_fraction_scaled(module.PI_HI, scale)
    corrected_hi = ceil_fraction_scaled(module.PI_HI, scale)
    common.require(defective_lo == defective_hi,
                   "defective PI endpoints no longer collapse")
    common.require(corrected_hi == defective_hi + 1,
                   "PI repair is not exactly one upper ulp")
    common.require(defective_hi % 4 == 1,
                   "defective PI endpoint residue changed")
    defective_half = half_scaled_bounds(defective_lo, defective_hi)
    corrected_half = half_scaled_bounds(defective_lo, corrected_hi)
    common.require(defective_half == corrected_half,
                   "missing PI upper ulp survives the half scaling")

    reference_lo, reference_hi = pi_reference_cell(common)
    defective_interval = module.Iv(defective_lo, defective_hi)
    corrected_interval = module.Iv(defective_lo, corrected_hi)
    common.require(module.PI_LO <= reference_lo <= reference_hi <= module.PI_HI,
                   "Machin bracket does not contain the hard pi cell")
    common.require(
        not (
            defective_interval.flo() < reference_lo
            and defective_interval.fhi() > reference_hi
        ),
        "defective PI singleton passed the hard-reference tripwire",
    )
    common.require(
        corrected_interval.flo() < reference_lo
        and corrected_interval.fhi() > reference_hi,
        "corrected PI interval failed the hard-reference tripwire",
    )
    return {
        "scale": scale_name,
        "scale_bits": scale.bit_length() - 1,
        "defective_width_ulps": defective_hi - defective_lo,
        "corrected_width_ulps": corrected_hi - defective_lo,
        "upper_repair_delta_ulps": corrected_hi - defective_hi,
        "defective_upper_parity": "odd",
        "defective_upper_mod4": defective_hi % 4,
        "defective_half_pi_integer_endpoints": list(defective_half),
        "corrected_half_pi_integer_endpoints": list(corrected_half),
        "half_pi_exactly_invariant": defective_half == corrected_half,
        "fixed_pi_delta": f"1/{scale}",
        "half_pi_delta": "0/1",
        "downstream_D_times_delta": "0/1",
        "hard_reference_decimal_places": 100,
        "defective_reference_tripwire_rejected": True,
        "corrected_reference_tripwire_passed": True,
    }


def fixed_pi_load_sites(repo: Path) -> list[dict]:
    sites: list[dict] = []
    for path in sorted((repo / "verification").glob("*.py")):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        parents = {
            child: parent
            for parent in ast.walk(tree)
            for child in ast.iter_child_nodes(parent)
        }
        candidates: list[ast.AST] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and node.id == "PI" and isinstance(
                node.ctx, ast.Load
            ):
                candidates.append(node)
            elif (
                isinstance(node, ast.Attribute)
                and node.attr == "PI"
                and isinstance(node.ctx, ast.Load)
                and isinstance(node.value, ast.Name)
                and node.value.id == "a1"
            ):
                candidates.append(node)
        for node in candidates:
            attribute = parents.get(node)
            call = parents.get(attribute) if attribute is not None else None
            safe_half_scale = (
                isinstance(attribute, ast.Attribute)
                and attribute.attr == "scal"
                and isinstance(call, ast.Call)
                and call.func is attribute
                and len(call.args) == 2
                and all(isinstance(arg, ast.Constant) for arg in call.args)
                and [arg.value for arg in call.args] == [1, 2]
            )
            sites.append({
                "path": str(path.relative_to(repo)),
                "line": node.lineno,
                "source": ast.get_source_segment(source, call or node),
                "first_operation": "scal(1,2)" if safe_half_scale else "other",
                "missing_upper_ulp_absorbed": safe_half_scale,
            })
    return sites


def capture_tail(a1, *, t: int):
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
        result = a1.section3_tail(M0, t=t, verbose=False, record=False)
    finally:
        sys.settrace(prior_trace)
    required = {"blo", "bhi", "clo", "chi"}
    if not required <= captured.keys():
        raise RuntimeError("tail envelope locals were not captured")
    return result, captured


def bridge_polynomials(a1, captured: dict[str, object], *, t: int, target: F):
    blo, bhi = captured["blo"], captured["bhi"]
    clo, chi = captured["clo"], captured["chi"]
    zero = [F(0)]

    def get(rows, index):
        return rows[index] if 0 <= index < len(rows) else zero

    x_lower = a1.padd(get(blo, t - 1), get(blo, t))
    y_upper = a1.padd(get(bhi, t - 2), get(bhi, t - 1))
    z_upper = a1.padd(get(bhi, t), get(bhi, t + 1))
    p_lower = a1.padd(
        a1.pmul(x_lower, x_lower), a1.pmul(y_upper, z_upper), -1
    )
    q_upper = a1.padd(
        a1.pmul(get(chi, t), get(chi, t)),
        a1.pmul(get(clo, t - 1), get(clo, t + 1)),
        -1,
    )
    q_lower = a1.padd(
        a1.pmul(get(clo, t), get(clo, t)),
        a1.pmul(get(chi, t - 1), get(chi, t + 1)),
        -1,
    )
    h_poly = a1.padd(
        p_lower, [2 * target * value for value in q_upper], -1
    )
    return {
        "blo": blo,
        "bhi": bhi,
        "clo": clo,
        "chi": chi,
        "Plo": p_lower,
        "Qlo": q_lower,
        "Qhi": q_upper,
        "H": h_poly,
    }


def flatten_families(families: dict[str, object]) -> dict[tuple[str, int, int], F]:
    flattened: dict[tuple[str, int, int], F] = {}
    for family, values in families.items():
        if family in {"blo", "bhi", "clo", "chi"}:
            for level, row in enumerate(values):
                for index, value in enumerate(row):
                    flattened[(family, level, index)] = value
        else:
            for index, value in enumerate(values):
                flattened[(family, -1, index)] = value
    return flattened


def sign(value: F) -> str:
    return "+" if value > 0 else "-" if value < 0 else "0"


def coefficient_direction_audit(common, a1, *, t: int, target: F) -> dict:
    scale = a1.D
    p = floor_fraction_scaled(a1.PI_LO, scale)
    defective = a1.Iv(p, p)
    corrected = a1.Iv(p, p + 1)
    upper_probe = a1.Iv(p, p + 2)
    lower_probe = a1.Iv(p - 2, p + 1)

    captures = {}
    results = {}
    for label, interval in (
        ("defective", defective),
        ("corrected", corrected),
        ("upper_probe", upper_probe),
        ("lower_probe", lower_probe),
    ):
        a1.PI = interval
        result, captured = capture_tail(a1, t=t)
        results[label] = result
        captures[label] = flatten_families(
            bridge_polynomials(a1, captured, t=t, target=target)
        )

    keys = captures["corrected"].keys()
    common.require(all(captures[label].keys() == keys for label in captures),
                   "coefficient family shapes changed under endpoint probes")
    records = []
    actual_changes = []
    h_values: list[tuple[tuple[str, int, int], F, F]] = []
    for key in sorted(keys):
        family, level, index = key
        defective_value = captures["defective"][key]
        corrected_value = captures["corrected"][key]
        upper_delta = captures["upper_probe"][key] - corrected_value
        lower_delta = captures["lower_probe"][key] - corrected_value
        actual_delta = corrected_value - defective_value
        if actual_delta:
            actual_changes.append(key)
        if family == "H":
            h_values.append((key, corrected_value, upper_delta))
        records.append({
            "family": family,
            "level": None if level < 0 else level,
            "coefficient_index": index,
            "proof_role": BOUND_ROLES[family],
            "corrected_sign": sign(corrected_value),
            "corrected_sha256":
                common.fraction_fingerprint(corrected_value)["sha256"],
            "actual_missing_upper_ulp_delta": "0/1" if not actual_delta
                else common.fraction_literal(actual_delta),
            "actual_missing_upper_ulp_reaches_coefficient": bool(actual_delta),
            "first_surviving_upper_probe_delta_sign": sign(upper_delta),
            "first_surviving_lower_probe_delta_sign": sign(lower_delta),
            "endpoint_direction_tag": (
                (["upper"] if upper_delta else [])
                + (["lower"] if lower_delta else [])
            ),
        })

    common.require(not actual_changes,
                   f"actual PI repair changed coefficients: {actual_changes[:8]}")
    common.require(all(value > 0 for _, value, _ in h_values),
                   "corrected H has a nonpositive coefficient")
    vulnerable_h = [
        (key, value) for key, value, upper_delta in h_values if upper_delta
    ]
    min_all_key, min_all = min(
        ((key, value) for key, value, _ in h_values),
        key=lambda item: item[1],
    )
    min_vulnerable_key, min_vulnerable = min(vulnerable_h, key=lambda item: item[1])
    simple_vulnerable_bound = F(1, 20) if t == 2 else F(1, 25_000)
    common.require(min_vulnerable > simple_vulnerable_bound,
                   "upper-reachable H margin lost its simple bound")
    if t == 3:
        common.require(min_all == F(1, 60_963_840),
                       "universal top t=3 coefficient changed")
    return {
        "t": t,
        "anchor_m0": M0,
        "target": common.fraction_literal(target),
        "coefficient_count": len(records),
        "actual_changed_coefficients": 0,
        "actual_downstream_sensitivity_bound": "D*delta=0 exactly",
        "upper_probe": (
            "widen corrected PI.hi by one more fixed-point ulp; "
            "this is the first widening that changes half_pi.hi"
        ),
        "lower_probe": (
            "widen corrected PI.lo downward by two fixed-point ulps; "
            "this changes half_pi.lo by one fixed-point ulp"
        ),
        "minimum_H_coefficient": {
            "index": min_all_key[2],
            **common.fraction_fingerprint(min_all),
            "decimal_truncated": common.decimal_truncation(min_all, 24),
        },
        "minimum_upper_reachable_H_coefficient": {
            "index": min_vulnerable_key[2],
            **common.fraction_fingerprint(min_vulnerable),
            "decimal_truncated": common.decimal_truncation(
                min_vulnerable, 24
            ),
            "strict_simple_lower_bound":
                common.fraction_literal(simple_vulnerable_bound),
        },
        "source_gate_corrected_passed": results["corrected"][0] is True,
        "coefficients": records,
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

    source_hashes = {}
    for relative, expected in EXPECTED_SOURCE_SHA256.items():
        digest = common.sha256_file(repo / relative)
        common.require(digest == expected, f"source hash mismatch: {relative}")
        source_hashes[relative] = digest

    primary = common.load_a1(repo)
    independent = load_module(
        repo / "verification" / "A1_verify_tail_independent.py",
        "a1_endpoint_independent",
    )
    constructors = {
        "primary": constructor_audit(common, primary, "D"),
        "independent": constructor_audit(common, independent, "SCALE"),
    }
    sites = fixed_pi_load_sites(repo)
    common.require(len(sites) == 3, f"fixed PI consumer count changed: {len(sites)}")
    common.require(all(site["missing_upper_ulp_absorbed"] for site in sites),
                   "a fixed PI consumer does not begin with scal(1,2)")

    c_star = None
    for r, numerator, denominator in common.reverse_bridge_stream(
        common.R_STAR, t=2
    ):
        if r == common.R_STAR:
            c_star = F(numerator, denominator)
    common.require(c_star is not None, "could not rebuild bridge target")
    t2 = coefficient_direction_audit(common, primary, t=2, target=c_star)
    t3 = coefficient_direction_audit(common, primary, t=3, target=F(4, 3))

    gates = {
        "source_hashes_pinned": source_hashes == EXPECTED_SOURCE_SHA256,
        "both_defective_constructors_reproduced": all(
            item["defective_width_ulps"] == 0 for item in constructors.values()
        ),
        "hard_reference_tripwire_fix_two": all(
            item["defective_reference_tripwire_rejected"]
            and item["corrected_reference_tripwire_passed"]
            for item in constructors.values()
        ),
        "all_fixed_PI_consumers_halve_first": all(
            site["missing_upper_ulp_absorbed"] for site in sites
        ),
        "half_pi_exactly_invariant": all(
            item["half_pi_exactly_invariant"] for item in constructors.values()
        ),
        "t2_coefficients_exactly_invariant":
            t2["actual_changed_coefficients"] == 0,
        "t3_coefficients_exactly_invariant":
            t3["actual_changed_coefficients"] == 0,
    }
    common.require(all(gates.values()), "one or more endpoint-audit gates failed")

    payload = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "status": "PASS",
        "passed": True,
        "claim": (
            "The inherited fixed PI constructors are invalid, but their missing "
            "upper ulp has zero theorem/result blast radius in every current "
            "consumer because scal(1,2) absorbs it exactly."
        ),
        "source_repo": repo.relative_to(common.RELEASE_ROOT).as_posix(),
        "source_hashes": source_hashes,
        "producer_sha256": common.sha256_file(Path(__file__).resolve()),
        "gates": gates,
        "constructors": constructors,
        "fixed_PI_load_sites": sites,
        "direct_exact_endpoint_uses": {
            "PI_LO": "unchanged; used as certified rational lower endpoint",
            "PI_HI": (
                "unchanged; positive upper-clock and injection computations "
                "consume the exact rational upper endpoint"
            ),
        },
        "coefficient_direction_audits": [t2, t3],
        "deposited_battery_conclusion": (
            "No full battery rerun is needed for this one-ulp defect: the AST "
            "inventory proves every current fixed-PI arithmetic path begins "
            "with the identical defective/corrected half_pi interval. A13 is "
            "direct; A35 and A39 are indirect through the same two constructors."
        ),
        "claim_boundary": (
            "The standalone PI interval constructors remain defective and must "
            "not be reused. The zero blast-radius proof applies only while the "
            "pinned consumer inventory remains unchanged."
        ),
    }
    common.write_json(args.output.resolve(), payload)
    print(
        "[PASS] defective PI upper endpoint rejected; corrected one-ulp "
        "interval passes the 100-decimal tripwire"
    )
    print(
        f"[PASS] {len(sites)} fixed-PI load sites all halve first; "
        "defective and corrected half_pi are identical"
    )
    print(
        f"[PASS] coefficient deltas: t=2 {t2['actual_changed_coefficients']}/"
        f"{t2['coefficient_count']}, t=3 {t3['actual_changed_coefficients']}/"
        f"{t3['coefficient_count']}; wrote {args.output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
