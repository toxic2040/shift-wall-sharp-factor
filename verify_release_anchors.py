#!/usr/bin/env python3
"""Fail-closed crosscheck of Paper 4's release anchors.

The full producer replay is orchestrated by verify_all.py.  This verifier
checks the resulting exact fields against the manuscript and release boundary.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.set_int_max_str_digits(1_000_000)

HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE
BRIDGE = RELEASE_ROOT / "verification" / "lower_stripes"
ALLT = RELEASE_ROOT / "verification" / "all_stripes"
A1 = (RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay" /
      "verification" / "A1_verify_tail_independent.py")
A1_NOTE = RELEASE_ROOT / "docs" / "FIXED_COLUMN_TAIL_T2.md"
TEX = RELEASE_ROOT / "paper" / "shift_wall_sharp_factor.tex"
MONOTONICITY_SCRIPT = BRIDGE / "verify_monotonicity_witness.py"
MONOTONICITY_RECORD = BRIDGE / "MONOTONICITY_WITNESS.json"

EXPECTED_FILE_HASHES = {
    BRIDGE / "BRIDGE_T2_MINIMUM.json":
        "ada93210b9b2b62f0c65e3fe69368a2a4de0f9e970da552a6431a47b40436640",
    BRIDGE / "BRIDGE_T3_BARRIER.json":
        "2913edd4aefe5d6cf18d7e57d1e993374ef43f75a92a9dbf7fbb7f8984cc9f14",
    BRIDGE / "SQ_T4_BRIDGE_TRANSFER.json":
        "813a4fd8ae6c2c9beec98c2796126de4e2225ea786b6de892438b1cbbc10d27e",
    ALLT / "join" / "ALLT_JOIN.json":
        "fcf6e9b6a58c800155cdc41a930b1f8fb373a4f546b43f29d7748cee28915070",
    ALLT / "join" / "verify_allt_join.py":
        "773152168c9c9d66b9baf2ebfbc284ac55ab754daa908aedd613c4392c996499",
    ALLT / "core" / "A_CORE_SWEEP_SUMMARY.json":
        "8f06ea626948f6cdd1171edbcc7a52d1f4fe992f946883a358277ae5b2f07095",
    ALLT / "join" / "verify_joint_core_t11.json":
        "684187d33c881923d4dd2538fc94ce7a34c2eaa962f21e1e921d10ed0297fa12",
    ALLT / "join" / "verify_joint_core_bidegree.json":
        "a4a79d659b822a8f3d87240510cbd644fd1e00f337950a119cf419015b375eb1",
    ALLT / "columns" / "H2_A_GATE.json":
        "c3fdef648550a2dade7f6fa6551135daf4ebf8420cd8652d4a1fd0dc8fe27b84",
    ALLT / "join" / "dlaw_leg1_prefix.py":
        "1ea01115e7abddb1512de5c62612676347bd639169d421e9b7517343f709c0e0",
    ALLT / "join" / "dlaw_leg1_prefix_anchor.json":
        "add4e1e03a3de48312d2ebb429b9f10d741a781789f009162950108f609ae47f",
    ALLT / "join" / "dlaw_leg2_bootstrap.py":
        "b99ad1b59f22b4a4a2e47aac3d6202129791ae87dfdc669cd8df42821c3755b2",
    ALLT / "join" / "DLAW_LEG2.json":
        "269ea38f1c17e3974b1c2c509f8f77f3677a4e175e0e7864fc52ed75c3c1e4d9",
    ALLT / "join" / "DLAW_LEG2_CONTROL.json":
        "b3b01b5d03ebf4ebe140365ecb80318a76c26947c49cf067893bec1010249903",
    ALLT / "join" / "verify_joint_core_sharp.py":
        "04aaf9457b8257d06a57f714741baa088fbe53b1f45b011cdc38339d4bd2e147",
    ALLT / "join" / "verify_joint_core_sharp.json":
        "b8019109d39e41c592e0bd85e26f9cdedc2edd4b704bdda3bcb4cdedb19afcee",
    A1:
        "ad7afc6ac798c7c2db565b2358325fcf9a68de3926b1596d24a16f358a1a6ec6",
    A1_NOTE:
        "8153e46b963def0943897071a1daaeefe4b4fc2d9f72d2d29890431a5c1d8f20",
    MONOTONICITY_SCRIPT:
        "1117b01c9e82fd53fee0b2e7d4c85127392fe0acd5c72d07a889bbc195922555",
    MONOTONICITY_RECORD:
        "5ef630ed42cd889c8c64372cebe4787a1bdaae27aadbb2b51e7898a9ddeb1c7e",
}

EXPECTED_COLUMNS = {
    5: ([5, 1000], 1001, 1000, 46,
        "d24df0232421bce4fef034648cd1808398cefc6554a6a6a8113ccf37b17c17cc"),
    6: ([6, 4006], 4007, 4006, 54,
        "bcfea50d53f8bcc20d8b6eab5e32cc9e3d86b6b2ca6b6dac19c7c034b6f7cff5"),
    7: ([7, 12398], 12399, 12398, 62,
        "4c7ec81b4dbfd66ae8a8fe6c5c23cc2fe610eb602b4c169cfcf52cad649641fa"),
    8: ([8, 19999], 20000, 19999, 70,
        "b9de8ff893babcc7e742b0a47732c34b18adeab744e08b7c93893c8fb71205b7"),
}

EXPECTED_CEILING_COLUMNS = {
    9: (68133, 68134, 78,
        "4d5499dcf7547a35c998149948df84a5e35ba40c9015ca8588cfad7e88e4276c"),
    10: (200000, 200001, 86,
         "e3cb2e78a0425e13d9617fee99e36495e6a949011e2c748288f810a472e0eaf6"),
}

failures: list[str] = []
passes = 0


def gate(name: str, condition: bool) -> None:
    global passes
    if condition:
        passes += 1
        print(f"[PASS] {name}")
    else:
        failures.append(name)
        print(f"[FAIL] {name}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> int:
    gate("all pinned result files are present",
         all(path.is_file() for path in EXPECTED_FILE_HASHES))
    for path, expected in EXPECTED_FILE_HASHES.items():
        gate(f"pinned digest: {path.name}", path.is_file() and digest(path) == expected)

    t2 = load(BRIDGE / "BRIDGE_T2_MINIMUM.json")
    gate("t=2 unique minimum record",
         t2["status"] == "PASS" and t2["passed"] is True
         and t2["minimizer_r"] == 1350
         and t2["prefix"]["equal_rows"] == [1350]
         and not t2["prefix"]["less_rows"])
    gate("t=2 target fingerprint",
         t2["target"]["sha256"]
         == "2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691")
    gate("sharp raw constant decimal",
         t2["raw_optimal_bridge_constant"]["decimal_truncated"]
         == "2.655018313913361199726406175676")

    t3 = load(BRIDGE / "BRIDGE_T3_BARRIER.json")
    gate("t=3 direct four-thirds barrier",
         t3["status"] == "PASS" and t3["passed"] is True
         and t3["c_star"]["exactly_below_barrier"] is True
         and not t3["prefix"]["non_strict_rows"])

    t4 = load(BRIDGE / "SQ_T4_BRIDGE_TRANSFER.json")
    gate("t=4 rho barrier and bridge transfer",
         t4["status"] == "PASS" and t4["passed"] is True
         and t4["sq_prefix"]["cell_count"] == 997
         and t4["sq_tail"]["retargeted_positive_coefficients"] == 39)

    monotonicity = load(MONOTONICITY_RECORD)
    gate("row-45082 exact monotonicity counterexample",
         monotonicity["schema"] == "shift-wall-monotonicity-witness-v1"
         and monotonicity["status"] == "PASS"
         and monotonicity["passed"] is True
         and monotonicity["r"] == 45082
         and monotonicity["sign_D"] == -1
         and monotonicity["strict_inequality_holds"] is True
         and all(monotonicity["denominator_minors_positive"].values())
         and monotonicity["definition_checks"]
             ["top_window_matches_full_definition"] is True
         and monotonicity["definition_checks"]
             ["recurrence_mutation_rejected"] is True
         and monotonicity["witness_bit_length"] == 5416154
         and monotonicity["witness_decimal_digits"] == 1630425
         and monotonicity["witness_bytes_length"] == 677020
         and monotonicity["witness_sha256"]
         == "0aed686f60e5eaf6056db6cb8e1f59aa99610dc26293b7f9f57fb53d2f8ed47e")

    core = load(ALLT / "core" / "A_CORE_SWEEP_SUMMARY.json")
    gate("finite core count, sign, and digest",
         core["status"] == "PASS" and core["total_cells"] == 123753
         and not core["failures"] and not core["degenerate"]
         and core["comparison_stream_sha256"]
         == "83797ef82eca71c6409a61f942bd81894d7352cf557c2bae2055bd40b9216402"
         and (core["global_minimum"]["r"], core["global_minimum"]["t"])
         == (503, 5))

    for t, (prefix, overlap, anchor, degree, coeff_hash) in EXPECTED_COLUMNS.items():
        rec = load(ALLT / "columns" / f"H2_COL_T{t}.json")
        gate(f"column t={t} exact prefix-tail join",
             rec["status"] == "PASS" and rec["passed"] is True
             and not rec["failures"] and all(rec["gates"].values())
             and rec["coverage"]["prefix"] == prefix
             and rec["coverage"]["overlap_cell"] == overlap
             and rec["tail"]["anchor_m0"] == anchor
             and rec["tail"]["retargeted_degree"] == degree
             and rec["tail"]["retargeted_sha256"] == coeff_hash)

    for t, (anchor, start, degree, coeff_hash) in EXPECTED_CEILING_COLUMNS.items():
        rec = load(ALLT / "columns" / f"H2_COL_T{t}.json")
        gate(f"column t={t} top-gate-tail join",
             rec["status"] == "TAIL_PLUS_TOP_GATE_PROVED"
             and rec["passed"] is True and not rec["failures"]
             and all(rec["gates"].values())
             and rec["tail"]["anchor_m0"] == anchor
             and rec["coverage"]["tail"] == f"r>={start}"
             and rec["tail"]["retargeted_degree"] == degree
             and rec["tail"]["retargeted_sha256"] == coeff_hash)

    joint = load(ALLT / "join" / "verify_joint_core_t11.json")
    bide = load(ALLT / "join" / "verify_joint_core_bidegree.json")
    gate("joint Sturm range t=11..16",
         joint["status"] == "PASS" and not joint["failures"]
         and all(joint["verdicts"][str(t)] == "sturm" for t in range(11, 17)))
    gate("uniform joint tail t>=17",
         bide["status"] == "PASS" and bide["T0"] == 17
         and bide["monomial_count"] == 105
         and bide["negative_monomials"] == 0
         and bide["constant_positive"] is True
         and bide["coefficient_sha256"]
         == "9726005502a811c97584bc22b98eeced6577d4b18a5c3c45df3eb178239c69b4")

    leg1 = load(ALLT / "join" / "dlaw_leg1_prefix_anchor.json")
    gate("D-LAW leg 1 exact ladder and anchor",
         leg1["schema"] == "dlaw-leg1-anchor-v1"
         and leg1["status"] == "PASS"
         and leg1["m_enforce"] == 60000 and leg1["m_pre"] == 262144
         and leg1["scale_bits"] == 192 and not leg1["failures"]
         and set(leg1["worst_margins"]) == {str(k) for k in range(5, 12)}
         and all(value[1] > 1 for value in leg1["worst_margins"].values()))

    leg2 = load(ALLT / "join" / "DLAW_LEG2.json")
    windows = leg2["windows"]
    gate("D-LAW leg 2 has 57 contiguous passing windows",
         leg2["schema"] == "dlaw-leg2-single-anchor-v2"
         and leg2["status"] == "PASS" and leg2["gate_constant"] == 216
         and leg2["m_lo"] == 262144 and leg2["h_step"] == "7/20"
         and leg2["window_count"] == len(windows) == 57
         and windows[0]["x"] == [0.0, 0.35]
         and windows[-1]["x"] == [19.6, 19.95]
         and all(item["window"] == i for i, item in enumerate(windows))
         and all(windows[i]["x"][1] == windows[i + 1]["x"][0]
                 for i in range(56))
         and all(item["envelope_positivity"] is True
                 and len(item["gates"]) == 7
                 and all(item["gates"].values()) for item in windows)
         and all(windows[i]["min_rel_margin"]
                 <= windows[i + 1]["min_rel_margin"] for i in range(56)))

    control = load(ALLT / "join" / "DLAW_LEG2_CONTROL.json")
    control_window = control["windows"][0]
    gate("D-LAW 4/5 negative control rejected",
         control["schema"] == "dlaw-leg2-single-anchor-v2"
         and control["status"] == "CONTROL_REJECTED"
         and control["gate_constant"] == 96
         and control["window_count"] == 1
         and control_window["envelope_positivity"] is True
         and {k for k, passed in control_window["gates"].items()
              if not passed} == {str(k) for k in range(6, 12)})

    sharp = load(ALLT / "join" / "verify_joint_core_sharp.json")
    gate("t=9 sharpened-corridor conditional chain is discharged",
         sharp["t9_conditional"]["all_ok_primary"] is True
         and "to be pinned" not in sharp["t9_conditional"]["d_law"]
         and all(piece.get("all_ok", True)
                 for piece in sharp["t9_conditional"]["pieces"])
         and sharp["t10"]["verdict"] == "sturm")

    join = load(ALLT / "join" / "ALLT_JOIN.json")
    gate("master join 22/22 unconditional PASS",
         join["status"] == "PASS" and join["gate_count"] == 22
         and len(join["gates"]) == 22
         and all(item["passed"] for item in join["gates"])
         and any(item["name"].startswith("t=9 second proof")
                 for item in join["gates"])
         and not join["conditional_legs"] and not join["failures"])
    gate("master exact floors",
         join["constants"]["joint_floor"] == "254783667121/190614029760"
         and join["constants"]["needed_core"]
         == "343105253568/254783667121"
         and join["constants"]["top_floor"]
         == "4100936855929/2058631521408")

    tex = TEX.read_text()
    literals = [
        "2.655018313913361199726406175676",
        "2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691",
        "123{,}753",
        "4100936855929}{2058631521408",
        "254783667121}{190614029760",
        "343105253568}{254783667121",
        "68133",
        "200000",
        "Fifty-seven contiguous windows",
        "certified decrement shell",
        "$22$-gate exact join",
        "2\\le t\\le j-2",
        "10.5281/zenodo.21866366",
        "9726005502a811c97584bc22b98eeced6577d4b18a5c3c45df3eb178239c69b4",
    ]
    gate("all load-bearing manuscript literals are present",
         all(literal in tex for literal in literals))

    if failures:
        print(f"\nRESULT: FAIL ({len(failures)} failing gates, {passes} passing)")
        for name in failures:
            print(f"  - {name}")
        return 1
    print(f"\nRESULT: PASS ({passes} gates)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
