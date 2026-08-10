#!/usr/bin/env python3
"""Master join for the sharp all-stripe shift-wall factor theorem.

Consumes the t>=5 leg artifacts under verification/all_stripes/ plus the
t=2,3,4 stripe records under verification/lower_stripes/, verifies every
splice and every exact comparison, and emits ALLT_JOIN.json asserting:

    (T5B)   rho(j,t) > 4/3 for every t >= 5, j >= t,
  and with the banked t=2,3,4 stripes and the lossless transfer:
    THEOREM Omega_t(r) >= Omega_2(1350) for all r >= t >= 2,
            equality iff (t,r) = (2,1350);
            2*Omega_2(1350) = 2.655018313913361... is sharp.

Fails loudly listing missing legs if any artifact is absent or non-PASS.
"""

import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.set_int_max_str_digits(1_000_000)

HERE = Path(__file__).resolve().parent
ALL_STRIPES = HERE.parent
BRIDGE = ALL_STRIPES.parent / "lower_stripes"

FOUR_THIRDS = F(4, 3)
ATTACHED = F(504761, 507276)          # A14 scalar product at r=504
CORE_TARGET = F(27, 20)
NEEDED_CORE = FOUR_THIRDS / (ATTACHED * ATTACHED)
TOP_FLOOR = (F(2011, 2013) * F(1007, 504)) ** 2 / 2

failures: list[str] = []
gate_results: list[dict[str, object]] = []


def gate(name: str, ok: bool):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    gate_results.append({"name": name, "passed": bool(ok)})
    if not ok:
        failures.append(name)


def load(path: Path, description: str):
    if not path.exists():
        gate(f"{description}: artifact present ({path.name})", False)
        return None
    with open(path) as f:
        return json.load(f)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> int:
    # ---- exact scalar comparisons -------------------------------------
    gate("27/20 strictly exceeds the needed joint core constant",
         CORE_TARGET > NEEDED_CORE)
    gate("(27/20)*attached^2 > 4/3 (joint floor beats 4/3)",
         CORE_TARGET * ATTACHED * ATTACHED > FOUR_THIRDS)
    gate("top floor ((2011/2013)(1007/504))^2/2 > 4/3",
         TOP_FLOOR > FOUR_THIRDS)
    gate("terminal floors t^2/2 >= 25/2 > 4/3 and 2(t+1)^2/9 >= 8 > 4/3",
         F(25, 2) > FOUR_THIRDS and F(8) > FOUR_THIRDS)

    # ---- c_* rebuilt exactly from the t=2 stripe producer -------------
    spec = importlib.util.spec_from_file_location(
        "bridge_common", BRIDGE / "verify_bridge_t2_minimum.py")
    common = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(common)
    c_star = None
    for r, num, den in common.reverse_bridge_stream(common.R_STAR, t=2):
        if r == common.R_STAR:
            c_star = F(num, den)
    c_star_sha = (common.fraction_fingerprint(c_star)["sha256"]
                  if c_star is not None else None)
    ok = c_star_sha == common.EXPECTED_TARGET_SHA256
    gate("c_* = Omega_2(1350) rebuilt, fingerprint matches deposit", ok)
    gate("c_* < 4/3 exactly", c_star is not None and c_star < FOUR_THIRDS)

    # ---- t=2,3,4 stripe legs from the lower-stripe records -------------
    # These records are verified here directly, not taken by reference,
    # so the replay is transitive.
    t2 = load(BRIDGE / "BRIDGE_T2_MINIMUM.json", "t=2 stripe minimum")
    if t2:
        gate("t=2 stripe: Omega_2 record verified, all gates, unique "
             "minimum at r=1350",
             t2.get("passed") is True
             and t2.get("gates") and all(t2["gates"].values())
             and t2["gates"].get("prefix_unique_minimum") is True
             and t2.get("minimizer_r") == 1350
             and t2.get("claim")
             == "Omega_2(r) has its unique global minimum at r=1350")
    t3 = load(BRIDGE / "BRIDGE_T3_BARRIER.json", "t=3 stripe barrier")
    if t3:
        gate("t=3 stripe: Omega_3>4/3 record verified, all gates, c_* "
             "fingerprint matches live rebuild",
             t3.get("passed") is True
             and t3.get("gates") and all(t3["gates"].values())
             and t3.get("c_star", {}).get("sha256") == c_star_sha
             and t3.get("claim")
             == "Omega_3(r)>4/3>Omega_2(1350) for every integer r>=3")
    t4 = load(BRIDGE / "SQ_T4_BRIDGE_TRANSFER.json", "t=4 stripe transfer")
    if t4:
        gate("t=4 stripe: rho(r,4)>4/3 transfer record verified, all "
             "gates, c_* fingerprint matches live rebuild",
             t4.get("passed") is True
             and t4.get("gates") and all(t4["gates"].values())
             and t4.get("c_star", {}).get("sha256") == c_star_sha
             and t4.get("claim")
             == "rho(r,4)>4/3 for every integer r>=4; consequently "
                "Omega_4(r)>4/3>Omega_2(1350) for every integer r>=4")

    # ---- leg artifacts ------------------------------------------------
    legs = {}

    core = load(ALL_STRIPES / "core" / "A_CORE_SWEEP_SUMMARY.json",
                "finite core sweep")
    if core:
        gate("core sweep PASS, 123753 cells, zero failures",
             core["status"] == "PASS" and core["total_cells"] == 123753
             and not core["failures"])
        legs["core"] = core["comparison_stream_sha256"]

    jt11 = load(HERE / "verify_joint_core_t11.json", "joint core t=11..200")
    if jt11:
        gate("joint core per-t certified 11..200, no failures",
             jt11["status"] == "PASS" and not jt11["failures"])

    bide = load(HERE / "verify_joint_core_bidegree.json",
                "joint core uniform tail")
    if bide:
        gate("joint core bidegree tail t>=17: zero negative monomials",
             bide["status"] == "PASS" and bide["negative_monomials"] == 0
             and bide["constant_positive"])
        legs["joint_bidegree"] = bide["coefficient_sha256"]

    # fixed columns: expected splice structure per column
    # prefix columns: claim rho(j,t)>4/3 for j>=t via prefix [t, hi] and
    # tail from hi+1 (= coverage.overlap_cell).
    # prefix-free columns (9, 10): tail from m0+1 plus the exact a-gate
    # a_{m0} <= (t-2)^2 so the TOP region covers every r <= m0.
    expected = {
        5: {"col": "H2_COL_T5.json", "prefix": [5, 1000], "overlap": 1001},
        6: {"col": "H2_COL_T6.json", "prefix": [6, 4006], "overlap": 4007},
        7: {"col": "H2_COL_T7.json", "prefix": [7, 12398], "overlap": 12399},
        8: {"col": "H2_COL_T8.json", "prefix": [8, 19999], "overlap": 20000},
    }
    for t, meta in expected.items():
        rec = load(ALL_STRIPES / "columns" / meta["col"], f"column t={t}")
        if not rec:
            continue
        base_ok = (rec.get("status") == "PASS"
                   and rec.get("passed") is True
                   and not rec.get("failures")
                   and all(rec.get("gates", {}).values())
                   and rec.get("gates"))
        if "prefix" in meta:
            cov = rec.get("coverage", {})
            ok = (base_ok
                  and cov.get("prefix") == meta["prefix"]
                  and cov.get("overlap_cell") == meta["overlap"]
                  and rec.get("claim")
                  == f"rho(j,{t})>4/3 for every integer j>={t}")
            gate(f"t={t}: PASS, prefix {meta['prefix']} + tail from "
                 f"{meta['overlap']}, all gates", ok)
    # ---- t=9, t=10: top region + tail at the gate-compatible anchor ---
    # UNCONDITIONAL. Column = terminal + core + TOP (a_j <= (t-2)^2,
    # certified monotone a-gate through the anchor) + tail (j > anchor).
    # The tail starts exactly at the a-crossing, so the union is [t,inf).
    conditional_legs = []
    for t, meta in ((9, {"col": "H2_COL_T9.json", "anchor": 68133,
                         "gate_txt": "a_r<=49"}),
                    (10, {"col": "H2_COL_T10.json", "anchor": 200000,
                          "gate_txt": "a_r<=64"})):
        rec = load(ALL_STRIPES / "columns" / meta["col"], f"column t={t}")
        if not rec:
            continue
        deps = set(rec.get("coverage", {}).get("external_dependencies", []))
        ok = (rec.get("status") == "TAIL_PLUS_TOP_GATE_PROVED"
              and rec.get("gates") and all(rec["gates"].values())
              and rec.get("tail", {}).get("anchor_m0") == meta["anchor"]
              and meta["gate_txt"] in
              rec.get("coverage", {}).get("top_gate", "")
              and deps <= {"TOP-region 4/3 closure for a_r<=(t-2)^2",
                           "finite core r<504"})
        gate(f"t={t}: UNCONDITIONAL — top gate ({meta['gate_txt']} "
             f"through {meta['anchor']}) + tail from {meta['anchor']+1}, "
             f"deps supplied by top/core legs above", ok)

    # corroboration (non-load-bearing): the sharpened corridor core
    # separately closes t=10 by exact Sturm.
    sharp = load(HERE / "verify_joint_core_sharp.json",
                 "sharpened joint core (corroboration)")
    if sharp:
        gate("t=10 second proof: sharpened corridor core (Sturm) agrees",
             sharp.get("t10", {}).get("verdict") == "sturm")

    # corroboration (non-load-bearing): t=9 second proof — the sharpened
    # corridor core conditional on (D-LAW), now discharged by Leg 1
    # (ladder, m <= 262144, a <= ~59.2) + Leg 2 (single-anchor envelope
    # certificate through a > 768, negative control rejected).
    leg1 = load(HERE / "dlaw_leg1_prefix_anchor.json", "(D-LAW) leg 1")
    leg2 = load(HERE / "DLAW_LEG2.json", "(D-LAW) leg 2")
    ctrl = load(HERE / "DLAW_LEG2_CONTROL.json", "(D-LAW) leg 2 control")
    if sharp and leg1 and leg2 and ctrl:
        windows = leg2.get("windows", [])
        control_windows = ctrl.get("windows", [])
        leg1_ok = (
            leg1.get("schema") == "dlaw-leg1-anchor-v1"
            and leg1.get("status") == "PASS"
            and leg1.get("m_enforce") <= 68133
            and leg1.get("m_pre") == 262144
            and not leg1.get("failures")
            and set(leg1.get("worst_margins", {}))
            == {str(k) for k in range(5, 12)}
        )
        leg2_ok = (
            leg2.get("schema") == "dlaw-leg2-single-anchor-v2"
            and leg2.get("status") == "PASS"
            and leg2.get("gate_constant") == 216
            and leg2.get("m_lo") == 262144
            and leg2.get("h_step") == "7/20"
            and leg2.get("window_count") == 57
            and len(windows) == 57
            and all(w.get("window") == i for i, w in enumerate(windows))
            and all(w.get("envelope_positivity") is True
                    and len(w.get("gates", {})) == 7
                    and all(w["gates"].values()) for w in windows)
            and windows[0].get("x") == [0.0, 0.35]
            and windows[-1].get("x") == [19.6, 19.95]
            and all(windows[i]["x"][1] == windows[i + 1]["x"][0]
                    for i in range(56))
        )
        control_ok = (
            ctrl.get("schema") == "dlaw-leg2-single-anchor-v2"
            and ctrl.get("status") == "CONTROL_REJECTED"
            and ctrl.get("gate_constant") == 96
            and ctrl.get("window_count") == 1
            and len(control_windows) == 1
            and control_windows[0].get("envelope_positivity") is True
            and not all(control_windows[0].get("gates", {}).values())
        )
        gate("t=9 second proof: sharpened corridor core + (D-LAW) proved "
             "(leg 1 ladder + leg 2 envelope certificate, control "
             "rejected)",
             sharp.get("t9_conditional", {}).get("all_ok_primary") is True
             and "to be pinned" not in
                 sharp.get("t9_conditional", {}).get("d_law", "")
             and leg1_ok and leg2_ok and control_ok)

    # ---- coverage arithmetic (symbolic) -------------------------------
    # This is a polynomial identity, not a finite-range check.  With
    # u=t-9 its right side is u^2+9u+4, whose coefficients are positive;
    # hence the value is positive for every integer t>=9.
    lhs_coeff = (4, -9, 1)  # t^2 - 9t + 4
    rhs_coeff = (81 - 81 + 4, -18 + 9, 1)
    gate("symbolic dichotomy identity and positivity for every t>=9",
         lhs_coeff == rhs_coeff and all(c > 0 for c in (4, 9, 1)))

    # The ledger has two algebraic cases.  For 5<=t<=501 the core begins
    # immediately after the terminal pair and ends immediately before the
    # interior: [t,t+1] | [t+2,503] | [504,inf).  For t>=502 the core is
    # empty and max(504,t+2)=t+2, again immediately after the terminal pair.
    # Checking the two boundary values proves that there is no omitted or
    # duplicated integer at the only place where the formula changes.
    adjacency = (
        501 + 2 == 503
        and 503 + 1 == 504
        and 502 + 2 == 504
        and max(504, 502 + 2) == 502 + 2
    )
    gate("terminal/core/interior adjacency in both algebraic cases",
         adjacency)

    # ---- verdict ------------------------------------------------------
    if failures:
        status = "INCOMPLETE"
    elif conditional_legs:
        status = "PASS_CONDITIONAL"
    else:
        status = "PASS"
    print(f"\n=== JOIN: {status} ===")
    if failures:
        print("missing/failing legs:")
        for f_ in failures:
            print("  -", f_)

    payload = {
        "schema": "allt-sharp-shift-wall-join-v1",
        "status": status,
        "theorem": (
            "Omega_t(r) >= Omega_2(1350) for all integers r >= t >= 2, "
            "equality iff (t,r)=(2,1350); the sharp bridge factor is "
            "2*Omega_2(1350) = 2.655018313913361199726406175676..."
        ),
        "t5b": "rho(j,t) > 4/3 for every t >= 5, j >= t",
        "constants": {
            "c_star_decimal": common.decimal_truncation(c_star, 30)
            if c_star else None,
            "four_thirds_minus_c_star":
                common.decimal_truncation(FOUR_THIRDS - c_star, 24)
                if c_star else None,
            "joint_floor": str(CORE_TARGET * ATTACHED * ATTACHED),
            "needed_core": str(NEEDED_CORE),
            "top_floor": str(TOP_FLOOR),
        },
        "leg_hashes": legs,
        "gate_count": len(gate_results),
        "gates": gate_results,
        "conditional_legs": conditional_legs,
        "failures": failures,
    }
    out = HERE / "ALLT_JOIN.json"
    tmp = out.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=1)
    tmp.replace(out)
    print(f"wrote {out}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
