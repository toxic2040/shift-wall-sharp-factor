# Certificate source map

All paths below are relative to the release root. Producer scripts and their
banked records sit together unless noted otherwise.

`verification/lower_stripes/` holds the `t=2,3,4` stripes, where the sharp
minimum lives. `verification/all_stripes/` holds the `t>=5` bulk and the join
that assembles every stripe into the theorem: `columns/` for the fixed columns
`t=5,...,10`, `core/` for the finite core, `join/` for the joint region, the
decrement law, and the master assembly.

| Claim | Producer | Banked record |
| --- | --- | --- |
| Unique minimum `Omega_2(1350)` | `verification/lower_stripes/verify_bridge_t2_minimum.py` | `BRIDGE_T2_MINIMUM.json` |
| `Omega_3(r)>4/3` | `verification/lower_stripes/verify_bridge_t3_barrier.py` | `BRIDGE_T3_BARRIER.json` |
| `rho(r,4)>4/3` and lossless transfer | `verification/lower_stripes/verify_sq_t4_transfer.py` | `SQ_T4_BRIDGE_TRANSFER.json` |
| Fixed-point endpoint audit | `verification/lower_stripes/audit_a1_pi_endpoint.py` | `A1_PI_ENDPOINT_AUDIT.json` |
| 123,753-cell finite core | `verification/all_stripes/core/A_core_sweep.py` | `A_core_sweep.jsonl`, `A_CORE_SWEEP_SUMMARY.json` |
| Independent finite-core spot check | `verification/all_stripes/core/A_crosscheck.py` | transcript gate in full replay |
| Exact prefixes `t=5,...,8` | `verification/all_stripes/columns/h2_prefix.py` | `H2_PREFIX_T5.json`, ..., `H2_PREFIX_T8.json` |
| Positive-coefficient tails `t=4,...,10` | `verification/all_stripes/columns/h2_tail_cert.py` | eleven `H2_TAIL_*.json` records |
| Level crossings and top anchors | `verification/all_stripes/columns/h2_a_gate.py` | `H2_A_GATE.json` |
| Complete columns `t=5,...,10` | `h2_join.py`, `h2_join_prefix_free.py` in the same directory | `H2_COL_T5.json`, ..., `H2_COL_T10.json` |
| Joint region `t=11,...,16` | `verification/all_stripes/join/verify_joint_core_t11.py` | `verify_joint_core_t11.json` |
| Uniform joint region `t>=17` | `verification/all_stripes/join/verify_joint_core_bidegree.py` | `verify_joint_core_bidegree.json` |
| D-LAW finite ladder | `verification/all_stripes/join/dlaw_leg1_prefix.py` | `dlaw_leg1_prefix_anchor.json` |
| D-LAW half-line envelope and control | `verification/all_stripes/join/dlaw_leg2_bootstrap.py` | `DLAW_LEG2.json`, `DLAW_LEG2_CONTROL.json` |
| Sharpened corridor second chain | `verification/all_stripes/join/verify_joint_core_sharp.py` | `verify_joint_core_sharp.json` |
| Complete theorem assembly | `verification/all_stripes/join/verify_allt_join.py` | `ALLT_JOIN.json` |
| Manuscript and release literals | `verify_release_anchors.py` | 40 fail-closed transcript gates |

The `t=2` fixed-column tail behind the shared envelope algebra is written up
in `docs/FIXED_COLUMN_TAIL_T2.md`; both of its verifiers ship under
`vendor/odd-cycles-square-tail-replay/verification/`.

The full runner first replays `vendor/odd-cycles-square-tail-replay/verify_all.py`
to check the transitive Paper 3 dependency. Paper 2 is not a proof dependency.

`SHA256SUMS` binds every released source, record, document, and PDF. The exact
target rational is additionally bound by the canonical fraction digest

```text
2fc6689fd42dd6eaa29088aeb0077e1120878fd2f6b54cf2e8dbd50e823a5691.
```
