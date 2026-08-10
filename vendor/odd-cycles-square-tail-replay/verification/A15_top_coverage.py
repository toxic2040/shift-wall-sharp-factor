#!/usr/bin/env python3
"""Exact top-end corollary and integer gap-free coverage audit for (SQ).

The analytic inputs are A10-TOP and A14's uniform joint-window theorem.  This
program derives the scalar top threshold, certifies the two a_r splice anchors
with A3's dyadic interval ladder, validates the frozen finite-certificate
manifests, and reads the coverage JSONL sequentially.  It does not rerun the
104-second upstream bulk computation.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
import time
import subprocess
from pathlib import Path
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE.parent / "sech_instrument_2026-08-04"
FAILS: list[str] = []


def gate(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}{'  ' + detail if detail else ''}")
    if not condition:
        FAILS.append(name)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest(root: Path, manifest_name: str) -> tuple[int, list[str]]:
    bad = []
    count = 0
    for literal in (root / manifest_name).read_text().splitlines():
        if not literal.strip():
            continue
        expected, relative = literal.split(maxsplit=1)
        relative = relative.lstrip("*")
        target = root / relative
        if not target.is_file() or sha256(target) != expected:
            bad.append(relative)
        count += 1
    return count, bad


def top_certificate() -> None:
    # A10-TOP gives
    #   r psi >= r eps + (1-eps)(2q-q^2/r),  q=(t-2)L1_lo.
    # For t>=sqrt(a)+2, q>=1.  Convexity of L and L1_lo<L1<=r/(r-1)
    # give q<r on the (SQ) domain.  The concave quadratic takes its minimum
    # over [1,r] at an endpoint.
    q, r = sp.symbols("q r", positive=True)
    top_quad = 2 * q - q**2 / r
    gate("top quadratic endpoint q=1", sp.factor(top_quad.subs(q, 1) - (2 - 1 / r)) == 0)
    gate("top quadratic endpoint q=r", sp.factor(top_quad.subs(q, r) - r) == 0)

    # A7/A8 give 0<eps<2/(4r-3).  The retained factor is increasing in r.
    r0 = 4
    factor = (1 - F(2, 4 * r0 - 3)) * (2 - F(1, r0))
    surplus = factor * factor - 2
    gate("top scalar clears sqrt(2) from r=4", surplus > 0,
         f"factor={factor}, squared surplus={surplus}")
    r3 = 3
    factor3 = (1 - F(2, 4 * r3 - 3)) * (2 - F(1, r3))
    gate("top scalar threshold control rejects r=3", factor3 * factor3 < 2)


def anchor_certificate() -> None:
    # Reuse the proved outward-rounded arithmetic and cancellation-free ladder
    # from A3.  C=0 means E_1 is unscaled.  The two anchor verdicts below are
    # exact integer comparisons against the common denominator 2^192.
    sys.path.insert(0, str(HERE))
    import A3_certify as a3

    a3.iv_setup(192)
    ladder = a3.Ladder(1, C=0)
    anchors = []
    for _ in range(2_000_000):
        m, *_ = ladder.step()
        if m in (1_999_999, 2_000_000):
            lo, hi = ladder.E[1]
            anchors.append((m + 1, lo, hi))

    one = 1 << a3.P
    by_r = {r: (lo, hi) for r, lo, hi in anchors}
    lo_left, hi_left = by_r[2_000_000]
    lo_right, hi_right = by_r[2_000_001]
    gate("left splice anchor a_2000000<121", 6 * hi_left < 121 * one)
    gate("right splice anchor a_2000001>49", 6 * lo_right > 49 * one)
    payload = "\n".join(":".join(map(str, row)) for row in anchors)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    gate("splice-anchor interval digest",
         digest == "f28e67c4395fc0e926b5734162a0ba099e15da7b705fc6091deafe8c709b6356",
         digest)
    # Width is reported only as an exact numerator over 2^192.
    gate("anchor enclosures nonempty and narrow",
         all(lo <= hi and 6 * (hi - lo) < one // 10**20 for _, lo, hi in anchors),
         f"width numerators={[6 * (hi - lo) for _, lo, hi in anchors]}/2^192")


JOBS = [
    # name, rlo, rhi, tcap; None means the entire row through r-2.
    ("fullt_a", 601, 1000, None),
    ("fullt_b", 1001, 1250, None),
    ("fullt_c", 1251, 1450, None),
    ("iv_t55", 4, 8000, 55),
    ("iv_t30", 4, 20000, 30),
    ("iv_t20", 4, 60000, 20),
    ("iv_t12", 4, 2_000_000, 12),
]


def stream_job(name: str, rlo: int, rhi: int, tcap: int | None) -> tuple[bool, int, int]:
    expected_r = rlo
    rows = cells = 0
    path = HERE / f"A3_cert_{name}.jsonl"
    with path.open() as source:
        for literal in source:
            record = json.loads(literal)
            r = record.get("r")
            want_hi = r - 2 if tcap is None else min(r - 2, tcap)
            if (r != expected_r or record.get("t_hi") != want_hi
                    or record.get("fail", 0) != 0
                    or record.get("indet", 0) != 0
                    or "error" in record
                    or record.get("pass", 0) != want_hi - 1):
                return False, rows, cells
            expected_r += 1
            rows += 1
            cells += want_hi - 1
    return expected_r == rhi + 1, rows, cells


def finite_packet_certificate() -> None:
    a5_count, a5_bad = validate_manifest(HERE, "A5_SHA256SUMS.txt")
    upstream_count, upstream_bad = validate_manifest(UPSTREAM, "SHA256SUMS")
    gate("A5 frozen manifest", a5_count == 21 and not a5_bad,
         f"entries={a5_count}, bad={a5_bad}")
    gate("r<=600 frozen manifest", upstream_count == 17 and not upstream_bad,
         f"entries={upstream_count}, bad={upstream_bad}")

    # R2: re-derive the bulk verdict from rho_bulk_rows.jsonl rather than
    # substring-matching RHO_GATE_RUN_2026-08-04.txt.  The row range, the cell
    # count and the rho margin are all computed here; the frozen log is no
    # longer load-bearing (its integrity is still covered by the UPSTREAM
    # manifest gate above).
    bulk_rows = [json.loads(line)
                 for line in (UPSTREAM / "rho_bulk_rows.jsonl").read_text().splitlines()]
    rs = [record["r"] for record in bulk_rows]
    contiguous = rs == list(range(4, 601))
    # every cell of row r is 2 <= t <= r-2, so the row contributes r-3 cells
    derived_cells = sum(r - 3 for r in rs)
    worst = min((record["rho_min_float"], record["r"]) for record in bulk_rows)
    bulk_ok = contiguous and derived_cells == 178_503 and worst[0] > 1.0
    gate("full rows r=4..600", bulk_ok,
         f"cells=sum(r-3)={derived_cells} derived; worst rho_min={worst[0]:.6f} "
         f"at r={worst[1]}")
    # Honest scope: the JSONL carries per-row float minima.  The claim that the
    # bulk cells clear 1 EXACTLY rests on the upstream exact gate, whose file
    # integrity is pinned by the UPSTREAM manifest, not on the floats here.
    gate("bulk exactness is upstream-anchored, not re-derived here", True,
         "float rows are a regression against the upstream exact gate")

    expected_cells = {
        "fullt_a": 319_000,
        "fullt_b": 280_625,
        "fullt_c": 269_500,
        "iv_t55": 430_407,
        "iv_t30": 579_507,
        "iv_t20": 1_139_772,
        "iv_t12": 21_999_912,
    }
    for name, rlo, rhi, tcap in JOBS:
        ok, rows, cells = stream_job(name, rlo, rhi, tcap)
        gate(f"finite JSONL {name}", ok and cells == expected_cells[name],
             f"rows={rows}, cells={cells}")

    # R2: by default RE-EXECUTE A1_verify_tail.py and read the four column
    # closures out of that run.  --trust-a1-log falls back to the frozen log
    # for fast iteration, and says so in the gate name so a reader can never
    # mistake an asserted gate for a re-derived one.
    wanted = ("ALL GATES PASS",
              "the t=2 column of (SQ) is CLOSED on r >= 4",
              "the t=3 column of (SQ) is CLOSED on r >= 5",
              "the t=4 column of (SQ) is CLOSED on r >= 6",
              "the t=5 column of (SQ) is CLOSED on r >= 7")
    if "--trust-a1-log" in sys.argv:
        a1_text = (HERE / "A5_A1_replay_2026-08-05.log").read_text()
        label = "fixed columns t=2..5 [TRUSTED FROZEN LOG, NOT RE-DERIVED]"
        detail = "rerun without --trust-a1-log to re-execute A1_verify_tail.py"
        rc = 0
    else:
        before = time.time()
        run = subprocess.run([sys.executable, "A1_verify_tail.py"], cwd=HERE,
                             text=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, check=False)
        a1_text, rc = run.stdout, run.returncode
        label = "fixed columns t=2..5 [RE-EXECUTED A1_verify_tail.py]"
        detail = f"exit {rc}, {time.time() - before:.1f}s"
    gate(label, rc == 0 and all(w in a1_text for w in wanted), detail)


def integer_partition_certificate() -> None:
    # The row intervals are adjacent, with no missing integer r.
    ranges = [(4, 600), (601, 1450), (1451, 8000), (8001, 20000),
              (20001, 60000), (60001, 2_000_000)]
    gate("finite row ranges are contiguous",
         ranges[0][0] == 4
         and all(ranges[j][1] + 1 == ranges[j + 1][0] for j in range(len(ranges) - 1))
         and ranges[-1][1] + 1 == 2_000_001)

    # a_r is strictly increasing by the positive E,D ladder.  Hence the left
    # anchor gives a_r<121 for every r<=2e6.  Every uncovered finite-packet
    # cell there has t>=13, so t>sqrt(a_r)+2 and A15-TOP applies.
    caps = [(1451, 8000, 55), (8001, 20000, 30),
            (20001, 60000, 20), (60001, 2_000_000, 12)]
    gate("every capped residual starts at t>=13",
         min(cap + 1 for _, _, cap in caps) == 13)
    gate("left anchor makes the top splice strict", 11 + 2 == 13)

    # On r>=2e6+1 the right anchor gives a_r>49.  For the complementary
    # joint window t<sqrt(a_r)+2, writing y=sqrt(a_r)>7 gives
    # a_r-5sqrt(a_r)-10 = y^2-5y-10 > 4, hence a_r>5t.
    y = sp.symbols("y", positive=True)
    gate("right anchor implies the A14 corridor condition",
         sp.expand((y**2 - 5 * y - 10).subs(y, 7)) == 4
         and sp.diff(y**2 - 5 * y - 10, y).subs(y, 7) > 0)

    # The partition is literal: t=2..5 is A1; every integer t>=6 is either
    # at/above the top threshold or strictly below it.  A14 applies in the
    # latter case because r>=2e6+1>504 and a_r>5t.
    gate("right half-line starts inside A14 range", 2_000_001 >= 504)


def main() -> int:
    print("A15 -- top-end theorem and gap-free integer coverage")
    top_certificate()
    anchor_certificate()
    finite_packet_certificate()
    integer_partition_certificate()
    print()
    if FAILS:
        print(f"*** {len(FAILS)} A15 GATE FAILURES ***")
        for failure in FAILS:
            print("  ", failure)
        return 1
    print("ALL A15 GATES PASS")
    print("THEOREM: (SQ) holds for every integer r>=4 and 2<=t<=r-2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
