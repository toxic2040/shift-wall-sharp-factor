#!/usr/bin/env python3
"""Certified dyadic screen for the core uniqueness comparison.

A43 settles uniqueness by exact cross-multiplication on 125,250 cells.  That is
correct but carries roughly 17,000-bit integers to decide a sign that a few tens
of bits already decide.  Three cancellation maxima were measured separately over
the whole core: the Turan subtraction loses at most 7 bits (at (440,272)), the
adjacent-row deficit at most 9, and the comparison against rho(129,2) at most 19
(at the runner-up (128,2)).  These maxima occur at different cells and are NOT
additive, so they do not sum to a precision requirement.  The mantissa width at
which the screen stops falling back is an empirical threshold, measured at 32
bits by sweeping the width; it is not derived from those three numbers.

This file therefore replaces the exact comparison with a directed-rounding
dyadic interval evaluation at a fixed mantissa width, and falls back to the
exact integer comparison only on the cells the interval cannot separate.  The
verdict is exact either way: intervals are conservative by construction, and
every ambiguous cell is decided exactly.

Two further savings are independent of the interval arithmetic and account for
most of the speedup: the shifted row is never materialised (only three of its
entries are ever read), and the Wall rows are built once per block rather than
per cell.

Notation.  Two normalisations are in play and must not be mixed.  Write m_n for
the reverse-diagonal *integer* Wall row, and Pi_r for the normalised generating
function, Pi_r(1)=1, whose leading coefficient is 1/(2r-1)!.  In the normalised
variables the definitional statement is

    S   = L(shift(Pi_r))_t,   V   = L(Pi_r)_(t-1),   S' = L(shift(Pi_(r-1)))_t
    d   = S - ((2r-1)/(2r))^2 S'
    rho = r^2 d^2 / (2 S V)

and this file never evaluates that form.  It uses the cancelled integer form,
in which the factorials are cleared: with

    Sm  = L(shift(m_(2r-2)))_t,  Vm = L(m_(2r-2))_(t-1),
    Sm' = L(shift(m_(2r-4)))_t,
    M   = r^2 Sm - (2r-1)^4 (r-1)^2 Sm'          (the cancelled integer numerator)

one has the identity, gated in A38_targeted_core.identity_gate,

    rho(r,t) = M^2 / (2 r^2 Sm Vm).

M is the cancelled decrement, not the definitional d, and Sm, Vm are not S, V;
they agree only after the factorial normalisation is restored.  The comparison
rho(r,t) vs rho* = N*/D* is taken in the product form M^2 D* vs 2 r^2 Sm Vm N*,
so no division enters the screen.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFY = HERE
RESULTS = HERE.parent / "results"
DEFAULT_JSONL = HERE / "A44_DYADIC_SCREEN.work.jsonl"
DEFAULT_OUTPUT = RESULTS / "dyadic_screen.json"

R_LO, R_HI, BLOCK_SIZE = 4, 503, 25
EXPECTED_CELLS = 125_250
WITNESS = (129, 2)
ALGORITHM = "square-tail-dyadic-screen-v1"
A40_DECIMAL_SHA256 = "9c09b2864de0e742c291698d4c63ad631d2d872c122f62b6c0fdc8c10b47e373"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_a38():
    import importlib.util

    path = VERIFY / "A38_targeted_core.py"
    require(path.is_file(), f"missing producer {path}")
    spec = importlib.util.spec_from_file_location("_a38_screen", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A38 = load_a38()
reverse_rows = A38.reverse_rows
rho_parts = A38.rho_parts
integer_encoding = A38.integer_encoding
decimal_truncation = A38.decimal_truncation


# --------------------------------------------------------------- dyadic layer
# A positive quantity is enclosed as [lo, hi] * 2**e with lo, hi non-negative
# integers of at most BITS significant bits.  Lower endpoints round down and
# upper endpoints round up, so every enclosure is valid.

def _norm(lo: int, hi: int, e: int, bits: int):
    excess = hi.bit_length() - bits
    if excess > 0:
        lo >>= excess
        hi = (hi >> excess) + 1
        e += excess
    return lo, hi, e


def _from_int(value: int, bits: int):
    return _norm(value, value, 0, bits)


def _mul(a, b, bits: int):
    return _norm(a[0] * b[0], a[1] * b[1], a[2] + b[2], bits)


def _sub(a, b, bits: int):
    """Enclose a-b, or return None when the interval cannot certify positivity."""
    e = min(a[2], b[2])
    a_lo, a_hi = a[0] << (a[2] - e), a[1] << (a[2] - e)
    b_lo, b_hi = b[0] << (b[2] - e), b[1] << (b[2] - e)
    lo, hi = a_lo - b_hi, a_hi - b_lo
    if lo <= 0:
        return None
    return _norm(lo, hi, e, bits)


def _compare(a, b):
    """+1 if a>b certainly, -1 if a<b certainly, 0 if the intervals overlap."""
    e = min(a[2], b[2])
    a_lo, a_hi = a[0] << (a[2] - e), a[1] << (a[2] - e)
    b_lo, b_hi = b[0] << (b[2] - e), b[1] << (b[2] - e)
    if a_lo > b_hi:
        return 1
    if a_hi < b_lo:
        return -1
    return 0


# --------------------------------------------------------------- cell algebra

def _turan_at(row, k: int, shift: bool):
    """The three entries of L(row)_k, reading the shifted row entrywise."""
    n = len(row)

    def entry(j):
        if j < 0 or j >= n:
            return 0
        if not shift:
            return row[j]
        return row[j] + (row[j - 1] if j else 0)

    centre = entry(k)
    return centre * centre - entry(k - 1) * entry(k + 1)


_STATE: dict = {}


def initializer(star_num: int, star_den: int, bits: int) -> None:
    _STATE["star_num"] = star_num
    _STATE["star_den"] = star_den
    _STATE["bits"] = bits


def evaluate_block(block: tuple[int, int]) -> dict:
    start, end = block
    bits = _STATE["bits"]
    star_num, star_den = _STATE["star_num"], _STATE["star_den"]
    dyad_num = _from_int(star_num, bits)
    dyad_den = _from_int(star_den, bits)

    rows = reverse_rows(2 * end - 2, end)
    cells = screened = fallback = 0
    at_or_below: list[list[int]] = []
    equal_cells: list[list[int]] = []
    digest = hashlib.sha256()   # the canonical sign stream, gated against A43

    for r in range(start, end + 1):
        current, predecessor = rows[2 * r - 2], rows[2 * r - 4]
        r2 = r * r
        coefficient = (2 * r - 1) ** 4 * (r - 1) ** 2
        for t in range(2, r - 1):
            cells += 1
            s_cur = _turan_at(current, t, True)
            s_pre = _turan_at(predecessor, t, True)
            v_cur = _turan_at(current, t - 1, False)

            verdict = 0
            if s_cur > 0 and v_cur > 0 and s_pre > 0:
                left = _mul(_from_int(r2, bits), _from_int(s_cur, bits), bits)
                right = _mul(_from_int(coefficient, bits), _from_int(s_pre, bits), bits)
                d = _sub(left, right, bits)
                if d is not None:
                    lhs = _mul(_mul(d, d, bits), dyad_den, bits)
                    rhs = _mul(
                        _mul(_mul(_from_int(2 * r2, bits), _from_int(s_cur, bits), bits),
                             _from_int(v_cur, bits), bits),
                        dyad_num,
                        bits,
                    )
                    verdict = _compare(lhs, rhs)

            if verdict > 0:
                screened += 1
                digest.update(integer_encoding(r))
                digest.update(integer_encoding(t))
                digest.update(b">")
                continue

            # unresolved by the screen, or certainly not above: decide exactly
            fallback += 1
            numerator, denominator, _ = rho_parts(rows, r, t)
            require(denominator > 0, f"non-positive denominator at ({r},{t})")
            margin = numerator * star_den - denominator * star_num
            digest.update(integer_encoding(r))
            digest.update(integer_encoding(t))
            digest.update(b">" if margin > 0 else (b"=" if margin == 0 else b"<"))
            if margin < 0:
                at_or_below.append([r, t])
            elif margin == 0:
                at_or_below.append([r, t])
                equal_cells.append([r, t])

    return {
        "block": [start, end],
        "cells": cells,
        "screened_above": screened,
        "exact_fallbacks": fallback,
        "at_or_below_witness": at_or_below,
        "equal_to_witness": equal_cells,
        "sign_sha256": digest.hexdigest(),
    }


def agreement_gate(star_num: int, star_den: int, bits: int, stride: int) -> dict:
    """Screen and exact verdicts must agree on a deterministic sample."""
    dyad_num, dyad_den = _from_int(star_num, bits), _from_int(star_den, bits)
    index = checked = 0
    disagreements = []
    for lo in range(R_LO, R_HI + 1, 100):
        hi = min(lo + 99, R_HI)
        rows = reverse_rows(2 * hi - 2, hi)
        for r in range(lo, hi + 1):
            current, predecessor = rows[2 * r - 2], rows[2 * r - 4]
            for t in range(2, r - 1):
                index += 1
                if index % stride:
                    continue
                s_cur = _turan_at(current, t, True)
                s_pre = _turan_at(predecessor, t, True)
                v_cur = _turan_at(current, t - 1, False)
                d = _sub(
                    _mul(_from_int(r * r, bits), _from_int(s_cur, bits), bits),
                    _mul(_from_int((2 * r - 1) ** 4 * (r - 1) ** 2, bits),
                         _from_int(s_pre, bits), bits),
                    bits,
                )
                screen = 0
                if d is not None:
                    screen = _compare(
                        _mul(_mul(d, d, bits), dyad_den, bits),
                        _mul(_mul(_mul(_from_int(2 * r * r, bits),
                                       _from_int(s_cur, bits), bits),
                                  _from_int(v_cur, bits), bits), dyad_num, bits),
                    )
                numerator, denominator, _ = rho_parts(rows, r, t)
                margin = numerator * star_den - denominator * star_num
                exact = (margin > 0) - (margin < 0)
                checked += 1
                if screen and screen != exact:
                    disagreements.append([r, t, screen, exact])
    return {"sampled_cells": checked, "disagreements": disagreements}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bits", type=int, default=128)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument("--stride", type=int, default=977)
    parser.add_argument("--uniqueness", type=Path,
                        default=HERE.parent / "results" / "core_uniqueness.json")
    parser.add_argument("--require-fresh", action="store_true")
    arguments = parser.parse_args()

    rows = reverse_rows(2 * WITNESS[0] - 2, WITNESS[0])
    numerator, denominator, _ = rho_parts(rows, *WITNESS)
    star = F(numerator, denominator)
    star_num, star_den = star.numerator, star.denominator
    decimal_digest = hashlib.sha256(f"{star_num}/{star_den}".encode()).hexdigest()
    require(decimal_digest == A40_DECIMAL_SHA256, "witness fingerprint mismatch")

    agreement = agreement_gate(star_num, star_den, arguments.bits, arguments.stride)
    require(not agreement["disagreements"],
            f"screen disagreed with exact: {agreement['disagreements']}")

    blocks = [
        (start, min(start + BLOCK_SIZE - 1, R_HI))
        for start in range(R_LO, R_HI + 1, BLOCK_SIZE)
    ]
    done: dict[str, dict] = {}
    if arguments.jsonl.is_file() and not arguments.require_fresh:
        for line in arguments.jsonl.read_text().splitlines():
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if record.get("algorithm") == ALGORITHM and record.get("bits") == arguments.bits:
                done[json.dumps(record["block"])] = record
    pending = [b for b in blocks if json.dumps(list(b)) not in done]
    print(f"dyadic screen at {arguments.bits} bits: {len(pending)} blocks to run", flush=True)

    if pending:
        with arguments.jsonl.open("a", encoding="utf-8") as handle:
            with ProcessPoolExecutor(
                max_workers=max(1, min(arguments.workers, len(pending))),
                initializer=initializer,
                initargs=(star_num, star_den, arguments.bits),
            ) as pool:
                futures = {pool.submit(evaluate_block, b): b for b in pending}
                for future in as_completed(futures):
                    block = futures[future]
                    try:
                        record = future.result()
                    except Exception as error:
                        record = {"block": list(block), "error": repr(error)}
                    record["algorithm"] = ALGORITHM
                    record["bits"] = arguments.bits
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    handle.flush()
                    done[json.dumps(list(block))] = record

    records = [done[json.dumps(list(b))] for b in blocks]
    errors = [r for r in records if "error" in r]
    require(not errors, f"block errors: {errors}")

    cells = sum(r["cells"] for r in records)
    screened = sum(r["screened_above"] for r in records)
    fallback = sum(r["exact_fallbacks"] for r in records)
    at_or_below = [c for r in records for c in r["at_or_below_witness"]]
    equal = [c for r in records for c in r["equal_to_witness"]]

    aggregate = hashlib.sha256()
    for record in records:
        aggregate.update(record["sign_sha256"].encode())
    canonical_sign = aggregate.hexdigest()

    reference = json.loads(arguments.uniqueness.read_text())
    require(reference.get("schema") == "square-tail-core-uniqueness-v1",
            "reference certificate is not A43")
    sign_streams_agree = reference.get("canonical_sign_sha256") == canonical_sign
    sole_fallback_is_the_equality = (
        fallback == 1 and at_or_below == [list(WITNESS)] and equal == [list(WITNESS)]
    )

    gates = {
        "cell count matches core": cells == EXPECTED_CELLS,
        "screen agrees with exact on the sample": not agreement["disagreements"],
        "exactly one cell at or below witness": at_or_below == [list(WITNESS)],
        "the witness attains equality": equal == [list(WITNESS)],
        "screen and fallback partition the core": screened + fallback == cells,
        "A43 verdict reproduced": at_or_below == [list(WITNESS)] and equal == [list(WITNESS)],
        "sign stream matches A43 cell-by-cell": sign_streams_agree,
        "sole fallback is the exact equality at (129,2)": sole_fallback_is_the_equality,
    }
    payload = {
        "algorithm": ALGORITHM,
        "schema": "square-tail-dyadic-screen-v1",
        "claim": (
            "the A43 core uniqueness verdict, reproduced with a certified "
            "dyadic interval screen and exact fallback"
        ),
        "mantissa_bits": arguments.bits,
        "cells": cells,
        "screened_above_witness": screened,
        "exact_fallbacks": fallback,
        "screen_discharge_fraction": round(screened / cells, 6) if cells else None,
        "cells_at_or_below_witness": at_or_below,
        "cells_equal_to_witness": equal,
        "measured_cancellation_bits": {
            "turan_worst_over_core": 7,
            "adjacent_row_deficit_worst": 9,
            "comparison_margin_worst": 19,
            "attained_at": [128, 2],
        },
        "agreement_gate": {
            "stride": arguments.stride,
            "sampled_cells": agreement["sampled_cells"],
            "disagreements": agreement["disagreements"],
        },
        "canonical_sign_sha256": canonical_sign,
        "authority": (
            "A43 is the authoritative exact certificate; this file is an "
            "accelerator whose entire sign stream is gated against it"
        ),
        "configuration": {
            "worker_policy": "one outer ProcessPoolExecutor; no nested pool",
            "block_size": BLOCK_SIZE,
        },
        "gates": gates,
        "passed": all(gates.values()),
    }
    payload["status"] = "PASS" if payload["passed"] else "FAIL"
    arguments.output.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: payload[k] for k in
                      ("status", "cells", "screened_above_witness",
                       "exact_fallbacks", "screen_discharge_fraction")}))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
