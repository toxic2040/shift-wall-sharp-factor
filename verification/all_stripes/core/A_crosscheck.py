#!/usr/bin/env python3
"""Spot crosscheck: sweep cells vs the independent top-window stream.

Compares the definition-level values of A_core_sweep.py (built from the
full wall row) against rho_stream() of the t=4 transfer producer
verification/lower_stripes/verify_sq_t4_transfer.py (a t-parameterized
top-window stream, and so a distinct implementation path).  That module
is imported read-only; nothing here writes to it.
"""
import importlib.util
import sys
from fractions import Fraction as Frac
from pathlib import Path

HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parents[2]
TRANSFER = (RELEASE_ROOT / "verification"
            / "lower_stripes"
            / "verify_sq_t4_transfer.py")

sys.set_int_max_str_digits(1_000_000)

spec = importlib.util.spec_from_file_location("sq_t4", TRANSFER)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

sweep_spec = importlib.util.spec_from_file_location(
    "core_sweep", HERE / "A_core_sweep.py")
sweep = importlib.util.module_from_spec(sweep_spec)
sweep_spec.loader.exec_module(sweep)


def sweep_cell(r, t):
    Fprev = Fcur = None
    for n, row in sweep.wall_stream(2 * r - 2):
        if n % 2:
            continue
        Fprev = Fcur
        Fcur = sweep.F_from_row(row, n // 2 + 1)
    Lg = sweep.Lop(sweep.shift_add(Fcur), t)
    Lv = sweep.Lop(Fcur, t - 1)
    Lg1 = sweep.Lop(sweep.shift_add(Fprev), t)
    M = r * r * Lg - (2 * r - 1) ** 4 * (r - 1) ** 2 * Lg1
    return Frac(M * M, 2 * r * r * Lg * Lv)


CELLS = [(7, 5), (20, 10), (60, 30), (100, 5), (200, 50), (503, 5),
         (503, 250), (503, 501)]

ok = True
for r, t in CELLS:
    stream_val = None
    for rr, num, den, _dec in mod.rho_stream(r, t=t):
        if rr == r:
            stream_val = Frac(num, den)
    mine = sweep_cell(r, t)
    match = stream_val == mine
    ok &= match
    print("(%3d,%3d)  %s  rho=%.9f" % (r, t, "MATCH" if match else
                                       "MISMATCH", float(mine)))
print("crosscheck:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
