#!/usr/bin/env python3
"""Exact characterization of the t=1 stripe on 4 <= r <= 3000.

The lower-stripe diagnostic in A36_bw_bridge.py evaluates the reversed
shift-wall ratio at q = r-1 on the finite window 4 <= r <= 69 from full Wall
rows.  At reverse index 1 the two Turan forms consume only the top four
coefficients of each Wall row, so the same exact ratio can be streamed far
past the full-row range.  This driver carries those top coefficients as
integers through the Wall recurrence M_n = M_{n-1} + n^2 u M_{n-2},
cross-checks the streamed values against the full-row computation on the
shared window, and then decides, in integer arithmetic only:

  * the failure set of omega_1(r) < 1 on 4 <= r <= RMAX is exactly r = 6..41;
  * the minimum is at r = 13 and equals the banked exact rational of the
    66-cell diagnostic digit for digit;
  * omega_1 is strictly increasing on 42 <= r <= RMAX.

All claims are finite-range: nothing here asserts the stripe beyond RMAX.
No floating-point value enters a verdict.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F


BANKED_MINIMUM = F(
    141759451874843777120613875493682110650895798,
    150947693570603679082853646116464639909331875,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def turan(a, k):
    def coeff(j):
        return a[j] if 0 <= j < len(a) else 0

    return coeff(k) ** 2 - coeff(k - 1) * coeff(k + 1)


def wall_rows(nmax):
    rows = [[1], [1, 1]]
    for n in range(2, nmax + 1):
        rows.append(add(rows[-1], [0, *[n * n * c for c in rows[-2]]]))
    return rows


def one_plus_u(a):
    return add(a, [0, *a])


WINDOW = 4


class TopWindow:
    """Top coefficients of one Wall row, indexed by absolute degree."""

    __slots__ = ("degree", "coefficients")

    def __init__(self, degree, coefficients):
        self.degree = degree
        self.coefficients = coefficients  # coefficients[i] is degree-i

    def coeff(self, j):
        if j < 0 or j > self.degree:
            return 0
        offset = self.degree - j
        require(
            offset < len(self.coefficients),
            f"streamed window read below tracked depth at degree {j}",
        )
        return self.coefficients[offset]


def step(previous: TopWindow, before: TopWindow, n: int, bump: int = 0) -> TopWindow:
    degree = max(previous.degree, before.degree + 1)
    weight = n * n + bump
    coefficients = [
        previous.coeff(j) + weight * before.coeff(j - 1)
        for j in range(degree, max(degree - WINDOW, -1), -1)
    ]
    return TopWindow(degree, coefficients)


def stripe_pair(m_even: TopWindow, m_odd: TopWindow, r: int):
    """Return (numerator, denominator) of omega_1(r) as positive integers."""
    g_top = m_even.coeff(r - 1)
    g_mid = m_even.coeff(r - 1) + m_even.coeff(r - 2)
    g_low = m_even.coeff(r - 2) + m_even.coeff(r - 3)
    turan_g = g_mid * g_mid - g_low * g_top
    turan_c = m_odd.coeff(r - 1) ** 2 - m_odd.coeff(r - 2) * m_odd.coeff(r)
    require(turan_c > 0, f"stripe denominator failed at r={r}")
    return (2 * r) ** 2 * turan_g, 2 * turan_c


def streamed_pairs(rmax, mutate_at=0):
    pairs = {}
    before = TopWindow(0, [1])
    previous = TopWindow(1, [1, 1])
    for n in range(2, 2 * rmax):
        bump = 1 if n == mutate_at else 0
        current = step(previous, before, n, bump)
        before, previous = previous, current
        if n % 2 == 1 and n >= 7:
            r = (n + 1) // 2
            pairs[r] = stripe_pair(before, previous, r)
    return pairs


def full_row_pairs(rmax):
    rows = wall_rows(2 * rmax - 1)
    pairs = {}
    for r in range(4, rmax + 1):
        q = r - 1
        numerator = (2 * r) ** 2 * turan(one_plus_u(rows[2 * r - 2]), q)
        denominator = 2 * turan(rows[2 * r - 1], q)
        require(denominator > 0, f"full-row denominator failed at r={r}")
        pairs[r] = (numerator, denominator)
    return pairs


def failure_rows(pairs):
    return [r for r in sorted(pairs) if pairs[r][0] < pairs[r][1]]


def gate_cross_check(pairs):
    reference = full_row_pairs(69)
    for r in sorted(reference):
        require(
            F(*pairs[r]) == F(*reference[r]),
            f"streamed ratio disagrees with full rows at r={r}",
        )
    print(f"[PASS] streamed/full-row cross-check: {len(reference)} cells agree exactly")


def gate_failure_window(pairs, rmax):
    failures = failure_rows(pairs)
    require(failures == list(range(6, 42)), "t=1 failure window changed")
    print(
        f"[PASS] t=1 failure set on 4<=r<={rmax}: exactly r=6..41 "
        f"({len(failures)} rows)"
    )


def gate_minimum(pairs):
    minimum_row = None
    for r in sorted(pairs):
        numerator, denominator = pairs[r]
        if minimum_row is None or (
            numerator * pairs[minimum_row][1]
            < pairs[minimum_row][0] * denominator
        ):
            minimum_row = r
    require(minimum_row == 13, "t=1 minimum row changed")
    minimum = F(*pairs[minimum_row])
    require(minimum == BANKED_MINIMUM, "t=1 minimum differs from the banked value")
    print(
        f"[PASS] exact minimum at r=13: "
        f"{minimum.numerator}/{minimum.denominator}"
    )


def gate_monotone_recovery(pairs, rmax):
    n42, d42 = pairs[42]
    require(n42 >= d42, "omega_1(42) is below one")
    for r in range(43, rmax + 1):
        n_prev, d_prev = pairs[r - 1]
        n_cur, d_cur = pairs[r]
        require(n_cur * d_prev > n_prev * d_cur,
                f"omega_1 is not strictly increasing at r={r}")
    print(f"[PASS] monotone recovery: omega_1 strictly increasing on 42<=r<={rmax}")


def gate_negative_controls():
    mutated = streamed_pairs(69, mutate_at=40)
    reference = full_row_pairs(69)
    detected = any(F(*mutated[r]) != F(*reference[r]) for r in sorted(reference))
    require(detected, "mutated streaming recurrence was not detected")

    corrupted = dict(streamed_pairs(69))
    n50, d50 = corrupted[50]
    corrupted[50] = (n50 - d50, d50)
    require(
        failure_rows(corrupted) != list(range(6, 42)),
        "corrupted failure row was not detected",
    )
    print("[PASS] negative controls: 2/2 rejected")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rmax", type=int, default=3000)
    args = parser.parse_args()
    require(args.rmax >= 70, "--rmax must be at least 70")

    pairs = streamed_pairs(args.rmax)
    gate_cross_check(pairs)
    gate_failure_window(pairs, args.rmax)
    gate_minimum(pairs)
    gate_monotone_recovery(pairs, args.rmax)
    gate_negative_controls()


if __name__ == "__main__":
    main()
