#!/usr/bin/env python3
"""Exact definition-level replay of the square-tail to shift-wall bridge.

The proof is symbolic in the manuscript.  This driver is a transcription
check: it constructs the Wall rows, all normalizers, the QDG inequality, the
fixed-t induction, and the final reversed-coefficient inequality using only
integers and Fraction arithmetic.  It also replays the finite lower-stripe
diagnostic quoted in the manuscript.  No floating-point value enters a verdict.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from math import factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def shift(a):
    return [F(0), *a]


def turan(a, k):
    def coeff(j):
        return a[j] if 0 <= j < len(a) else F(0)

    return coeff(k) ** 2 - coeff(k - 1) * coeff(k + 1)


def wall_rows(nmax):
    rows = [[F(1)], [F(1), F(1)]]
    for n in range(2, nmax + 1):
        rows.append(add(rows[-1], [F(0), *[F(n * n) * c for c in rows[-2]]]))
    return rows


def normalized_reverse(row):
    lead = row[-1]
    return [c / lead for c in reversed(row)]


def pi_row(rows, r):
    return [c / factorial(2 * r - 1) for c in reversed(rows[2 * r - 2])]


def p_row(rows, r):
    return normalized_reverse(rows[2 * r - 1])


def q_row(rows, r):
    return normalized_reverse(rows[2 * r])


def x_value(rows, r):
    return rows[2 * r - 2][-1] / rows[2 * r - 1][-1]


def alpha(r):
    numerator = 1
    denominator = 1
    for k in range(1, r):
        numerator *= 2 * k
        denominator *= 2 * k + 1
    return F(numerator, denominator)


def sqrt_ge_sum(a, b, c):
    """Decide sqrt(a) >= sqrt(b)+sqrt(c), exactly, for nonnegative rationals."""
    require(min(a, b, c) >= 0, "negative radicand in >= comparison")
    if a < b + c:
        return False
    return (a - b - c) ** 2 >= 4 * b * c


def sqrt_le_sum(a, b, c):
    """Decide sqrt(a) <= sqrt(b)+sqrt(c), exactly, for nonnegative rationals."""
    require(min(a, b, c) >= 0, "negative radicand in <= comparison")
    if a <= b + c:
        return True
    return (a - b - c) ** 2 <= 4 * b * c


def one_plus_u(a):
    return add(a, shift(a))


def g_row(rows, r):
    return one_plus_u(rows[2 * r - 2])


def gate_lower_stripe(rows):
    values = []
    for r in range(4, 70):
        q = r - 1
        numerator = F((2 * r) ** 2) * turan(g_row(rows, r), q)
        denominator = 2 * turan(rows[2 * r - 1], q)
        require(denominator > 0, f"lower-stripe denominator failed at r={r}")
        values.append((numerator / denominator, r))

    failures = [r for ratio, r in values if ratio < 1]
    minimum_ratio, minimum_row = min(values)
    require(failures == list(range(6, 42)), "lower-stripe failure window changed")
    require(minimum_row == 13, "lower-stripe minimum row changed")
    require(
        all(ratio >= 1 for ratio, r in values if r <= 5 or r >= 42),
        "lower-stripe complement contains a failure",
    )
    print(
        "[PASS] lower-stripe diagnostic: 66 cells, failures r=6..41, "
        f"exact minimum at r={minimum_row} "
        f"({minimum_ratio.numerator}/{minimum_ratio.denominator})"
    )


def gate_bridge(rmax):
    rows = wall_rows(max(2 * rmax, 138))
    recurrence_cells = 0
    normalization_cells = 0
    qdg_cells = 0
    triangle_cells = 0
    induction_cells = 0
    wall_cells = 0

    for r in range(2, rmax + 1):
        p = p_row(rows, r)
        p_prev = p_row(rows, r - 1)
        q_prev = q_row(rows, r - 1)
        x = x_value(rows, r)
        require(p == add(p_prev, [F(0), *[x * c for c in q_prev]]),
                f"P/Q recurrence failed at r={r}")
        recurrence_cells += 1

        pi = pi_row(rows, r)
        c0 = pi[0]
        require(q_prev == [value / c0 for value in pi],
                f"Q/Pi normalization failed at r={r}")
        require(x / c0 == alpha(r), f"x/c0 identity failed at r={r}")

        for t in range(2, r + 1):
            s = turan(one_plus_u(pi), t)
            v = turan(pi, t - 1)
            gamma = F(2 * r - 1, 2 * r)
            pi_prev = pi_row(rows, r - 1)
            s_prev = turan(one_plus_u(pi_prev), t)
            d = s - gamma ** 2 * s_prev

            math_s = turan(one_plus_u(q_prev), t)
            math_v = turan(q_prev, t - 1)
            a2 = F(r * r) * x * x * math_s
            b2 = F((r - 1) ** 2) * x_value(rows, r - 1) ** 2 * turan(
                one_plus_u(q_row(rows, r - 2)), t
            )
            capital_b = x * x * math_v
            delta = a2 - b2

            require(a2 == F(r * r) * alpha(r) ** 2 * s,
                    f"A normalization failed at (r,t)=({r},{t})")
            require(capital_b == alpha(r) ** 2 * v,
                    f"B normalization failed at (r,t)=({r},{t})")
            require(delta == F(r * r) * alpha(r) ** 2 * d,
                    f"Delta normalization failed at (r,t)=({r},{t})")
            require(d > 0, f"decrement is not positive at (r,t)=({r},{t})")
            require(delta * delta >= 2 * a2 * capital_b,
                    f"square-tail normalization failed at (r,t)=({r},{t})")
            require(sqrt_ge_sum(a2, b2, capital_b / 2),
                    f"QDG failed at (r,t)=({r},{t})")
            normalization_cells += 1
            qdg_cells += 1

            lp = turan(p, t)
            lp_prev = turan(p_prev, t)
            require(sqrt_le_sum(lp, lp_prev, x * x * math_v),
                    f"Turan triangle failed at (r,t)=({r},{t})")
            require(lp <= 2 * a2,
                    f"fixed-t induction bound failed at (r,t)=({r},{t})")
            triangle_cells += 1
            induction_cells += 1

            q = r - t
            c_row = rows[2 * r - 1]
            g = g_row(rows, r)
            o0 = c_row[-1]
            require(lp == turan(c_row, q) / (o0 * o0),
                    f"C reversal failed at (r,t)=({r},{t})")
            require(x * x * math_s == turan(g, q) / (o0 * o0),
                    f"G reversal failed at (r,t)=({r},{t})")
            require(F((2 * r) ** 2) * turan(g, q) >= 2 * turan(c_row, q),
                    f"shift-wall inequality failed at (r,q)=({r},{q})")
            wall_cells += 1

    # Negative controls target the two most error-prone transcription points.
    r = 4
    bad_x = x_value(rows, r) + 1
    require(
        p_row(rows, r)
        != add(p_row(rows, r - 1), [F(0), *[bad_x * c for c in q_row(rows, r - 1)]]),
        "mutated P/Q coefficient was not detected",
    )
    q = 1
    bad_g = list(g_row(rows, r))
    bad_g[q] += 1
    require(
        x_value(rows, r) ** 2
        * turan(one_plus_u(q_row(rows, r - 1)), r - q)
        != turan(bad_g, q) / rows[2 * r - 1][-1] ** 2,
        "mutated reversal identity was not detected",
    )

    print(f"[PASS] P/Q recurrence: {recurrence_cells} rows")
    print(f"[PASS] normalizers and square-tail form: {normalization_cells} cells")
    print(f"[PASS] QDG radical comparison: {qdg_cells} cells")
    print(f"[PASS] proper-position triangle replay: {triangle_cells} cells")
    print(f"[PASS] fixed-t induction bound: {induction_cells} cells")
    print(f"[PASS] reversed shift-wall inequality: {wall_cells} cells")
    print("[PASS] negative controls: 2/2 rejected")
    gate_lower_stripe(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rmax", type=int, default=24)
    args = parser.parse_args()
    require(args.rmax >= 4, "--rmax must be at least 4")
    gate_bridge(args.rmax)


if __name__ == "__main__":
    main()
