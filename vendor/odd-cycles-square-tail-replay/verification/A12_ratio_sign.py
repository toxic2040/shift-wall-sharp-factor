#!/usr/bin/env python3
"""Exact gates for the scaled-ratio sign theorem.

The manuscript proof uses the total-nonnegative odd coefficient triangle. This
program derives both connection formulae, checks their literal normalizations,
proves the load-bearing connection-ratio sign symbolically, and replays the
resulting row inequalities with exact integers.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import math

import sympy as sp


def odd_harmonic(d: int) -> F:
    return sum((F(1, 2 * r + 1) for r in range(d)), F(0))


def derivative_connection(m: int, j: int) -> F:
    d = m - j
    return F(math.factorial(2 * m + 1), math.factorial(2 * j + 1)) * odd_harmonic(d) / d


def low_connection(m: int, j: int) -> F:
    d = m - j
    return F(
        math.factorial(2 * m + 1) * (2 * m + 2 * j + 3),
        2 * math.factorial(2 * j + 1) * (2 * d - 1) * (2 * d + 1),
    )


def A_rows(nmax: int, y: sp.Symbol):
    rows = [sp.Integer(1), y]
    for n in range(1, nmax):
        rows.append(sp.expand(y * rows[-1] + n * n * rows[-2]))
    return rows


def P_rows(mmax: int, u: sp.Symbol):
    rows = [sp.Integer(1), u + 5]
    for m in range(1, mmax):
        rows.append(sp.expand((u + 8 * m * m + 12 * m + 5) * rows[-1]
                              - (2 * m * (2 * m + 1)) ** 2 * rows[-2]))
    return rows


def gate_connection_identities() -> None:
    y, u = sp.symbols("y u")
    A = A_rows(19, y)
    P = P_rows(9, u)

    # d^2/dy^2 is multiplication by atanh(x)^2 in the bivariate EGF.
    # Its coefficient is derived from the convolution, not inserted.
    for d in range(1, 10):
        convolution = sum(F(1, (2 * r + 1) * (2 * (d - 1 - r) + 1))
                          for r in range(d))
        assert convolution == odd_harmonic(d) / d
    print("[PASS] [x^(2d)] atanh(x)^2 = H_d^odd/d derived for d=1..9")

    checked = 0
    for m in range(1, 9):
        derived = sum(sp.Rational(derivative_connection(m, j).numerator,
                                  derivative_connection(m, j).denominator) * A[2 * j + 1]
                      for j in range(m))
        assert sp.expand(sp.diff(A[2 * m + 1], y, 2) - derived) == 0
        checked += 1
    print(f"[PASS] derivative lowering A_(2m+1)''=sum d_(m,j)A_(2j+1), m=1..8 ({checked} rows)")

    theta = lambda f: sp.expand(u * sp.diff(f, u))
    for m in range(1, 9):
        lowered = sum(sp.Rational(low_connection(m, j).numerator,
                                  low_connection(m, j).denominator) * P[j]
                      for j in range(m))
        assert sp.expand((m * P[m] - theta(P[m])) - lowered) == 0
    print("[PASS] LOW literal (m-theta)P_m=sum ell_(m,j)P_j, m=1..8")

    # A one-factor mutation can preserve sampled ratio signs but must not pass
    # the derived connection identity.
    m = 3
    mutated = sum(sp.Rational(math.factorial(2 * m + 1), math.factorial(2 * j + 1))
                  * A[2 * j + 1] for j in range(m))
    assert sp.expand(sp.diff(A[2 * m + 1], y, 2) - mutated) != 0
    print("[PASS] mutation control: deleting H_d^odd/d breaks derivative lowering at m=3")


def gate_connection_ratio() -> None:
    d, j, H = sp.symbols("d j H", integer=True, positive=True)
    m = d + j

    def ratio(dd, hh):
        return 2 * (4 * dd**2 - 1) * hh / (dd * (4 * m - 2 * dd + 3))

    difference = sp.factor(ratio(d + 1, H + 1 / (2 * d + 1)) - ratio(d, H))
    numerator, denominator = sp.together(difference).as_numer_denom()
    polynomial = sp.Poly(sp.expand(numerator), H, d, j)
    assert all(value > 0 for value in polynomial.coeffs())
    assert denominator == d * (d + 1) * (2 * d + 4 * j + 1) * (2 * d + 4 * j + 3)
    digest = hashlib.sha256(sp.srepr(sp.expand(numerator)).encode()).hexdigest()
    print("[PASS] d_(m,j)/ell_(m,j) strictly decreases with j; "
          f"positive monomials={len(polynomial.terms())}, sha256={digest}")

    # Pin the two determinant orientations used in the proof.
    b0, b1, b2 = sp.symbols("b0 b1 b2", positive=True)
    r0, r1, r2 = sp.symbols("r0 r1 r2")
    k = sp.symbols("k")
    det2 = sp.Matrix([[b0, b1], [b0 * r0, b1 * r1]]).det()
    assert sp.expand(det2 - b0 * b1 * (r1 - r0)) == 0
    det3 = sp.Matrix([
        [b0, b1, b2],
        [k * b0, (k + 1) * b1, (k + 2) * b2],
        [r0 * b0, r1 * b1, r2 * b2],
    ]).det()
    assert sp.expand(det3 - b0 * b1 * b2 * (r2 - 2 * r1 + r0)) == 0
    print("[PASS] determinant orientations: first and second ratio differences")


def gate_exact_rows(mmax: int = 250) -> None:
    # Coefficients of P_m, exact integers.  Column scaling b_(m,k)=(2k+1)!p_(m,k).
    old = [1] + [0] * mmax
    row = [5, 1] + [0] * (mmax - 1)
    cells = convex_cells = 0
    min_gap = None
    min_convex = None
    gap_cell = convex_cell = None
    digest = hashlib.sha256()
    for m in range(1, mmax + 1):
        if m >= 2:
            X = [None] + [F(math.factorial(2 * k + 1) * row[k],
                              math.factorial(2 * k - 1) * row[k - 1])
                          for k in range(1, m + 1)]
            for k in range(2, m + 1):
                gap = X[k - 1] - X[k]
                assert gap >= 0
                cells += 1
                if min_gap is None or gap < min_gap:
                    min_gap, gap_cell = gap, (m, k)
                digest.update(f"g,{m},{k},{gap.numerator}/{gap.denominator}\n".encode())
            for k in range(2, m):
                curvature = X[k - 1] - 2 * X[k] + X[k + 1]
                assert curvature >= 0
                convex_cells += 1
                if min_convex is None or curvature < min_convex:
                    min_convex, convex_cell = curvature, (m, k)
                # The telescoped corridor is an exact consequence, checked
                # here independently at every finite replay cell.
                assert X[k] >= X[1] - (k - 1) * (X[1] - X[2])
                digest.update(f"c,{m},{k},{curvature.numerator}/{curvature.denominator}\n".encode())
        if m == mmax:
            break
        diagonal = 8 * m * m + 12 * m + 5
        back = (2 * m * (2 * m + 1)) ** 2
        new = [0] * (mmax + 1)
        for k in range(m + 2):
            new[k] = (row[k - 1] if k else 0) + diagonal * row[k] - back * old[k]
        old, row = row, new
    print(f"[PASS] exact row replay m<= {mmax}: decreasing cells={cells}, "
          f"convex cells={convex_cells}, min_gap_cell={gap_cell}, "
          f"min_convex_cell={convex_cell}, sha256={digest.hexdigest()}")


def main() -> int:
    gate_connection_identities()
    gate_connection_ratio()
    gate_exact_rows()
    print("ALL A12 SIGN GATES PASS")
    print("THEOREM: X_k is decreasing and convex; every decrement is bounded by X_1-X_2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
