#!/usr/bin/env python3
"""Theorem-backed algebra checks for the CONN2 connection ratio.

The manuscript proves the all-row connection theorem. This program checks its
load-bearing factorial cancellation, odd-harmonic closed form,
adjacent-difference identity, and positive-coefficient certificate symbolically.
A finite Fraction replay guards the endpoints and the original factorial
indexing.
"""

from fractions import Fraction as F
from math import factorial as fac

import sympy


MS = (6, 8, 12, 20, 30, 45, 60)


def a(m, j):
    return F(
        fac(2 * m + 1) * (2 * m + 2 * j + 3),
        2 * fac(2 * j + 1) * (2 * m - 2 * j - 1) * (2 * m - 2 * j + 1),
    )


def g(m, i):
    return (m - 1 - i) * a(m, i) + sum(
        a(m, j) * a(j, i) for j in range(i + 1, m)
    )


def odd_harmonic(d):
    """H(d)=sum_{r=0}^{d-2} 1/(2r+1)."""
    return sum((F(1, 2 * r + 1) for r in range(d - 1)), F(0))


def G_closed(m, j):
    d = m - j
    eta = 2 * m + 2 * j + 3
    H = odd_harmonic(d)
    return F(d - 1) + F(4 * d * d - 1, 4 * eta * d * (d * d - 1)) * (
        (eta * eta - 2 * d * d + 1) * H
        + F((d - 1) ** 2 * (eta * eta + 2 * d + 1), 2 * d - 1)
    )


def certificate_polynomials(n, j):
    A0 = 36920 + 58918*n + 36550*n**2 + 11056*n**3 + 1640*n**4 + 96*n**5
    A1 = 38200 + 50040*n + 23704*n**2 + 4800*n**3 + 352*n**4
    A2 = 13520 + 14208*n + 4800*n**2 + 512*n**3
    A3 = 1600 + 1280*n + 256*n**2
    B0 = 44336 + 76132*n + 51022*n**2 + 16690*n**3 + 2672*n**4 + 168*n**5
    B1 = 36544 + 48048*n + 23584*n**2 + 5064*n**3 + 400*n**4
    B2 = 10496 + 9600*n + 3024*n**2 + 320*n**3
    B3 = 1024 + 512*n + 64*n**2
    return A0 + j*A1 + j**2*A2 + j**3*A3, B0 + j*B1 + j**2*B2 + j**3*B3


def gate_symbolic_certificate():
    n, j, h, x = sympy.symbols("n j h x", nonnegative=True)
    d = n + 3
    eta = 2*n + 4*j + 9

    # One summand of sum_l a(m,l)a(l,j)/a(m,j), with x=l-j.
    raw = ((4*d**2 - 1) * (eta + 2*x) * (eta - 2*d + 2*x)
           / (2*eta*(4*(d-x)**2 - 1)*(4*x**2 - 1)))
    common = (4*d**2 - 1) / (16*d*eta)
    c1 = (eta - 1)*(2*d - eta + 1)/(d + 1)
    c2 = (eta + 1)*(2*d + eta + 1)/(d + 1)
    c3 = (eta - 1)*(2*d + eta - 1)/(d - 1)
    c4 = (eta + 1)*(2*d - eta - 1)/(d - 1)
    apart = common * (
        c1/(2*x + 1) + c2/(-2*d + 2*x - 1)
        - c3/(-2*d + 2*x + 1) - c4/(2*x - 1)
    )
    assert sympy.cancel(raw - apart) == 0

    # Summing x=1..d-1: the first two odd-harmonic strings equal
    # H-1+1/(2d-1), and the last two equal H.
    s1 = h - 1 + 1/(2*d - 1)
    summed = common * ((c1 - c2)*s1 + (c3 - c4)*h)
    claimed_sum = ((4*d**2 - 1)/(4*eta*d*(d**2 - 1))
                   * ((eta**2 - 2*d**2 + 1)*h
                      + (d - 1)**2*(eta**2 + 2*d + 1)/(2*d - 1)))
    assert sympy.cancel(summed - claimed_sum) == 0

    def G(dd, ee, hh):
        return (dd - 1 + (4*dd**2 - 1)/(4*ee*dd*(dd**2 - 1))
                * ((ee**2 - 2*dd**2 + 1)*hh
                   + (dd - 1)**2*(ee**2 + 2*dd + 1)/(2*dd - 1)))

    A, B = certificate_polynomials(n, j)
    den = 4*eta*(eta + 2)*d*(d - 1)*(d - 2)*(d + 1)
    difference = G(d - 1, eta + 2, h - 1/(2*d - 3)) - G(d, eta, h)
    certificate = A/den*(h - sympy.Rational(4, 3)) + B/(3*den)
    assert sympy.cancel(difference - certificate) == 0

    a_coeffs = sympy.Poly(A, j, n).coeffs()
    b_coeffs = sympy.Poly(B, j, n).coeffs()
    assert len(a_coeffs) == 18 and len(b_coeffs) == 18
    assert all(value > 0 for value in a_coeffs + b_coeffs)
    print("[THEOREM-BACKED PASS] symbolic G sum and CONN2+ identity; "
          "36/36 polynomial coefficients positive")


def gate_determinant_orientations():
    k = sympy.symbols("k")
    r0, r1, r2, r3 = sympy.symbols("r0 r1 r2 r3")
    d3 = sympy.Matrix([
        [1, 1, 1],
        [k, k + 1, k + 2],
        [r0, r1, r2],
    ]).det()
    assert sympy.expand(d3 - (r2 - 2*r1 + r0)) == 0

    d4 = sympy.Matrix([
        [1, 1, 1, 1],
        [k, k + 1, k + 2, k + 3],
        [k**2, (k + 1)**2, (k + 2)**2, (k + 3)**2],
        [r0, r1, r2, r3],
    ]).det()
    assert sympy.expand(d4 - 2*(r3 - 3*r2 + 3*r1 - r0)) == 0

    m = sympy.symbols("m")
    row_change = sympy.Matrix([
        [1, 0, 0],
        [m, -1, 0],
        [m*(m - 1), -(2*m - 1), 1],
    ])
    assert row_change.det() == -1
    # [m,i,j,m-1] has three inversions relative to [i,j,m-1,m].
    assert sum(left > right for pos, left in enumerate((3, 0, 1, 2))
               for right in (3, 0, 1, 2)[pos + 1:]) == 3
    print("[THEOREM-BACKED PASS] 3x3/4x4 determinant factors and row signs")


def gate_factorial_replay():
    checked_g = checked_d = 0
    for m in range(2, 31):
        ratios = []
        for j in range(m - 1):
            direct = g(m, j) / a(m, j)
            assert direct == G_closed(m, j)
            ratios.append(direct)
            checked_g += 1
        for j in range(m - 2):
            d = m - j
            eta = 2*m + 2*j + 3
            n = d - 3
            A, B = certificate_polynomials(n, j)
            den = 4*eta*(eta + 2)*d*(d - 1)*(d - 2)*(d + 1)
            certificate = F(A, den)*(odd_harmonic(d) - F(4, 3)) + F(B, 3*den)
            assert ratios[j + 1] - ratios[j] == certificate > 0
            checked_d += 1
    print(f"[THEOREM-BACKED PASS] factorial/index replay: {checked_g} G values, "
          f"{checked_d} strict adjacent differences, m=2..30")


def report_selected_rows():
    print("A8-CONN2 selected-row regression")
    print(f"  {'m':>4} {'pairs':>7} {'viol':>5} {'min D':>12} {'argmin j':>9}")
    total = bad_total = 0
    for m in MS:
        ratio = {i: G_closed(m, i) for i in range(m - 1)}
        pairs = [(i, j) for i in range(m - 2) for j in range(i + 1, m - 1)]
        bad = [(i, j) for i, j in pairs if ratio[j] <= ratio[i]]
        diff = {i: ratio[i + 1] - ratio[i] for i in range(m - 2)}
        lo = min(diff, key=diff.get)
        total += len(pairs)
        bad_total += len(bad)
        print(f"  {m:>4} {len(pairs):>7} {len(bad):>5} "
              f"{float(diff[lo]):>12.8f} {lo:>9}")
    assert bad_total == 0
    print(f"  [PASS] strict CONN2 on {total} selected-row pairs; theorem range is all m")


def main():
    gate_symbolic_certificate()
    gate_determinant_orientations()
    gate_factorial_replay()
    report_selected_rows()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
