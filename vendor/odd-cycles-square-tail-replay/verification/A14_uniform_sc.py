#!/usr/bin/env python3
"""Exact certificate for the uniform joint-window scalar comparison.

This closes A10's remaining (LSC) obligation on

    r >= 504,  t >= 6,  a >= 5t,  t < sqrt(a)+2.

The current row enters only through A12--A13's FIVE-CORRIDOR.  The core
estimate is reduced to a rational inequality in two nonnegative variables.
For t>=8 its cleared polynomial has positive monomial coefficients; t=6,7
are handled by exact Sturm sequences.  No sampled row or floating-point value
enters a verdict.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib

import sympy as sp


FAILS: list[str] = []


def gate(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}{'  ' + detail if detail else ''}")
    if not condition:
        FAILS.append(name)


def sign_variations(values: list[sp.Expr]) -> int:
    signs = [sp.sign(value) for value in values if value != 0]
    assert all(value in (-1, 1) for value in signs)
    return sum(signs[j] != signs[j - 1] for j in range(1, len(signs)))


def primitive_sturm_digest(sequence: list[sp.Expr], variable: sp.Symbol) -> str:
    rows = []
    for expression in sequence:
        polynomial = sp.Poly(expression, variable, domain=sp.QQ)
        _, primitive = polynomial.primitive()
        if primitive.LC() < 0:
            primitive = -primitive
        rows.append(",".join(str(value) for value in primitive.all_coeffs()))
    return hashlib.sha256("\n".join(rows).encode()).hexdigest()


def core_certificate() -> None:
    # Local row variables.  U=xy and W=z/x; strict Newton gives 0<U,W<1.
    x, U, W = sp.symbols("x U W", positive=True)
    t = sp.symbols("t", integer=True, positive=True)
    v0 = 1 - U
    v1 = x * (1 - U * W)
    v2 = x**2 * (1 - W)
    s = v0 + v1 + v2
    Arow = v1 + 2 * v2
    C1 = 2 * (t - 1) * s + Arow

    gate("three-term decomposition of s", sp.expand(s - (v0 + v1 + v2)) == 0)
    gate("weighted decomposition of C1",
         sp.expand(C1 - (2 * (t - 1) * v0 + (2 * t - 1) * v1 + 2 * t * v2)) == 0)
    gate("C1 lower-weight remainder is v2",
         sp.expand(C1 - (2 * (t - 1) * v0 + (2 * t - 1) * (v1 + v2)) - v2) == 0)
    gate("exact s/V quotient",
         sp.cancel(s / v0 - (1 + x + (x * U + x**2) * (1 - W) / (1 - U))) == 0)

    # Closed corridor parameter.  If Y=X_(t-1), then
    #   Y>=a-5(t-2), X_t>=Y-5, U>=alpha(1-5/Y), W<=beta.
    # The hypothesis a>=5t gives Y>=10.
    Q, T, S = sp.symbols("Q T S", nonnegative=True)
    tt = T + 6
    Y = Q + 10
    aa_upper = Y + 5 * (tt - 2)
    alpha = (tt - 1) * (2 * tt - 1) / (tt * (2 * tt + 1))
    beta = tt * (2 * tt + 1) / ((tt + 1) * (2 * tt + 3))
    x0 = (Y - 5) / (2 * tt * (2 * tt + 1))
    u0 = alpha * (1 - 5 / Y)
    ratio0 = (1 - beta) / (1 - u0)
    R0 = sp.factor(1 + x0 + (x0 * u0 + x0**2) * ratio0)

    # For R=s/V>=1, ((2t-1)R-1)^2/R is increasing.  The following
    # rational inequality therefore proves C1^2 >= (101/100) 2 a s V.
    R = sp.symbols("R", positive=True)
    w = 2 * tt - 1
    derivative = sp.factor(sp.diff((w * R - 1)**2 / R, R))
    gate("core quotient increases with R",
         sp.factor(derivative - (w * R - 1) * (w * R + 1) / R**2) == 0)

    constant = sp.Rational(101, 100)
    core = sp.together((w * R0 - 1)**2 - 2 * constant * aa_upper * R0)
    numerator, denominator = core.as_numer_denom()
    polynomial = sp.Poly(numerator, Q, T).primitive()[1]
    gate("core denominator positive form",
         sp.factor(denominator) ==
         400 * (T + 6)**2 * (T + 7)**2 * (2 * T + 13)**2 * (2 * T + 15)**2
         * (4 * Q * T + 23 * Q + 10 * T**2 + 145 * T + 505)**2)

    # t>=8: T=S+2, Q,S>=0.
    shifted = sp.Poly(sp.expand(polynomial.as_expr().subs(T, S + 2)), Q, S)
    payload = "\n".join(f"{monomial}:{coefficient}"
                          for monomial, coefficient in shifted.terms())
    digest = hashlib.sha256(payload.encode()).hexdigest()
    gate("CORE101 t>=8 coefficientwise positivity",
         len(shifted.terms()) == 74 and all(value > 0 for value in shifted.coeffs()),
         f"degree={shifted.degree_list()}, terms={len(shifted.terms())}")
    gate("CORE101 coefficient digest",
         digest == "4be33f93a2cc465b1389d7f89325da73ad6d4f5a3ec4e72eb220fe7c93998bb1",
         digest)

    # The two endpoint columns are exact half-line statements, not samples.
    expected = {
        0: "8f47094b06d2fb32c8588c98a23138e34a246967b6b29374644fa8d9b007a5cf",
        1: "a2d378c626345c1b4e887e9adb1a4de27ae5b9ca92025b829e494c97ab9ffd64",
    }
    for Tv in (0, 1):
        univariate = sp.Poly(polynomial.as_expr().subs(T, Tv), Q, domain=sp.QQ)
        sturm = sp.sturm(univariate.as_expr(), Q)
        at_zero = [sp.Poly(row, Q).eval(0) for row in sturm]
        at_infinity = [sp.Poly(row, Q).LC() for row in sturm]
        v0_count = sign_variations(at_zero)
        vinf_count = sign_variations(at_infinity)
        sdigest = primitive_sturm_digest(sturm, Q)
        gate(f"CORE101 t={Tv + 6} has no root on Q>=0",
             univariate.eval(0) > 0 and v0_count == vinf_count,
             f"Sturm variations={v0_count}->{vinf_count}, length={len(sturm)}")
        gate(f"CORE101 t={Tv + 6} Sturm digest", sdigest == expected[Tv], sdigest)

    # Mutation: 103/100 exceeds the true lower bound of this corridor
    # reduction.  A conclusion-only or sign-only gate would miss this.
    mutation = sp.together((w * R0 - 1)**2
                           - 2 * sp.Rational(103, 100) * aa_upper * R0)
    mutation_value = mutation.as_numer_denom()[0].subs({T: 0, Q: 95})
    gate("CORE101 mutation control rejects 103/100", mutation_value < 0,
         f"exact numerator={mutation_value}")


def correction_certificate() -> None:
    x, U, W, t = sp.symbols("x U W t", positive=True)
    s = (1 - U) + x * (1 - U * W) + x**2 * (1 - W)
    Arow = x * (1 - U * W) + 2 * x**2 * (1 - W)
    D0 = x**2 + U + 2 * U * W * x
    C1 = 2 * s * (t - 1) + Arow
    C2 = s * (t - 1)**2 + Arow * (t - 1) + D0
    gate("quadratic correction identity",
         sp.expand((t - 1) * C1 - C2 - (s * (t - 1)**2 - D0)) == 0)

    # Degree-aware Newton implies U<=(t-1)/t, W<=t/(t+1), and hence
    # UW<=(t-1)/(t+1).  Every x coefficient is then nonnegative for t>=3.
    difference = sp.Poly(sp.expand(s * (t - 1)**2 - D0), x)
    bounds = [
        sp.factor(difference.coeff_monomial(x**0).subs(U, (t - 1) / t)),
        sp.factor(difference.coeff_monomial(x**1).subs(U * W, (t - 1) / (t + 1))),
        sp.factor(difference.coeff_monomial(x**2).subs(W, t / (t + 1))),
    ]
    expected = [
        (t - 2) * (t - 1) / t,
        2 * (t - 2) * (t - 1) / (t + 1),
        t * (t - 3) / (t + 1),
    ]
    gate("C2 <= (t-1)C1 coefficient bounds",
         all(sp.factor(left - right) == 0 for left, right in zip(bounds, expected)))

    h, r = sp.symbols("h r", positive=True)
    gate("finite-r correction keeps its correlation",
         sp.factor(C1 / h - C2 / (r * h**2)
                   - C1 / h * (1 - (t - 1) / (r * h))
                   - ((t - 1) * C1 - C2) / (r * h**2)) == 0)


def scalar_factor_certificate() -> None:
    # A7 gives beta=b(1+A)-A, eps=beta(2-beta), and R_(r-1)>=4r-3.
    # Since 0<beta<b<1, eps<2/(4r-3).
    b, A = sp.symbols("b A", positive=True)
    beta = b * (1 + A) - A
    gate("normalizer comparison beta-b=A(b-1)", sp.expand(beta - b - A * (b - 1)) == 0)

    # In the joint window t<sqrt(a)+2 and a>=1,
    # 1-(t-1)/(r sqrt(a)) > 1-2/r.  Both retained factors increase in r.
    r0 = 504
    factor = (1 - F(2, 4 * r0 - 3)) * (1 - F(2, r0))
    surplus = F(101, 100) * factor * factor - 1
    gate("uniform finite-r factor at r=504", surplus > 0,
         f"factor={factor}, squared surplus={surplus}")
    gate("r=503 is correctly rejected",
         F(101, 100) * (1 - F(2, 4 * 503 - 3))**2 * (1 - F(2, 503))**2 < 1)


def main() -> int:
    print("A14 -- exact uniform (SC) certificate")
    core_certificate()
    correction_certificate()
    scalar_factor_certificate()
    print()
    if FAILS:
        print(f"*** {len(FAILS)} A14 GATE FAILURES ***")
        for failure in FAILS:
            print("  ", failure)
        return 1
    print("ALL A14 GATES PASS")
    print("THEOREM: (LSC), hence (SQ), holds for r>=504, t>=6, "
          "a>=5t, t<sqrt(a)+2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
