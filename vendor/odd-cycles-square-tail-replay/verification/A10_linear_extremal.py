#!/usr/bin/env python3
"""Exact certificate for the linear-lowering extremal theorem.

For a PF row F of degree m and an adjacent-row lowering profile L satisfying

    L_0 = 0,  L_m = r,
    Delta^2 L_k >= 0,
    Delta^3 L_k <= 0,

put ell=L_1, mu_k=1-L_k/r, and mu^0_k=1-k*ell/r.  The theorem certified here
is

    L((1+u)(mu*F))_t <= L((1+u)(mu^0*F))_t,   3 <= t <= m-1.

The proof has two pieces.  The manuscript proves the pure defect row
log-concave by hand.  The mixed polarization has three coefficients in the
positive scale x=F_t/F_{t-1}; this program proves all three nonnegative on the
exact admissible polytope.  It uses degree-aware Newton bounds and a
Bernstein-basis certificate.  Every operation below is symbolic over Q;
sampled rows play no role.
"""

from __future__ import annotations

import hashlib
import itertools
import math

import sympy as sp


FAILS: list[str] = []


def gate(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}{'  ' + detail if detail else ''}")
    if not condition:
        FAILS.append(name)


def bernstein_coefficients(poly: sp.Expr, variables: list[sp.Symbol]):
    """Power-basis to tensor-product Bernstein coefficients on [0,1]^d."""
    p = sp.Poly(sp.expand(poly), *variables)
    degrees = p.degree_list()
    power = p.as_dict()
    out = {}
    for index in itertools.product(*(range(d + 1) for d in degrees)):
        value = sp.S.Zero
        for monomial in itertools.product(*(range(i + 1) for i in index)):
            coefficient = power.get(monomial, sp.S.Zero)
            if coefficient == 0:
                continue
            factor = sp.S.One
            for i, j, degree in zip(index, monomial, degrees):
                factor *= sp.Rational(math.comb(i, j), math.comb(degree, j))
            value += coefficient * factor
        out[index] = sp.factor(value)
    return degrees, out


def nonnegative_kn_polynomial(value: sp.Expr, K: sp.Symbol, N: sp.Symbol) -> bool:
    """K,N are nonnegative integers; nonnegative monomial coefficients suffice."""
    return all(c >= 0 for c in sp.Poly(sp.expand(value), K, N).coeffs())


def certificate_digest(values: dict) -> str:
    payload = "\n".join(f"{key}:{sp.srepr(values[key])}" for key in sorted(values))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def main() -> int:
    # k=t-2 >= 1; n=m-k >= 3.  K=k-1 and N=n-3 are therefore nonnegative.
    K, N = sp.symbols("K N", integer=True, nonnegative=True)
    k, n, t, m = K + 1, N + 3, K + 3, K + N + 4

    # Unit-cube parameters.  D is the used endpoint-slope budget; P is the
    # curvature share of the local slope; T and PH are the next two curvature
    # ratios; E locates the accumulated defect between its sharp two bounds.
    D, P, T, PH, E = sp.symbols("D P T PH E", nonnegative=True)
    cube = [D, P, T, PH, E]

    x, U, W = sp.symbols("x U W", positive=True)
    a, p, h, q, v = sp.symbols("a p h q v", nonnegative=True)
    d = a + p

    # Local normalized actual profile at indices k,...,k+3.  The normalization
    # is mu_k=1.  h is the defect L_k-kL_1 in the same scale.
    mu = [1, 1 - d, 1 - 2 * d - q, 1 - 3 * d - 2 * q - v]
    defect = [h, h + p, h + 2 * p + q, h + 3 * p + 2 * q + v]
    f = [U / x, 1, x, W * x**2]

    actual = [mu[j] * f[j] for j in range(4)]
    delta = [defect[j] * f[j] for j in range(4)]
    # Coefficients of (1+u)(mu*F) and (1+u)(defect*F) at t-1,t,t+1.
    Gm, G0, Gp = actual[0] + actual[1], actual[1] + actual[2], actual[2] + actual[3]
    Hm, H0, Hp = delta[0] + delta[1], delta[1] + delta[2], delta[2] + delta[3]
    cross = sp.expand(2 * G0 * H0 - Gm * Hp - Hm * Gp)
    cross_poly = sp.Poly(cross, x)
    gate("mixed polarization has degree two in x", cross_poly.degree() == 2)
    C = [sp.factor(cross_poly.coeff_monomial(x**j)) for j in range(3)]

    # Each coefficient is decreasing in the only row-ratio variable it uses.
    # The derivatives are certified on a weaker box (d<=1/n), so substituting
    # the degree-aware Newton maxima is safe on the sharper endpoint polytope.
    Umax = sp.factor((t - 1) * (n - 1) / (t * n))
    Wmax = sp.factor(t * (n - 2) / ((t + 1) * (n - 1)))
    Z = sp.symbols("Z", positive=True)
    Zmax = sp.factor(Umax * Wmax)

    weak_d = D / n
    weak_p = P * weak_d
    weak_q = T * weak_p / k
    weak_v = PH * weak_q
    weak_h = weak_p * (k - 1) * (1 + E) / 2
    weak_a = weak_d - weak_p
    weak_sub = {a: weak_a, p: weak_p, h: weak_h, q: weak_q, v: weak_v}

    derivative_expressions = [
        -sp.diff(C[0], U),
        -sp.diff(C[1].subs(U * W, Z), Z),
        -sp.diff(C[2], W),
    ]
    derivative_counts = []
    for j, expression in enumerate(derivative_expressions):
        reduced = sp.cancel(expression.subs(weak_sub) / weak_p)
        numerator, denominator = sp.fraction(reduced)
        degrees, coeffs = bernstein_coefficients(numerator, cube)
        good = all(nonnegative_kn_polynomial(c, K, N) for c in coeffs.values())
        derivative_counts.append(len(coeffs))
        gate(f"row-ratio monotonicity C{j}", good,
             f"Bernstein degree={degrees}, coefficients={len(coeffs)}")
        gate(f"row-ratio monotonicity C{j} denominator positive",
             nonnegative_kn_polynomial(denominator, K, N))

    # Strong endpoint parameterization.  Convexity gives p<=d.  Decreasing
    # curvature gives q<=p/k and v<=q.  The endpoint mu_m=0 gives
    # n*d+(n-1)*q+(n-2)*v <= 1.  These are represented exactly by five unit
    # parameters; no relaxation remains in the slope budget.
    alpha = P * T / k
    beta = alpha * PH
    endpoint_denominator = sp.factor(n + (n - 1) * alpha + (n - 2) * beta)
    strong_d = D / endpoint_denominator
    strong_p = P * strong_d
    strong_q = alpha * strong_d
    strong_v = beta * strong_d
    strong_h = strong_p * (k - 1) * (1 + E) / 2
    strong_a = strong_d - strong_p
    strong_sub = {a: strong_a, p: strong_p, h: strong_h,
                  q: strong_q, v: strong_v}

    bounded = [C[0].subs(U, Umax), C[1].subs(U * W, Zmax), C[2].subs(W, Wmax)]
    all_coefficients = {}
    expected_degrees = [(1, 1, 2, 1, 1), (1, 1, 2, 2, 1), (1, 1, 2, 2, 1)]
    expected_counts = [48, 72, 72]
    for j, expression in enumerate(bounded):
        reduced = sp.cancel(expression.subs(strong_sub) / strong_p)
        numerator, denominator = sp.fraction(reduced)
        degrees, coeffs = bernstein_coefficients(numerator, cube)
        good = all(nonnegative_kn_polynomial(c, K, N) for c in coeffs.values())
        strict = all(c != 0 for c in coeffs.values())
        gate(f"mixed coefficient C{j} nonnegative", good and strict,
             f"Bernstein degree={degrees}, coefficients={len(coeffs)}")
        gate(f"mixed coefficient C{j} certificate shape",
             tuple(degrees) == expected_degrees[j] and len(coeffs) == expected_counts[j])
        gate(f"mixed coefficient C{j} denominator positive",
             nonnegative_kn_polynomial(denominator, K, N))
        for index, value in coeffs.items():
            all_coefficients[(j, index)] = value

    digest = certificate_digest(all_coefficients)
    print(f"  certificate SHA-256: {digest}")
    gate("Bernstein certificate digest",
         digest == "6a0bd9ec9d03f693df0b8dd5babd6bd4771014daefc98fe852480dd035e0efc4")

    # The quadratic correction for the linear profile is increasing in ell on
    # 0<=ell<=r/m.  It is enough to prove m*C1-2*C2 >=0.  Degree-aware Newton
    # again makes all three x-coefficients positive.
    s = (1 - U) + x * (1 - U * W) + x**2 * (1 - W)
    Arow = x * (1 - U * W) + 2 * x**2 * (1 - W)
    D0 = x**2 + U + 2 * U * W * x
    C1 = sp.expand(2 * s * (t - 1) + Arow)
    C2 = sp.expand(s * (t - 1)**2 + Arow * (t - 1) + D0)
    lam = sp.symbols("lam", nonnegative=True)
    y = U / x
    z = W * x
    mu0 = lambda index: 1 - index * lam
    sigma_linear = ((x * mu0(t) + mu0(t - 1))**2
                    - (mu0(t - 1) + y * mu0(t - 2))
                    * x * (z * mu0(t + 1) + mu0(t)))
    gate("linear-profile reduction",
         sp.expand(s - sigma_linear - C1 * lam + C2 * lam**2) == 0)
    monotone = sp.Poly(sp.expand(m * C1 - 2 * C2), x)
    monotone_coeffs = []
    for power in range(3):
        coefficient = monotone.coeff_monomial(x**power)
        if power == 0:
            lower = coefficient.subs(U, Umax)
        elif power == 1:
            lower = coefficient.subs(U * W, Zmax)
        else:
            lower = coefficient.subs(W, Wmax)
        lower = sp.factor(lower)
        monotone_coeffs.append(lower)
        gate(f"linear correction monotonicity x^{power}",
             nonnegative_kn_polynomial(sp.together(lower).as_numer_denom()[0], K, N))

    # The still coarser top comparison sigma_linear <= (1-(t-2)lam)^2 s.
    # Its derivative on 0<=lam<=1/m again has three positive Newton-bounded
    # coefficients.  This is the theorem gate behind (TOP) in the note.
    ktop = t - 2
    top_linear = 2 * s + Arow
    top_quadratic = s * (2 * t - 3) + Arow * (t - 1) + D0
    top_monotone = sp.Poly(sp.expand(m * top_linear - top_quadratic), x)
    for power in range(3):
        coefficient = top_monotone.coeff_monomial(x**power)
        if power == 0:
            lower = coefficient.subs(U, Umax)
        elif power == 1:
            lower = coefficient.subs(U * W, Zmax)
        else:
            lower = coefficient.subs(W, Wmax)
        gate(f"top comparison monotonicity x^{power}",
             nonnegative_kn_polynomial(sp.together(lower).as_numer_denom()[0], K, N))
    gate("top comparison identity",
         sp.expand((s - sigma_linear)
                   - (s - (1 - ktop * lam)**2 * s)
                   - top_linear * lam + top_quadratic * lam**2) == 0)

    print()
    if FAILS:
        print(f"*** {len(FAILS)} GATE FAILURES ***")
        for failure in FAILS:
            print(f"  {failure}")
        return 1
    print("ALL GATES PASS")
    print("THEOREM: the linear lowering profile maximizes the predecessor Turan determinant")
    print(f"Bernstein coefficients: mixed={sum(expected_counts)}, monotonicity={sum(derivative_counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
