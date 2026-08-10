#!/usr/bin/env python3
"""Independent exact certificate for a fixed-column tail.

The analytic state is enclosed by rational fixed-point intervals.  Pi is
enclosed by Machin's formula, logarithms by the atanh series, and every
polynomial gate is an exact Fraction comparison.  Decimal output is
orientation only.

This verifier is deliberately independent of A1_verify_tail.py.  The exact
wall-to-ladder identity is replayed separately by A1_reduce_check.py.

Usage:
    python3 A1_verify_tail_independent.py [t] [M0] [--controls]

The default remains the historical t=2, M0=312 certificate.  For a supplied
column t, M0 is the certified tail anchor and the conclusion starts at
r=M0+1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from fractions import Fraction as F
from pathlib import Path


PREC = 256
SCALE = 1 << PREC
M0 = 312
T = 2
KMAX = T + 1


class Iv:
    """Closed interval with endpoints in SCALE**-1 Z."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo: int, hi: int):
        assert lo <= hi
        self.lo = lo
        self.hi = hi

    @staticmethod
    def exact(p: int, q: int = 1) -> "Iv":
        n = p * SCALE
        assert q > 0 and n % q == 0
        return Iv(n // q, n // q)

    @staticmethod
    def frac(p: int, q: int = 1) -> "Iv":
        assert q > 0
        n = p * SCALE
        return Iv(n // q, -((-n) // q))

    def __add__(self, other: "Iv") -> "Iv":
        return Iv(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other: "Iv") -> "Iv":
        return Iv(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self) -> "Iv":
        return Iv(-self.hi, -self.lo)

    def __mul__(self, other: "Iv") -> "Iv":
        vals = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Iv(min(vals) // SCALE, -((-max(vals)) // SCALE))

    def inv(self) -> "Iv":
        assert self.lo > 0
        ss = SCALE * SCALE
        return Iv(ss // self.hi, -((-ss) // self.lo))

    def __truediv__(self, other: "Iv") -> "Iv":
        return self * other.inv()

    def scal(self, p: int, q: int = 1) -> "Iv":
        assert q > 0
        lo, hi = self.lo * p, self.hi * p
        if p < 0:
            lo, hi = hi, lo
        return Iv(lo // q, -((-hi) // q))

    def flo(self) -> F:
        return F(self.lo, SCALE)

    def fhi(self) -> F:
        return F(self.hi, SCALE)

    def width(self) -> F:
        return F(self.hi - self.lo, SCALE)


ZERO = Iv(0, 0)
ONE = Iv.exact(1)


def pi_bounds(terms: int = 60) -> tuple[F, F]:
    """Machin formula with alternating-series remainder brackets."""

    def atan_inv(n: int) -> tuple[F, F]:
        total = F(0)
        lo = hi = F(0)
        for k in range(terms):
            term = F(1, (2 * k + 1) * n ** (2 * k + 1))
            total += term if k % 2 == 0 else -term
            nxt = F(1, (2 * k + 3) * n ** (2 * k + 3))
            if k % 2 == 0:
                lo, hi = total - nxt, total
            else:
                lo, hi = total, total + nxt
        return lo, hi

    a5l, a5h = atan_inv(5)
    a239l, a239h = atan_inv(239)
    return 16 * a5l - 4 * a239h, 16 * a5h - 4 * a239l


PI_LO, PI_HI = pi_bounds()
PI = Iv(int(PI_LO * SCALE), -int(-PI_HI * SCALE))


def log_bounds(x: F, terms: int = 200) -> tuple[F, F]:
    """Exact atanh-series enclosure of log(x), x > 0."""

    assert x > 0
    exponent = 0
    while x >= 2:
        x /= 2
        exponent += 1
    while x < 1:
        x *= 2
        exponent -= 1

    def atanh(z: F) -> tuple[F, F]:
        total = F(0)
        power = z
        for k in range(terms):
            total += power / (2 * k + 1)
            power *= z * z
        remainder = power / ((2 * terms + 1) * (1 - z * z))
        return total, total + remainder

    l2l, l2h = atanh(F(1, 3))
    l2l, l2h = 2 * l2l, 2 * l2h
    z = (x - 1) / (x + 1)
    xl, xh = atanh(z)
    xl, xh = 2 * xl, 2 * xh
    if exponent >= 0:
        return xl + exponent * l2l, xh + exponent * l2h
    return xl + exponent * l2h, xh + exponent * l2l


def padd(a: list[F], b: list[F], sign: int = 1) -> list[F]:
    n = max(len(a), len(b))
    return [
        (a[i] if i < len(a) else F(0))
        + sign * (b[i] if i < len(b) else F(0))
        for i in range(n)
    ]


def pmul(a: list[F], b: list[F]) -> list[F]:
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def pshift(a: list[F], base: F) -> list[F]:
    """Coefficients of a(base+x), ascending in x."""

    out = [F(0)] * len(a)
    for i, coeff in enumerate(a):
        for j in range(i + 1):
            out[j] += coeff * math.comb(i, j) * base ** (i - j)
    return out


def pplus(a: list[F]) -> list[F]:
    """Coefficientwise positive part; dominates a on x >= 0."""

    return [max(c, F(0)) for c in a]


def pint_discrete_upper(a: list[F], max_step: F) -> list[F]:
    """Integral of a plus max_step*a, for a with nonnegative coefficients."""

    assert all(c >= 0 for c in a)
    integral = [F(0)] + [c / (i + 1) for i, c in enumerate(a)]
    return padd(integral, [max_step * c for c in a])


def peval(a: list[F], x: F) -> F:
    return sum(c * x**i for i, c in enumerate(a))


def interval_ladder(limit: int, kmax: int | None = None):
    """Outward-rounded h,T,B,C ladder through m=limit."""

    if kmax is None:
        kmax = KMAX

    h = ONE
    tm = ONE
    b = [ONE] + [ZERO] * kmax
    c = [ONE, ONE] + [ZERO] * (kmax - 1)
    for m in range(1, limit + 1):
        h = h.scal((2 * m) ** 2, (2 * m - 1) ** 2)
        q = h.inv()
        v = h.scal(1, (2 * m + 1) ** 2)
        bp, cp = b[:], c[:]
        tm = tm + q
        b = [tm] + [bp[k] + q * cp[k] for k in range(1, kmax + 1)]
        c = [ONE] + [cp[k] + v * b[k - 1] for k in range(1, kmax + 1)]
        yield m, h, q, tm, b, c


def exact_ladder(limit: int, kmax: int | None = None):
    if kmax is None:
        kmax = KMAX
    h = F(1)
    tm = F(1)
    b = [F(1)] + [F(0)] * kmax
    c = [F(1), F(1)] + [F(0)] * (kmax - 1)
    for m in range(1, limit + 1):
        h *= F((2 * m) ** 2, (2 * m - 1) ** 2)
        q = 1 / h
        v = h / (2 * m + 1) ** 2
        bp, cp = b[:], c[:]
        tm += q
        b = [tm] + [bp[k] + q * cp[k] for k in range(1, kmax + 1)]
        c = [F(1)] + [cp[k] + v * b[k - 1] for k in range(1, kmax + 1)]
        yield m, h, q, tm, b, c


def continuum_polys(a: list[F], b: list[F], kmax: int) -> tuple[list[list[F]], list[list[F]]]:
    """Polynomial solution beta_k'=gamma_k, gamma_k'=beta_{k-1}."""

    beta: list[list[F]] = []
    gamma: list[list[F]] = []
    for k in range(kmax + 1):
        poly = [F(0)] * (2 * k + 2)
        for j in range(k + 1):
            d = k - j
            poly[2 * d] += a[j] / math.factorial(2 * d)
            poly[2 * d + 1] += b[j] / math.factorial(2 * d + 1)
        while len(poly) > 1 and poly[-1] == 0:
            poly.pop()
        beta.append(poly)
        gamma.append([i * coeff for i, coeff in enumerate(poly)][1:] or [F(0)])
    return beta, gamma


def all_nonnegative(polys: list[list[F]]) -> bool:
    return all(all(c >= 0 for c in poly) for poly in polys)


def build_certificate(
    m0: int = M0, t: int = T, injection_scale: int = 1
) -> dict:
    assert m0 >= 2
    kmax = t + 1
    last = None
    for last in interval_ladder(m0, kmax):
        pass
    assert last is not None
    m, h0, q0, tm0, b0, c0 = last
    assert m == m0

    half_pi = PI.scal(1, 2)
    tau0 = half_pi * tm0
    beta0 = [half_pi * z for z in b0]

    powers = [ONE]
    for _ in range(1, 2 * kmax + 3):
        powers.append(powers[-1] * tau0)

    const_a = [ZERO]
    const_b = [ONE]
    for k in range(1, kmax + 1):
        prior_gamma = ZERO
        for j in range(k):
            d = k - j
            prior_gamma = prior_gamma + const_a[j] * powers[2 * d - 1].scal(
                1, math.factorial(2 * d - 1)
            )
            prior_gamma = prior_gamma + const_b[j] * powers[2 * d].scal(
                1, math.factorial(2 * d)
            )
        bk = c0[k] - prior_gamma

        prior_beta = ZERO
        for j in range(k):
            d = k - j
            prior_beta = prior_beta + const_a[j] * powers[2 * d].scal(
                1, math.factorial(2 * d)
            )
            prior_beta = prior_beta + const_b[j] * powers[2 * d + 1].scal(
                1, math.factorial(2 * d + 1)
            )
        ak = beta0[k] - prior_beta - bk * tau0
        const_a.append(ak)
        const_b.append(bk)

    alo = [z.flo() for z in const_a]
    ahi = [z.fhi() for z in const_a]
    blo = [z.flo() for z in const_b]
    bhi = [z.fhi() for z in const_b]
    alo[0] = ahi[0] = F(0)
    blo[0] = bhi[0] = F(1)
    t0l, t0h = tau0.flo(), tau0.fhi()

    beta_l, gamma_l = continuum_polys(alo, blo, kmax)
    beta_h, gamma_h = continuum_polys(ahi, bhi, kmax)
    beta_l = [pshift(poly, t0l) for poly in beta_l]
    beta_h = [pshift(poly, t0l) for poly in beta_h]
    gamma_l = [pshift(poly, t0l) for poly in gamma_l]
    gamma_h = [pshift(poly, t0l) for poly in gamma_h]
    beta_l[0] = beta_h[0] = [t0l, F(1)]
    gamma_l[0] = gamma_h[0] = [F(1)]

    log43_lo, log43_hi = log_bounds(F(4, 3))
    logm_lo, logm_hi = log_bounds(F(m0))
    tau_upper_const = PI_HI / 2 + log43_hi / 2 + F(3, 8 * m0)

    beta_abs_upper = [pplus(poly) for poly in continuum_polys(ahi, bhi, kmax)[0]]
    gamma_abs_upper = [pplus(poly) for poly in continuum_polys(ahi, bhi, kmax)[1]]
    tail_degree = max(
        len(poly) - 1 for poly in beta_abs_upper[:kmax] + gamma_abs_upper[:kmax]
    )
    assert 2 * logm_lo >= tail_degree

    def log_tail_moment(i: int) -> F:
        return sum(math.perm(i, p) * logm_hi ** (i - p) for p in range(i + 1)) / m0

    def tail_sum(poly: list[F]) -> F:
        total = F(0)
        for n, coeff in enumerate(poly):
            for i in range(n + 1):
                total += (
                    coeff
                    * math.comb(n, i)
                    * tau_upper_const ** (n - i)
                    * F(1, 2**i)
                    * log_tail_moment(i)
                )
        return total / 4

    d_beta = [F(0)] * (kmax + 1)
    d_gamma = [F(0)] * (kmax + 1)
    for k in range(1, kmax + 1):
        d_beta[k] = injection_scale * tail_sum(beta_abs_upper[k - 1])
        d_gamma[k] = injection_scale * tail_sum(
            padd(beta_abs_upper[k - 1], gamma_abs_upper[k - 1])
        )

    max_step = F(1, 2 * m0)
    err_beta = [[F(0)]]
    err_gamma = [[F(0)]]
    for k in range(1, kmax + 1):
        eg = [d_gamma[k]]
        if k >= 2:
            eg = padd(eg, pint_discrete_upper(err_beta[k - 1], max_step))
        eb = padd([d_beta[k]], pint_discrete_upper(eg, max_step))
        err_gamma.append(eg)
        err_beta.append(eb)

    beta_lo = [padd(beta_l[k], err_beta[k], -1) for k in range(kmax + 1)]
    beta_hi = [padd(beta_h[k], err_beta[k]) for k in range(kmax + 1)]
    gamma_lo = [padd(gamma_l[k], err_gamma[k], -1) for k in range(kmax + 1)]
    gamma_hi = [padd(gamma_h[k], err_gamma[k]) for k in range(kmax + 1)]
    beta_lo[0] = beta_hi[0] = beta_l[0]
    gamma_lo[0] = gamma_hi[0] = [F(1)]

    beta_prev_lo = [
        padd(beta_lo[k], [max_step * z for z in gamma_hi[k]], -1)
        for k in range(kmax + 1)
    ]
    beta_prev_hi = beta_hi
    gamma_prev_lo = [
        padd(
            gamma_lo[k],
            [max_step * z for z in (beta_hi[k - 1] if k >= 1 else [F(0)])],
            -1,
        )
        for k in range(kmax + 1)
    ]
    gamma_prev_hi = gamma_hi
    beta_prev_lo[0] = padd(beta_lo[0], [max_step], -1)
    gamma_prev_lo[0] = [F(1)]

    def get(rows: list[list[F]], i: int) -> list[F]:
        return rows[i] if 0 <= i < len(rows) else [F(0)]

    xlo = padd(get(beta_lo, t - 1), get(beta_lo, t))
    xhi = padd(get(beta_hi, t - 1), get(beta_hi, t))
    ylo = padd(get(beta_lo, t - 2), get(beta_lo, t - 1))
    zlo = padd(get(beta_lo, t), get(beta_lo, t + 1))
    zhi = padd(get(beta_hi, t), get(beta_hi, t + 1))
    xprev_lo = padd(get(beta_prev_lo, t - 1), get(beta_prev_lo, t))
    yprev_hi = padd(get(beta_prev_hi, t - 2), get(beta_prev_hi, t - 1))

    p_upper = padd(pmul(xhi, xhi), pmul(ylo, zlo), -1)
    q_upper = padd(
        pmul(get(beta_hi, t - 1), get(beta_hi, t - 1)),
        pmul(get(beta_lo, t - 2), get(beta_lo, t)),
        -1,
    )

    cx_lo = padd(get(gamma_prev_lo, t - 1), get(gamma_prev_lo, t))
    cy_hi = padd(get(gamma_prev_hi, t - 2), get(gamma_prev_hi, t - 1))
    cz_hi = padd(get(gamma_prev_hi, t), get(gamma_prev_hi, t + 1))
    g_lower = padd(
        padd(pmul(cx_lo, padd(xlo, xprev_lo)), pmul(cy_hi, zhi), -1),
        pmul(yprev_hi, cz_hi),
        -1,
    )

    anchor_width = tau0.width()
    anchor_factor = F(1) - 4 * anchor_width
    v_index = [
        anchor_factor * m0,
        anchor_factor * 2 * m0,
        anchor_factor * 2 * m0,
    ]
    correction_factor = F(1001, 1000)
    u_lower = padd(
        pmul(v_index, g_lower),
        [correction_factor * c for c in p_upper],
        -1,
    )
    certificate_poly = padd(
        pmul(u_lower, u_lower),
        pmul([F(8)], pmul(p_upper, pmul(q_upper, pmul(v_index, v_index)))),
        -1,
    )

    gates = {
        "pi bracket": PI_LO < PI_HI,
        "anchor width": 0 <= anchor_width < F(1, 10**60),
        "anchor factor": 0 < anchor_factor < 1,
        "continuum lower": all_nonnegative(beta_lo[1:] + gamma_lo[1:]),
        "predecessor lower": all_nonnegative(beta_prev_lo + gamma_prev_lo),
        "P and Q upper": all_nonnegative([p_upper, q_upper]),
        "index polynomial": all_nonnegative([v_index]),
        "numerator lower": all_nonnegative([u_lower]),
        "certificate": all_nonnegative([certificate_poly]),
        "certificate strict": all(c > 0 for c in certificate_poly),
    }

    anchor_numerator = peval(u_lower, F(0)) ** 2
    anchor_denominator = (
        8 * peval(p_upper, F(0)) * peval(q_upper, F(0)) * peval(v_index, F(0)) ** 2
    )
    anchor_ratio = anchor_numerator / anchor_denominator
    return {
        "ok": all(gates.values()),
        "gates": gates,
        "tau0": tau0,
        "anchor_width": anchor_width,
        "anchor_factor": anchor_factor,
        "constants_a": const_a,
        "constants_b": const_b,
        "d_beta": d_beta,
        "d_gamma": d_gamma,
        "tail_degree": tail_degree,
        "certificate_poly": certificate_poly,
        "anchor_ratio": anchor_ratio,
    }


def check_interval_ladder(t: int = T) -> bool:
    kmax = t + 1
    exact = {row[0]: row for row in exact_ladder(30, kmax)}
    for row in interval_ladder(30, kmax):
        m, h, q, tm, b, c = row
        if m not in (1, 10, 30):
            continue
        _, he, qe, te, be, ce = exact[m]
        pairs = [(h, he), (q, qe), (tm, te)]
        pairs += list(zip(b, be)) + list(zip(c, ce))
        if not all(iv.flo() <= value <= iv.fhi() for iv, value in pairs):
            return False
    return True


def fraction_fingerprint(value: F) -> dict[str, str | int]:
    literal = f"{value.numerator}/{value.denominator}"
    return {
        "sha256": hashlib.sha256(literal.encode()).hexdigest(),
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
    }


def write_summary(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as sink:
        json.dump(payload, sink, indent=2, sort_keys=True)
        sink.write("\n")
        sink.flush()
        os.fsync(sink.fileno())
    os.replace(temporary, path)


def main() -> int:
    # Higher fixed columns produce exact certificate coefficients with tens of
    # thousands of decimal digits.  This explicit local limit affects only
    # deterministic serialization after the rational gates have run.
    sys.set_int_max_str_digits(1_000_000)
    parser = argparse.ArgumentParser()
    parser.add_argument("t", type=int, nargs="?", default=T)
    parser.add_argument("m0", type=int, nargs="?", default=M0)
    parser.add_argument("--controls", action="store_true")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    # Preserve the historical no-argument replay, including its two controls;
    # higher-column audits can omit those extra full certificate builds.
    controls = args.controls or len(sys.argv) == 1
    t, m0 = args.t, args.m0
    if t < 2 or m0 < 2:
        raise SystemExit("require t>=2 and M0>=2")

    print(f"A1 independent fixed-t={t} analytic tail certificate")
    print(f"precision: {PREC} bits; anchor m0={m0}; R0={m0 + 1}")
    print(f"pi enclosure width: {float(PI_HI - PI_LO):.3e}")
    ladder_ok = check_interval_ladder(t)
    print(f"interval ladder brackets exact checkpoints: {ladder_ok}")

    result = build_certificate(m0=m0, t=t)
    for name, passed in result["gates"].items():
        print(f"{name:24s} {'PASS' if passed else 'FAIL'}")
    print(f"tau_m0 interval width: {float(result['anchor_width']):.3e}")
    print(f"tail log-moment max degree: {result['tail_degree']}")
    print("injection budgets beta:", " ".join(f"{float(x):.3e}" for x in result["d_beta"][1:]))
    print("injection budgets gamma:", " ".join(f"{float(x):.3e}" for x in result["d_gamma"][1:]))
    print(f"certificate polynomial degree: {len(result['certificate_poly']) - 1}")
    print(f"anchor certificate ratio: {float(result['anchor_ratio']):.12f}")

    controls_ok = True
    if controls:
        inflated = build_certificate(m0=m0, t=t, injection_scale=2000)
        negative_control = not inflated["gates"]["certificate"]
        print(f"negative control, injection budgets x2000 rejected: {negative_control}")

        prior = build_certificate(m0=m0 - 1, t=t)
        prior_rejected = not prior["gates"]["certificate"]
        print(f"boundary control, anchor m0={m0 - 1} rejected: {prior_rejected}")
        controls_ok = negative_control and prior_rejected

    ok = ladder_ok and result["ok"] and controls_ok
    coefficient_literal = "\n".join(
        f"{coefficient.numerator}/{coefficient.denominator}"
        for coefficient in result["certificate_poly"]
    )
    payload = {
        "schema": "square-tail-fixed-column-v1",
        "version": 1,
        "passed": ok,
        "status": "PASS" if ok else "FAIL",
        "claim": f"rho(r,{t})>=1 for every integer r>={m0 + 1}",
        "configuration": {
            "t": t,
            "anchor_m0": m0,
            "interval_bits": PREC,
            "workers": 1,
            "controls": controls,
        },
        "gates": {**result["gates"], "interval_ladder": ladder_ok,
                  "optional_controls": controls_ok},
        "tail_log_moment_max_degree": result["tail_degree"],
        "certificate_polynomial": {
            "degree": len(result["certificate_poly"]) - 1,
            "coefficient_count": len(result["certificate_poly"]),
            "sha256": hashlib.sha256(coefficient_literal.encode()).hexdigest(),
        },
        "anchor_ratio": fraction_fingerprint(result["anchor_ratio"]),
        "anchor_width": fraction_fingerprint(result["anchor_width"]),
    }
    if args.json is not None:
        write_summary(args.json, payload)
        print(f"summary: {args.json}")
    if ok:
        print("ALL EXACT GATES PASS")
        print(f"THEOREM CERTIFIED: rho(r,{t}) >= 1 for every integer r >= {m0 + 1}")
        return 0
    print("GATE FAILURE")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
