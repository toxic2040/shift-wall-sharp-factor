#!/usr/bin/env python3
"""Exact certificate for the first scaled-ratio decrement Delta_m < 5.

The tail reuses the proved A1 Wallis-clock INJ/ENV envelopes.  Only the final
functional is new.  All interval endpoints and polynomial coefficients are
Fractions; the finite leg is exact integer arithmetic.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1_000_000)
import A1_verify_tail as a1


def envelope(M0: int, K: int = 2):
    """A1 section 3a--3d, returned at levels b_0..b_K."""
    last = None
    for rec in a1.ladder_interval(M0, K=K):
        last = rec
    m, _, _, T0, B0, C0, _, _ = last
    assert m == M0
    half_pi = a1.PI.scal(1, 2)
    tau0 = half_pi * T0
    bI = [half_pi * value for value in B0]
    cI = C0

    powers = [a1.Iv.exact(1)]
    for _ in range(1, 2 * K + 3):
        powers.append(powers[-1] * tau0)

    cA = [a1.ZERO]
    cB = [a1.Iv.exact(1)]
    for k in range(1, K + 1):
        accum = a1.ZERO
        for j in range(k):
            d = k - j
            accum = accum + cA[j] * powers[2 * d - 1].scal(1, math.factorial(2 * d - 1))
            accum = accum + cB[j] * powers[2 * d].scal(1, math.factorial(2 * d))
        Bk = cI[k] - accum
        accum = a1.ZERO
        for j in range(k):
            d = k - j
            accum = accum + cA[j] * powers[2 * d].scal(1, math.factorial(2 * d))
            accum = accum + cB[j] * powers[2 * d + 1].scal(1, math.factorial(2 * d + 1))
        Ak = bI[k] - accum - Bk * tau0
        cA.append(Ak)
        cB.append(Bk)

    Alo, Ahi = [v.flo() for v in cA], [v.fhi() for v in cA]
    Blo, Bhi = [v.flo() for v in cB], [v.fhi() for v in cB]
    Alo[0] = Ahi[0] = F(0)
    Blo[0] = Bhi[0] = F(1)
    t0l = tau0.flo()
    bl, cl = a1.bC_polys(Alo, Blo, K)
    bh, ch = a1.bC_polys(Ahi, Bhi, K)
    BL = [a1.pshift(poly, t0l) for poly in bl]
    BH = [a1.pshift(poly, t0l) for poly in bh]
    CL = [a1.pshift(poly, t0l) for poly in cl]
    CH = [a1.pshift(poly, t0l) for poly in ch]
    BL[0] = BH[0] = [t0l, F(1)]
    CL[0] = CH[0] = [F(1)]

    _, log43_hi = a1.log_bounds(F(4, 3))
    c0const = a1.PI_HI / 2 + log43_hi / 2 + F(3, 8 * M0)
    _, logM_hi = a1.log_bounds(F(M0))

    def J(power: int) -> F:
        return sum(math.perm(power, p) * logM_hi ** (power - p)
                   for p in range(power + 1)) / M0

    def tailsum(poly: list[F]) -> F:
        total = F(0)
        for n, coefficient in enumerate(a1.pplus(poly)):
            for i in range(n + 1):
                total += (coefficient * math.comb(n, i) * c0const ** (n - i)
                          * F(1, 2**i) * J(i))
        return total / 4

    bplus = [a1.pplus(poly) for poly in bh]
    cplus = [a1.pplus(poly) for poly in ch]
    Db = [F(0)] * (K + 1)
    Dc = [F(0)] * (K + 1)
    for k in range(1, K + 1):
        Db[k] = tailsum(bplus[k - 1])
        Dc[k] = tailsum(a1.padd(bplus[k - 1], cplus[k - 1]))

    dmax = F(1, 2 * M0)
    Eb, Ec = [[F(0)]], [[F(0)]]
    for k in range(1, K + 1):
        Eck = (a1.padd([Dc[k]], a1.pint(Eb[k - 1], dmax))
               if k >= 2 else [Dc[1]])
        Ebk = a1.padd([Db[k]], a1.pint(Eck, dmax))
        Ec.append(Eck)
        Eb.append(Ebk)

    lower = [a1.padd(BL[k], Eb[k], -1) for k in range(K + 1)]
    upper = [a1.padd(BH[k], Eb[k], +1) for k in range(K + 1)]
    lower[0], upper[0] = BL[0], BH[0]
    return lower, upper, tau0, Db, Dc


def delta_polynomial(M0: int):
    lower, upper, tau0, Db, Dc = envelope(M0)
    b0 = lower[0]
    positive = a1.padd(
        a1.pmul([F(5)], a1.pmul(b0, lower[1])),
        a1.pmul([F(20)], a1.pmul(b0, lower[2])),
    )
    certified = a1.padd(positive, a1.pmul([F(6)], a1.pmul(upper[1], upper[1])), -1)
    return certified, lower, upper, tau0, Db, Dc


def finite_gate(mmax: int = 1000):
    old = [1, 0, 0]
    row = [5, 1, 0]
    minimum = None
    witness = None
    digest = hashlib.sha256()
    checked = 0
    for m in range(1, mmax + 1):
        if m >= 2:
            p0, p1, p2 = row
            numerator = 5 * p0 * p1 - 6 * p1 * p1 + 20 * p0 * p2
            denominator = p0 * p1
            assert numerator > 0 and denominator > 0
            slack = F(numerator, denominator)
            if minimum is None or slack < minimum:
                minimum, witness = slack, m
            digest.update(f"{m},{numerator}/{denominator}\n".encode())
            checked += 1
        if m == mmax:
            break
        diagonal = 8 * m * m + 12 * m + 5
        back = (2 * m * (2 * m + 1)) ** 2
        new = [(row[k - 1] if k else 0) + diagonal * row[k] - back * old[k]
               for k in range(3)]
        old, row = row, new
    print(f"[PASS] finite Delta_m<5: m=2..{mmax}, cells={checked}, "
          f"minimum slack={float(minimum):.12f} at m={witness}, "
          f"sha256={digest.hexdigest()}")


def tail_gate(M0: int = 1000):
    polynomial, lower, upper, tau0, Db, Dc = delta_polynomial(M0)
    assert all(all(c >= 0 for c in row) for row in lower[1:])
    assert all(c >= 0 for c in polynomial)
    assert polynomial[0] > 0
    payload = "\n".join(f"{i}:{c.numerator}/{c.denominator}"
                        for i, c in enumerate(polynomial))
    digest = hashlib.sha256(payload.encode()).hexdigest()
    print(f"[PASS] tail Delta_m<5 for every m>={M0}: degree={len(polynomial)-1}, "
          f"coefficients={len(polynomial)}, sha256={digest}")
    print(f"       tau0 width={float(tau0.wid()):.3e}, "
          f"Db={[float(v) for v in Db[1:]]}, Dc={[float(v) for v in Dc[1:]]}")

    # Certificate-depth control.  The theorem remains true at m=500, but this
    # particular coefficientwise tail enclosure is not strong enough there.
    shallow, *_ = delta_polynomial(500)
    negative = [(i, c) for i, c in enumerate(shallow) if c < 0]
    assert negative and negative[0][0] == 4
    print("[PASS] depth control: anchor 500 rejected by a negative degree-4 coefficient")

    # Algebra literal: after b_k=(pi/2)B_k and E_k=B_k/T, this is exactly
    # b0*b1*(5-(6E1-20E2/E1)).
    import sympy as sp
    b0, b1, b2 = sp.symbols("b0 b1 b2", positive=True)
    E1, E2 = b1 / b0, b2 / b0
    assert sp.expand(b0 * b1 * (5 - (6 * E1 - 20 * E2 / E1))
                     - (5 * b0 * b1 - 6 * b1**2 + 20 * b0 * b2)) == 0
    print("[PASS] cap functional derived from the literal X_1-X_2 normalization")


def main() -> int:
    finite_gate()
    tail_gate()
    print("ALL A13 DELTA5 GATES PASS")
    print("THEOREM: 0 <= X_1-X_2 < 5 on every canonical row")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
