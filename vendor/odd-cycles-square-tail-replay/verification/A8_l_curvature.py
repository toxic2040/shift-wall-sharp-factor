#!/usr/bin/env python3
"""A8 exact gates for the adjacent-row L-curvature reduction.

G1 and G3--G8 check algebraic statements in exact arithmetic. G2 is a finite
replay of the manuscript's all-row profile-curvature theorem; that proof does
not depend on the enumeration bound here.

No project modules are imported.  The recurrence is rebuilt from the A7
(E,D) definitions using fractions.Fraction throughout.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib


def ladder(m_max: int, k_max: int):
    """Return the exact E,D,R arrays through row m_max and degree k_max."""
    e = [[F(0) for _ in range(k_max + 1)] for _ in range(m_max + 1)]
    d = [[F(0) for _ in range(k_max + 1)] for _ in range(m_max + 1)]
    rseq = [F(1) for _ in range(m_max + 1)]
    e[0][0] = F(1)
    if k_max:
        d[0][1] = F(1)
    for m in range(1, m_max + 1):
        rseq[m] = 1 + F(2 * m, 2 * m - 1) ** 2 * rseq[m - 1]
        inv_r = 1 / rseq[m]
        v_m = rseq[m] / F((2 * m + 1) ** 2)
        for k in range(k_max + 1):
            e[m][k] = e[m - 1][k] + d[m - 1][k] * inv_r
        for k in range(k_max + 1):
            d[m][k] = d[m - 1][k] * (1 - inv_r)
            if k:
                d[m][k] += v_m * e[m][k - 1]
    return e, d, rseq


def l_row(e, m: int):
    """L_k(r), k=0..m, with r=m+1 and L_0=0."""
    r = m + 1
    return [F(0)] + [F(r) * (e[m][k] - e[m - 1][k]) / e[m][k]
                      for k in range(1, m + 1)]


def gate_exact_decomposition():
    """Gate Q, Theta, H and the exact finite-r rho reduction."""
    e, _, rseq = ladder(40, 40)
    cells = 0
    for m in range(4, 41):
        r = m + 1
        ll = l_row(e, m)
        for t in range(3, m):
            x = e[m][t] / e[m][t - 1]
            y = e[m][t - 2] / e[m][t - 1]
            z = e[m][t + 1] / e[m][t]
            v = 1 - x * y
            s = v + x * (1 - y * z) + x * (x - z)
            p = 2 * v + x * (1 - y * z)
            aa = 2 * s - p
            bb = x * y * (1 + z)
            cc = x * z * (1 + y)
            d1 = ll[t] - ll[t - 1]
            q0 = ll[t] - 2 * ll[t - 1] + ll[t - 2]
            q1 = ll[t + 1] - 2 * ll[t] + ll[t - 1]

            coef = {
                0: x * (2 * x + 1 - y),
                -1: 2 + x - x * z,
                1: -x * z * (1 + y),
                -2: -x * y * (1 + z),
            }
            q_direct = sum(coef[o] * ll[t + o] for o in coef)
            q_form_1 = 2 * s * ll[t - 1] + aa * d1 - bb * q0 - cc * q1
            q_form_2 = (2 * s * ll[t - 2] + (2 * s + aa) * d1
                        - (2 * s + bb) * q0 - cc * q1)
            assert q_direct == q_form_1 == q_form_2

            # Abelization about k*L1.  Full ULC (checked below) forces every
            # prefix of these weights, although xy<1 alone does not.
            qall = [ll[k + 2] - 2 * ll[k + 1] + ll[k] for k in range(t)]
            weights = [2 * s * (t - 2 - k) + aa for k in range(t - 2)]
            weights.extend([aa - bb, -cc])
            assert len(weights) == len(qall) == t
            assert q_direct == (2 * s * t - p) * ll[1] + sum(
                weights[k] * qall[k] for k in range(t)
            )

            # Newton/ULC coefficient budget Theta >= 1/(t-1).
            u = x * y
            w = z / x
            num = 2 * s + aa
            den = 2 * s + bb + cc
            assert num == 2 * (1 - u) + 3 * x * (1 - u * w) + 4 * x * x * (1 - w)
            assert den == 2 - u + 2 * x + x * x * (2 - w)
            assert u <= F(t - 1, t)
            assert w <= F(t, t + 1)
            c0 = 2 * t - 4 - (2 * t - 3) * u
            c1 = 3 * (t - 1) * (1 - u * w) - 2
            c2 = 4 * t - 6 - (4 * t - 5) * w
            assert (t - 1) * num - den == c0 + c1 * x + c2 * x * x
            assert c0 >= F(t - 3, t)
            assert c1 >= F(4 * (t - 2), t + 1)
            assert c2 >= F(3 * (t - 2), t + 1)
            assert num * (t - 1) >= den
            prefixes = []
            running = F(0)
            for value in weights:
                running += value
                prefixes.append(running)
            assert all(value >= 0 for value in prefixes)
            total_formula = ((t - 1) * (2 * s + aa) - (2 * s + bb + cc)
                             + s * (t - 2) * (t - 3))
            assert prefixes[-1] == total_formula

            # Exact quadratic correction H and its sign-exposed expansion.
            ell = ll[t - 1]
            d0 = x * (x + y * (1 + 2 * z))
            h_direct = ((x * ll[t] + ll[t - 1]) ** 2
                        - (ll[t - 1] + y * ll[t - 2])
                        * x * (z * ll[t + 1] + ll[t]))
            h_expand = (s * ell * ell + aa * ell * d1 + d0 * d1 * d1
                        - x * y * ((1 + z) * ell + (1 + 2 * z) * d1) * q0
                        - x * z * ((1 + y) * ell - y * d1) * q1
                        - x * y * z * q0 * q1)
            assert h_direct == h_expand

            mu = [1 - value / r for value in ll]
            sigma = ((x * mu[t] + mu[t - 1]) ** 2
                     - (mu[t - 1] + y * mu[t - 2])
                     * (x * z * mu[t + 1] + x * mu[t]))
            assert s - sigma == q_direct / r - h_direct / F(r * r)

            alpha = F(1, 4 * r * (r - 1))
            beta = (1 + alpha) / rseq[m] - alpha
            eps = beta * (2 - beta)
            rho_a7 = F(r * r) * (1 - (1 - eps) * sigma / s) ** 2 * s / (2 * v)
            rho_a8 = (r * eps * s + (1 - eps) * (q_direct - h_direct / r)) ** 2 / (2 * s * v)
            assert rho_a7 == rho_a8
            cells += 1
    print(f"[PASS] G1 exact Q/Theta/H/rho identities: {cells} Fraction cells")


def gate_l3_theorem_replay():
    """Exact finite regression check of the proved L3 theorem."""
    e, _, _ = ladder(120, 120)
    q_count = 0
    for m in range(2, 121):
        ll = l_row(e, m)
        qq = [ll[k + 2] - 2 * ll[k + 1] + ll[k] for k in range(m - 1)]
        assert all(q >= 0 for q in qq)
        assert all(qq[k + 1] <= qq[k] for k in range(len(qq) - 1))
        assert all(ll[k] >= k * ll[1] for k in range(1, m + 1))
        q_count += len(qq)
    print(f"[THEOREM-BACKED PASS] G2 L3 through m=120: "
          f"{q_count} exact curvatures; global proof in manuscript")


def gate_ar1_finite_seeds():
    """The only computer-assisted signs in the global AR-1 proof."""
    h = tsum = rcur = d1 = psum = F(1)
    e1 = F(0)
    min_defect = None
    c1000 = None
    for m in range(1, 1002):
        h *= F(2 * m, 2 * m - 1) ** 2
        tsum += 1 / h
        rcur = h * tsum
        c = d1 / rcur
        e1 += c
        defect = 6 * F((m + 1) ** 2) * c * c - e1
        assert defect > F(19, 40)
        if min_defect is None or defect < min_defect[0]:
            min_defect = (defect, m)
        d1 = d1 * (1 - 1 / rcur) + rcur / F((2 * m + 1) ** 2)
        psum += h * tsum * tsum / F((2 * m + 1) ** 2)
        assert psum == tsum * d1
        if m == 1000:
            g = h / F(2 * m + 1)
            c1000 = psum - g * g * tsum ** 3 / 3
            assert c1000 < -F(1, 5)

    def digest(value: F):
        raw = f"{value.numerator}/{value.denominator}".encode()
        return hashlib.sha256(raw).hexdigest()

    c_gap = c1000 + F(1, 5)
    h_gap = min_defect[0] - F(19, 40)
    assert digest(c_gap) == "5efc5131d94c95041aabe27ad8cedf222e15a48feab1dbf9c8818d743c58f8ca"
    assert digest(h_gap) == "78da37f6a33b4c9a0cc26ab354cc52e7ad46b8962c5a343d73a3d2af85612cf3"
    print("[PASS] G3 AR-1 finite seeds: C_1000<-1/5; "
          f"min H_m at m={min_defect[1]}; canonical hashes match")


def gate_ar1_cleared_algebra():
    """Exact checks of the affine-envelope identity and Z expansion."""
    points = [
        (2, F(11, 2), F(7)),
        (7, F(137, 5), F(49, 2)),
        (19, F(1001, 9), F(133, 2)),
        (1001, F(1234567, 97), F(7007, 2)),
    ]
    for m, rr, h in points:
        dc = (2 * m + 1) ** 2 + 4 * (m + 1) ** 2 * rr
        nn = (2 * m + 3) * (2 * m + 1) ** 2 + (3 * m + 2) * rr
        alpha = F((m + 2) * (rr - 1) * (2 * m + 1) ** 2, (m + 1) * dc)
        beta = F((m + 2) * rr, dc)
        bound = F(m + 1, rr) * (
            F((2 * m - 1) ** 2 * (rr - 1) ** 2, 48 * m ** 4)
            - F(h, 5 * (rr - 1))
        )
        delta = beta - (1 - alpha) * bound
        xx = (48 * m ** 4 * (m + 2) * rr ** 2 - 4 * m ** 3 * rr * dc
              - nn * (2 * m - 1) ** 2 * (rr - 1) ** 2)
        lhs = 240 * m ** 4 * rr * (rr - 1) * dc * (delta - F(1, 12 * m))
        rhs = 5 * (rr - 1) * xx + 48 * m ** 4 * nn * h
        assert lhs == rhs

        c4 = -60*m**3 + 20*m**2 + 25*m - 10
        c3 = 80*m**4 + 180*m**3 + 60*m**2 - 85*m + 15
        c2 = 240*m**5 + 320*m**4 - 360*m**3 - 300*m**2 + 105*m + 15
        c1 = 432*m**6 - 112*m**5 - 640*m**4 + 320*m**3 + 340*m**2 - 55*m - 35
        c0 = (1152*m**8 + 2880*m**7 + 2016*m**6 + 592*m**5 + 240*m**4
              - 80*m**3 - 120*m**2 + 10*m + 15)
        zz = c4*rr**4 + c3*rr**3 + c2*rr**2 + c1*rr + c0
        assert zz == 5 * (rr - 1) * xx + 144 * m ** 5 * nn
        assert rhs >= zz  # the sample points obey h >= 3m
    print(f"[PASS] G4 AR-1 cleared algebra: {len(points)} exact rational points")


def poly_mul(left, right):
    out = [F(0)] * (len(left) + len(right) - 1)
    for i, aa in enumerate(left):
        for j, bb in enumerate(right):
            out[i + j] += aa * bb
    return out


def gate_t6_shape_certificate():
    """Positive-coefficient certificate for the deep t=6 scalar limit."""
    s = [F(23, 78), F(5, 1638), F(1, 94640)]
    aa = [F(0), F(5, 1638), F(1, 47320)]
    nn = [11 * s[i] + aa[i] for i in range(3)]
    poly = poly_mul(nn, nn)
    for i, value in enumerate(s):
        poly[i + 1] -= F(23, 39) * value
    expected = [
        F(64009, 6084),
        F(1357, 21294),
        F(7739, 17886960),
        F(33, 8612240),
        F(1, 52998400),
    ]
    assert poly == expected
    assert all(value > 0 for value in poly)
    print("[PASS] G5 t=6 deep-shape certificate: five positive coefficients")


def series_mul(left, right, degree):
    out = [F(0)] * (degree + 1)
    for i, aa in enumerate(left):
        if aa == 0:
            continue
        for j, bb in enumerate(right[:degree + 1 - i]):
            if bb:
                out[i + j] += aa * bb
    return out


def gate_atanh_triangle():
    """Gate the explicit atanh coefficient triangle used in the L3 theorem."""
    m_max = 12
    degree = 2 * m_max + 1
    atanh = [F(0)] * (degree + 1)
    inv_sqrt = [F(0)] * (degree + 1)
    for j in range(m_max + 1):
        atanh[2 * j + 1] = F(1, 2 * j + 1)
        # [x^(2j)] (1-x^2)^(-1/2) = binom(2j,j)/4^j.
        import math
        inv_sqrt[2 * j] = F(math.comb(2 * j, j), 4 ** j)

    e, _, _ = ladder(m_max, m_max)
    powers = {1: atanh}
    for exponent in range(3, 2 * m_max + 2, 2):
        powers[exponent] = series_mul(powers[exponent - 2],
                                     series_mul(atanh, atanh, degree), degree)

    qcoef = [[F(0)] * (m_max + 1) for _ in range(m_max + 1)]
    for m in range(m_max + 1):
        for k in range(m + 1):
            product = series_mul(powers[2 * k + 1], inv_sqrt, degree)
            qcoef[m][k] = product[2 * m + 1] / math.factorial(2 * k + 1)
            assert e[m][k] == qcoef[m][k] / qcoef[m][0]

    # L/r = 1-(s_m/s_(m-1)) h_k, h_k=q_(m-1,k)/q_(m,k).
    for m in range(2, m_max + 1):
        ll = l_row(e, m)
        sm, sp = qcoef[m][0], qcoef[m - 1][0]
        for k in range(1, m):
            hk = qcoef[m - 1][k] / qcoef[m][k]
            assert ll[k] / (m + 1) == 1 - sm * hk / sp
    print(f"[PASS] G6 atanh/Riordan triangle: all coefficients through m={m_max}")


def gate_quadratic_path_model():
    """Gate the exact i^2 hard-core-path representation of L3."""
    m_max = 12

    # Monic odd-ladder polynomials P_m in ascending powers of u.
    polys = [[1], [5, 1]]
    for m in range(1, m_max):
        bb = 8 * m * m + 12 * m + 5
        cc = (2 * m * (2 * m + 1)) ** 2
        cur, prev = polys[m], polys[m - 1]
        nxt = [0] * (len(cur) + 1)
        for k, value in enumerate(cur):
            nxt[k] += bb * value
            nxt[k + 1] += value
        for k, value in enumerate(prev):
            nxt[k] -= cc * value
        polys.append(nxt)

    # Z[M][p]: p nonadjacent sites in {1,...,M}, activity i^2 at site i.
    ztab = [[0] * (m_max + 1) for _ in range(2 * m_max + 1)]
    ztab[0][0] = 1
    for mm in range(1, 2 * m_max + 1):
        ztab[mm][0] = 1
        for p in range(1, m_max + 1):
            ztab[mm][p] = ztab[mm - 1][p]
            if mm == 1:
                if p == 1:
                    ztab[mm][p] += 1
            else:
                ztab[mm][p] += mm * mm * ztab[mm - 2][p - 1]

    e, _, _ = ladder(m_max, m_max)
    for m in range(1, m_max + 1):
        mm = 2 * m
        for k in range(m + 1):
            assert polys[m][k] == ztab[mm][m - k]
        for k in range(1, m):
            p = m - k
            occupancy = F(mm * mm * ztab[mm - 2][p - 1], ztab[mm][p])
            h_ratio = F((mm + 1) * mm * polys[m - 1][k], polys[m][k])
            assert h_ratio == F(mm + 1, mm) * occupancy
            l_over_r = 1 - F(polys[m][0],
                              polys[m - 1][0] * mm * (mm + 1)) * h_ratio
            ll = l_row(e, m)
            assert ll[k] / (m + 1) == l_over_r
    print(f"[PASS] G7 quadratic hard-core path model: all states through m={m_max}")


def gate_lowering_formula():
    """Gate the exact positive lowering connection for the odd ladder."""
    import math

    m_max = 12
    polys = [[1], [5, 1]]
    for m in range(1, m_max):
        bb = 8 * m * m + 12 * m + 5
        cc = (2 * m * (2 * m + 1)) ** 2
        nxt = [0] * (len(polys[m]) + 1)
        for k, value in enumerate(polys[m]):
            nxt[k] += bb * value
            nxt[k + 1] += value
        for k, value in enumerate(polys[m - 1]):
            nxt[k] -= cc * value
        polys.append(nxt)

    for m in range(1, m_max + 1):
        lhs = [F((m - k) * value) for k, value in enumerate(polys[m])]
        rhs = [F(0)] * (m + 1)
        for j in range(m):
            coeff = F(
                math.factorial(2 * m + 1) * (2 * m + 2 * j + 3),
                2 * math.factorial(2 * j + 1)
                * (2 * m - 2 * j - 1) * (2 * m - 2 * j + 1),
            )
            assert coeff > 0
            for k, value in enumerate(polys[j]):
                rhs[k] += coeff * value
        assert lhs == rhs
    print(f"[PASS] G8 positive lowering formula: all rows through m={m_max}")


def main():
    gate_exact_decomposition()
    gate_l3_theorem_replay()
    gate_ar1_finite_seeds()
    gate_ar1_cleared_algebra()
    gate_t6_shape_certificate()
    gate_atanh_triangle()
    gate_quadratic_path_model()
    gate_lowering_formula()
    print("ALL A8 GATES PASS (G2 REPLAYS THE PROVED L3 THEOREM)")


if __name__ == "__main__":
    main()
