#!/usr/bin/env python3
"""Uniform tail of CORE-27/20: all integer t >= 17 in one certificate.

Same cleared polynomial P_t(Y) as verify_joint_core_t11.py. Substitute
t = S + 17, Y = Q + (t-2)(t-7) = Q + (S+15)(S+10), expand exactly over
the integers, and check that EVERY monomial coefficient in (Q, S) is
nonnegative (with a strictly positive constant term). That proves
P_t(Y) > 0 for all real Q >= 0, S >= 0 — i.e. for every integer t >= 17
and every real Y >= (t-2)(t-7) — and the per-t leg (Sturm, t=11..16)
covers the rest. Together: CORE-27/20 for every t >= 11 on the joint
domain a > (t-2)^2.
"""

import json
import sys

sys.set_int_max_str_digits(1_000_000)

T0 = 17


def mul(a, b):
    out = {}
    for (i1, j1), x in a.items():
        for (i2, j2), y in b.items():
            k = (i1 + i2, j1 + j2)
            out[k] = out.get(k, 0) + x * y
    return {k: v for k, v in out.items() if v}


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def neg(p):
    return {k: -v for k, v in p.items()}


def scal(p, c):
    return {k: c * v for k, v in p.items()}


def main():
    # variables: key = (deg_Y, deg_t) initially
    Y = {(1, 0): 1}
    t = {(0, 1): 1}
    one = {(0, 0): 1}

    def C(c):
        return {(0, 0): c} if c else {}

    A1 = mul(add(t, C(-1)), add(scal(t, 2), C(-1)))       # (t-1)(2t-1)
    A2 = mul(t, add(scal(t, 2), C(1)))                    # t(2t+1)
    B1 = A2
    B2 = mul(add(t, C(1)), add(scal(t, 2), C(3)))         # (t+1)(2t+3)
    G = scal(A2, 2)                                       # 2t(2t+1)
    Wt = add(mul(add(A2, neg(A1)), Y), scal(A1, 5))       # (A2-A1)Y + 5A1
    Ym5 = add(Y, C(-5))

    r0n = add(
        mul(Wt, mul(B2, mul(G, G))),
        mul(mul(Wt, Ym5), mul(B2, G)),
        mul(mul(mul(Ym5, Ym5), add(B2, neg(B1))),
            add(mul(A1, G), mul(A2, Y))),
    )
    r0d = mul(Wt, mul(mul(G, G), B2))
    lin = add(mul(r0n, add(scal(t, 2), C(-1))), neg(r0d))
    a_hi = add(Y, scal(add(t, C(-2)), 5))                 # Y + 5(t-2)
    P = add(scal(mul(lin, lin), 20), neg(scal(mul(a_hi, mul(r0n, r0d)), 54)))

    # substitute t = S + T0: reinterpret key (i, j) via t^j -> (S+T0)^j
    from math import comb

    def subst_t(p):
        out = {}
        for (i, j), v in p.items():
            for k in range(j + 1):
                key = (i, k)  # (deg_Y, deg_S)
                out[key] = out.get(key, 0) + v * comb(j, k) * T0 ** (j - k)
        return {k: v for k, v in out.items() if v}

    P2 = subst_t(P)  # keys now (deg_Y, deg_S)

    # substitute Y = Q + (S+15)(S+10) = Q + S^2 + 25 S + 150
    base = {(0, 2): 1, (0, 1): 25, (0, 0): 150, (1, 0): 1}  # (deg_Q, deg_S)
    powers = [{(0, 0): 1}]
    max_y = max(i for i, _ in P2)
    for _ in range(max_y):
        powers.append(mul(powers[-1], base))
    out = {}
    for (i, j), v in P2.items():
        for (qi, sj), w in powers[i].items():
            key = (qi, sj + j)
            out[key] = out.get(key, 0) + v * w
    out = {k: v for k, v in out.items() if v}

    negs = {k: v for k, v in out.items() if v < 0}
    const = out.get((0, 0), 0)
    n_terms = len(out)
    print(f"monomials: {n_terms}; negative: {len(negs)}; "
          f"constant term positive: {const > 0}")
    if negs:
        worst = sorted(negs.items())[:10]
        print("negative monomials (Q_deg, S_deg):", worst)
    ok = not negs and const > 0
    print(f"UNIFORM TAIL t>={T0}: {'PASS' if ok else 'FAIL'}")

    import hashlib
    digest = hashlib.sha256()
    for k in sorted(out):
        digest.update(f"{k[0]},{k[1]}:{out[k]}\n".encode())
    summary = {
        "schema": "joint-core-27-20-bidegree-v1",
        "claim": (
            "cleared CORE-27/20 polynomial, t=S+17, Y=Q+(S+15)(S+10): "
            "all monomial coefficients in (Q,S) nonnegative, constant "
            "positive => CORE-27/20 for every integer t>=17, "
            "Y >= (t-2)(t-7)"
        ),
        "T0": T0,
        "monomial_count": n_terms,
        "negative_monomials": len(negs),
        "constant_positive": const > 0,
        "coefficient_sha256": digest.hexdigest(),
        "status": "PASS" if ok else "FAIL",
    }
    with open(__file__.replace(".py", ".json"), "w") as f:
        json.dump(summary, f, indent=1)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
