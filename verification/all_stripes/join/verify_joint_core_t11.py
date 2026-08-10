#!/usr/bin/env python3
"""Exact CORE-27/20 on the joint domain, t >= 11.

Statement being certified (the A14 uniform-corridor chain of Paper 3, see
vendor/odd-cycles-square-tail-replay/verification/A14_uniform_sc.py,
retargeted here to the core constant 27/20):

    For every integer t >= 11 and every real Y >= Ymin(t) = (t-2)(t-7),
    with a = Y + 5(t-2):
        [(2t-1) R0 - 1]^2 / (2 a R0)  >=  27/20,
    where R0 = 1 + x0 + (x0 u0 + x0^2)(1-beta)/(1-u0),
          x0 = (Y-5)/G, u0 = alpha (Y-5)/Y,
          G = 2t(2t+1),
          alpha = (t-1)(2t-1)/(t(2t+1)), beta = t(2t+1)/((t+1)(2t+3)).

Since A14 proves rho >= Core * (1-2/(4r-3))^2 (1-2/r)^2 on the joint
window (r>=504, t>=6, a>=5t, t<sqrt(a)+2), and the scalar product at
r=504 is (504761/507276)^2, CORE-27/20 gives
    rho >= (27/20)(504761/507276)^2 = 1.33667... > 4/3
on every joint cell with t>=11 and a>(t-2)^2 (which forces
Y >= a-5(t-2) > (t-2)(t-7) and a >= 5t+4 > 5t for t>=9).

Integer clearing: with A1=(t-1)(2t-1), A2=t(2t+1)=B1, B2=(t+1)(2t+3),
Wt=(A2-A1)Y+5*A1:
    R0N' = B2*Wt*(G^2 + G(Y-5)) + (B2-B1)(Y-5)^2 (A1*G + A2*Y)
    R0D' = G^2*B2*Wt
    P_t(Y) = 20[(2t-1)R0N' - R0D']^2 - 54(Y+5(t-2)) R0N' R0D'  >= 0.

Per-t: exact integer polynomial in Y, positivity on [Ymin, inf) by
shifted-coefficient signs with Sturm fallback. Uniform tail: bivariate
expansion t = S + T0, Y = Ymin(t) + Q, all monomials >= 0.
"""

import json
import sys
from fractions import Fraction as F

sys.set_int_max_str_digits(1_000_000)


# ---------------- univariate integer polynomial helpers (dense lists) ----


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def padd(a, b, s=1):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + s * (b[i] if i < len(b) else 0)
            for i in range(n)]


def pshift(p, c):
    """p(Y) -> p(Q + c) exactly (c integer)."""
    out = [0] * len(p)
    for coef in reversed(p):
        prev = out[:]
        out = [0] * len(p)
        out[0] = coef
        for i in range(len(prev)):
            out[i] += c * prev[i]
            if i + 1 < len(out):
                out[i + 1] += prev[i]
    return out


def sturm_no_root_ge0(p):
    """True if p has no root in [0, inf) (p as integer list, p(0)!=0)."""
    def dv(a):
        return [i * a[i] for i in range(1, len(a))]

    def prem(a, b):
        a = [F(x) for x in a]
        b = [F(x) for x in b]
        while len(a) >= len(b) and any(a):
            if a[-1] == 0:
                a.pop()
                continue
            q = a[-1] / b[-1]
            d = len(a) - len(b)
            for i in range(len(b)):
                a[d + i] -= q * b[i]
            a.pop()
        while a and a[-1] == 0:
            a.pop()
        return [-x for x in a]

    chain = [[F(x) for x in p], [F(x) for x in dv(p)]]
    while chain[-1] and any(chain[-1]):
        r = prem(chain[-2], chain[-1])
        if not r:
            break
        chain.append(r)

    def sgn_at0(q):
        for c in q:
            if c:
                return 1 if c > 0 else -1
        return 0

    def sgn_inf(q):
        return 1 if q[-1] > 0 else -1

    def variations(signs):
        s = [x for x in signs if x]
        return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])

    v0 = variations([sgn_at0(q) for q in chain])
    vi = variations([sgn_inf(q) for q in chain])
    return v0 - vi == 0


# ---------------- the certificate polynomial ---------------------------


def cert_poly(t):
    A1 = (t - 1) * (2 * t - 1)
    A2 = t * (2 * t + 1)
    B1 = A2
    B2 = (t + 1) * (2 * t + 3)
    G = 2 * t * (2 * t + 1)
    Wt = [5 * A1, A2 - A1]                      # (A2-A1) Y + 5 A1
    Ym5 = [-5, 1]
    r0n = padd(
        pmul(Wt, [B2 * G * G]),
        pmul(pmul(Wt, Ym5), [B2 * G]),
    )
    r0n = padd(
        r0n,
        pmul(pmul(pmul(Ym5, Ym5), [(B2 - B1)]), [A1 * G, A2]),
    )
    r0d = pmul(Wt, [G * G * B2])
    lin = padd(pmul(r0n, [2 * t - 1]), r0d, -1)
    p = padd(
        pmul(pmul(lin, lin), [20]),
        pmul([54], pmul([5 * (t - 2), 1], pmul(r0n, r0d))),
        -1,
    )
    return p


def check_t(t):
    p = cert_poly(t)
    ymin = (t - 2) * (t - 7)
    shifted = pshift(p, ymin)
    all_pos = all(c > 0 for c in shifted)
    verdict = "coeff-positive" if all_pos else None
    if not all_pos:
        if shifted[0] > 0 and shifted[-1] > 0 and sturm_no_root_ge0(shifted):
            verdict = "sturm"
    # float margin scout at a few Y (diagnostic only)
    def corefloat(Y):
        G = 2 * t * (2 * t + 1)
        al = (t - 1) * (2 * t - 1) / (t * (2 * t + 1))
        be = t * (2 * t + 1) / ((t + 1) * (2 * t + 3))
        x0 = (Y - 5) / G
        u0 = al * (Y - 5) / Y
        R0 = 1 + x0 + (x0 * u0 + x0 * x0) * (1 - be) / (1 - u0)
        a = Y + 5 * (t - 2)
        return ((2 * t - 1) * R0 - 1) ** 2 / (2 * a * R0)

    ys = [ymin + k for k in
          [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]]
    mn = min(corefloat(y) for y in ys if y > 5)
    return verdict, float(mn), shifted


def main():
    results = {}
    bad = []
    for t in range(11, 201):
        verdict, mn, _ = check_t(t)
        results[t] = (verdict, mn)
        if verdict is None:
            bad.append(t)
    ok = [t for t, (v, _) in results.items() if v]
    first_coeff = min(
        (t for t, (v, _) in results.items() if v == "coeff-positive"),
        default=None,
    )
    print(f"t=11..200: certified {len(ok)}/190; failures={bad}")
    print(f"first coefficient-positive t: {first_coeff}")
    for t in (11, 12, 13, 15, 20, 50, 100, 200):
        v, mn = results[t]
        print(f"  t={t}: {v}; float core min over probe grid ~ {mn:.4f}")
    summary = {
        "schema": "joint-core-27-20-t11-v1",
        "claim": (
            "[(2t-1)R0-1]^2 >= (27/20) 2 a R0 for every integer 11<=t<=200 "
            "and real Y >= (t-2)(t-7), a=Y+5(t-2)"
        ),
        "per_t_range": [11, 200],
        "verdicts": {str(t): results[t][0] for t in results},
        "failures": bad,
        "first_coeff_positive_t": first_coeff,
        "uniform_tail": "see verify_joint_core_bidegree.py",
        "status": "PASS" if not bad else "FAIL",
    }
    with open(__file__.replace(".py", ".json"), "w") as f:
        json.dump(summary, f, indent=1)
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
