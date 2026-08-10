#!/usr/bin/env python3
"""Sharpened corridor core for t=9 and t=10 on their joint domains.

Sharpened chain (A14 section 1 with three kept terms):
  - C1/V = (2t-1)R - 1 + x^2(1-W)/V   (exact identity, term kept)
  - L1_lo kept exactly: the (CORR) bracket keeps the factor (1+57/(20a))
    in the squared object (the C2-term at r >= 68134 is retained as the
    explicit deduction (t-1)/(r sqrt(a)) with r at the region edge)
  - corridor constant 987/200 > pi^2/2 (A33 Delta-theorem) instead of 5
  - window slack delta on the two local decrements d_t, d_{t+1} and
    row-sum slack Sigma on a - Y:
        x0 = (Y-delta)/G, u0 = alpha(1-delta/Y), a <= Y + Sigma.

Certified statement per (t, delta, Sigma, target):
  For every real Y >= (t-2)^2 - Sigma:
    [(2t-1)R0 - 1 + tau0]^2 * (1 + 57/(20(Y+Sigma)))
        >= target * 2 (Y+Sigma) R0,
  where R0 = 1 + x0 + (x0 u0 + x0^2)(1-beta)/(1-u0),
        tau0 = x0^2 (1-beta)/(1-u0).

Runs:
  t=10 UNCONDITIONAL: delta = Sigma-per-step = 987/200 (banked corridor),
       Sigma = (987/200)(t-2). Target 1077/800 = 1.34625 (covers 4/3
       divided by the region-edge scalars with margin).
  t=9 CONDITIONAL on (D-LAW) d_k <= (6/5) sqrt(a)/k for 4<=k<=10:
       delta = (6/5)sqrt(a)/9 at the dip -- implemented as exact
       rational majorants on a Y-grid decomposition, since delta then
       depends on a. Here: piecewise certificate over Y-intervals with
       interval-constant rational delta/Sigma majorants.
"""

import json
import sys
from fractions import Fraction as F

sys.set_int_max_str_digits(1_000_000)

CORRIDOR = F(987, 200)  # > pi^2/2, banked A33 Delta-theorem


def core_poly(t, delta, sigma, target):
    """Cleared certificate polynomial in Y (list of Fractions).

    Positivity for Y in the stated range proves the sharpened core
    inequality at `target` with window slack `delta`, row-sum slack
    `sigma`.
    """
    al = F((t - 1) * (2 * t - 1), t * (2 * t + 1))
    be = F(t * (2 * t + 1), (t + 1) * (2 * t + 3))
    G = 2 * t * (2 * t + 1)

    # rational function pieces over the common denominator DEN(Y):
    #   x0 = (Y - delta)/G
    #   u0 = al (Y - delta)/Y      1 - u0 = ((1-al)Y + al*delta)/Y
    #   E := (1-be)/(1-u0) = (1-be) Y / W,  W = (1-al)Y + al*delta
    #   R0 = 1 + x0 + x0(u0 + x0) E
    #   tau0 = x0^2 E
    # Multiply through by G^2 * W:
    #   R0N = G^2 W + G(Y-delta) W + (1-be)(Y-delta)^2 (al G + Y) ... note
    #   u0 + x0 = (Y-delta)(al G + Y)/(G Y) and x0(u0+x0)E =
    #   (Y-delta)^2 (al G + Y)(1-be) / (G^2 W)
    #   tau0 = (Y-delta)^2 (1-be) Y / (G^2 W)
    one = [F(1)]
    Yp = [F(0), F(1)]
    Ymd = [-delta, F(1)]
    W = [al * delta, 1 - al]

    def pmul(a, b):
        out = [F(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        out[i + j] += x * y
        return out

    def padd(a, b, s=1):
        n = max(len(a), len(b))
        return [
            (a[i] if i < len(a) else F(0))
            + s * (b[i] if i < len(b) else F(0))
            for i in range(n)
        ]

    def pscal(a, c):
        return [c * x for x in a]

    G2W = pscal(W, F(G * G))
    r0n = padd(
        padd(G2W, pscal(pmul(Ymd, W), F(G))),
        pscal(pmul(pmul(Ymd, Ymd), padd(pscal(one, al * G), Yp)), 1 - be),
    )
    taun = pscal(pmul(pmul(Ymd, Ymd), Yp), 1 - be)   # tau0 numerator
    r0d = G2W
    # bracket = (2t-1) R0 - 1 + tau0  ->  ((2t-1) r0n - r0d + taun)/r0d
    brk = padd(padd(pscal(r0n, F(2 * t - 1)), r0d, -1), taun)
    a_of_Y = padd(Yp, pscal(one, sigma))
    # inequality (cleared, x 20a x r0d^2):
    #   brk^2 (20 a + 57) >= target * 40 a^2 r0d ... more precisely:
    #   brk^2/r0d^2 * (1 + 57/(20a)) >= target * 2 a * (r0n/r0d)
    #   -> brk^2 (20a + 57) * ??? ; multiply both sides by 20 a r0d^2:
    #   brk^2 (20a + 57) >= target * 40 a^2 r0n r0d
    lhs = pmul(pmul(brk, brk), padd(pscal(a_of_Y, F(20)), pscal(one, F(57))))
    rhs = pscal(pmul(pmul(a_of_Y, a_of_Y), pmul(r0n, r0d)), 40 * target)
    return padd(lhs, rhs, -1)


def pshift(p, c):
    out = [F(0)] * len(p)
    for coef in reversed(p):
        prev = out[:]
        out = [F(0)] * len(p)
        out[0] = coef
        for i in range(len(prev)):
            out[i] += c * prev[i]
            if i + 1 < len(out):
                out[i + 1] += prev[i]
    return out


def certify_interval(t, delta, sigma, target, ylo, yhi=None):
    """Positivity of the cleared polynomial on [ylo, yhi] (or [ylo,inf))."""
    p = core_poly(t, delta, sigma, target)
    sh = pshift(p, F(ylo))
    if yhi is None:
        if all(c > 0 for c in sh):
            return "coeff-positive"
        return sturm_ok(sh) and "sturm" or None
    # bounded interval: substitute Y = ylo + (yhi-ylo) z/(1+z), z>=0 maps
    # onto [ylo, yhi); simpler: check on [ylo, inf) restricted via shift +
    # Sturm on the compact part is overkill -- use Sturm on the shifted
    # polynomial for no root in [0, yhi-ylo] via two shifts.
    sh2 = pshift(p, F(yhi))
    if all(c > 0 for c in sh) or sturm_ok(sh):
        return "whole-halfline"
    # sign-constant on [ylo,yhi]: no roots in [0, yhi-ylo] for sh and
    # positive at both ends
    val_lo = sh[0]
    val_hi = sh2[0]
    if val_lo > 0 and val_hi > 0 and no_root_in(sh, F(yhi - ylo)):
        return "interval-sturm"
    return None


def sturm_chain(p):
    def dv(a):
        return [i * a[i] for i in range(1, len(a))]

    def prem(a, b):
        a = a[:]
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

    chain = [p[:], dv(p)]
    while chain[-1] and any(chain[-1]):
        r = prem(chain[-2], chain[-1])
        if not r:
            break
        chain.append(r)
    return chain


def _variations(vals):
    s = [1 if v > 0 else -1 for v in vals if v != 0]
    return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])


def sturm_ok(p):
    """No root of p in [0, inf), p(0) > 0."""
    if p[0] <= 0:
        return False
    chain = sturm_chain(p)

    def at0(q):
        for c in q:
            if c:
                return c
        return F(0)

    v0 = _variations([at0(q) for q in chain])
    vi = _variations([q[-1] for q in chain])
    return v0 - vi == 0


def no_root_in(p, hi):
    chain = sturm_chain(p)

    def ev(q, x):
        acc = F(0)
        for c in reversed(q):
            acc = acc * x + c
        return acc

    def at0(q):
        for c in q:
            if c:
                return c
        return F(0)

    v0 = _variations([at0(q) for q in chain])
    vh = _variations([ev(q, hi) for q in chain])
    return v0 - vh == 0


def main():
    results = {}

    # ---- t=10 UNCONDITIONAL ------------------------------------------
    t = 10
    delta = CORRIDOR
    sigma = CORRIDOR * (t - 2)
    target = F(1077, 800)  # 1.34625 > (4/3)*edge-scalar allowance
    ylo = (t - 2) ** 2 - int(sigma) - 1  # conservative floor (a > 64)
    verdict = certify_interval(t, delta, sigma, target, ylo)
    print(f"t=10 unconditional: delta=987/200, Sigma=987/25, "
          f"target=1077/800, Y>={ylo}: {verdict}")
    results["t10"] = {
        "delta": str(delta), "sigma": str(sigma), "target": str(target),
        "y_floor": ylo, "verdict": verdict,
    }

    # ---- t=9 CONDITIONAL on (D-LAW) ----------------------------------
    # (D-LAW): d_k <= (6/5) sqrt(a)/k for 4 <= k <= 10, plus banked
    # d_k < 987/200 for k <= 3. Piecewise in a: on a-interval [alo, ahi]:
    #   delta = (6/5) sqrt(ahi)/9      (window decrements k=9,10)
    #   Sigma = 3*(987/200) + (6/5) sqrt(ahi) (1/4+...+1/8)
    # rational sqrt majorant: isqrt+1.
    t = 9
    target9 = F(267, 200)  # 1.335 > needed 1.333431 at the region edge
    H48 = F(1, 4) + F(1, 5) + F(1, 6) + F(1, 7) + F(1, 8)
    corr_sigma9 = CORRIDOR * 7  # A14 counting, t-2 = 7 decrements
    pieces = []
    import math
    all_ok = True

    # Piece L: a in [49, A_LO] with the UNCONDITIONAL corridor only.
    # Piece H: a in [A_HI, inf) likewise (halfline).
    # Conditional (D-LAW) slices only on the dip [A_LO, A_HI].
    A_LO, A_HI = 49, 768
    vH = certify_interval(t, CORRIDOR, corr_sigma9, target9,
                          A_HI - int(corr_sigma9) - 1)
    pieces.append({"a": [A_HI, "inf"], "mode": "unconditional-corridor",
                   "verdict": vH})
    all_ok &= bool(vH)
    print(f"t=9 | a >= {A_HI} unconditional halfline: {vH}")

    # Conditional dip slices, both window-index readings (kwin=9 primary,
    # kwin=8 the conservative convention-proof variant), 16 slices.
    for kwin in (9, 8):
        ok_k = True
        sl = []
        edges = [A_LO + (A_HI - A_LO) * i // 24 for i in range(25)]
        for alo, ahi in zip(edges, edges[1:]):
            rt = F(math.isqrt(ahi * 4) + 1, 2)   # sqrt(ahi) <= rt, half-step
            delta9 = F(6, 5) * rt / kwin
            sigma9 = 3 * CORRIDOR + F(6, 5) * rt * H48
            ylo9 = max(10, alo - int(sigma9) - 1)
            v = certify_interval(t, delta9, sigma9, target9, ylo9, ahi)
            sl.append({"a": [alo, ahi], "delta": str(delta9),
                       "sigma": str(sigma9), "verdict": v})
            if not v:
                ok_k = False
        print(f"t=9 | dip [{A_LO},{A_HI}] conditional, kwin={kwin}: "
              f"{'ALL PASS' if ok_k else 'FAIL in ' + str([s['a'] for s in sl if not s['verdict']])}")
        pieces.append({"a": [A_LO, A_HI], "mode": f"conditional-kwin{kwin}",
                       "slices": sl, "all_ok": ok_k})
        if kwin == 9:
            all_ok &= ok_k

    results["t9_conditional"] = {
        "d_law": ("(k-1)(X_(k-1)-X_k) <= (6/5) sqrt(a) for k=5..11 "
                  "in the deposited A12 indexing; the equivalent shifted "
                  "/k convention is consumed here; d_k < 987/200 for the "
                  "first three decrements"),
        "target": str(target9),
        "needed_at_region_edge": "1.333431 (r >= 68134)",
        "pieces": pieces,
        "all_ok_primary": all_ok,
    }

    with open(__file__.replace(".py", ".json"), "w") as f:
        json.dump(results, f, indent=1)
    print("wrote", __file__.replace(".py", ".json"))


if __name__ == "__main__":
    main()
