#!/usr/bin/env python3
"""(D-LAW) Leg 2 — single-anchor envelope certificate over the dip window.

Extends Leg 1's certified anchor (m = 262144, 192-bit outward enclosures
of the deposited A1 ladder objects) through the whole t=9 conditional
window a in [49, 768], proving

    (D-LAW)  (k-1)(X_{m,k-1} - X_{m,k}) <= (6/5) sqrt(a_{m+1}),
             k = 5..11 (deposited A12 indexing),

for every row with tau_m >= tau(262144) until a > 768 (beyond which the
sharpened-core certificate is unconditional).

Mechanism (single anchor, global envelopes):
  1. Anchor: Leg 1's 192-bit interval enclosures at m = 262144
     (beta = (pi/2) B; gamma = the RAW C-ladder per the verified
     dictionary), outward-rounded to 96-bit dyadics.
  2. Deposited constants/envelope construction (verbatim algebra from
     A1_verify_tail_independent.build_certificate, anchor-parameterized):
     two-sided polynomial envelopes beta_lo/hi_k(x), x = tau - tau0,
     valid for ALL x >= 0 — the continuum polys (3.2) are the exact
     polynomial closure of the flow beta_k' = gamma_k, gamma_k' =
     beta_{k-1}, so no re-anchoring is needed (the deposited t=2
     certificate covers an infinite tail from one anchor the same
     way).  Injection budgets are computed in the anchor-SHIFTED
     variable at m_lo = 262144 (conservative for every row).
  3. Gates on x in [0, x_end], x_end = first H-multiple with
     a_lo(x_end) > 768 (a is certified increasing, so rows beyond
     x_end are outside the law's window; the anchor row is x = 0):
     for k = 5..11, with L(x) = 2(k-1)(2k-1) bhi_{k-1}^2
                              - 2k(2k+1) blo_k blo_{k-2}:
       Q_k(x) := 216 blo_1 (blo_{k-2} blo_{k-1})^2
                 - 25 (k-1)^2 max(L,0)^2 bhi_0  >= 0
     certified by outward interval subdivision (depth-capped
     bisection), in H-width windows for progress reporting, plus
     lower-envelope nonnegativity per window (squared lower
     envelopes are lower bounds only when nonnegative).

Soundness note (recorded for the writeup): step 2 cites the deposited
INJ/ENV lemma with the Leg-1 anchor; the lemma's algebra is
anchor-agnostic, and the pi bracket / anchor-width hypotheses are
gated.  Envelope-coefficient positivity (a convenience gate of the
rho-certificates) is NOT required here and is not asserted.
"""

import importlib.util
import json
import math
import sys
import time
from fractions import Fraction as F
from pathlib import Path

sys.set_int_max_str_digits(10_000_000)

HERE = Path(__file__).resolve().parent
RELEASE_ROOT = HERE.parents[2]
A1_PATH = (RELEASE_ROOT / "vendor" / "odd-cycles-square-tail-replay"
           / "verification" / "A1_verify_tail_independent.py")
ANCHOR_JSON = HERE / "dlaw_leg1_prefix_anchor.json"
OUT = HERE / "DLAW_LEG2.json"

KMAX = 12
H_STEP = F(7, 20)          # tau step 0.35
A_STOP = 768
M_LO0 = 262_144
DY = 1 << 96               # dyadic rounding grid


def load_a1():
    spec = importlib.util.spec_from_file_location("a1leg2", A1_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def outround(lo: F, hi: F):
    return (F(math.floor(lo * DY), DY), F(-((-hi * DY).__floor__()), DY))


def ieval(poly, xlo: F, xhi: F):
    """Outward interval Horner of poly over [xlo, xhi], x >= 0."""
    lo, hi = F(0), F(0)
    for c in reversed(poly):
        # multiply [lo,hi] by [xlo,xhi] (x >= 0)
        cands = (lo * xlo, lo * xhi, hi * xlo, hi * xhi)
        lo, hi = min(cands), max(cands)
        lo, hi = lo + c, hi + c
    return lo, hi


def certify_nonneg(poly, xlo, xhi, depth=0, max_depth=14):
    lo, _ = ieval(poly, xlo, xhi)
    if lo >= 0:
        return True
    if depth >= max_depth:
        return False
    mid = (xlo + xhi) / 2
    return (certify_nonneg(poly, xlo, mid, depth + 1, max_depth)
            and certify_nonneg(poly, mid, xhi, depth + 1, max_depth))


def build_stage(a1, tau0, betas, gammas, m_lo):
    """Deposited envelope algebra from a synthetic interval anchor.

    tau0 = (lo, hi) Fractions; betas/gammas = per-k (lo, hi) at anchor.
    Returns beta_lo/beta_hi polynomial lists (x = tau - tau0_lo).
    """
    fac = math.factorial
    t0l, t0h = tau0
    powers_lo = [F(1)]
    powers_hi = [F(1)]
    for _ in range(1, 2 * KMAX + 3):
        powers_lo.append(powers_lo[-1] * t0l)
        powers_hi.append(powers_hi[-1] * t0h)

    alo = [F(0)]
    ahi = [F(0)]
    blo_c = [F(1)]
    bhi_c = [F(1)]
    for k in range(1, KMAX + 1):
        pg_lo, pg_hi = F(0), F(0)
        for j in range(k):
            d = k - j
            pg_lo += (alo[j] * powers_lo[2 * d - 1]
                      + blo_c[j] * powers_lo[2 * d]) / fac(2 * d - 1) \
                if False else F(0)
        # exact port of the deposited constant recursion, two-sided:
        pg_lo = sum(
            min(alo[j] * powers_lo[2 * (k - j) - 1],
                alo[j] * powers_hi[2 * (k - j) - 1],
                ahi[j] * powers_lo[2 * (k - j) - 1],
                ahi[j] * powers_hi[2 * (k - j) - 1]) / fac(2 * (k - j) - 1)
            + min(blo_c[j] * powers_lo[2 * (k - j)],
                  blo_c[j] * powers_hi[2 * (k - j)],
                  bhi_c[j] * powers_lo[2 * (k - j)],
                  bhi_c[j] * powers_hi[2 * (k - j)]) / fac(2 * (k - j))
            for j in range(k)
        )
        pg_hi = sum(
            max(alo[j] * powers_lo[2 * (k - j) - 1],
                alo[j] * powers_hi[2 * (k - j) - 1],
                ahi[j] * powers_lo[2 * (k - j) - 1],
                ahi[j] * powers_hi[2 * (k - j) - 1]) / fac(2 * (k - j) - 1)
            + max(blo_c[j] * powers_lo[2 * (k - j)],
                  blo_c[j] * powers_hi[2 * (k - j)],
                  bhi_c[j] * powers_lo[2 * (k - j)],
                  bhi_c[j] * powers_hi[2 * (k - j)]) / fac(2 * (k - j))
            for j in range(k)
        )
        bk_lo = gammas[k][0] - pg_hi
        bk_hi = gammas[k][1] - pg_lo

        pb_lo = sum(
            min(alo[j] * powers_lo[2 * (k - j)],
                alo[j] * powers_hi[2 * (k - j)],
                ahi[j] * powers_lo[2 * (k - j)],
                ahi[j] * powers_hi[2 * (k - j)]) / fac(2 * (k - j))
            + min(blo_c[j] * powers_lo[2 * (k - j) + 1],
                  blo_c[j] * powers_hi[2 * (k - j) + 1],
                  bhi_c[j] * powers_lo[2 * (k - j) + 1],
                  bhi_c[j] * powers_hi[2 * (k - j) + 1]) / fac(2 * (k - j) + 1)
            for j in range(k)
        )
        pb_hi = sum(
            max(alo[j] * powers_lo[2 * (k - j)],
                alo[j] * powers_hi[2 * (k - j)],
                ahi[j] * powers_lo[2 * (k - j)],
                ahi[j] * powers_hi[2 * (k - j)]) / fac(2 * (k - j))
            + max(blo_c[j] * powers_lo[2 * (k - j) + 1],
                  blo_c[j] * powers_hi[2 * (k - j) + 1],
                  bhi_c[j] * powers_lo[2 * (k - j) + 1],
                  bhi_c[j] * powers_hi[2 * (k - j) + 1]) / fac(2 * (k - j) + 1)
            for j in range(k)
        )
        ak_lo = betas[k][0] - pb_hi - max(bk_lo * t0l, bk_lo * t0h,
                                          bk_hi * t0l, bk_hi * t0h)
        ak_hi = betas[k][1] - pb_lo - min(bk_lo * t0l, bk_lo * t0h,
                                          bk_hi * t0l, bk_hi * t0h)
        alo.append(ak_lo)
        ahi.append(ak_hi)
        blo_c.append(bk_lo)
        bhi_c.append(bk_hi)

    beta_l, _ = a1.continuum_polys(alo, blo_c, KMAX)
    beta_h, _ = a1.continuum_polys(ahi, bhi_c, KMAX)
    beta_l = [a1.pshift(p, t0l) for p in beta_l]
    beta_h = [a1.pshift(p, t0l) for p in beta_h]
    beta_l[0] = [t0l, F(1)]
    beta_h[0] = [t0h, F(1)]

    # Injection budgets — deposited (3.3)/(3.5) machinery, instantiated
    # in the ANCHOR-SHIFTED variable.  Two corrections to the first
    # blind port, both against docs/FIXED_COLUMN_TAIL_T2.md:
    #
    # 1. The deposited bound (3.4) is GLOBAL, not windowed:
    #    tau_j <= c0 + (1/2) log j for every j, with
    #    c0 = pi/2 + log(4/3)/2 + 3/(8 m_lo); log(4/3)/2 is the
    #    constant of integration of log((j+3/4)/(3/4)).  The earlier
    #    tau_upper = t0h + 20 double-counted the log growth against
    #    the moments (~1e37 inflation).
    # 2. The deposited abs-majorant in RAW tau powers is sound but
    #    catastrophically loose at high k (cancellation ~2.5e4 at
    #    k=11, tau0 ~ 7.9: abs-poly ~ 400 vs value 0.016), swamping
    #    beta_10/beta_11.  Majorize in the shifted variable instead —
    #    Taylor coefficients about the anchor carry no cancellation:
    #      x_j = tau_j - t0l <= s0 + (1/2) log(j/m_lo),  j > m_lo,
    #      s0 = max(0, c0 + logm_hi/2 - t0l)  (outward; ~0.04 stage 0),
    #    with the exact shifted moments over the tail
    #      sum_{j>m} ((1/2) log(j/m))^i / j^2
    #        <= i!/(2^i m) + 2 (i/8)^i / m^2
    #    (integral of the unimodal integrand plus twice its peak
    #    value (i/4)^i/(m^2 e^i) <= (i/8)^i/m^2).  Coefficientwise
    #    brackets survive pshift (positive weights), so
    #    max(|lo|,|hi|) per shifted coefficient majorizes |P|.
    log43_hi = a1.log_bounds(F(4, 3))[1]
    logm_hi = a1.log_bounds(F(m_lo))[1]
    c0 = a1.PI_HI / 2 + log43_hi / 2 + F(3, 8 * m_lo)
    if m_lo == M_LO0:
        # Exact only at stage 0 (m known there; later m_lo lags the true
        # row).  Catches dictionary/scale errors against deposited (3.4).
        assert t0h <= c0 + logm_hi / 2, "anchor violates deposited (3.4)"
    # round up to the dyadic grid: conservative, bounds bit-growth
    s0 = F(-((-max(F(0), c0 + logm_hi / 2 - t0l) * DY).__floor__()), DY)

    def abs_shifted(lo_polys, hi_polys):
        out = []
        for pl, ph in zip(lo_polys, hi_polys):
            pl = a1.pshift(pl, t0l)
            ph = a1.pshift(ph, t0l)
            n = max(len(pl), len(ph))
            pl = pl + [F(0)] * (n - len(pl))
            ph = ph + [F(0)] * (n - len(ph))
            out.append([max(abs(a), abs(b)) for a, b in zip(pl, ph)])
        return out

    cont_lo_b, cont_lo_g = a1.continuum_polys(alo, blo_c, KMAX)
    cont_hi_b, cont_hi_g = a1.continuum_polys(ahi, bhi_c, KMAX)
    beta_abs = abs_shifted(cont_lo_b, cont_hi_b)
    gamma_abs = abs_shifted(cont_lo_g, cont_hi_g)

    deg_max = max(len(p) for p in beta_abs + gamma_abs)
    moments = [F(math.factorial(i), 2 ** i * m_lo)
               + 2 * F(i, 8) ** i / (F(m_lo) ** 2)
               for i in range(deg_max)]

    def tail_sum(poly):
        total = F(0)
        for n, coeff in enumerate(poly):
            for i in range(n + 1):
                total += (coeff * math.comb(n, i) * s0 ** (n - i)
                          * moments[i])
        return F(-((-total / 4 * DY).__floor__()), DY)  # dyadic, up

    max_step = F(1, 2 * m_lo)
    err_beta = [[F(0)]]
    err_gamma = [[F(0)]]
    for k in range(1, KMAX + 1):
        eg = [tail_sum(a1.padd(beta_abs[k - 1], gamma_abs[k - 1]))]
        if k >= 2:
            eg = a1.padd(eg, a1.pint_discrete_upper(err_beta[k - 1],
                                                    max_step))
        eb = a1.padd([tail_sum(beta_abs[k - 1])],
                     a1.pint_discrete_upper(eg, max_step))
        err_gamma.append(eg)
        err_beta.append(eb)

    beta_lo = [a1.padd(beta_l[k], err_beta[k], -1) for k in range(KMAX + 1)]
    beta_hi = [a1.padd(beta_h[k], err_beta[k]) for k in range(KMAX + 1)]
    beta_lo[0] = [t0l, F(1)]
    beta_hi[0] = [t0h, F(1)]
    return beta_lo, beta_hi


def main():
    # Negative control (--negative-control): retarget the gate from 6/5
    # to 4/5, strictly below the measured saturation (k-1)d_k/sqrt(a)
    # ~ 0.88 — the certificate MUST fail, or the harness proves nothing.
    # Gate constant: target^2 * 6 * 25 -> 216 at 6/5, 96 at 4/5.
    control = "--negative-control" in sys.argv[1:]
    gate_const = 96 if control else 216
    out_path = OUT.with_name("DLAW_LEG2_CONTROL.json") if control else OUT
    t0c = time.time()
    a1 = load_a1()
    # mandatory pi repair (contract), template pattern
    scale = a1.SCALE
    lo = a1.PI_LO.numerator * scale // a1.PI_LO.denominator
    bad_hi = a1.PI_HI.numerator * scale // a1.PI_HI.denominator
    good_hi = -((-a1.PI_HI.numerator * scale) // a1.PI_HI.denominator)
    assert lo == bad_hi and good_hi == bad_hi + 1
    a1.PI = a1.Iv(lo, good_hi)
    half_pi_lo = a1.PI.scal(1, 2).flo()
    half_pi_hi = a1.PI.scal(1, 2).fhi()

    with open(ANCHOR_JSON) as f:
        anc = json.load(f)
    assert anc["status"] == "PASS", "leg1 anchor not certified"
    S = F(1 << anc["scale_bits"])
    b_lo = [F(int(x, 16)) / S for x in anc["anchor_hex"]["b_lo"]]
    b_hi = [F(int(x, 16)) / S for x in anc["anchor_hex"]["b_hi"]]
    c_lo = [F(int(x, 16)) / S for x in anc["anchor_hex"]["c_lo"]]
    c_hi = [F(int(x, 16)) / S for x in anc["anchor_hex"]["c_hi"]]
    # dictionary: beta_k = (pi/2) B_k, gamma_k = (pi/2) C_k * ... —
    # deposited: beta0 = half_pi * b0, and gammas enter only through the
    # constants recursion with the same scaling (verbatim template).
    betas = [outround(half_pi_lo * b_lo[k], half_pi_hi * b_hi[k])
             for k in range(KMAX + 1)]
    # deposited convention: the gamma-system is the RAW C-ladder (only
    # beta carries the pi/2 scale; beta' = gamma matches b' ~ q c).
    gammas = [outround(c_lo[k], c_hi[k]) for k in range(KMAX + 1)]
    tm_lo = F(int(anc["anchor_hex"]["tm"][0], 16)) / S
    tm_hi = F(int(anc["anchor_hex"]["tm"][1], 16)) / S
    tau0 = outround(half_pi_lo * tm_lo, half_pi_hi * tm_hi)
    m_lo = M_LO0

    # SINGLE-ANCHOR GLOBAL CERTIFICATE.  The continuum polys (3.2) are
    # the exact solution of the flow (polynomial closure beta_k' =
    # gamma_k, gamma_k' = beta_{k-1}, degree 2k+1), so one envelope
    # construction at the Leg-1 anchor is valid for every x >= 0 —
    # the deposited t=2 certificate itself covers an infinite tail
    # from one anchor.  Error polys grow like D x^d/d! (~30 at
    # x=19.8, D~1e-6) while values grow like tau^{2k+1}/(2k+1)!
    # (~5.8e10 at tau=27.7): worst relative error is AT the anchor
    # (~0.4% at k=11).  The originally architected 56-stage
    # re-anchoring is unnecessary — and unsound in practice: the
    # anchor re-inversion amplifies interval widths by the
    # partial-sum magnitudes (~1e2-1e3 per index), fatal from
    # stage 1 (measured: stage-1 k=5 margin -1.9e8).
    beta_lo, beta_hi = build_stage(a1, tau0, betas, gammas, m_lo)

    def a_lo_at(x):
        return (6 * ieval(beta_lo[1], x, x)[0]
                / ieval(beta_hi[0], x, x)[1])

    # certification range: first H-multiple with a_lo > A_STOP
    x_end = F(0)
    while a_lo_at(x_end) <= A_STOP:
        x_end += H_STEP
        if x_end > 40:
            raise RuntimeError("a_lo never exceeded A_STOP")

    windows = []
    all_ok = True
    x = F(0)
    win = 0
    while x < x_end and all_ok:
        w_hi = min(x + H_STEP, x_end)
        # Soundness of the cleared gate polynomial Q requires the lower
        # envelopes it squares (indices 1 and k-2..k, k=5..11) to be
        # nonnegative on the gate range: a negative lower envelope would
        # make its square an UPPER, not lower, bound of B^2.
        env_pos = all(certify_nonneg(beta_lo[j], x, w_hi)
                      for j in (1, 3, 4, 5, 6, 7, 8, 9, 10, 11))
        gate_ok = env_pos
        margins = {}
        for k in range(5, 12):
            Lpoly = a1.padd(
                [2 * (k - 1) * (2 * k - 1) * c
                 for c in a1.pmul(beta_hi[k - 1], beta_hi[k - 1])],
                [2 * k * (2 * k + 1) * c
                 for c in a1.pmul(beta_lo[k], beta_lo[k - 2])],
                -1,
            )
            # Q = 216 blo1 (blo_{k-2} blo_{k-1})^2 - 25 (k-1)^2 L^2 bhi0
            prod = a1.pmul(beta_lo[k - 2], beta_lo[k - 1])
            rhs = a1.pmul(beta_lo[1], a1.pmul(prod, prod))
            lhs = a1.pmul([F(25 * (k - 1) ** 2)],
                          a1.pmul(a1.pmul(Lpoly, Lpoly), beta_hi[0]))
            Q = a1.padd([gate_const * c for c in rhs], lhs, -1)
            ok = certify_nonneg(Q, x, w_hi)
            if not ok:
                # L might be negative there (gate trivially true); check
                ok = ieval(Lpoly, x, w_hi)[1] <= 0
            gate_ok &= ok
            qlo = ieval(Q, x, w_hi)[0]
            rlo = ieval(rhs, x, w_hi)[0]
            margins[k] = (bool(ok),
                          float(qlo / (216 * rlo)) if rlo > 0 else 0.0)
        windows.append({
            "window": win, "x": [float(x), float(w_hi)],
            "a_lo_start": float(a_lo_at(x)),
            "envelope_positivity": env_pos,
            "gates": {str(k): margins[k][0] for k in margins},
            "min_rel_margin": min(margins[k][1] for k in margins),
        })
        print(f"window {win}: x=[{float(x):.2f},{float(w_hi):.2f}] "
              f"a_lo={windows[-1]['a_lo_start']:.1f} gates "
              f"{'ALL OK' if gate_ok else 'FAIL ' + str([k for k in margins if not margins[k][0]])} "
              f"min_margin={windows[-1]['min_rel_margin']:.4f} "
              f"[{time.time()-t0c:.0f}s]", flush=True)
        all_ok &= gate_ok
        x = w_hi
        win += 1

    done = bool(windows) and all_ok and a_lo_at(x_end) > A_STOP
    if control:
        status = "CONTROL_REJECTED" if not all_ok else "CONTROL_NOT_REJECTED"
        exit_ok = not all_ok
    else:
        status = "PASS" if done else "FAIL"
        exit_ok = done
    payload = {
        "schema": "dlaw-leg2-single-anchor-v2",
        "status": status,
        "gate_constant": gate_const,
        "law": anc["law"],
        "window": f"tau(262144) through a > {A_STOP}",
        "x_end": float(x_end),
        "m_lo": m_lo,
        "h_step": str(H_STEP),
        "window_count": len(windows),
        "windows": windows,
        "soundness_note": (
            "single-anchor envelopes cite the deposited INJ/ENV algebra "
            "with the Leg-1 anchor; injection budgets are the "
            "deposited tail_sum with the (3.4) global constant "
            "c0 = pi/2 + log(4/3)/2 + 3/(8 m_lo) (infinite-tail moments "
            "from m_lo soundly cover the full range).  Gamma-from-beta-"
            "derivative review DONE: "
            "coefficientwise brackets survive differentiation and "
            "d/dx E^beta_k = E^gamma_k + d_* (E^gamma_k)' >= E^gamma_k "
            "since error polys have nonnegative coefficients.  Lower-"
            "envelope nonnegativity certified per window (env_pos gate) "
            "so squared lower envelopes are valid lower bounds."
        ),
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=1)
    print(f"[{status}] wrote {out_path}")
    return 0 if exit_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
