#!/usr/bin/env python3
"""A2 -- cancellation-free streaming engine for rho(r,t), plus exact helpers.

WHY THIS EXISTS.  engine.py forms  d = S_r - w_r S_{r-1}  as a literal
subtraction.  Those two are equal to relative order 1/(r log r), so in doubles
the subtraction destroys ~log10(r log r) digits: ~9 at r=1e8, ~11 at r=1e9.
This module rewrites d so that every cancellation is O(t) instead of O(r log r).

TWO EXACT REWRITES (both proved in A2_JOINT_REGIME.md section 2):

 (K)  kappa_r = ((2r-1)/(2r-2)) * (1 - 1/R_{r-1})
      [identical to CONTEXT.md (R3) kappa_r = (2r-1)(2r-2)R_{r-2}/(R_{r-1}(2r-3)^2),
       using R_{r-1} = 1 + ((2r-2)/(2r-3))^2 R_{r-2}]
      hence  w_r = ((2r-1)/(2r))^2 kappa_r^2 = ((1+A)(1-b))^2,
             A = 1/(4r(r-1)),  b = 1/R_{r-1},
      and    eps_r := 1 - w_r = -y(2+y),  y = A - b(1+A).      [no cancellation]

 (G)  the successive normalized rows differ by a POSITIVE vector:
      with E_k(m) = f_k(2m), O_k(m) = f_k(2m+1), D_k(m) = O_k(m) - E_k(m),
        R_1 = 5,  R_m = 1 + (2m/(2m-1))^2 R_{m-1},   x_m = R_m/(2m+1)^2
        E_k(m) = E_k(m-1) + D_k(m-1)/R_m
        D_k(m) = D_k(m-1)(1 - 1/R_m) + x_m E_{k-1}(m)
      all terms positive, so E, D, R are computed with no cancellation at all.
      g_k := E_k(r-1) - E_k(r-2) = D_k(r-2)/R_{r-1}  is available exactly.

Then  S_r - S_{r-1} = L(p+q)_t - L(p)_t  with p = (1+u)E(r-2), q = (1+u)g,
expanded bilinearly, and  d = (S_r - S_{r-1}) + eps_r S_{r-1}.

Validated in A2_verify_joint.py against exact_rho from engine.py.
"""
from fractions import Fraction as F
import math

__all__ = ["StableEngine", "exact_row", "exact_a"]


def _shift1(f, zero):
    out = [zero] * (len(f) + 1)
    for i, c in enumerate(f):
        out[i] = out[i] + c
        out[i + 1] = out[i + 1] + c
    return out


class StableEngine:
    """Stream r = 5..RMAX yielding (r, {t: rho}, a_r) with a_r = 6 f_1(2r-2).

    num: float (default), or mpmath.mpf, or fractions.Fraction (exact but slow).
    """

    def __init__(self, KMAX, num=float):
        self.KMAX = KMAX
        self.num = num

    def run(self, RMAX, ts, callback, rmin=5):
        num, K = self.num, self.KMAX
        one, zero = num(1), num(0)
        tmax = max(ts)
        # m = 0 initial data
        E = [one] + [zero] * K                      # E_k(0) = delta_{k0}
        D = [zero, one] + [zero] * (K - 1)          # D_k(0) = delta_{k1}
        R = num(5)                                  # R_1
        Eprev = list(E)                             # E(m-1) for the r = m+1 report
        Dprev = list(D)
        for m in range(1, RMAX):
            Eprev, Dprev = E, D
            Rm = num(5) if m == 1 else one + (num(2 * m) / num(2 * m - 1)) ** 2 * R
            R = Rm
            invR = one / Rm
            xm = Rm / num((2 * m + 1) ** 2)
            E = [Eprev[k] + Dprev[k] * invR for k in range(K + 1)]
            D = [Dprev[k] * (one - invR) + (xm * E[k - 1] if k >= 1 else zero)
                 for k in range(K + 1)]
            # rows for r: Pi_r  <-> E(r-1),  Pi_{r-1} <-> E(r-2).  m = r-1.
            r = m + 1
            if r < rmin or r < 4:
                continue
            g = [Dprev[k] * invR for k in range(K + 1)]      # = E(m) - E(m-1)
            out = self._row(r, Eprev, g, Rm, ts, tmax)
            callback(r, out, num(6) * E[1])

    def _row(self, r, Ep, g, Rm, ts, tmax):
        num = self.num
        one, zero = num(1), num(0)
        p = _shift1(Ep, zero)                # (1+u) Pihat_{r-1}
        q = _shift1(g, zero)                 # (1+u) (Pihat_r - Pihat_{r-1})
        f = [Ep[k] + g[k] for k in range(len(Ep))]   # Pihat_r row
        A = one / num(4 * r * (r - 1))
        b = one / Rm
        y = A - b * (one + A)
        eps = -y * (num(2) + y)              # = 1 - w_r, cancellation-free
        P = lambda v, i: v[i] if 0 <= i < len(v) else zero
        out = {}
        for t in ts:
            if t > r - 2 or t + 1 > self.KMAX:
                continue
            Sm = P(p, t) ** 2 - P(p, t - 1) * P(p, t + 1)              # S_{r-1}
            dS = (2 * P(p, t) * P(q, t) - P(p, t - 1) * P(q, t + 1)
                  - P(q, t - 1) * P(p, t + 1)) \
                 + (P(q, t) ** 2 - P(q, t - 1) * P(q, t + 1))          # S_r - S_{r-1}
            S = Sm + dS
            V = P(f, t - 1) ** 2 - P(f, t - 2) * P(f, t)
            d = dS + eps * Sm
            out[t] = num(r * r) * d * d / (num(2) * S * V)
        return out


# ------------------------------------------------------------ exact utilities

def exact_row(r, K):
    """Exact [f_0..f_K] for Pihat_r = Pi_r normalised to constant term 1."""
    eng = StableEngine(K, num=F)
    got = {}
    eng.run(r + 1, [2], lambda rr, row, a: None, rmin=10 ** 18)
    # simpler: replay the recurrence directly
    E = [F(1)] + [F(0)] * K
    D = [F(0), F(1)] + [F(0)] * (K - 1)
    R = F(5)
    for m in range(1, r):
        Rm = F(5) if m == 1 else 1 + F(2 * m, 2 * m - 1) ** 2 * R
        R = Rm
        xm = Rm / F((2 * m + 1) ** 2)
        Ep, Dp = E, D
        E = [Ep[k] + Dp[k] / Rm for k in range(K + 1)]
        D = [Dp[k] * (1 - 1 / Rm) + (xm * E[k - 1] if k >= 1 else F(0))
             for k in range(K + 1)]
    return E                                   # = f_k(2r-2)


def exact_a(r):
    """a_r := 6 f_1(2r-2), exact Fraction."""
    return 6 * exact_row(r, 1)[1]
