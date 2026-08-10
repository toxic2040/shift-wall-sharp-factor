#!/usr/bin/env python3
"""Shared engines for the (SQ) tail slice.  Standard library + optional mpmath.

Two engines, both for  rho(r,t) = r^2 d^2 / (2 S V):

  exact_rho(r, t)      exact Fraction, builds M_{2r-2} outright.  Ground truth.
  Engine(KMAX)         normalized top-diagonal recurrence, O(1) memory per step,
                       reaches r ~ 1e6+.  Floats by default; pass mp=True for
                       mpmath at the current mp.dps.

The definitions are those in the accompanying manuscript.  The normalized
engine agrees with exact rho to 12 significant figures at r=10,100,129,300;
release decisions use the exact engine, not that floating-point regression.
"""
from fractions import Fraction as F
import math

__all__ = ["wall_M", "Pi_exact", "exact_rho", "Engine", "Lop", "shift1"]


# ---------------------------------------------------------------- exact side

def wall_M(N):
    """M_n(z) coefficient lists, n = 0..N.  M_n = M_{n-1} + n^2 z M_{n-2}."""
    M = [[1], [1, 1]]
    for n in range(2, N + 1):
        a = M[n - 1][:]
        b = M[n - 2]
        row = a + [0] * max(0, len(b) + 1 - len(a))
        for i, c in enumerate(b):
            row[i + 1] += n * n * c
        M.append(row)
    return M


def Lop(f, t):
    """L(f)_t = f_t^2 - f_{t-1} f_{t+1}, zero outside the row."""
    g = lambda i: f[i] if 0 <= i < len(f) else 0 * f[0]
    return g(t) * g(t) - g(t - 1) * g(t + 1)


def shift1(f):
    """(1+u) * f"""
    out = [0 * f[0]] * (len(f) + 1)
    out = list(out)
    for i, c in enumerate(f):
        out[i] = out[i] + c
        out[i + 1] = out[i + 1] + c
    return out


def Pi_exact(r, M=None):
    """Pi_r coefficient list, exact Fractions.  [u^k] = m_{2r-2,r-1-k}/(2r-1)!."""
    if M is None:
        M = wall_M(2 * r)
    m = M[2 * r - 2]
    fa = math.factorial(2 * r - 1)
    return [F(m[r - 1 - k], fa) for k in range(0, r)]


def exact_rho(r, t, M=None):
    """Exact rho(r,t) as a Fraction."""
    if M is None:
        M = wall_M(2 * r)
    P, Pm = Pi_exact(r, M), Pi_exact(r - 1, M)
    S = Lop(shift1(P), t)
    Sm = Lop(shift1(Pm), t)
    V = Lop(P, t - 1)
    d = S - F((2 * r - 1) ** 2, (2 * r) ** 2) * Sm
    return F(r * r) * d * d / (2 * S * V)


# ----------------------------------------------------------- normalized side

class Engine:
    """Streaming normalized top-diagonal recurrence.

        A_k(2m)   = A_k(2m-1)   + (2m)^2   A_k(2m-2)
        A_k(2m+1) = A_{k-1}(2m) + (2m+1)^2 A_k(2m-1)
        f_k = A_k/A_0 :
        f_k(2m)   = s f_k(2m-1) + (1-s) f_k(2m-2),  s = 1/R_m
        f_k(2m+1) = u f_{k-1}(2m) + f_k(2m-1),      u = R_m/(2m+1)^2
        R_m = A_0(2m)/A_0(2m-1) = 1 + (2m/(2m-1))^2 R_{m-1},  R_1 = 5

    step() advances one r and yields (r, rho_row) for the requested t values.
    Keeps R_{r-2}, R_{r-1} and the two rows it needs; O(KMAX) memory.
    """

    def __init__(self, KMAX, mp=False):
        self.KMAX = KMAX
        if mp:
            from mpmath import mpf
            self.one, self.zero, self.num = mpf(1), mpf(0), mpf
        else:
            self.one, self.zero, self.num = 1.0, 0.0, float
        self.KMAX = KMAX

    def run(self, RMAX, ts, callback, rmin=5):
        """Iterate r = rmin..RMAX; call callback(r, {t: rho}) for t in ts, t<=r-2.

        rmin defaults to 5 for compatibility with archived diagnostics.  The
        r=4 row is computable by passing rmin=4.  The release theorem does not
        use this streaming sweep: A38_targeted_core.py rebuilds every interior
        cell on 4<=r<=503 in exact integer arithmetic.
        """
        one, zero, num = self.one, self.zero, self.num
        K = self.KMAX
        f = {0: [one] + [zero] * K, 1: [one] + [zero] * K}
        f[1][1] = one
        Rprev = {}
        for n in range(2, 2 * RMAX + 1):
            if n % 2 == 0:
                m = n // 2
                Rm = num(5) if m == 1 else one + (num(2 * m) / (2 * m - 1)) ** 2 * Rprev[m - 1]
                Rprev[m] = Rm
                s = one / Rm
                f[n] = [s * f[n - 1][k] + (one - s) * f[n - 2][k] for k in range(K + 1)]
            else:
                m = (n - 1) // 2
                u = (Rprev[m] if m >= 1 else one) / num((2 * m + 1) ** 2)
                f[n] = [(u * f[n - 1][k - 1] if k >= 1 else zero) + f[n - 2][k]
                        for k in range(K + 1)]
            if n - 3 in f:
                del f[n - 3]
            if m - 3 in Rprev:
                del Rprev[m - 3]
            if n % 2 == 0:
                r = n // 2 + 1
                if r >= rmin and (2 * r - 4) in f:
                    callback(r, self._row(r, f, Rprev, ts))

    def _row(self, r, f, R, ts):
        one, num = self.one, self.num
        fr, frm = f[2 * r - 2], f[2 * r - 4]
        kap = num((2 * r - 1) * (2 * r - 2)) * R[r - 2] / (R[r - 1] * num((2 * r - 3) ** 2))
        shr, shm = shift1(fr), shift1(frm)
        w = (num(2 * r - 1) / (2 * r)) ** 2 * kap ** 2
        out = {}
        for t in ts:
            if t > r - 2 or t + 1 > self.KMAX:
                continue
            S, V, Sm = Lop(shr, t), Lop(fr, t - 1), Lop(shm, t)
            d = S - w * Sm
            out[t] = num(r * r) * d * d / (2 * S * V)
        return out


if __name__ == "__main__":
    # self-check: normalized engine vs exact rationals
    M = wall_M(700)
    targets = {10: None, 100: None, 129: None, 300: None}
    for r in targets:
        targets[r] = float(exact_rho(r, 2, M))
    got = {}
    Engine(8).run(301, [2], lambda r, row: got.__setitem__(r, row.get(2)))
    print("SELF-CHECK  normalized engine vs exact rho(r,2)")
    ok = True
    for r in sorted(targets):
        rel = abs(got[r] - targets[r]) / targets[r]
        ok &= rel < 1e-11
        print(f"  r={r:4d}  exact={targets[r]:.12f}  engine={got[r]:.12f}  rel={rel:.2e}")
    print("PASS" if ok else "*** FAIL ***")
    raise SystemExit(0 if ok else 1)
