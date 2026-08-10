#!/usr/bin/env python3
"""A1_verify_tail.py -- exact/rigorous verifier for the fixed-t tail bound.

    THEOREM (A1).  rho(r,2) >= 1 for every r >= R_0(2) = 313, and the same
    machinery certifies R_0(3)=116, R_0(4)=250, R_0(5)=700.
    plus  COROLLARY.  rho(r,t) > 1 for every r >= t+2, t = 2..5  (columns closed).

Everything that gates is exact integer or exact rational arithmetic:
  * the wall/ladder identities are Fraction-exact,
  * the analytic tail runs on fixed-point RATIONAL INTERVALS (endpoints are
    integers over 2^PREC; every operation rounds outward), so each printed
    bound is a true rational inequality,
  * pi and log are enclosed by exact rational series with proven remainders.
No float is load-bearing anywhere below; floating-point values are reporting
only.

Proof skeleton, matching the fixed-column argument in the manuscript:

 (I)   IDENTITY.  With m = r-1, h_m = ((2m)!!/(2m-1)!!)^2, q_m = 1/h_m,
       v_m = h_m/(2m+1)^2, T_m = 1 + sum_{j<=m} 1/h_j, and the ladder
           B_k(m) = B_k(m-1) + q_m C_k(m-1),  C_k(m) = C_k(m-1) + v_m B_{k-1}(m)
           B_0(m) = T_m,  C_0 = 1,  B_k(0)=0 (k>=1), C_1(0)=1, C_k(0)=0 (k>=2)
       one has  f_k(2m) = B_k(m)/T_m  and, with
           P(m) = (B_1+B_2)^2 - (B_0+B_1)(B_2+B_3),  Q(m) = B_1^2 - B_0 B_2,
           g_r  = ((2r-1)^2/((2r)(2r-2)))^2,
           rho(r,2) = r^2 (P(m) - g_r P(m-1))^2 / (2 P(m) Q(m)).
 (II)  WALLIS.  pi(m+1/4) < h_m < pi(m+1/2)   (exact monotone form).
 (III) CONTINUUM.  tau := (pi/2) T_m, b_k := (pi/2) B_k(m).  The ladder is
       symplectic Euler, step dtau_m = (pi/2)q_m, for  b_k' = C_k, C_k' = b_{k-1}
       with source damping eta_m = (2 h_m/(pi(2m+1)))^2 in (1-1/(2m+1), 1).
 (IV)  ENVELOPE.  Freezing the six connection constants at m = M0 gives
       explicit two-sided polynomial envelopes for b_k, C_k valid for all
       m >= M0, with injection budgets D^b_k, D^C_k bounded in closed form.
 (V)   POSITIVITY.  rho >= 1 follows from  (G - corr)^2 >= 8 P Q  because
       (r dtau_m)^2 > 1/4.  Substituting the envelopes turns this into a
       single polynomial in s = tau - tau_{M0} >= 0 whose coefficients are
       certified nonnegative.

exit 0 iff every gate passes.
"""
import sys
import math
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from engine import wall_M, Pi_exact, exact_rho, Lop, shift1

PREC = 256                      # fixed-point scale: values are ints / 2^PREC
D = 1 << PREC
FAILS = []
KMAX = 3                        # t = 2 needs B_0..B_3


def gate(name, cond, extra=""):
    if not cond:
        FAILS.append(f"{name} {extra}")
        print(f"    FAIL  {name} {extra}")
    return bool(cond)


# =============================================================== SECTION 0
# fixed-point rational intervals.  value(x) in [x.lo/D, x.hi/D], both exact.

class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    @staticmethod
    def frac(p, q=1):
        n, d = p * D, q
        return Iv(n // d, -((-n) // d))

    @staticmethod
    def exact(p, q=1):
        n, d = p * D, q
        assert n % d == 0, "not exact in fixed point"
        return Iv(n // d, n // d)

    def __add__(s, o):
        return Iv(s.lo + o.lo, s.hi + o.hi)

    def __sub__(s, o):
        return Iv(s.lo - o.hi, s.hi - o.lo)

    def __neg__(s):
        return Iv(-s.hi, -s.lo)

    def __mul__(s, o):
        a = (s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi)
        lo, hi = min(a), max(a)
        return Iv(lo // D, -((-hi) // D))

    def inv(s):
        assert s.lo > 0, "reciprocal of interval containing 0"
        return Iv((D * D) // s.hi, -((-(D * D)) // s.lo))

    def __truediv__(s, o):
        return s * o.inv()

    def scal(s, p, q=1):
        """multiply by the exact rational p/q  (q > 0)."""
        a, b = s.lo * p, s.hi * p
        if p < 0:
            a, b = b, a
        return Iv(a // q, -((-b) // q))

    def flo(s):
        return F(s.lo, D)

    def fhi(s):
        return F(s.hi, D)

    def __repr__(s):
        return f"[{float(s.lo)/D:.12f},{float(s.hi)/D:.12f}]"

    def wid(s):
        return F(s.hi - s.lo, D)


ZERO = Iv(0, 0)
ONE = Iv.exact(1)


# ------------------------------------------------ exact rational pi and log
def pi_bounds(K=60):
    """Machin: pi = 16 atan(1/5) - 4 atan(1/239); alternating -> brackets."""
    def atan_inv(n, K):
        s_lo = F(0)
        s_hi = F(0)
        term = F(1, n)
        acc = F(0)
        for k in range(K):
            t = F(1, (2 * k + 1) * n ** (2 * k + 1))
            acc += t if k % 2 == 0 else -t
            nxt = F(1, (2 * k + 3) * n ** (2 * k + 3))
            if k % 2 == 0:
                s_hi, s_lo = acc, acc - nxt
            else:
                s_lo, s_hi = acc, acc + nxt
        return s_lo, s_hi
    a5l, a5h = atan_inv(5, K)
    a2l, a2h = atan_inv(239, K)
    return 16 * a5l - 4 * a2h, 16 * a5h - 4 * a2l


PI_LO, PI_HI = pi_bounds()
PI = Iv(int(PI_LO * D), -int(-PI_HI * D))


def log_bounds(x):
    """rational two-sided bounds for log(x), x = Fraction > 0."""
    assert x > 0
    e = 0
    while x >= 2:
        x /= 2
        e += 1
    while x < 1:
        x *= 2
        e -= 1
    # log 2 = 2 atanh(1/3)
    def atanh_bounds(z, K=200):
        s = F(0)
        zz = z
        for k in range(K):
            s += zz / (2 * k + 1)
            zz *= z * z
        rem = zz / ((2 * K + 1) * (1 - z * z))
        return s, s + rem
    l2l, l2h = atanh_bounds(F(1, 3))
    l2l, l2h = 2 * l2l, 2 * l2h
    z = (x - 1) / (x + 1)
    al, ah = atanh_bounds(z)
    lo = 2 * al + (e * l2l if e >= 0 else e * l2h)
    hi = 2 * ah + (e * l2h if e >= 0 else e * l2l)
    return lo, hi


# ----------------------------------------------- polynomials with Fraction c
def pmul(a, b):
    o = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                o[i + j] += x * y
    return o


def padd(a, b, s=1):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else F(0)) + s * (b[i] if i < len(b) else F(0))
            for i in range(n)]


def pev(a, x):
    return sum(c * x ** i for i, c in enumerate(a))


def pshift(a, t):
    """coefficients of a(t+x) in x, exact."""
    o = [F(0)] * len(a)
    for i, c in enumerate(a):
        if not c:
            continue
        for j in range(i + 1):
            o[j] += c * math.comb(i, j) * t ** (i - j)
    return o


def pplus(a):
    """dominating polynomial with nonneg coefficients"""
    return [c if c > 0 else F(0) for c in a]


def pint(a, dmax):
    """upper bound for sum_j dtau_j a(s_j) : integral + dmax * a(s).
    a must have nonneg coefficients."""
    out = [F(0)] + [c / (i + 1) for i, c in enumerate(a)]
    return padd(out, [dmax * c for c in a])


# =============================================================== SECTION 1
def section1_identities():
    print("[1] exact identities, Wallis sandwich, pinning control")
    MM = 45
    Mw = wall_M(2 * MM + 4)
    # ladder in exact Fractions (small m only)
    W = F(1)
    h, T, B, C = [F(1)], [F(1)], [[F(1), F(0), F(0), F(0)]], [[F(1), F(1), F(0), F(0)]]
    for m in range(1, MM + 2):
        W *= F(2 * m, 2 * m - 1)
        hm = W * W
        q, v = 1 / hm, hm / F((2 * m + 1) ** 2)
        h.append(hm)
        T.append(T[-1] + q)
        Bm = [T[-1]] + [B[-1][k] + q * C[-1][k] for k in (1, 2, 3)]
        Cm = [F(1)] + [C[-1][k] + v * Bm[k - 1] for k in (1, 2, 3)]
        B.append(Bm)
        C.append(Cm)

    def PQ(m):
        t, b1, b2, b3 = B[m][0], B[m][1], B[m][2], B[m][3]
        return (b1 + b2) ** 2 - (t + b1) * (b2 + b3), b1 * b1 - t * b2

    ok = True
    for r in range(4, MM + 1):
        m = r - 1
        Pn, Qn = PQ(m)
        Pp, _ = PQ(m - 1)
        g = F((2 * r - 1) ** 2, (2 * r) * (2 * r - 2)) ** 2
        rho = F(r * r) * (Pn - g * Pp) ** 2 / (2 * Pn * Qn)
        ok &= gate("identity", rho == exact_rho(r, 2, Mw), f"r={r}")
        ok &= gate("rho>1", rho > 1, f"r={r} rho={float(rho):.9f}")
        # f_k = B_k/T_m
        Pi = Pi_exact(r, Mw)
        for k in range(4):
            ok &= gate("f_k", Pi[k] / Pi[0] == B[m][k] / T[m], f"r={r} k={k}")
    print(f"    identity + f_k + rho>1 on 4<=r<={MM}: {'ok' if ok else 'FAILED'}")

    # NEGATIVE CONTROL: g_r is load-bearing -- replacing it by 1 breaks it
    bad = 0
    for r in range(4, 20):
        m = r - 1
        Pn, Qn = PQ(m)
        Pp, _ = PQ(m - 1)
        if F(r * r) * (Pn - Pp) ** 2 / (2 * Pn * Qn) != exact_rho(r, 2, Mw):
            bad += 1
    gate("NC g_r load-bearing", bad == 16, f"mismatches={bad}/16")

    # Wallis sandwich, exact monotone form
    okw = True
    for m in range(2, 3000):
        okw &= F((2 * m) ** 2 * (4 * m - 3), (2 * m - 1) ** 2 * (4 * m + 1)) < 1
        okw &= F((2 * m) ** 2 * (2 * m - 1), (2 * m - 1) ** 2 * (2 * m + 1)) > 1
    gate("Wallis monotone", okw)
    okw2 = all(h[m] / F(4 * m + 1, 4) > PI_HI and h[m] / F(2 * m + 1, 2) < PI_LO
               for m in range(1, MM))
    gate("Wallis sandwich vs certified pi", okw2)
    print(f"    Wallis  pi(m+1/4) < h_m < pi(m+1/2): ok")

    # pinning negative control (r,t) = (2,1)
    Mw3 = wall_M(6)
    P2, P1 = Pi_exact(2, Mw3), Pi_exact(1, Mw3)
    S, Sm, V = Lop(shift1(P2), 1), Lop(shift1(P1), 1), Lop(P2, 0)
    d = S - F(9, 16) * Sm
    gate("NC (2,1) value", d * d - S * V / 2 == F(-4351, 20736))
    gate("NC (2,1) rho<1", F(4) * d * d / (2 * S * V) < 1)
    print(f"    (2,1) control  d^2-SV/2 = -4351/20736, rho = "
          f"{float(F(4)*d*d/(2*S*V)):.9f} < 1: ok")
    return B, C, T, h


# =============================================================== SECTION 2
def ladder_interval(M, K=3):
    """rigorous interval ladder, levels k = 0..K."""
    h = Iv.exact(1)
    T = Iv.exact(1)
    B = [Iv.exact(1)] + [ZERO] * K
    Cc = [Iv.exact(1), Iv.exact(1)] + [ZERO] * (K - 1)
    for m in range(1, M + 1):
        h = h.scal((2 * m) ** 2, (2 * m - 1) ** 2)
        q = h.inv()
        v = h.scal(1, (2 * m + 1) ** 2)
        Bp, Cp = B[:], Cc[:]
        T = T + q
        B = [T] + [Bp[k] + q * Cp[k] for k in range(1, K + 1)]
        Cc = [Iv.exact(1)] + [Cp[k] + v * B[k - 1] for k in range(1, K + 1)]
        yield m, h, q, T, B, Cc, Bp, Cp


def PQ_iv(B, t):
    g = lambda i: B[i] if 0 <= i < len(B) else ZERO
    P = (g(t - 1) + g(t)) * (g(t - 1) + g(t)) - (g(t - 2) + g(t - 1)) * (g(t) + g(t + 1))
    Q = g(t - 1) * g(t - 1) - g(t - 2) * g(t)
    return P, Q


def section2_finite(RCERT, t=2):
    print(f"[2] finite certificate: rho(r,{t}) > 1 for {t+2} <= r <= {RCERT} "
          f"(rigorous rational intervals, PREC={PREC})")
    Pprev = None
    worst = worstr = None
    ok = True
    n = 0
    for m, h, q, T, B, Cc, Bp, Cp in ladder_interval(RCERT - 1, K=t + 1):
        P, Q = PQ_iv(B, t)
        r = m + 1
        if Pprev is not None and r >= t + 2:
            g = Iv.frac((2 * r - 1) ** 2, (2 * r) * (2 * r - 2))
            g = g * g
            dP = P - g * Pprev
            rho = (dP * dP).scal(r * r) / ((P * Q).scal(2))
            n += 1
            if rho.lo <= D:
                ok &= gate(f"finite rho(.,{t})>1", False, f"r={r} rho={rho}")
            if worst is None or rho.lo < worst.lo:
                worst, worstr = rho, r
        Pprev = P
    gate(f"finite certificate t={t}", ok)
    print(f"    {n} cells certified;  min lower bound at r={worstr}: "
          f"{float(worst.lo)/D:.15f}  (width {float(worst.wid()):.2e})")
    # cross-check the interval ladder against exact rho
    Mw = wall_M(2 * 60)
    for r in (t + 2, 10, 40, 55):
        Pprev = None
        for m, h, q, T, B, Cc, Bp, Cp in ladder_interval(r - 1, K=t + 1):
            P, Q = PQ_iv(B, t)
            if m == r - 1:
                g = Iv.frac((2 * r - 1) ** 2, (2 * r) * (2 * r - 2))
                g = g * g
                dP = P - g * Pprev
                rho = (dP * dP).scal(r * r) / ((P * Q).scal(2))
                ex = exact_rho(r, t, Mw)
                gate("interval brackets exact", rho.flo() <= ex <= rho.fhi(),
                     f"r={r} t={t}")
            Pprev = P
    print(f"    interval enclosures bracket exact_rho at r = {t+2},10,40,55: ok")
    return worstr, worst


# =============================================================== SECTION 3
def bC_polys(A, B, K):
    """continuum solution: b_k(tau), C_k(tau) for k=0..K, exact rationals.

        b_k = sum_{j<=k} [ A_j tau^{2(k-j)}/(2(k-j))! + B_j tau^{2(k-j)+1}/(2(k-j)+1)! ]
        C_k = b_k'                        with A_0 = 0, B_0 = 1.
    """
    b, c = [], []
    for k in range(K + 1):
        pb = [F(0)] * (2 * k + 2)
        for j in range(k + 1):
            d = k - j
            pb[2 * d] += A[j] / math.factorial(2 * d)
            pb[2 * d + 1] += B[j] / math.factorial(2 * d + 1)
        while len(pb) > 1 and pb[-1] == 0:
            pb.pop()
        b.append(pb)
        c.append([i * x for i, x in enumerate(pb)][1:] or [F(0)])
    return b, c


def PQt(b, t):
    """P_t = L((1+u)b)_t = (b_{t-1}+b_t)^2 - (b_{t-2}+b_{t-1})(b_t+b_{t+1})
       Q_t = L(b)_{t-1}  = b_{t-1}^2 - b_{t-2} b_t          (b_{-1} = 0)."""
    g = lambda i: b[i] if 0 <= i < len(b) else [F(0)]
    P = padd(pmul(padd(g(t - 1), g(t)), padd(g(t - 1), g(t))),
             pmul(padd(g(t - 2), g(t - 1)), padd(g(t), g(t + 1))), -1)
    Q = padd(pmul(g(t - 1), g(t - 1)), pmul(g(t - 2), g(t)), -1)
    return P, Q


def section_constant(TS=(2, 3, 4, 5, 6, 7, 8)):
    """EXACT gate: leading coefficient of the continuum profile, every t.

    deg b_k = 2k+1 with leading coeff 1/(2k+1)!, independent of the
    connection constants, so the leading coefficients of P_t, Q_t, P_t'
    are constant-free and  Rho_t(tau)/tau^2 -> an absolute rational.
    """
    print("[0] exact leading coefficient of the continuum profile, all fixed t")
    print(f"    {'t':>3} {'[P_t]':>14} {'[Q_t]':>14} {'c_t = Rho/tau^2':>18} "
          f"{'rho/L^2 = c_t/4':>18}")
    import random
    rnd = random.Random(20260804)
    for t in TS:
        lb = lambda k: F(1, math.factorial(2 * k + 1)) if k >= 0 else F(0)
        cP = lb(t) ** 2 - lb(t - 1) * lb(t + 1)
        cQ = lb(t - 1) ** 2 - lb(t - 2) * lb(t)
        ct = (4 * t + 2) ** 2 * cP / (8 * cQ)
        gate(f"[P_{t}] closed form",
             cP == F(8 * t + 6, math.factorial(2 * t + 1) ** 2 * (2 * t + 2) * (2 * t + 3)))
        gate(f"[Q_{t}] closed form",
             cQ == F(8 * t - 2, math.factorial(2 * t - 1) ** 2 * (2 * t) * (2 * t + 1)))
        gate(f"c_{t} closed form",
             ct == F((2 * t + 1) * (4 * t + 3), 8 * t * (t + 1) * (2 * t + 3) * (4 * t - 1)))
        # the constants really drop out: 3 arbitrary rational constant sets
        for _ in range(3):
            A = [F(0)] + [F(rnd.randint(-40, 40), rnd.randint(1, 17)) for _ in range(t + 1)]
            Bc = [F(1)] + [F(rnd.randint(-40, 40), rnd.randint(1, 17)) for _ in range(t + 1)]
            bb, cc = bC_polys(A, Bc, t + 1)
            P, Q = PQt(bb, t)
            gate(f"constant-free t={t}",
                 len(P) == 4 * t + 3 and P[4 * t + 2] == cP and
                 len(Q) == 4 * t - 1 and Q[4 * t - 2] == cQ)
        print(f"    {t:>3} {str(cP)[:14]:>14} {str(cQ)[:14]:>14} {str(ct):>18} {str(ct/4):>18}")
    print("    closed forms (exact, all t):")
    print("      [P_t] = (8t+6)/(((2t+1)!)^2 (2t+2)(2t+3))")
    print("      [Q_t] = (8t-2)/(((2t-1)!)^2 (2t)(2t+1))")
    print("      rho(r,t)/L^2 -> (2t+1)(4t+3) / (32 t(t+1)(2t+3)(4t-1))"
          "   ->  1/(32 t^2)")
    print(f"      t=2:  {F(55,9408)} = {float(F(55,9408)):.12f}"
          f"   [the chat-only constant, DERIVED]")
    return F(55, 9408)


def section3_tail(M0, t=2, verbose=True, record=True, perturb=None, dscale=1):
    """Rigorous tail certificate for column t, anchored at m = M0.

    Variable: x := tau - t0l >= 0, where t0l is an EXACT rational lower bound
    for tau_{M0}.  Every polynomial below has exact rational coefficients and
    every inequality asserted holds for every real x >= 0, hence at every
    tau_m with m >= M0.  Returns (ok, tau0, constants).
    """
    K = t + 1                                   # levels b_0..b_{t+1} needed

    def g(name, cond, extra=""):
        if record:
            return gate(name, cond, extra)
        return bool(cond)

    # ---- 3a: certified state at M0
    last = None
    for rec in ladder_interval(M0, K=K):
        last = rec
    m, h0, q0, T0, B0, C0, _, _ = last
    assert m == M0
    half_pi = PI.scal(1, 2)
    tau0 = half_pi * T0
    bI = [half_pi * z for z in B0]
    cI = C0

    # ---- 3b: frozen connection constants (intervals), exact sequential inversion
    tt = [Iv.exact(1)]
    for i in range(1, 2 * K + 3):
        tt.append(tt[-1] * tau0)                # tt[i] = tau0^i

    cA = [ZERO]                                 # A_0 = 0
    cB = [Iv.exact(1)]                          # B_0 = 1
    for k in range(1, K + 1):
        acc = ZERO
        for j in range(k):
            d = k - j
            if 2 * d - 1 >= 0:
                acc = acc + cA[j] * tt[2 * d - 1].scal(1, math.factorial(2 * d - 1))
            acc = acc + cB[j] * tt[2 * d].scal(1, math.factorial(2 * d))
        Bk = cI[k] - acc
        acc2 = ZERO
        for j in range(k):
            d = k - j
            acc2 = acc2 + cA[j] * tt[2 * d].scal(1, math.factorial(2 * d))
            acc2 = acc2 + cB[j] * tt[2 * d + 1].scal(1, math.factorial(2 * d + 1))
        Ak = bI[k] - acc2 - Bk * tau0
        cA.append(Ak)
        cB.append(Bk)

    def two(iv):
        return (iv.flo(), iv.fhi())

    Alo = [two(z)[0] for z in cA]
    Ahi = [two(z)[1] for z in cA]
    Blo = [two(z)[0] for z in cB]
    Bhi = [two(z)[1] for z in cB]
    Alo[0] = Ahi[0] = F(0)
    Blo[0] = Bhi[0] = F(1)
    t0l, t0h = two(tau0)
    if perturb:                                 # negative-control hook
        d = F(perturb)
        for j in range(1, K + 1):
            Alo[j] += d; Ahi[j] += d
            Blo[j] += d; Bhi[j] += d

    # coefficientwise lower / upper polynomials in tau (valid for tau >= 0),
    # then EXACT shift tau = t0l + x.
    bl, cl = bC_polys(Alo, Blo, K)
    bh, ch = bC_polys(Ahi, Bhi, K)
    BL = [pshift(p, t0l) for p in bl]
    BH = [pshift(p, t0l) for p in bh]
    CL = [pshift(p, t0l) for p in cl]
    CH = [pshift(p, t0l) for p in ch]
    BL[0] = BH[0] = [t0l, F(1)]
    CL[0] = CH[0] = [F(1)]

    # ---- 3c: injection budgets  D^b_k, D^C_k   (Lemma INJ)
    l43l, l43h = log_bounds(F(4, 3))
    c0const = PI_HI / 2 + l43h / 2 + F(3, 8 * M0)      # tau_j <= c0const + log(j)/2
    lM_lo, lM_hi = log_bounds(F(M0))

    def J(i):
        s = F(0)
        for pp in range(i + 1):
            s += math.perm(i, pp) * lM_hi ** (i - pp)
        return s / M0

    def tailsum(poly):
        f = pplus(poly)
        tot = F(0)
        for n, cn in enumerate(f):
            if cn == 0:
                continue
            for i in range(n + 1):
                tot += cn * math.comb(n, i) * c0const ** (n - i) * F(1, 2 ** i) * J(i)
        return tot / 4

    bplus = [pplus(p) for p in bh]
    cplus = [pplus(p) for p in ch]
    Db = [F(0)] * (K + 1)
    Dc = [F(0)] * (K + 1)
    for k in range(1, K + 1):
        Db[k] = dscale * tailsum(bplus[k - 1])
        Dc[k] = dscale * tailsum(padd(bplus[k - 1], cplus[k - 1]))

    dmax = F(1, 2 * M0)

    # ---- 3d: error envelopes in x  (Lemma ENV)
    Eb = [[F(0)]] * 1
    Ec = [[F(0)]] * 1
    Eb = [[F(0)]]
    Ec = [[F(0)]]
    for k in range(1, K + 1):
        Eck = padd([Dc[k]], pint(Eb[k - 1], dmax)) if k >= 2 else [Dc[1]]
        Ebk = padd([Db[k]], pint(Eck, dmax))
        Ec.append(Eck)
        Eb.append(Ebk)

    blo = [padd(BL[k], Eb[k], -1) for k in range(K + 1)]
    bhi = [padd(BH[k], Eb[k], +1) for k in range(K + 1)]
    clo = [padd(CL[k], Ec[k], -1) for k in range(K + 1)]
    chi = [padd(CH[k], Ec[k], +1) for k in range(K + 1)]
    blo[0], bhi[0] = BL[0], BH[0]
    clo[0], chi[0] = CL[0], CH[0]

    # ---- 3e: index m-1 versions (tau -> tau - dtau, 0 < dtau <= dmax)
    bmlo = [padd(blo[k], [dmax * z for z in chi[k]], -1) for k in range(K + 1)]
    bmhi = bhi
    cmlo = [padd(clo[k], [dmax * z for z in (bhi[k - 1] if k >= 1 else [F(0)])], -1)
            for k in range(K + 1)]
    cmhi = chi
    bmlo[0] = padd(blo[0], [dmax], -1)
    cmlo[0] = [F(1)]

    # ---- 3f: P_t, Q_t, G_t, corr    (b_{-1} = 0)
    def gl(a, i):
        return a[i] if 0 <= i < len(a) else [F(0)]

    Xlo = padd(gl(blo, t - 1), gl(blo, t))          # b_{t-1}+b_t
    Xhi = padd(gl(bhi, t - 1), gl(bhi, t))
    Ylo = padd(gl(blo, t - 2), gl(blo, t - 1))      # b_{t-2}+b_{t-1}
    Yhi = padd(gl(bhi, t - 2), gl(bhi, t - 1))
    Zlo = padd(gl(blo, t), gl(blo, t + 1))          # b_t+b_{t+1}
    Zhi = padd(gl(bhi, t), gl(bhi, t + 1))
    Xmlo = padd(gl(bmlo, t - 1), gl(bmlo, t))
    Ymhi = padd(gl(bmhi, t - 2), gl(bmhi, t - 1))
    Phi = padd(pmul(Xhi, Xhi), pmul(Ylo, Zlo), -1)
    Qhi = padd(pmul(gl(bhi, t - 1), gl(bhi, t - 1)), pmul(gl(blo, t - 2), gl(blo, t)), -1)

    CXlo = padd(gl(cmlo, t - 1), gl(cmlo, t))       # (C_{t-1}+C_t)(m-1) lower
    CYhi = padd(gl(cmhi, t - 2), gl(cmhi, t - 1))   # (C_{t-2}+C_{t-1})(m-1) upper
    CZhi = padd(gl(cmhi, t), gl(cmhi, t + 1))       # (C_t+C_{t+1})(m-1) upper
    # G = (C_{t-1}+C_t)(m-1) (X(m)+X(m-1)) - (C_{t-2}+C_{t-1})(m-1) Z(m)
    #     - (b_{t-2}+b_{t-1})(m-1) (C_t+C_{t+1})(m-1)
    Glo = padd(padd(pmul(CXlo, padd(Xlo, Xmlo)), pmul(CYhi, Zhi), -1),
               pmul(Ymhi, CZhi), -1)
    # corr <= (1001/1000) P(m)/m  and  m >= M0 (1+2x+2x^2) (1-4w)   [Lemma MG]
    #   s = tau_m - tau_{M0} >= x - w  with  w = width(tau0),  so
    #   m >= M0 e^{2s} >= M0 e^{2x} e^{-2w} >= M0 (1+2x+2x^2)(1-2w) >= M0 W(x)(1-4w).
    # The (1-4w) factor MUST be carried: dropping it asserts a larger lower
    # bound on m than the lemma licenses, which shrinks corr and makes the
    # certificate easier -- unsafe direction.  (Defect found by the replay lane
    # in the first version of this file; repaired here.  w ~ 5e-75, so R_0 is
    # unchanged, but the omission was a genuine soundness gap.)
    wfac = 1 - 4 * tau0.wid()
    V = [F(M0) * wfac, F(2 * M0) * wfac, F(2 * M0) * wfac]
    Ucl = padd(pmul(V, Glo), [F(1001, 1000) * c for c in Phi], -1)
    Fpoly = padd(pmul(Ucl, Ucl), pmul([F(8)], pmul(Phi, pmul(Qhi, pmul(V, V)))), -1)

    # ---- 3g: gates
    ok = True
    w = tau0.wid()
    ok &= g(f"t={t} M0={M0} tau0 width < 1e-6", w < F(1, 10 ** 6))
    ok &= g(f"t={t} M0={M0} b_lo nonneg", all(all(c >= 0 for c in p) for p in blo[1:]))
    ok &= g(f"t={t} M0={M0} c_lo nonneg", all(all(c >= 0 for c in p) for p in clo[1:]))
    ok &= g(f"t={t} M0={M0} Ylo Zlo nonneg",
            all(c >= 0 for c in Ylo) and all(c >= 0 for c in Zlo))
    ok &= g(f"t={t} M0={M0} Phi Qhi nonneg",
            all(c >= 0 for c in Phi) and all(c >= 0 for c in Qhi))
    ok &= g(
        f"t={t} M0={M0} V*G-corr oriented",
        all(c >= 0 for c in Ucl) and Ucl[0] > 0,
    )
    negs = [(i, c) for i, c in enumerate(Fpoly) if c < 0]
    ok &= g(f"t={t} M0={M0} F coefficients nonneg", not negs,
            f"negative at {[i for i, _ in negs][:8]}")
    if verbose:
        num = pev(Ucl, F(0)) ** 2
        den = 8 * pev(Phi, F(0)) * pev(Qhi, F(0)) * pev(V, F(0)) ** 2
        print(f"    t={t}  M0={M0}   tau0 = {tau0}   width {float(w):.2e}")
        print(f"          certified lower bound on rho at the anchor: "
              f"{float(num/den):.9f}  (margin {float(num/den)-1:+.9f})")
        print(f"          frozen constants at m=M0 (A_j, B_j):")
        for j in range(1, K + 1):
            print(f"            j={j}:  A={float(cA[j].flo()):+.12f}   "
                  f"B={float(cB[j].flo()):+.12f}")
        print(f"          injection budgets  D^b = {['%.3e' % float(z) for z in Db[1:]]}")
        print(f"                             D^C = {['%.3e' % float(z) for z in Dc[1:]]}")
        print(f"          deg F = {len(Fpoly)-1}; all {len(Fpoly)} coefficients >= 0: "
              f"{'YES' if not negs else 'NO'}")
    return ok, tau0, (cA, cB)


def section4_controls(M0=1000):
    """Negative controls for the tail machinery -- each MUST be rejected."""
    print("[4] negative controls on the tail certificate")
    bad, _, _ = section3_tail(6, t=2, verbose=False, record=False)
    gate("NC anchor below crossing rejected", not bad)
    print(f"    anchor M0=6 (tau0 ~ 3.0, inside the region where Rho<1): "
          f"{'REJECTED (correct)' if not bad else 'ACCEPTED  *** WRONG ***'}")
    print("    slack probe: shift EVERY A_j, B_j by a constant d, M0=1000")
    for d in ("-1/100000", "-1/1000", "-1/100", "-1/50", "-1/20", "-1/10", "-1/4"):
        okd, _, _ = section3_tail(M0, t=2, verbose=False, record=False, perturb=d)
        print(f"      d={d:>10}: {'accepted' if okd else 'REJECTED'}")
        if d == "-1/4":
            gate("NC large constant shift rejected", not okd)
        if d == "-1/100000":
            gate("PC tiny constant shift accepted", okd)
    print("    slack probe: inflate the injection budgets by a factor")
    for sc in (1, 10, 50, 200, 2000):
        oks, _, _ = section3_tail(M0, t=2, verbose=False, record=False, dscale=sc)
        print(f"      x{sc:<5}: {'accepted' if oks else 'REJECTED'}")
        if sc == 2000:
            gate("NC injection budgets inflated x2000 rejected", not oks)
        if sc == 1:
            gate("PC unscaled injection budgets accepted", oks)
    return True


def threshold(t, anchors, lo_floor=4):
    results = {}
    for M in anchors:
        ok, _, _ = section3_tail(M, t=t, verbose=False, record=False)
        results[M] = ok
        print(f"    t={t}  M0={M:8d}  certificate {'PASSES' if ok else 'fails'}")
    passing = sorted([M for M, v in results.items() if v])
    if not passing:
        return None
    lo = max([M for M in anchors if M < min(passing)] + [lo_floor])
    hi = min(passing)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        okm, _, _ = section3_tail(mid, t=t, verbose=False, record=False)
        if okm:
            hi = mid
        else:
            lo = mid
    print(f"    t={t}  smallest certified anchor M0 = {hi}  ->  R_0({t}) = {hi+1}")
    return hi


# =============================================================== main
def main():
    RCERT = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    M0 = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    TMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    print("=" * 78)
    print("A1_verify_tail.py   (SQ) fixed-t tail bound")
    print("=" * 78)
    print(f"certified pi enclosure width = {float(PI_HI-PI_LO):.3e}")
    print()
    section_constant()
    print()
    section1_identities()
    print()
    wr, wv = section2_finite(RCERT, t=2)
    print()
    print("[3] analytic tail certificate, t = 2")
    R2 = threshold(2, (130, 200, 250, 300, 320, 400, 500, 1000, 10000))
    print()
    print(f"    detail at the adopted anchor M0={M0}, t=2:")
    section3_tail(M0, t=2, verbose=True, record=True)
    print()
    section4_controls(M0)
    print()
    print("[5] the same machinery at higher fixed t")
    cols = {2: (R2, wr, wv, RCERT)}
    for t in range(3, TMAX + 1):
        Rt = threshold(t, (60, 200, 1000, 5000, 30000))
        if Rt is None:
            print(f"    t={t}: no anchor certified in the probed range")
            continue
        box = max(Rt + 1, t + 2)
        wtr, wtv = section2_finite(box, t=t)
        cols[t] = (Rt, wtr, wtv, box)
    print()
    print("=" * 78)
    if FAILS:
        print(f"*** {len(FAILS)} GATE FAILURES ***")
        for f in FAILS[:25]:
            print("   ", f)
        return 1
    print("ALL GATES PASS")
    print()
    print("  THEOREM (fixed-t tail).  For each t below, rho(r,t) >= 1 for all r >= R_0(t):")
    print(f"    {'t':>3} {'R_0(t)':>8} {'anchor M0':>10} {'finite box certified':>24} "
          f"{'min in box':>12} {'at r':>7}")
    for t in sorted(cols):
        Rt, wtr, wtv, box = cols[t]
        print(f"    {t:>3} {Rt+1:>8} {Rt:>10} {'%d <= r <= %d' % (t+2, box):>24} "
              f"{float(wtv.lo)/D:>12.9f} {wtr:>7}")
    print()
    for t in sorted(cols):
        Rt, wtr, wtv, box = cols[t]
        if box >= Rt + 1:
            print(f"  => the t={t} column of (SQ) is CLOSED on r >= {t+2}"
                  f"  (tail from {Rt+1}, finite box up to {box}).")
    print()
    print("  NOT covered: the joint regime t ~ c log r.  rho(r,t)/L^2 -> "
          "(2t+1)(4t+3)/(32 t(t+1)(2t+3)(4t-1)) ~ 1/(32 t^2),")
    print("  so the leading term alone needs L >= sqrt(32) t = 5.657 t, while the")
    print("  measured argmins sit at L/t ~ 4.2.  A fixed-t result does not reach there.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
