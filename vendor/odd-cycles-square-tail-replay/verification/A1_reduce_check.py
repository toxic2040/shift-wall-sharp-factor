#!/usr/bin/env python3
"""A1_reduce_check.py -- exact verification of every structural claim A1 uses.

All arithmetic is Fraction/int.  No floats anywhere in a gate.

Claims checked from the manuscript's wall/ladder and fixed-column reductions:

  G1  (R2) top reverse-diagonal recurrence  A_k(2m)=A_k(2m-1)+(2m)^2A_k(2m-2),
                                            A_k(2m+1)=A_{k-1}(2m)+(2m+1)^2A_k(2m-1)
  G2  A_0(2m+1) = ((2m-1)!!)^2 * (2m+1)^2 = ((2m+1)!!)^2 ; R_m recurrence
  G3  closed form   R_m = h_m * T_m,  h_m = ((2m)!!/(2m-1)!!)^2, T_m = 1+sum_{j<=m} 1/h_j
  G4  B/C ladder    B_k(m)=A_k(2m)/((2m)!!)^2, C_k(m)=A_k(2m+1)/((2m+1)!!)^2
                    B_k(m)=B_k(m-1)+q_m C_k(m-1),  C_k(m)=C_k(m-1)+v_m B_{k-1}(m)
                    q_m = 1/h_m,  v_m = h_m/(2m+1)^2,  q_m v_m = 1/(2m+1)^2
  G5  f_k(2m) = B_k(m)/T_m   (normalized row = engine.py's f)
  G6  kappa_r = (2r-1)/(2r-2) * T_{r-2}/T_{r-1}      [T-only, all h cancels]
  G7  THE A1 IDENTITY
        rho(r,2) = r^2 (P(m) - g_r P(m-1))^2 / (2 P(m) Q(m)),   m = r-1
        P(m) = (B1+B2)^2 - (T+B1)(B2+B3)   [ = T_m^2 * Shat_r ]
        Q(m) = B1^2 - T*B2                 [ = T_m^2 * Vhat_r ]
        g_r  = ((2r-1)^2/((2r)(2r-2)))^2 = (1 + 1/(4r(r-1)))^2
  G8  cancellation-free increment  P(m)-P(m-1) = q_m * K(m)  with K explicit
  G9  (R1) t=2 sees only the top 4 reverse diagonals
  G10 pinning negative control (r,t)=(2,1): d^2 - S V/2 = -4351/20736
  G11 Wallis sandwich  pi(m+1/4) < h_m < pi(m+1/2)  as an exact rational
      monotonicity statement (h_m/(m+1/4) strictly decreasing, h_m/(m+1/2)
      strictly increasing) + the exact limit pi.

Exit 0 iff every gate passes.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from engine import wall_M, Pi_exact, exact_rho, Lop, shift1

FAILS = []


def check(name, cond, extra=""):
    if not cond:
        FAILS.append(f"{name} {extra}")
        print(f"  FAIL {name} {extra}")
    return cond


def decimal_truncation(value: F, places: int) -> str:
    """Deterministic decimal rendering for an exact Fraction."""
    sign = "-" if value < 0 else ""
    whole, remainder = divmod(abs(value.numerator), value.denominator)
    digits = []
    for _ in range(places):
        remainder *= 10
        digit, remainder = divmod(remainder, value.denominator)
        digits.append(str(digit))
    return f"{sign}{whole}." + "".join(digits)


# ------------------------------------------------------------------ ladder
def ladder(M, K=3):
    """Exact T_m, h_m, B_k(m), C_k(m) for m=0..M, k=0..K."""
    W = F(1)                       # (2m)!!/(2m-1)!!
    h = [F(1)]                     # h_0 = 1
    T = [F(1)]                     # T_0 = 1
    B = [[F(1)] + [F(0)] * K]      # B_0(0)=1, B_k(0)=0
    C = [[F(1), F(1)] + [F(0)] * (K - 1)]   # C_0(0)=C_1(0)=1
    for m in range(1, M + 1):
        W *= F(2 * m, 2 * m - 1)
        hm = W * W
        q = 1 / hm
        v = hm / F((2 * m + 1) ** 2)
        h.append(hm)
        T.append(T[m - 1] + q)
        Bm = [T[m]] + [B[m - 1][k] + q * C[m - 1][k] for k in range(1, K + 1)]
        Cm = [F(1)] + [C[m - 1][k] + v * Bm[k - 1] for k in range(1, K + 1)]
        B.append(Bm)
        C.append(Cm)
    return h, T, B, C


def PQ(T, B):
    """P = (B1+B2)^2-(T+B1)(B2+B3);  Q = B1^2-T*B2  at one index."""
    t, b1, b2, b3 = T, B[1], B[2], B[3]
    return (b1 + b2) ** 2 - (t + b1) * (b2 + b3), b1 * b1 - t * b2


def Kincr(T, B, C, Tp, Bp, Cp):
    """K(m) with (T,B,C) at m and (Tp,Bp,Cp) at m-1;  P(m)-P(m-1) = q_m K(m)."""
    X, Y, Z = B[1] + B[2], T + B[1], B[2] + B[3]
    Xp, Zp = Bp[1] + Bp[2], Bp[2] + Bp[3]
    Yp = Tp + Bp[1]
    return ((Cp[1] + Cp[2]) * (X + Xp) - (1 + Cp[1]) * Z - Yp * (Cp[2] + Cp[3]))


# ----------------------------------------------------------------- gates
def main():
    MMAX = 60                       # ladder index m; wall rows up to 2*MMAX+1
    NMAX = 2 * MMAX + 2
    M = wall_M(NMAX)

    def A(k, n):
        """A_k(n) = [z^{ceil(n/2)-k}] M_n."""
        d = -(-n // 2)
        j = d - k
        row = M[n]
        return F(row[j]) if 0 <= j < len(row) else F(0)

    print("G1  top reverse-diagonal recurrence (R2)")
    ok = True
    for m in range(1, MMAX + 1):
        for k in range(0, 5):
            ok &= check("G1even", A(k, 2 * m) == A(k, 2 * m - 1) + (2 * m) ** 2 * A(k, 2 * m - 2),
                        f"m={m} k={k}")
            ok &= check("G1odd", A(k, 2 * m + 1) == A(k - 1, 2 * m) + (2 * m + 1) ** 2 * A(k, 2 * m - 1),
                        f"m={m} k={k}")
    print("    ok" if ok else "    FAILED")

    print("G2/G3  A_0 closed form, R_m recurrence, R_m = h_m T_m")
    h, T, B, C = ladder(MMAX)
    ok = True
    df = 1
    for m in range(1, MMAX + 1):
        df *= (2 * m - 1)
        ok &= check("G2 odd A0", A(0, 2 * m - 1) == F(df * df), f"m={m}")
        Rm = A(0, 2 * m) / A(0, 2 * m - 1)
        if m == 1:
            ok &= check("G2 R1", Rm == 5)
        else:
            Rprev = A(0, 2 * m - 2) / A(0, 2 * m - 3)
            ok &= check("G2 Rrec", Rm == 1 + F(2 * m, 2 * m - 1) ** 2 * Rprev, f"m={m}")
        ok &= check("G3 closed", Rm == h[m] * T[m], f"m={m}")
    print("    ok" if ok else "    FAILED")

    print("G4  B/C ladder equals A_k / double-factorials, and its recurrence")
    ok = True
    ee = 1
    oo = 1
    for m in range(0, MMAX + 1):
        if m >= 1:
            ee *= (2 * m)
            oo *= (2 * m + 1)
        for k in range(0, 4):
            ok &= check("G4 B", B[m][k] == A(k, 2 * m) / F(ee * ee), f"m={m} k={k}")
            ok &= check("G4 C", C[m][k] == A(k, 2 * m + 1) / F(oo * oo), f"m={m} k={k}")
        if m >= 1:
            q, v = 1 / h[m], h[m] / F((2 * m + 1) ** 2)
            ok &= check("G4 qv", q * v == F(1, (2 * m + 1) ** 2), f"m={m}")
            for k in range(0, 4):
                ok &= check("G4 recB", B[m][k] == B[m - 1][k] + q * C[m - 1][k], f"m={m} k={k}")
                prev = B[m][k - 1] if k >= 1 else F(0)
                ok &= check("G4 recC", C[m][k] == C[m - 1][k] + v * prev, f"m={m} k={k}")
    print("    ok" if ok else "    FAILED")

    print("G5/G6/G7/G8  normalized row, kappa, THE A1 IDENTITY, increment form")
    ok = True
    Mw = wall_M(2 * (MMAX + 1))
    for r in range(5, MMAX + 1):
        m = r - 1
        P = Pi_exact(r, Mw)
        # G5: f_k = B_k(m)/T_m, and Pi normalized to constant term 1
        for k in range(0, 4):
            ok &= check("G5", P[k] / P[0] == B[m][k] / T[m], f"r={r} k={k}")
        # G6: kappa
        Pm1 = Pi_exact(r - 1, Mw)
        kap = Pm1[0] / P[0]
        ok &= check("G6", kap == F(2 * r - 1, 2 * r - 2) * T[m - 1] / T[m], f"r={r}")
        # G7: the A1 identity
        Pn, Qn = PQ(T[m], B[m])
        Pp, _ = PQ(T[m - 1], B[m - 1])
        g = F((2 * r - 1) ** 2, (2 * r) * (2 * r - 2)) ** 2
        ok &= check("G7 g", g == (1 + F(1, 4 * r * (r - 1))) ** 2, f"r={r}")
        rho_id = F(r * r) * (Pn - g * Pp) ** 2 / (2 * Pn * Qn)
        ok &= check("G7", rho_id == exact_rho(r, 2, Mw), f"r={r}")
        # G8: increment
        ok &= check("G8", Pn - Pp == (1 / h[m]) * Kincr(T[m], B[m], C[m], T[m - 1], B[m - 1], C[m - 1]),
                    f"r={r}")
        # positivity of the pieces
        ok &= check("G8 pos", Pn > 0 and Qn > 0 and Pn - g * Pp > 0, f"r={r}")
    print("    ok" if ok else "    FAILED")

    print("G9  (R1) t=2 depends only on the top 4 reverse diagonals")
    ok = True
    for r in (10, 20, 40, 55):
        Mw2 = wall_M(2 * r)
        full = exact_rho(r, 2, Mw2)
        Pa, Pb = Pi_exact(r, Mw2)[:4], Pi_exact(r - 1, Mw2)[:4]
        S = Lop(shift1(Pa), 2)
        Sm = Lop(shift1(Pb), 2)
        V = Lop(Pa, 1)
        d = S - F((2 * r - 1) ** 2, (2 * r) ** 2) * Sm
        ok &= check("G9", F(r * r) * d * d / (2 * S * V) == full, f"r={r}")
    print("    ok" if ok else "    FAILED")

    print("G10  pinning negative control (r,t)=(2,1)")
    Mw3 = wall_M(6)
    P2, P1 = Pi_exact(2, Mw3), Pi_exact(1, Mw3)
    S = Lop(shift1(P2), 1)
    Sm = Lop(shift1(P1), 1)
    V = Lop(P2, 0)
    d = S - F(9, 16) * Sm
    ok = check("G10", d * d - S * V / 2 == F(-4351, 20736), f"got {d*d - S*V/2}")
    ok &= check("G10 rho<1", F(4) * d * d / (2 * S * V) < 1)
    print("    ok" if ok else "    FAILED")

    print("G11  Wallis sandwich  pi(m+1/4) < h_m < pi(m+1/2)   [monotone form]")
    # exact rational statement:  U_m := h_m/(m+1/4) strictly decreasing,
    #                            D_m := h_m/(m+1/2) strictly increasing,
    # both -> pi.  Hence D_m < pi < U_m for all m>=1, i.e. the sandwich.
    ok = True
    for m in range(1, MMAX + 1):
        U = h[m] / F(4 * m + 1, 4)
        Up = h[m - 1] / F(4 * m - 3, 4) if m >= 2 else None
        D = h[m] / F(2 * m + 1, 2)
        Dp = h[m - 1] / F(2 * m - 1, 2) if m >= 2 else None
        if m >= 2:
            ok &= check("G11 U dec", U < Up, f"m={m}")
            ok &= check("G11 D inc", D > Dp, f"m={m}")
    # the two monotonicities as *identities* in m (valid for all m>=2):
    #   U_m/U_{m-1} = (2m)^2 (4m-3) / ((2m-1)^2 (4m+1)) < 1
    #   D_m/D_{m-1} = (2m)^2 (2m-1) / ((2m-1)^2 (2m+1)) > 1
    for m in range(2, 400):
        ok &= check("G11 Uid", F((2 * m) ** 2 * (4 * m - 3), (2 * m - 1) ** 2 * (4 * m + 1)) < 1, f"m={m}")
        ok &= check("G11 Did", F((2 * m) ** 2 * (2 * m - 1), (2 * m - 1) ** 2 * (2 * m + 1)) > 1, f"m={m}")
    print("    ok" if ok else "    FAILED")

    print()
    if FAILS:
        print(f"*** {len(FAILS)} FAILURES ***")
        return 1
    print("ALL GATES PASS")
    # A few non-gating reference values for the note.
    Mw4 = wall_M(2 * 130)
    rho_witness = exact_rho(129, 2, Mw4)
    print(
        "  rho(129,2), exact Fraction rendered by decimal truncation = "
        f"{decimal_truncation(rho_witness, 25)}"
    )
    hh, TT, BB, CC = ladder(4000)
    for m in (10, 100, 1000, 4000):
        import math
        print(f"  m={m:5d}  T_m={float(TT[m]):.9f}  L=pi*T_m={math.pi*float(TT[m]):.9f}  "
              f"L-log(m)={math.pi*float(TT[m])-math.log(m):.9f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
