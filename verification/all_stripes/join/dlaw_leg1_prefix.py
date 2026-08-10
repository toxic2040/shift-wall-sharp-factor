#!/usr/bin/env python3
"""(D-LAW) Leg 1 — exact outward-interval prefix on the deposited ladder.

Statement being certified on the prefix (deposited A12 indexing,
X_{m,k} = 2k(2k+1) B_k(m)/B_{k-1}(m), a_{m+1} = X_{m,1}):

    For every m in [M_ENFORCE, M_PRE] and every k in {5,...,11}:
        (k-1) * (X_{m,k-1} - X_{m,k}) <= (6/5) sqrt(a_{m+1}).

Gate per (m,k), all outward fixed-point integer arithmetic at 192 bits:
    L := (k-1) * ( 2(k-1)(2k-1) B_{k-1}^2 - 2k(2k+1) B_k B_{k-2} )
    (upper endpoint; if <= 0 the gate is trivially true since d_k >= 0
    is the banked A12-SIGN theorem)
    gate: 25 * Lhi^2 * B0lo <= 36 * 6 * B1lo * (B_{k-2} B_{k-1})_lo^2
    which is the squared, cross-multiplied form of the display above.
    (a = 6 B1/B0.)

Ladder: A1_verify_tail_independent.interval_ladder recurrence, verbatim
in raw ints: h *= (2m)^2/(2m-1)^2, q = 1/h, v = h/(2m+1)^2,
tm += q, b_k += q c_k, c_k += v b_{k-1}; B_0 = tm.
Writes progress JSONL (flush) and a final exact anchor dump for Leg 2.
"""

import json
import sys
import time

sys.set_int_max_str_digits(10_000_000)

S = 1 << 192
M_PRE = 1 << 18          # 262144
M_ENFORCE = 60_000
KMAX = 12                # b_0..b_12 tracked; X up to k=12 available

OUT = __file__.replace(".py", "_progress.jsonl")
ANCHOR = __file__.replace(".py", "_anchor.json")


def main():
    t0 = time.time()
    # intervals as (lo, hi) scaled by S
    hlo = hhi = S
    tmlo = tmhi = S
    blo = [S] + [0] * KMAX
    bhi = [S] + [0] * KMAX
    clo = [S, S] + [0] * (KMAX - 1)
    chi = [S, S] + [0] * (KMAX - 1)

    worst = {k: None for k in range(5, 12)}
    failures = []
    prog = open(OUT, "a")

    for m in range(1, M_PRE + 1):
        n2, d2 = (2 * m) ** 2, (2 * m - 1) ** 2
        hlo = hlo * n2 // d2
        hhi = -((-hhi * n2) // d2)
        qlo = S * S // hhi
        qhi = -((-S * S) // hlo)
        w2 = (2 * m + 1) ** 2
        vlo = hlo // w2
        vhi = -((-hhi) // w2)
        tmlo += qlo
        tmhi += qhi
        nbl = [tmlo] + [0] * KMAX
        nbh = [tmhi] + [0] * KMAX
        for k in range(1, KMAX + 1):
            nbl[k] = blo[k] + qlo * clo[k] // S
            nbh[k] = bhi[k] + (-((-qhi * chi[k]) // S))
        ncl = [S] + [0] * KMAX
        nch = [S] + [0] * KMAX
        for k in range(1, KMAX + 1):
            ncl[k] = clo[k] + vlo * nbl[k - 1] // S
            nch[k] = chi[k] + (-((-vhi * nbh[k - 1]) // S))
        blo, bhi, clo, chi = nbl, nbh, ncl, nch
        # B_0 = tm
        blo[0], bhi[0] = tmlo, tmhi

        if m >= M_ENFORCE:
            for k in range(5, 12):
                # L upper endpoint (scaled by S^2 after products)
                lhi = (k - 1) * (
                    2 * (k - 1) * (2 * k - 1) * bhi[k - 1] * bhi[k - 1]
                    - 2 * k * (2 * k + 1) * blo[k] * blo[k - 2]
                )
                if lhi <= 0:
                    continue
                lhs = 25 * lhi * lhi * bhi[0]          # upper, ~ S^5
                rhs = 216 * blo[1] * (blo[k - 2] * blo[k - 1]) ** 2
                # scales: lhi^2 ~ S^4, x bhi[0] -> S^5;
                # blo[1] ~ S x (S^2)^2 -> S^5. Same scale, no shift.
                if lhs > rhs:
                    failures.append((m, k))
                    if len(failures) > 20:
                        break
                margin = rhs / lhs if lhs > 0 else float("inf")
                if worst[k] is None or margin < worst[k][1]:
                    worst[k] = (m, margin)
            if failures and len(failures) > 20:
                break

        if m % 20_000 == 0 or m == M_PRE:
            a_lo = 6 * blo[1] / bhi[0]
            prog.write(json.dumps({
                "m": m, "a_lo": round(a_lo, 6),
                "elapsed": round(time.time() - t0, 1),
                "failures": len(failures),
                "worst": {str(k): (worst[k][0], round(worst[k][1], 4))
                          for k in worst if worst[k]},
            }) + "\n")
            prog.flush()

    status = "PASS" if not failures else "FAIL"
    anchor = {
        "schema": "dlaw-leg1-anchor-v1",
        "status": status,
        "m_pre": M_PRE,
        "m_enforce": M_ENFORCE,
        "scale_bits": 192,
        "law": "(k-1)(X_{k-1}-X_k) <= (6/5) sqrt(a_{m+1}), k=5..11, "
               "deposited A12 indexing",
        "failures": failures[:20],
        "worst_margins": {str(k): (worst[k][0], worst[k][1])
                          for k in worst if worst[k]},
        "anchor_hex": {
            "h": [hex(hlo), hex(hhi)],
            "tm": [hex(tmlo), hex(tmhi)],
            "b_lo": [hex(x) for x in blo],
            "b_hi": [hex(x) for x in bhi],
            "c_lo": [hex(x) for x in clo],
            "c_hi": [hex(x) for x in chi],
        },
    }
    tmp = ANCHOR + ".tmp"
    with open(tmp, "w") as f:
        json.dump(anchor, f, indent=1)
    import os
    os.replace(tmp, ANCHOR)
    print(f"[{status}] leg1 to m={M_PRE}; failures={len(failures)}; "
          f"{round(time.time()-t0,1)}s; anchor -> {ANCHOR}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
