# Fixed-column tail for \(\rho(r,2)\)

**Scope.** The single column \(t=2\). The statement below is a fixed-column
bound; it is not a joint \((r,t)\) tail bound, and nothing here asserts that
\(t=2\) is the binding column at every \(r\).

**Why it is in this release.** This is the write-up of the exact analytic
certificate built by
`vendor/odd-cycles-square-tail-replay/verification/A1_verify_tail_independent.py`.
Paper 4 consumes that certificate in two places: the fixed-column tails of
`verification/all_stripes/columns/h2_tail_cert.py` reuse its builder at the
target \(4/3\), and the D-LAW envelope certificate of
`verification/all_stripes/join/dlaw_leg2_bootstrap.py` instantiates the
injection-budget and envelope algebra of sections 3.1--3.3 below in an
anchor-shifted variable. The sections are numbered as the verifier's
equations reference them.

## Theorem

For every integer \(r\ge 313\),

\[
\boxed{\rho(r,2)\ge 1.}
\]

Thus an explicit fixed-column threshold is

\[
\boxed{R_0(2)=313.}
\]

The proof gives strict inequality; the non-strict form is the one used
downstream. All certificate arithmetic is integer or rational. Decimal values
below are diagnostics and carry no claim.

## 1. Exact reduction

Put \(m=r-1\),

\[
h_m=\left(\frac{(2m)!!}{(2m-1)!!}\right)^2,
\qquad q_m=h_m^{-1},
\qquad v_m=\frac{h_m}{(2m+1)^2},
\]

and define the triangular ladder

\[
\begin{aligned}
T_m&=T_{m-1}+q_m, & T_0&=1,\\
B_k(m)&=B_k(m-1)+q_mC_k(m-1),\\
C_k(m)&=C_k(m-1)+v_mB_{k-1}(m).
\end{aligned}
\tag{1.1}
\]

Here \(B_0=T\), \(C_0=1\), \(B_k(0)=0\) for \(k\ge1\), \(C_1(0)=1\), and
\(C_k(0)=0\) for \(k\ge2\). These are the top four reverse diagonals of the
wall recurrence after division by the relevant double factorials.

Set

\[
\begin{aligned}
P_m&=(B_1+B_2)^2-(B_0+B_1)(B_2+B_3),\\
Q_m&=B_1^2-B_0B_2,\\
g_r&=\left(\frac{(2r-1)^2}{(2r)(2r-2)}\right)^2
    =\left(1+\frac1{4r(r-1)}\right)^2.
\end{aligned}
\tag{1.2}
\]

Exact coefficient reversal gives

\[
\boxed{
\rho(r,2)=\frac{r^2(P_m-g_rP_{m-1})^2}{2P_mQ_m}.
}
\tag{1.3}
\]

The factors \(P_m,Q_m\) are positive: up to positive row scalings, they are
the interior Turán surpluses of \((1+u)\Pi_r\) and \(\Pi_r\). The
Euler--Bernoulli bridge factors \(\Pi_r\) into linear factors with positive
Bernoulli parameters, so strict Newton inequalities apply.

The identities (1.1)--(1.3), including the cross-row factor \(g_r\), are
checked by
`vendor/odd-cycles-square-tail-replay/verification/A1_reduce_check.py`
against the original wall recurrence using only `Fraction` arithmetic.

## 2. Wallis bounds and the scaled ladder

The Wallis ratio obeys

\[
\pi\left(m+\frac14\right)<h_m<
\pi\left(m+\frac12\right).
\tag{2.1}
\]

Indeed, \(h_m/(m+1/4)\) is strictly decreasing, \(h_m/(m+1/2)\) is strictly
increasing, and both tend to \(\pi\). Their successive quotients reduce to
rational inequalities:

\[
\frac{(2m)^2(4m-3)}{(2m-1)^2(4m+1)}<1,
\qquad
\frac{(2m)^2(2m-1)}{(2m-1)^2(2m+1)}>1.
\]

Define

\[
\tau_m=\frac\pi2T_m,\qquad
a_m=\tau_m-\tau_{m-1}=\frac{\pi}{2h_m},
\qquad b_k(m)=\frac\pi2B_k(m),
\]

and

\[
\eta_m=\left(\frac{2h_m}{\pi(2m+1)}\right)^2.
\]

Then (1.1) becomes

\[
\begin{aligned}
b_k(m)&=b_k(m-1)+a_mC_k(m-1),\\
C_k(m)&=C_k(m-1)+a_m\eta_m b_{k-1}(m).
\end{aligned}
\tag{2.2}
\]

Equation (2.1) implies

\[
\frac1{2m+1}<a_m<\frac1{2m},
\qquad 0<1-\eta_m<a_m,
\qquad r a_m>\frac12.
\tag{2.3}
\]

For the middle inequality, the lower Wallis bound gives

\[
\sqrt{\eta_m}>1-\frac1{4m+2},
\]

so \(1-\eta_m<1/(2m+1)<a_m\).

Let \(P,Q\) now denote the expressions (1.2) formed from the \(b_k\)'s. The
common scaling by \(\pi/2\) cancels from (1.3). If

\[
X=b_1+b_2,\qquad Y=b_0+b_1,\qquad Z=b_2+b_3,
\]

direct expansion of \(P_m-P_{m-1}\) using (2.2) gives

\[
P_m-P_{m-1}=a_mG_m,
\tag{2.4}
\]

where

\[
\begin{aligned}
G_m={}&(C_1^-+C_2^-)(X_m+X^-)
       -(C_0^-+C_1^-)Z_m
       -Y^-(C_2^-+C_3^-).
\end{aligned}
\tag{2.5}
\]

A superscript minus means index \(m-1\). Consequently

\[
\rho(r,2)=
\frac{(ra_m)^2\bigl(G_m-c_m\bigr)^2}{2P_mQ_m},
\qquad
c_m=\frac{(g_r-1)P_{m-1}}{a_m}.
\tag{2.6}
\]

## 3. Polynomial comparison system

Fix \(m_0=312\). At \(\tau_{m_0}\), match the discrete state to the
polynomial solution

\[
\beta_k'=\gamma_k,\qquad
\gamma_k'=\beta_{k-1},\qquad
\beta_0=\tau,\quad\gamma_0=1.
\tag{3.1}
\]

Its general form is

\[
\beta_k(\tau)=
\sum_{j=0}^k\left[
 A_j\frac{\tau^{2(k-j)}}{(2(k-j))!}
 +B_j\frac{\tau^{2(k-j)+1}}{(2(k-j)+1)!}
\right],
\tag{3.2}
\]

with \(A_0=0,B_0=1\), and \(\gamma_k=\beta_k'\). The six constants
\(A_j,B_j\), \(1\le j\le3\), are obtained by triangular inversion of the
enclosed state at \(m_0\).

The verifier encloses \(\pi\) using Machin's formula

\[
\pi=16\arctan(1/5)-4\arctan(1/239)
\]

with alternating-series remainders. It propagates the ladder in fixed-point
rational intervals with denominator \(2^{256}\). Thus the constants in (3.2),
rather than rounded approximations to them, are enclosed.

### 3.1 Local error

On the final certified envelopes, every \(\beta_k,\gamma_k\) used below is
nonnegative. Taylor's integral formula and (2.2) then give the one-step
defects

\[
\begin{aligned}
0\le R^\beta_{k,j}
  &\le a_j^2\beta_{k-1}(\tau_j),\\
|R^\gamma_{k,j}|
  &\le a_j^2\bigl(\beta_{k-1}(\tau_j)
                       +\gamma_{k-1}(\tau_j)\bigr).
\end{aligned}
\tag{3.3}
\]

For the second line, the right-rectangle error contributes at most
\(a_j^2\gamma_{k-1}(\tau_j)\), while the damping error contributes at most
\(a_j(1-\eta_j)\beta_{k-1}(\tau_j)<a_j^2\beta_{k-1}(\tau_j)\). This is a
closed bootstrap: the budgets are first constructed under the nonnegativity
hypothesis, and the resulting coefficientwise lower envelopes prove that
hypothesis on the whole tail.

### 3.2 Closed tail budgets

For \(j\ge m_0\), (2.1) gives \(a_j^2<1/(4j^2)\). Convex midpoint integration
gives

\[
\sum_{k=1}^j\frac1{k+1/4}
\le \log\frac{j+3/4}{3/4}
\le \log j+\log(4/3)+\frac{3}{4j}.
\]

Hence

\[
\tau_j\le c_0+\frac12\log j,
\qquad
c_0=\frac\pi2+\frac12\log(4/3)+\frac{3}{8m_0}.
\tag{3.4}
\]

This bound is global in \(j\), not windowed; the \(\tfrac12\log(4/3)\) term is
the constant of integration of \(\log((j+3/4)/(3/4))\).

For \(0\le i\le5\), the function \(x^{-2}(\log x)^i\) is decreasing on
\([m_0,\infty)\). Its exact integral is

\[
\int_{m_0}^{\infty}\frac{(\log x)^i}{x^2}\,dx
=\frac1{m_0}\sum_{p=0}^i
  \frac{i!}{(i-p)!}(\log m_0)^{i-p}.
\tag{3.5}
\]

The verifier encloses both logarithms by the positive atanh series and uses
(3.4)--(3.5) coefficientwise. Applied to (3.3), this produces global
injection budgets \(D_k^\beta,D_k^\gamma\). For orientation only, they are

\[
\begin{aligned}
D^\beta&\approx(4.076\!\times10^{-3},
                 1.839\!\times10^{-2},
                 2.899\!\times10^{-2}),\\
D^\gamma&\approx(4.878\!\times10^{-3},
                 2.886\!\times10^{-2},
                 5.420\!\times10^{-2}).
\end{aligned}
\]

The exact values are the rational expressions constructed in
`A1_verify_tail_independent.py`.

### 3.3 Global envelopes

Let \(x=\tau-t_0^-\), where \([t_0^-,t_0^+]\) encloses \(\tau_{m_0}\), and put
\(d_*=1/(2m_0)\). For a polynomial \(E\) with nonnegative coefficients define

\[
\mathcal I E(x)=\int_0^xE(s)\,ds+d_*E(x).
\tag{3.6}
\]

If the mesh increments are at most \(d_*\), monotonicity gives

\[
\sum_j a_jE(x_j)\le\mathcal I E(x).
\tag{3.7}
\]

For right endpoints, subtract the integral on each interval and bound the
excess by \(d_*[E(x_j)-E(x_{j-1})]\); the excess telescopes. Left endpoints
are bounded directly by the integral.

Starting from the total budgets, define triangular error polynomials

\[
\begin{aligned}
E^\gamma_1&=D^\gamma_1,\\
E^\gamma_k&=D^\gamma_k+\mathcal I E^\beta_{k-1}\quad(k\ge2),\\
E^\beta_k&=D^\beta_k+\mathcal I E^\gamma_k.
\end{aligned}
\tag{3.8}
\]

Induction in \(k\) and in the mesh index, using (3.3) and (3.7), proves

\[
|b_k(m)-\beta_k(\tau_m)|\le E^\beta_k(x),
\qquad
|C_k(m)-\gamma_k(\tau_m)|\le E^\gamma_k(x).
\tag{3.9}
\]

The interval connection constants in (3.2), shifted by \(t_0^-\), and (3.9)
give explicit polynomial lower and upper envelopes for every factor in (2.5).
The exact gate verifies coefficientwise nonnegativity of all current-index
lower envelopes and all predecessor lower envelopes.

## 4. The index and correction bounds

Let \(w=t_0^+-t_0^-\). Since \(a_j<1/(2j)\),

\[
\tau_m-\tau_{m_0}<\frac12\log(m/m_0).
\]

Since \(\tau_m-\tau_{m_0}\ge x-w\),

\[
m\ge V(x):=m_0(1-4w)(1+2x+2x^2).
\tag{4.1}
\]

Here \(e^{2x}\ge1+2x+2x^2\), and \(e^{-2w}\ge1-2w\ge1-4w\). The factor
\(1-4w\) is retained in the exact certificate; omitting it would use
\(t_0^-\) as though it were the exact anchor.

Let \(\overline P,\overline Q\) be the polynomial upper bounds formed from the
envelopes, and let \(\underline G\) be the lower bound obtained by using lower
envelopes in the positive product of (2.5) and upper envelopes in its two
negative products. The predecessor sign gates make these product bounds valid.

The correction in (2.6) satisfies

\[
c_m<\frac{P_m}{m}.
\tag{4.2}
\]

To see this, first note that the certificate below gives \(\underline G\ge0\),
so (2.4) yields \(P_{m-1}\le P_m\). Also \(a_m>1/(2m+1)\), and direct algebra
gives

\[
m(g_r-1)(2m+1)
=\frac{2m+1}{2(m+1)}
 +\frac{2m+1}{16m(m+1)^2}<1.
\]

Combining (4.1)--(4.2), with a small rational slack, gives

\[
G_m-c_m\ge
\frac{U(x)}{V(x)},
\qquad
U(x)=V(x)\underline G(x)-\frac{1001}{1000}\overline P(x).
\tag{4.3}
\]

## 5. Exact polynomial certificate

The verifier constructs

\[
\mathcal F(x)=
U(x)^2-8V(x)^2\overline P(x)\overline Q(x).
\tag{5.1}
\]

At \(m_0=312\):

- \(U\) has nonnegative rational coefficients;
- \(\overline P\) and \(\overline Q\) have nonnegative rational coefficients;
- \(\mathcal F\) has degree \(22\), and all \(23\) rational coefficients are
  strictly positive;
- the anchor interval width is less than \(10^{-60}\), and the exact factor
  \(1-4w\) in (4.1) is positive;
- the maximum logarithmic moment degree is \(5\), and the exact enclosure
  verifies \(2\log m_0\ge5\).

Therefore (5.1) is nonnegative for every real \(x\ge0\). Since \(U\ge0\),
(4.3) gives

\[
(G_m-c_m)^2
\ge 8\overline P(x)\overline Q(x)
\ge 8P_mQ_m.
\tag{5.2}
\]

Finally, (2.3), (2.6), and (5.2) yield

\[
\rho(r,2)
\ge \frac{(ra_m)^2\,8P_mQ_m}{2P_mQ_m}
>\frac14\cdot4=1
\]

for every \(m\ge312\), equivalently every \(r\ge313\). This proves the
theorem.

## 6. Replay

Both verifiers ship with this release. From the release root:

```bash
cd vendor/odd-cycles-square-tail-replay/verification
python3 A1_reduce_check.py
python3 A1_verify_tail_independent.py
```

The first checks the wall reduction, the cross-row normalization, and a
pinning negative control. The second checks the analytic certificate,
brackets the interval ladder against exact `Fraction` states at three
checkpoints, and rejects a negative control in which every injection budget is
inflated by a factor of \(2000\).

The successful tail replay ends with:

```text
certificate polynomial degree: 22
negative control, injection budgets x2000 rejected: True
boundary control, anchor m0=311 rejected: True
ALL EXACT GATES PASS
THEOREM CERTIFIED: rho(r,2) >= 1 for every integer r >= 313
```

`SHA256SUMS` binds both verifiers, so no separate source pin is needed here.

## 7. What this does not cover

The bound is for fixed \(t=2\). It says nothing about the joint regime
\(t\asymp\log r\), and no joint-tail claim is made from it.
