# Open problems after the sharp bridge theorem

The global shift-wall constant is closed. The questions below concern finer
structure; none is needed for the theorem in the manuscript.

## 1. Exact minima of the higher stripes

Determine, for each fixed `t>=3`,

```text
inf_(r>=t) Omega_t(r)
```

and, for `t>=4`,

```text
inf_(r>=t) rho(r,t).
```

The current proof deliberately uses `4/3` as a separator. It neither asserts
that `4/3` is approached nor identifies the minimizer of any higher stripe.
Numerical scouting can propose targets, but a durable answer needs an exact
prefix/half-line join or an all-parameter analytic proof.

## 2. The profile behind the proved decrement shell

In the deposited ladder indexing, define

```text
X_(m,k) = 2k(2k+1) B_k(m) / B_(k-1)(m),
d_(m,k) = X_(m,k-1) - X_(m,k).
```

The following bound is now an exact computer-assisted theorem:

```text
d_(m,k) <= (6/5) sqrt(a_(m+1)) / (k-1)
```

for `k=5,...,11`, `m>=68133`, and `a_(m+1)<=768`, where `a` is the increasing
Wallis-level parameter.

Leg 1 is a 192-bit outward interval ladder through `m=262144`. Leg 2 builds
one global envelope at that anchor and certifies 57 contiguous rational
windows until the lower level bound exceeds 768. All seven gates and every
lower-envelope positivity gate pass. Retargeting the same certificate to
`4/5` is rejected in its first window.

The open problem is the sharper profile suggested by the data:

```text
d_(m,k) ~ sqrt(a_(m+1)) g((k-1)/sqrt(a_(m+1))).
```

Determine `g`, prove or refute that the saturation constant tends to 1, and
derive a uniform-in-`t` deficit transfer. The proved `6/5` shell supplies
slack around this profile; it is no longer itself an open obligation.

## 3. A sharper joint mechanism

The joint proof discards part of the available two-sided corridor information.
Its row-free core cannot reach the target required for columns `t=9,10`.
Fixed-column tails give the load-bearing coverage, while D-LAW and the
sharpened corridor now provide second certificate chains. A single joint
argument that explains both columns uniformly may expose the true large-`t`
deficit profile.

This is a request for a new proof mechanism, not a missing leg in the present
theorem.

## 4. Formal replay

A proof-assistant program can be staged without redefining the theorem:

1. formalize the wall recurrence, reversal identities, and lossless transfer;
2. import the exact finite comparison streams with checked parsers;
3. formalize coefficientwise-positive polynomial certificates;
4. check the Sturm root counts and the 105-monomial joint certificate;
5. assemble the coverage partition and equality statement.

Until that program is complete, the correct label remains “exact
computer-assisted theorem,” not “machine-verified theorem.”
