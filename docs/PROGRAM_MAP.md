# Program map

The four papers answer different questions about one recurrence and its
determinantal structure. Only the arrows below are logical dependencies.

```text
wall / odd-cycle recurrence
|
+-- Paper 1: sharp multiplicative Delannoy TP2
|   `-- produces the shift-wall boundary question
|
+-- Paper 2: exact grounded-path DPP algorithms
|   `-- exposes a related determinant and sampling mechanism
|
`-- Paper 3: odd cycles, sech spectra, and the square-tail inequality
    |-- proves the cross-row quotient rho >= 1
    |-- proves the bridge with structural factor 2
    `-- supplies the lossless transfer ingredients
        `-- Paper 4: the sharp shift-wall constant
            |-- isolates the unique extremal cell (1350,2)
            |-- replaces 2 by 2.655018313913361... exactly
            `-- proves a quantitative decrement shell behind the t=9 corridor
```

Paper 2 is a conceptual sibling, not a hidden premise of Paper 4. Paper 1 is
the application context that generated the boundary problem. Paper 3 is the
load-bearing mathematical input.

## What each paper contributes

| Work | Main public object | Result used here | Logical role |
| --- | --- | --- | --- |
| [Sharp Delannoy TP2](https://doi.org/10.5281/zenodo.21778524) | continued Delannoy multiplication table | identifies the shift-wall bridge as a boundary leaf | motivation and application |
| [Grounded-path DPP algorithms](https://doi.org/10.5281/zenodo.21778518) | exact rational Green-kernel algorithms | none required for the proof | structural sibling |
| [Odd cycles and square tails](https://doi.org/10.5281/zenodo.21853982) | odd-cycle polynomials, Jacobi structure, square-tail quotient | positivity, terminal formulas, top and corridor machinery, original bridge | direct premise |
| [Sharp shift-wall constant](https://doi.org/10.5281/zenodo.21866366) | bridge ratio `Omega_t(r)` and decrement shell | closes the exact global infimum and proves D-LAW on its stated window | present paper |

## The public narrative

The recurrence first appears as a Delannoy boundary object. Its reversed rows
are also odd-cycle enumerators and Gram determinants. Paper 3 proves enough
cross-row curvature to make the bridge true with factor 2. Paper 4 asks the
remaining extremal question: how much larger can that factor be while staying
valid at every coefficient and every row?

The answer is controlled by the second reverse stripe:

```text
best raw factor = 2 Omega_2(1350)
                = 2.655018313913361199726406175676...
```

All other stripes lie above the rational separator `4/3` after normalization,
so the equality cell is unique. A second certificate route for the ninth
stripe also proves a quantitative decrement shell; the remaining frontier is
the limiting profile behind that shell, not the displayed `6/5` inequality.

## Suggested reading order

1. Paper 4, Sections 1 and 3: the theorem and lossless transfer.
2. Paper 4, Section 6: the five-region `4/3` barrier.
3. Paper 3, bridge and coverage sections: the machinery reused by Paper 4.
4. Paper 1 for the larger fixed-product application.
5. Paper 2 for the determinant/algorithm branch.

The certificate map in `docs/SOURCE_MAP.md` is for auditing, not a prerequisite
for reading the proof.

## Post-release structural branch

The four releases above are closed. Follow-up work now splits into two scoped
surfaces:

- [green-path-dpp PR #1](https://github.com/toxic2040/green-path-dpp/pull/1)
  develops the two-boundary Green identity into exact checkpoint-design
  algorithms for Paper 2's next version;
- [odd-cycles-square-tail-replay issue #1](https://github.com/toxic2040/odd-cycles-square-tail-replay/issues/1)
  tracks a structural sequel around the Gram dictionary, reflected zero mode,
  Hermite entry-cell boundary, and the open all-order terminal-positivity
  transport.

The sequel begins from Paper 3's Jacobi and Gram carrier. Paper 2 supplies a
direct algorithmic offshoot, not a proof of the nonlinear sign theorem. It
does not revise Paper 1's sharp TP2 threshold or Paper 4's sharp bridge
constant. Finite-dimensional projector identities, spectral interpretation,
and the open positivity-transport theorem remain separate claim classes.
