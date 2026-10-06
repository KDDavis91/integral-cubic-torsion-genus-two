# Integral cubic torsion and Frobenius in genus two

**Kyle Davis — 6 October 2026**

This repository accompanies the preprint **“Integral cubic torsion and Frobenius in genus two.”**

**Zenodo DOI:** [10.5281/zenodo.23196462](https://doi.org/10.5281/zenodo.23196462)

## Summary

Let \(T\) be a symplectic \(\mathbf Z_2\)-lattice of rank four and let \(L\) be its closed-surface cubic Lie quotient. For a nondegenerate symmetric monodromy pairing \(B\), the paper computes the two-torsion in

\[
\operatorname{coker}(\rho_3(\tau_B)-1)
\]

and its normaliser action.

For primitive \(B\), the two-torsion is canonically identified with the free cubic Lie module of the rank-two quotient by the vanishing-cycle Lagrangian, reduced modulo two. For even \(B\), the Frobenius-fixed dimension depends on integral extension data and can take the values

\[
2,3,4,5,6,8,10.
\]

The paper also gives explicit genus-two curves over \(\mathbf Q_5\) with identical metric reduction graphs, graph Frobenius actions and ordinary component groups, but different cubic invariants.

## Preprint

The DOI-backed version of record for the first public preprint is available on Zenodo:

- https://doi.org/10.5281/zenodo.23196462

A copy of the PDF is included in this repository release.

## Verification

The `verification/` directory contains reproducibility material for the exact lattice and finite computations used in the paper.

Requirements:

- Python 3.10 or later
- standard library only
- run without Python's `-O` option, because assertions are part of the checks

Run:

```bash
python3 verification/verify_cubic_torsion.py --output replay.json
```

The verifier produces a deterministic JSON certificate which can be compared byte-for-byte with the supplied reference certificate.

The verifier supports the displayed algebraic and arithmetic calculations; the proofs and geometric arguments are contained in the preprint itself.

## Citation

Please cite the Zenodo preprint:

> Kyle Davis, *Integral cubic torsion and Frobenius in genus two*, 2026. DOI: 10.5281/zenodo.23196462.

Machine-readable citation metadata is provided in [`CITATION.cff`](CITATION.cff).

## Version

The initial release corresponds to the public version dated **6 October 2026**.
