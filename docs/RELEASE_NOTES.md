# v2026.10.20

Initial public snapshot of the October 19 COSTAR implementation, stored results, and mathematical companion.

- Paper-facing MATLAB sources and six benchmark methods.
- Deterministic SPSC with full-space certification, tangent inertia, and negative-direction recovery.
- NNLS feasible reference descent, exact-HVP tCG, restoration, and precision coordination.
- Six figure entry points and stored results, with inherited Figure 1 explicitly identified.
- A 118-page Chinese derivation companion.
- MIT license and SHA-256 manifests.

## Publication checks

MATLAB R2022b, single computational thread:

- All six stored-result replot entries completed without optimization reruns.
- Current geometry, SPSC, stopping, quartic-cache, and physical-endpoint checks passed.
- Sixteen dense inertia comparisons passed.
- Structured curvature checks covered full-space, tangent, negative-curvature, and empty-tangent branches.
- Eighty-eight numerical source files match the October 19 input snapshot byte for byte.

These publication checks validate packaging and the numerical mechanisms; they do not replace the frozen experiment records or relabel historical results as new simulations.

## Release assets

| File | Contents |
|---|---|
| `COSTAR_results_20261020.zip` | Complete MAT archives from the newer code package and compact Figure 1 plotting arrays; 207 MAT files, with metadata and inventory |
| `COSTAR_fig01_full_archive_20260907.zip` | Optional full historical Figure 1 archive, including solver histories |
| `COSTAR_Chinese_Companion.pdf` | Chinese extended mathematical reading companion |
| `SHA256SUMS.txt` | SHA-256 checksums for the three data/document assets |

The two ZIPs preserve repository-relative paths. The first ZIP is sufficient to replot all six figures. The full historical Figure 1 archive is optional for examining its raw histories.
