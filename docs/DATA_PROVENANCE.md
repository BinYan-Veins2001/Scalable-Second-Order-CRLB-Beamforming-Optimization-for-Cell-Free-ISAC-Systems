# Published data and reproduction

## Versioned files

The initial public release is `v2026.10.20`. Core numerical sources come from `COSTAR_PaperFigures_20261019`; the release manifest records the numerical-source signature, file hashes, asset sizes, and SHA-256 checksums. Publication preparation does not change the optimization/configuration/baseline/model source bytes.

The repository contains source files, CSV summaries, experiment provenance, and exported figures. The large MAT files are distributed separately to keep ordinary clones small. `download_paper_results` verifies the ZIP checksum before extraction. The companion PDF is stored both under `docs/` and as a release attachment.

## Archive contents

`COSTAR_results_20261020.zip` extracts relative to the repository root:

- `results/paper_fig01/data.mat`: compact arrays used to replot the inherited dimension study.
- `results/paper_fig02/` through `results/paper_fig06/`: stored figure archives and raw experiment records.
- `results/shared_cases/`: shared physical instances and recorded solver runs.
- `benchmark/validation/`: inputs, endpoints, timing repetitions, and summary tables for the old/new alignment study.
- `audit/` and `tests/`: numerical mechanism and provenance records accompanying the published snapshot.

The optional `COSTAR_fig01_full_archive_20260907.zip` retains the full historical dimension-study archive, including per-solver histories. It extracts to `data/fig15_data.mat`. To regenerate the compact arrays, run `import_legacy_fig01(path_to_that_file)`. This import loads saved plot arrays and does not execute a dimension sweep.

## Figure-by-figure interpretation

| Figure | Source and measurement |
|---|---|
| 1 | Existing dimension-study results from `COSTAR_TWC_Simulations_ExactHVP_PTR_20260907/data/fig15_data.mat`; no new dimension sweep during publication preparation. The compact archive retains per-realization CRLB, runtime, feasibility, and KKT arrays, specification, seeds, and original metadata; heavy `*_runs` histories remain in the optional full archive. |
| 2 | Newer COSTAR trajectory. The first-order plot is the full KKT residual. The horizontal-curvature trace is an estimate, not the algorithm's lower-bound certificate. |
| 3 | Seven methods at their attained endpoints. Archived baseline endpoints are reused only after shared-input/model checks; reuse verification and raw archive remain supplied. Negative horizontal Rayleigh values identify directions; positive estimates do not independently certify the whole tangent space. |
| 4 | Ten realizations each for the default, 20-dB SINR, and 32-AP settings. All 30 satisfy the displayed full-KKT/dense-tangent stopping ratio. All terminate through full-space SPSC, with independent dense tangent checks. |
| 5 | Eleven qualified controlled local Newton perturbation samples. Fitted order 2.097, quadratic-ratio CV 12.8%. Quality-based inclusion is fixed independently of observed contraction. The reference-condition table records LICQ but does not verify strict complementarity or strict horizontal positive definiteness, so this is a controlled residual-contraction measurement. |
| 6 | AP and user sweeps, three seeds and three warm repeated paired timings. There are 24 distinct returned states; one optimizer budget exit is retained. The timed SPSC calls all take the full-space success branch. Tangent certification and negative-direction recovery are checked separately by mechanism tests. |

## Timing definitions

Paper end-to-end timing includes the common initialization and complete solver execution, including terminal refinement and SPSC. Independent post-solve dense audits are outside that timed region. The alignment experiment separately reports optimizer time with common initialization and independent audits excluded. These are different timing measurements.

Figure 6 compares the same final state, multipliers, normal basis, and threshold. Dense timing includes Hessian construction, tangent basis, projection, and spectral computation. Kernel timing does not include finding that state. A median of paired speedup ratios is distinct from the ratio of median times.

## Replotting and new computation

After extracting the result ZIP, `main_replot_paper_figures(1:6)` reads the supplied archives without solving optimization problems. `main_run_paper_figures` performs new computations for Figures 2-6 by default. The cache checks numerical and experiment signatures; changed numerical sources cannot be silently represented as the saved snapshot.

`main_fig01_network_dimension_scaling` deliberately requests a full new sweep. It is not needed to inspect the inherited Figure 1 results.

Historical alignment preparation references the separately maintained older solver; saved old/new comparison records are supplied even when that older source is not installed. Historical tests and compatibility functions are not the recommended paper-facing entry points. Use `test_latest_alignment` and `test_paperfig_efficiency` for the current numerical mechanisms.

## Default physical problem

16 APs, 4 transmit and 4 receive antennas per AP, 8 users, and a 10-dB user SINR requirement. Paper stopping tolerances are `epsilon_KKT = epsilon_H = 1e-5`; solver feasibility and plotting feasibility use their separately documented thresholds. Budget and unresolved numerical exits are retained as statuses, not relabeled as certified returns.
