# COSTAR: Scalable Second-Order CRLB Beamforming Optimization for Cell-Free ISAC Systems

MATLAB research code, stored experiment results, publication figures, and a detailed Chinese mathematical companion for **COSTAR** (CRLB-Oriented Second-order Tangent-space Algorithm with Restoration).

[中文运行说明](README_CN.md) · [Chinese companion, 118 pages](docs/COSTAR_Chinese_Companion.pdf) · [Download numerical results](https://github.com/BinYan-Veins2001/Scalable-Second-Order-CRLB-Beamforming-Optimization-for-Cell-Free-ISAC-Systems/releases/tag/v2026.10.20) · [MIT license](LICENSE)

COSTAR minimizes localization CRLB under user SINR and per-AP power constraints in cell-free ISAC systems. It combines an NNLS feasible reference decrease, tangent-space trust-region steps, exact Hessian-vector products, feasibility restoration, and a deterministic structured positive-semidefiniteness check (SPSC). SPSC uses user-block factorizations and matrix inertia to certify tangential curvature or recover a negative-curvature direction. The default solver does not invoke randomized block Krylov searches.

## Quick start

Tested with **MATLAB R2022b**. The COSTAR implementation uses MATLAB and Java numerical facilities and does not require CVX. The convex-subproblem baselines require **CVX**, with **SDPT3** selected by default; install those dependencies separately. Neither MATLAB nor third-party solvers are redistributed.

Open this repository as the MATLAB working directory:

```matlab
addpath(genpath(pwd), '-begin');
result = main_costar_paper_demo(2026);  % default: 16 APs, 4 Tx/Rx antennas, 8 users
```

To inspect the published results without running simulations:

```matlab
download_paper_results;                 % downloads the Release ZIP and verifies SHA-256
report = main_replot_paper_figures(1:6); % replot only; does not launch optimization
```

Alternatively, download `COSTAR_results_20261020.zip` from the release and extract it **into the repository root**, preserving `results/`, `benchmark/`, `audit/`, and `tests/` paths. The downloader refuses to replace existing paper MAT archives unless called with `download_paper_results(true)`.

To compute new results:

```matlab
main_run_paper_figures;  % Figures 2-6; can take substantial time
main_fig02_costar_convergence;
main_fig04_joint_stationarity;
main_fig06_structured_curvature_efficiency;
```

Figure 1 is an explicit, expensive sweep and is excluded from the default runner. To request a new sweep deliberately:

```matlab
main_run_paper_figures(struct('figures',1,'run_dimension_scaling',true));
```

On Windows, `run_paper_figures.ps1` launches independent single-computational-thread MATLAB processes. Set `-Matlab` to your installed executable when necessary.

## Paper figures and stored results

| Figure | Entry point | Measurement |
|---|---|---|
| 1 | `main_fig01_network_dimension_scaling` | CRLB and end-to-end time versus AP, antenna, and user dimensions |
| 2 | `main_fig02_costar_convergence` | CRLB, full KKT residual, and horizontal curvature trajectory |
| 3 | `main_fig03_endpoint_curvature` | Endpoint horizontal-curvature estimates for seven methods |
| 4 | `main_fig04_joint_stationarity` | Joint first-order and tangential second-order stopping accuracy |
| 5 | `main_fig05_local_quadratic` | Controlled local Newton tangential-residual contraction |
| 6 | `main_fig06_structured_curvature_efficiency` | Structured curvature evaluation versus dense tangent-space spectra |

The removed endpoint-multiplier plot is not part of this numbering. The latest manuscript includes Figures 1-5; Figure 6 is the additional structured-curvature efficiency experiment.

The stored October 19 code/results are preserved rather than presented as newly run simulations. Figure 1 is inherited from the archived dimension study; the release supplies compact plotting arrays and an optional full historical archive. Figures 2-6 use the newer experiment package. Some Figure 3 baseline endpoints are verified archived runs, with reuse checks retained. See [data provenance](docs/DATA_PROVENANCE.md) and [recorded experiment details](docs/实验结果与验证.md).

Figure 4 reports 30/30 joint stopping successes across its three settings, checked against dense tangent-space spectra. Figure 5 measures controlled residual contraction; its separate regularity record does not verify every sufficient condition of the conditional local-distance theorem. Figure 6 times the same final states, multipliers, normals, and curvature threshold; all stored timing calls use the full-space SPSC branch, with tangent and negative-direction branches checked separately. Per-run solver status, including iteration-budget exits, remains in the archives.

## Algorithms and repository structure

The six baselines are ZF-BF, SO-BF, PFO-ALM, PTR-ALM, MRFP-SCA, and SDR-SCA. Iterative methods use the shared feasible stacked-ZF initialization.

```text
optimization/        COSTAR, geometry, NNLS, restoration, HVPs, SPSC
config/              Physical and numerical defaults
baselines/           Fixed-direction baseline wrappers
channels/, scenario/, crlb/, utils/
paper_figures/       Figure 1-6 experiment and plotting functions
sim_twc/             Shared simulation helpers and historical compatibility functions
results/             CSV summaries and provenance; MAT archives come from the release
figures/             Editable MATLAB figures, PNG, PDF, EPS
benchmark/           Endpoint audits and recorded alignment comparisons
tests/               Numerical mechanism and diagnostic tests
docs/                Algorithm notes, experiment notes, and Chinese companion
```

Recommended numerical checks:

```matlab
test_latest_alignment;       % geometry, SPSC, stopping, reference curve, physical endpoint
test_paperfig_efficiency;    % structured-curvature branches and timing-input consistency
```

Historical `sim_twc` entries and older diagnostic tests are retained for traceability. The paper-facing entry points are those listed above. Re-running the historical old/new solver comparison requires its separate older solver source; its saved inputs and outputs are included in the release.

## Mathematical companion

The [118-page Chinese PDF](docs/COSTAR_Chinese_Companion.pdf) expands the current paper in chapter order. It derives NNLS reference points, exact HVPs, inertia identities, negative-curvature direction recovery, feasible descent, KKT convergence, finite-accuracy second-order stationarity, and conditional local convergence. Each chapter introduces its logic and explains the conditions used in its proofs. This is supplementary reading, not the English manuscript.

## Version and citation

Initial public release: **2026-10-20**, using the **2026-10-19** MATLAB package and mathematical companion. Core numerical source bytes are unchanged in this publication package; publication additions concern documentation, data distribution, and inherited Figure 1 import. [CITATION.cff](CITATION.cff) supplies software citation metadata.

Accompanying manuscript: *Scalable Second-Order CRLB Beamforming Optimization for Cell-Free ISAC Systems*, Bin Yan, Kaijie Zeng, Zheng Wang, and Amin Sakzad. Publication identifiers can be added when available.

## License

[MIT](LICENSE). Third-party MATLAB/CVX/solver dependencies retain their own licenses and must be installed separately.
