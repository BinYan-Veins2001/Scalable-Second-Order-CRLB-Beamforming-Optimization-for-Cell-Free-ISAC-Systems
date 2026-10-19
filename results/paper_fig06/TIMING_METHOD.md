# Computational efficiency of structured curvature evaluation

Returned endpoint of each fresh COSTAR trajectory, with the same rule for successful and iteration-budget exits; solver_success and actual_terminal_spsc are recorded for every case. This benchmark evaluates curvature computation at fixed feasible iterates, independently of whether the optimization has reached its stationarity tolerances. Both methods receive the same W, nonnegative NNLS lambda, orthonormal working-normal Q, and scalar sigma selected by SPSC. SPSC recomputes all block factorizations and its full-space/tangent tests, including direction recovery when needed. Dense tangent EVD independently assembles the analytic ambient Hessian, completes Q by full QR, projects to its orthogonal complement, and computes the complete eigensystem. Both methods are warmed once, then measured in interleaved order in a single-computational-thread MATLAB process. Per-realization median of repeated calls; curve median across declared seeds. Problem construction, COSTAR optimization, geometry/NNLS construction, independent HVP audits and all graphics are excluded. The saved lower bound is not an eigenvalue estimate. No older experiment timings are imported.

MATLAB: 9.13.0.2049777 (R2022b). Computational threads: 1.

Seeds: [2026 2027 2028]. Warm timed repeats: 3.

Point-trial rows: 27. Unique trajectory jobs: 24 (the common M/U intersection is reused across panels).

Each trace contributes its returned endpoint. Kernel branch counts appear in `curve_summary.csv`; all checkpoints and raw repeated timings are saved under `jobs/`. An ambient certificate bypasses tangent solves by design; its timing is the actual SPSC workflow.
