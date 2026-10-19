# Fig. 6 provenance and branch coverage

Initial experiment provenance backfill: timed dense-Hessian/EVD and SPSC kernel implementations and physical problem-generation files were verified unchanged. The compute wrapper was subsequently amended to use all returned endpoints, preserve unsuccessful optimizer exits, and validate complete cache fingerprints. Older successful samples were retained only after the saved final-iteration proof. No timed calls or raw timing arrays were replaced during this provenance step.

Unique returned endpoints: 24; optimizer-certified endpoints: 23; iteration-budget exits: 1.

Timed SPSC branches: 24 ambient certificates, 0 tangent certificates, 0 negative-direction recoveries. The plotted curves therefore measure the ambient certification path of SPSC. Tangent certification and negative-direction recovery are separately verified by `tests/test_paperfig_efficiency.m`; their runtime is not represented by these curves.

Benchmark source SHA-256: `75d2addf127697c1a581768994d4e311780b574c5d30d73d3c23e654ffe83d27`. Per-file hashes are saved in `source_manifest.json`.
