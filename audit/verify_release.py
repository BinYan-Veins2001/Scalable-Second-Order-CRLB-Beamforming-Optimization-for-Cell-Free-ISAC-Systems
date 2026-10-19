"""Inspect the completed local release without launching any simulation."""
from pathlib import Path
import hashlib
import json
import struct

from pypdf import PdfReader

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "audit/source_manifest.json").read_text(encoding="utf-8"))
source = Path(manifest["source_root"])
checks = []
for item in manifest["files"]:
    original = source / item["path"]
    checks.append(original.is_file() and hashlib.sha256(original.read_bytes()).hexdigest() == item["sha256"])

unchanged_components = {}
for component in ("baselines", "channels", "scenario", "crlb"):
    items = [i for i in manifest["files"] if i["path"].startswith(component + "/")]
    unchanged_components[component] = all(
        hashlib.sha256((root / i["path"]).read_bytes()).hexdigest() == i["sha256"] for i in items
    )

stems = {
    2: "fig02_costar_convergence",
    3: "fig03_endpoint_curvature",
    4: "fig04_joint_stationarity",
    5: "fig05_local_quadratic",
    6: "fig06_structured_curvature_efficiency",
}
figures = []
for number, stem in stems.items():
    paths = {suffix: root / "figures" / (stem + "." + suffix) for suffix in ("fig", "png", "pdf", "eps")}
    present = all(p.is_file() and p.stat().st_size > 1024 for p in paths.values())
    pages = len(PdfReader(str(paths["pdf"])).pages) if paths["pdf"].is_file() else None
    dimensions = None
    if paths["png"].is_file():
        header = paths["png"].read_bytes()[:24]
        if header[:8] == b"\x89PNG\r\n\x1a\n":
            dimensions = struct.unpack(">II", header[16:24])
    figures.append({
        "figure": number, "all_four_formats": present, "pdf_pages": pages,
        "png_dimensions": dimensions,
        "archive_exists": (root / "results" / f"paper_fig{number:02d}" / "data.mat").is_file(),
    })

report = {
    "root": str(root), "original_source_files_unchanged": all(checks),
    "unchanged_components": unchanged_components,
    "figure1_excluded_from_this_run": not (root / "results/paper_fig01/data.mat").exists(),
    "figures": figures,
    "all_requested_figure_artifacts_exist": all(
        i["all_four_formats"] and i["pdf_pages"] == 1 and i["archive_exists"] for i in figures
    ),
}
(root / "audit/final_release.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False, indent=2))
if not (report["original_source_files_unchanged"] and all(unchanged_components.values())
        and report["figure1_excluded_from_this_run"] and report["all_requested_figure_artifacts_exist"]):
    raise SystemExit(1)
