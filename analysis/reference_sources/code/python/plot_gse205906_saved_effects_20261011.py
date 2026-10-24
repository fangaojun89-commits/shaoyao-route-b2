"""Figure 6: descriptive display of the retained October GSE205906 results.

Reads only the pinned public saved-effect projection; no counts, refits or SSH.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SOURCE = "results/diagnostics/GSE205906_successor_v6_approval_and_readonly_observation_20261004/ROOT_CURRENT_SAVED_EFFECT_REPORT_READONLY_V2_OBSERVED_20261006.json"
SOURCE_SHA = "b0ee89154185b361f6fdabd815fa5fbae70af6492fef3d576dc5bc5422297c45"
BASE = "$.transport_result.current_result.effect_report"
METHODS = ("gene_edger", "gene_voom", "gene_deseq")
LABELS = ("edgeR (primary)", "limma-voom", "DESeq2")
CONTRASTS = ("local_disease", "intervention")
CONTRAST_LABELS = ("SNI minus sham", "PLM minus SNI")
COLORS = ("#0072B2", "#D55E00")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError(message)


def table(path: Path, rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    source = root / SOURCE
    require(digest(source) == SOURCE_SHA, "Saved projection version changed")
    obj = json.loads(source.read_bytes())
    effect = obj["transport_result"]["current_result"]["effect_report"]
    require(effect["status"] == "DESCRIPTIVE_EFFECT_REPORT_ONLY_NOT_SCIENTIFIC_ACCEPTANCE", "Not the reviewed reporting branch")
    require(effect["q_rule"] == "STRICT_Q_LT_0_05_WITHIN_METHOD_COLLECTION_CONTRAST", "Changed q convention")
    require(effect["original_library_n"] == 9 and not effect["verified_animal_independence"], "Changed sample-unit scope")
    require(effect["gate"]["pass_n"] is None, "Unexpected E2 pass count")
    require(not effect["scientific_acceptance"] and not effect["independent_replication"], "Unexpected scientific scope")

    genes = effect["gene_distributions"]
    by_gene = {(row["table"], row["contrast"]): row for row in genes}
    require(len(genes) == len(by_gene) == 6, "Incomplete or duplicate gene summaries")
    require(set(by_gene) == {(m, c) for m in METHODS for c in CONTRASTS}, "Different method/contrast set")
    activity = effect["activity"]
    require(len(activity["columns"]) == len(set(activity["columns"])), "Duplicate activity fields")
    require(activity["row_n"] == len(activity["rows"]) == 190, "Incomplete activity projection")
    require(all(len(r) == len(activity["columns"]) for r in activity["rows"]), "Malformed activity row")
    all_activity = [dict(zip(activity["columns"], row)) for row in activity["rows"]]
    rows = [r for r in all_activity if r["collection"] == "PROGENy"]
    by_activity = {(r["program_id"], r["contrast"]): r for r in rows}
    programs = [r["program_id"] for r in rows if r["contrast"] == CONTRASTS[0]]
    require(len(programs) == len(set(programs)) == 14 and len(rows) == len(by_activity) == 28, "Not all 14 PROGENy programs")
    require(set(by_activity) == {(p, c) for p in programs for c in CONTRASTS}, "Incomplete paired activity display")
    for row in genes:
        require(row["tested_n"] == 14217 and 0 <= row["q_lt_0_05_n"] <= row["tested_n"], "Invalid gene denominator")
        width = row["ci"]["width"]
        require(width["finite_n"] == row["ci"]["available_n"] == row["tested_n"], "Missing interval-width observations")
        require(0 <= width["quartile_25"] <= width["median"] <= width["quartile_75"], "Invalid width quantiles")
    for row in rows:
        require(all(math.isfinite(row[k]) for k in ("effect", "ci_low", "ci_high", "q_within_collection_contrast")), "Missing plotted activity value")
        require(row["ci_low"] <= row["effect"] <= row["ci_high"], "Invalid activity interval")
        require(0 <= row["q_within_collection_contrast"] <= 1, "Invalid activity q")

    out.mkdir(parents=True, exist_ok=False)
    inventory: list[dict] = []

    def claim(label: str, value, location: str) -> None:
        inventory.append({"Claim": label, "Reported value": value, "Source file": SOURCE,
                          "Source location/object": BASE + location, "Verification status": "VERIFIED",
                          "Discrepancy": "None", "Required action": "None"})

    gene_export = []
    for i, row in enumerate(genes):
        width = row["ci"]["width"]
        exported = {k: row[k] for k in ("table", "method", "contrast", "tested_n", "q_lt_0_05_n")}
        exported.update({"width_q25": width["quartile_25"], "width_median": width["median"], "width_q75": width["quartile_75"], "ci_definition": row["ci_definition"]})
        gene_export.append(exported)
        for field in ("tested_n", "q_lt_0_05_n"):
            claim(f"Panels A/B {row['table']} {row['contrast']} {field}", row[field], f".gene_distributions[{i}].{field}")
        for field in ("quartile_25", "median", "quartile_75", "finite_n"):
            claim(f"Panel B {row['table']} {row['contrast']} interval width {field}", width[field], f".gene_distributions[{i}].ci.width.{field}")
    table(out / "GENE_FIGURE_SOURCE.tsv", gene_export)
    table(out / "PROGENY_FIGURE_SOURCE.tsv", rows)
    for i, row in enumerate(all_activity):
        if row["collection"] != "PROGENy":
            continue
        for field in ("effect", "ci_low", "ci_high", "q_within_collection_contrast", "n_omissions", "n_estimable", "minimum_estimable_omissions"):
            claim(f"Panel C {row['program_id']} {row['contrast']} {field}", row[field], f".activity.rows[{i}][{activity['columns'].index(field)}]")
    claim("Figure original library count", effect["original_library_n"], ".original_library_n")
    claim("Panel C complete unique PROGENy programs", len(programs), ".activity.rows[*]; unique program_id where collection=PROGENy")
    claim("Panel C paired PROGENy effect rows", len(rows), ".activity.rows[*]; collection=PROGENy")
    claim("Caption complete activity appendix rows", activity["row_n"], ".activity.row_n")
    claim("Caption nominal interval confidence level", effect["ci_caveats"]["confidence_level"], ".ci_caveats.confidence_level")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.labelsize": 8.5,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "GSE205906_saved_effects_20261011"})
    fig = plt.figure(figsize=(7.2, 8.4))
    axa = fig.add_axes((0.16, 0.715, 0.29, 0.15))
    axb = fig.add_axes((0.67, 0.715, 0.28, 0.15))
    axc = fig.add_axes((0.16, 0.19, 0.59, 0.435))
    fig.suptitle("GSE205906 | same-data method sensitivity and activity effects", x=0.5, y=0.976, fontsize=10.8, fontweight="bold")
    fig.text(0.5, 0.949, "October Ensembl101 branch; 9 submitted libraries; unresolved animal / pooling / dose histories", ha="center", fontsize=8)
    handles = [Line2D([], [], color=color, marker=marker, markersize=5, label=label)
               for color, marker, label in zip(COLORS, ("o", "D"), CONTRAST_LABELS)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.935), ncol=2, frameon=False, fontsize=8.5)
    axa.set_title("A  Gene findings: strict BH q < 0.05", loc="left", pad=9, fontsize=9, fontweight="bold")
    axb.set_title("B  Gene-interval widths", loc="left", pad=9, fontsize=9, fontweight="bold")
    axc.set_title("C  All 14 PROGENy programs: inferred score differences", loc="left", pad=12, fontsize=9, fontweight="bold")
    plotted = {"gene_counts": [], "width_quantiles": [], "activity": []}
    for ci, contrast in enumerate(CONTRASTS):
        y = np.arange(3) + (-0.16 if ci == 0 else 0.16)
        hits = [by_gene[m, contrast]["q_lt_0_05_n"] for m in METHODS]
        bars = axa.barh(y, hits, height=0.25, color=COLORS[ci], alpha=0.85)
        require([b.get_width() for b in bars] == hits, "Bar values differ")
        for yi, n in zip(y, hits):
            axa.text(n + 8, yi, str(n), va="center", fontsize=8, color=COLORS[ci])
        plotted["gene_counts"].append({"contrast": contrast, "values": hits})
        widths = [by_gene[m, contrast]["ci"]["width"] for m in METHODS]
        mid = np.array([w["median"] for w in widths])
        lo = np.array([w["quartile_25"] for w in widths])
        hi = np.array([w["quartile_75"] for w in widths])
        result = axb.errorbar(mid, y, xerr=[mid - lo, hi - mid], fmt="o" if ci == 0 else "D", color=COLORS[ci], markersize=4.0, capsize=2, linewidth=1)
        require(np.array_equal(result.lines[0].get_xdata(), mid), "Width markers differ")
        segments = result.lines[2][0].get_segments()
        require(all(np.allclose(s[:, 0], [a, b], atol=1e-14, rtol=0) for s, a, b in zip(segments, lo, hi)), "Width IQR endpoints differ")
        plotted["width_quantiles"].append({"contrast": contrast, "q25": lo.tolist(), "median": mid.tolist(), "q75": hi.tolist()})
    for ax in (axa, axb):
        ax.set_yticks(np.arange(3), LABELS, fontsize=8)
        ax.set_ylim(2.55, -0.55)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color="#E4E7EB", linewidth=0.6)
    axa.set_xlim(0, 510)
    axa.set_xticks([0, 200, 400])
    axa.set_xlabel("Genes; 14,217 tested per method / contrast")
    axb.set_xlim(0, 1.25)
    axb.set_xticks([0, 0.5, 1.0])
    axb.set_xlabel("Median and IQR across all tested genes\nWidth on log2 scale; not an effect CI")
    axc.axvline(0, color="#727B86", linewidth=0.8, linestyle=(0, (3, 3)))
    axc.set_yticks(np.arange(14), programs)
    axc.set_ylim(13.6, -0.6)
    lower = min(r["ci_low"] for r in rows)
    upper = max(r["ci_high"] for r in rows)
    axc.set_xlim(math.floor((lower - 0.1) * 4) / 4, math.ceil((upper + 0.1) * 4) / 4)
    axc.set_xticks([-1.0, -0.5, 0, 0.5, 1.0])
    axc.set_xlabel("ULM activity-score difference; stored nominal 95% interval")
    axc.grid(axis="x", color="#E4E7EB", linewidth=0.6)
    axc.tick_params(axis="y", length=0)
    axc.text(1.035, 1.025, "SNI–sham q", transform=axc.transAxes, ha="center", fontsize=8)
    axc.text(1.235, 1.025, "PLM–SNI q", transform=axc.transAxes, ha="center", fontsize=8)
    for ci, contrast in enumerate(CONTRASTS):
        for i, program in enumerate(programs):
            row = by_activity[program, contrast]
            x, lo, hi, q = (row[k] for k in ("effect", "ci_low", "ci_high", "q_within_collection_contrast"))
            y = i + (-0.15 if ci == 0 else 0.15)
            item = axc.errorbar([x], [y], xerr=[[x - lo], [hi - x]], fmt="o" if ci == 0 else "D", color=COLORS[ci], markerfacecolor=COLORS[ci] if q < 0.05 else "white", markersize=4.2, linewidth=1, capsize=2)
            require(float(item.lines[0].get_xdata()[0]) == x, "Activity point differs")
            seg = item.lines[2][0].get_segments()[0]
            require(np.allclose(seg[:, 0], [lo, hi], atol=1e-14, rtol=0), "Activity interval differs")
            axc.text(1.035 if ci == 0 else 1.235, i, f"{q:.3f}", transform=axc.get_yaxis_transform(), ha="center", va="center", fontsize=8, color=COLORS[ci])
            plotted["activity"].append({"program_id": program, "contrast": contrast, "effect": x, "ci_low": lo, "ci_high": hi, "q": q, "q_label": f"{q:.3f}", "filled": q < 0.05})
    for ax in (axa, axb, axc):
        ax.set_axisbelow(True)
        for edge in ("top", "right"):
            ax.spines[edge].set_visible(False)
    footers = ["Open markers: q >= 0.05; filled marker: q < 0.05. q is adjusted within collection and contrast.",
               "Intervals are nominal, not simultaneous or coverage-validated; panel B summarizes widths, not average-effect precision.",
               "The edgeR gene interval is a project approximate Fisher-t construction, not inversion of glmQLFTest.",
               "Both contrasts reuse SNI libraries: opposite signs are not independent reversal evidence.",
               "Descriptive molecular findings do not establish target engagement, whole-formula rescue, PSS efficacy or E2 admission."]
    for i, line in enumerate(footers):
        fig.text(0.06, 0.108 - i * 0.018, line, fontsize=7.5, va="top")
    for suffix in ("pdf", "svg", "png"):
        options = {"dpi": 300} if suffix == "png" else {}
        if suffix == "pdf":
            options["metadata"] = {"CreationDate": None, "ModDate": None, "Title": "Figure 6: GSE205906 October saved effects"}
        elif suffix == "svg":
            options["metadata"] = {"Date": None}
        fig.savefig(out / f"Figure_6_GSE205906_saved_effects.{suffix}", **options)
    plt.close(fig)
    (out / "ARTIST_VALUES.json").write_text(json.dumps(plotted, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    table(out / "NUMERIC_CLAIM_INVENTORY.tsv", inventory)
    caption = ("### Figure 6. GSE205906 October branch: method sensitivity and inferred activity effects\n\n"
               "The separately versioned Rnor6/Ensembl101 branch uses the same nine submitted libraries, not an independent cohort. "
               "SNI denotes spared nerve injury; PLM is the submitter-labelled constituent-combination group. "
               "Panel A displays all three gene methods and both contrasts using strict BH q<0.05 within method and contrast; "
               "14,217 genes were tested in every family. edgeR is primary; limma-voom and DESeq2 are same-data sensitivities, "
               "not independent validation or competing methods from which to choose the most discoveries. "
               "Panel B displays the median and interquartile range of stored nominal 95% gene-interval widths across all tested genes "
               "on the log2 scale. These are across-gene width distributions, not intervals for an average gene effect, "
               "coverage validation, or equivalence assessments. The edgeR interval uses the project approximate Fisher-t "
               "construction and is not inversion of glmQLFTest. Panel C displays every recorded PROGENy program in original order, "
               "with score differences and stored nominal 95% moderated-t intervals from ULM scores followed by limma. "
               "All paired effects and q values are displayed, with no significance or direction selection. Filled markers denote "
               "strict q<0.05; open markers denote q>=0.05 within each collection and contrast. Printed q values are rounded, "
               "but marker status uses the original values. Both contrasts reuse SNI libraries; opposite directions and differences "
               "in threshold crossing do not test a between-contrast difference or demonstrate reversal. Intervals are not "
               "simultaneous, selection-adjusted or empirically coverage-validated. Animal linkage, pooling, batch and actual-dose "
               "histories remain unresolved. Inferred score differences do not establish target engagement, a whole-formula "
               "effect, a PSS mechanism or clinical efficacy. The complete activity appendix retains all 190 PROGENy/DoRothEA "
               "records, including failed minimum omission-estimability rows; this figure does not replace it or reopen E2.\n")
    (out / "CAPTION.md").write_text(caption, encoding="utf-8")
    require(digest(source) == SOURCE_SHA, "Source changed while rendering")
    products = {p.name: {"bytes": p.stat().st_size, "sha256": digest(p)} for p in sorted(out.iterdir()) if p.is_file()}
    receipt = {"status": "SCOPED_DISPLAY_SOURCE_AND_ARTIST_CONSISTENCY_PASS", "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
               "source": SOURCE, "source_sha256": SOURCE_SHA, "renderer_sha256": digest(Path(__file__)), "command": sys.argv,
               "working_directory": str(Path.cwd()), "python": sys.version, "platform": platform.platform(),
               "matplotlib": matplotlib.__version__, "numpy": np.__version__, "randomness": "NONE",
               "gene_summary_rows": len(genes), "progeny_rows": len(rows), "all_activity_rows_retained_in_source": len(all_activity),
               "numeric_inventory_rows": len(inventory), "outputs": products, "source_unchanged": True,
               "model_run": False, "raw_data_read": False, "new_statistics_or_q": False,
               "visual_review_completed": False, "scientific_acceptance": False, "whole_manuscript_QA": False}
    (out / "FIGURE_QA.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "output": str(out), "inventory_rows": len(inventory), "bytes": sum(p.stat().st_size for p in out.iterdir())}))


if __name__ == "__main__":
    main()
