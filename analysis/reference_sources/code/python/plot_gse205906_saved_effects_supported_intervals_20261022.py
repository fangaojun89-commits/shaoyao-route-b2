"""Presentation correction; retain historical source, omit unsupported edgeR CI."""
from pathlib import Path
import hashlib

base = Path(__file__).with_name('plot_gse205906_saved_effects_20261011.py')
raw = base.read_bytes()
assert hashlib.sha256(raw).hexdigest() == 'fad8c7f87bbf7c6042e2843231aad780e15fa58e3d5a54cfc48122c4ed8cf6ff'
text = raw.decode('utf-8')

def replace(old, new):
    global text
    assert text.count(old) == 1, 'CORRECTION_SOURCE_SLOT_CHANGED'
    text = text.replace(old, new)

replace('axc = fig.add_axes((0.16, 0.19, 0.59, 0.435))',
        'axc = fig.add_axes((0.16, 0.19, 0.59, 0.410))')
replace('        width = row["ci"]["width"]\n        require(width["finite_n"]',
        '        if row["table"] == "gene_edger":\n            continue\n        width = row["ci"]["width"]\n        require(width["finite_n"]')
replace('        exported.update({"width_q25": width["quartile_25"], "width_median": width["median"], "width_q75": width["quartile_75"], "ci_definition": row["ci_definition"]})',
        '        exported.update({"width_q25": width["quartile_25"] if row["table"] != "gene_edger" else "NA", "width_median": width["median"] if row["table"] != "gene_edger" else "NA", "width_q75": width["quartile_75"] if row["table"] != "gene_edger" else "NA", "ci_definition": row["ci_definition"] if row["table"] != "gene_edger" else "NOT_REPORTED_UNSUPPORTED_HISTORICAL_INTERVAL"})')
replace('        for field in ("quartile_25", "median", "quartile_75", "finite_n"):\n            claim(',
        '        for field in (() if row["table"] == "gene_edger" else ("quartile_25", "median", "quartile_75", "finite_n")):\n            claim(')
replace('        widths = [by_gene[m, contrast]["ci"]["width"] for m in METHODS]',
        '        widths = [by_gene[m, contrast]["ci"]["width"] for m in METHODS[1:]]\n        y = y[1:]')
replace('        plotted["width_quantiles"].append({"contrast": contrast, "q25": lo.tolist(), "median": mid.tolist(), "q75": hi.tolist()})',
        '        plotted["width_quantiles"].append({"contrast": contrast, "methods": list(METHODS[1:]), "edgeR_interval": "NOT_REPORTED", "q25": lo.tolist(), "median": mid.tolist(), "q75": hi.tolist()})')
replace('    axb.set_xlim(0, 1.25)',
        '    axb.text(0.02, 0, "Interval unavailable", va="center", fontsize=7)\n    axb.set_xlim(0, 1.25)')
replace('The edgeR gene interval is a project approximate Fisher-t construction, not inversion of glmQLFTest.',
        'Unsupported historical edgeR intervals are not displayed; no replacement interval is inferred.')
replace('Panel B displays the median and interquartile range of stored nominal 95% gene-interval widths across all tested genes ',
        'Panel B displays the median and interquartile range of stored nominal 95% gene-interval widths for limma-voom and DESeq2 across all tested genes ')
replace('coverage validation, or equivalence assessments. The edgeR interval uses the project approximate Fisher-t "\n               "construction and is not inversion of glmQLFTest. Panel C displays',
        'coverage validation, or equivalence assessments. Unsupported historical edgeR intervals are omitted; "\n               "no replacement interval is inferred. Panel C displays')
namespace = {'__name__': '__main__', '__file__': __file__}
exec(compile(text, str(base), 'exec'), namespace)
