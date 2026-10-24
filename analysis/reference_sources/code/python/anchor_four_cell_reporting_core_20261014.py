"""Pure four-cell metric node extracted from the reviewed source (not a model fit).

The numerical loop is retained verbatim. Historical table labels are display
metadata only; actual reads use the explicitly supplied resolved input paths.
"""
import itertools
import statistics

SPECS = ("primary26_current", "primary26_original", "all27_current", "all27_original")


def compute(groups, historical_labels, metrics):
    def direction(row, group):
        return (1 if row["direction"] == "Up" else -1) if group == "Reactome_CAMERA" else metrics.sign(float(row["estimate"]))

    families, pairwise, curve, missing, features = [], [], [], [], []
    for group, paths in groups.items():
        maps = {}
        for spec, path in zip(SPECS, paths):
            data, selected = metrics.validate_table(path, group)
            maps[spec] = selected
            q = [float(r["q_value_BH"]) for r in selected.values()]
            families.append({"layer": group, "specification": spec, "tested": len(selected),
                             "q_le_005": sum(v <= .05 for v in q), "minimum_q": min(q), "native_table": path.as_posix()})
        common = sorted(set.intersection(*(set(m) for m in maps.values())))
        union = sorted(set.union(*(set(m) for m in maps.values())))
        assert len(common) > 1
        ranks = {s: metrics.average_ranks({k: float(m[k]["p_value"]) for k in common}) for s, m in maps.items()}
        for key in union:
            available = [s for s in SPECS if key in maps[s]]
            signs = {s: direction(maps[s][key], group) for s in available}
            primary = signs.get(SPECS[0])
            common_ranks = []
            for spec in SPECS:
                if key not in maps[spec]:
                    missing.append({"layer": group, "feature_id": key, "missing_specification": spec,
                                    "status": "NOT_ESTIMATED_NO_SOURCE_ROW"})
                    continue
                r = maps[spec][key]
                rank = ranks[spec].get(key)
                if rank is not None:
                    common_ranks.append(rank)
                curve.append({"layer": group, "feature_id": key, "specification": spec,
                              "native_estimate": r.get("estimate", "NA_NO_CAMERA_COEFFICIENT"),
                              "native_ci_low": r.get("ci_low", "NA_NO_CAMERA_CI"),
                              "native_ci_high": r.get("ci_high", "NA_NO_CAMERA_CI"),
                              "native_p": r["p_value"], "native_q": r["q_value_BH"], "native_direction": signs[spec],
                              "fixed_four_cell_common_n": len(common), "average_p_rank_common": rank if rank is not None else "NA_NOT_COMMON"})
            features.append({"layer": group, "feature_id": key, "available_specs": len(available),
                             "missing_specs": 4 - len(available),
                             "native_q_threshold_frequency": sum(float(maps[s][key]["q_value_BH"]) <= .05 for s in available) / len(available),
                             "native_sign_agreement_primary": sum(v == primary for v in signs.values()) / len(available) if primary is not None else "NA_NO_PRIMARY",
                             "median_fixed_common_p_rank": statistics.median(common_ranks) if common_ranks else "NA_NOT_COMMON",
                             "minimum_fixed_common_p_rank": min(common_ranks) if common_ranks else "NA_NOT_COMMON",
                             "maximum_fixed_common_p_rank": max(common_ranks) if common_ranks else "NA_NOT_COMMON",
                             "pooled_effect": "NA_NATIVE_SCALES_DIFFER", "gate_selection_frequency": "NA_EMPTY_S1_FAMILY"})
        for a, b in itertools.combinations(SPECS, 2):
            pa = {k: float(maps[a][k]["p_value"]) for k in common}
            pb = {k: float(maps[b][k]["p_value"]) for k in common}
            ta, tb = metrics.top20(pa), metrics.top20(pb)
            pairwise.append({"layer": group, "left": a, "right": b, "four_cell_common_n": len(common),
                             "Kendall_tau_b_p_ranks": metrics.tau_b(list(pa.values()), list(pb.values())),
                             "native_direction_agreement": sum(direction(maps[a][k], group) == direction(maps[b][k], group) for k in common) / len(common),
                             "top20_left_n": len(ta) if ta is not None else "NA_FEWER_THAN20",
                             "top20_right_n": len(tb) if tb is not None else "NA_FEWER_THAN20",
                             "top20_intersection_n": len(ta & tb) if ta is not None else "NA_FEWER_THAN20",
                             "top20_union_n": len(ta | tb) if ta is not None else "NA_FEWER_THAN20",
                             "top20_Jaccard": len(ta & tb) / len(ta | tb) if ta is not None else "NA_FEWER_THAN20"})

    for row in families:
        row["native_table"] = historical_labels[row["layer"]][row["specification"]]
    return {"FAMILY_SUMMARY.tsv": families, "PAIRWISE_METRICS.tsv": pairwise,
            "CURVE_LONG.tsv": curve, "MISSING_COVERAGE.tsv": missing,
            "FEATURE_METRICS.tsv": features}

