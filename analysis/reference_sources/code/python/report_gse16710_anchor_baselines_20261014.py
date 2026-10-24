"""Generate a bounded source-linked report; no scientific model is run here."""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--verification", default=None)
    args = parser.parse_args()
    root, run, out = map(Path, (args.project, args.run, args.out))
    out.mkdir(parents=True, exist_ok=True)
    target = out / "SCIENCE_REPORT.md"
    if target.exists():
        raise FileExistsError(target)
    s = json.loads((run / "RUN_SUMMARY.json").read_text(encoding="utf-8"))
    assert s["status"] == "COMPUTED_NUMERIC_CHECKS_PASS"
    before, after = tsv(run / "INPUT_MANIFEST.tsv"), tsv(run / "INPUT_MANIFEST_AFTER.tsv")
    assert before == after
    for item in before:
        path = root / item["relative_path"]
        assert path.stat().st_size == int(item["bytes"]) and digest(path) == item["sha256"]
    for item in s["outputs"]:
        path = run / item["file"]
        assert path.stat().st_size == item["bytes"] and digest(path) == item["sha256"]
    collections = s["collections"]
    assert sum(v["n"] for v in collections.values()) == s["features"] == 1771
    verified = args.verification is not None
    verification = Path(args.verification) if verified else None
    if verified:
        assert verification.is_file()
        verified_receipt = json.loads(verification.read_text(encoding="utf-8"))
        assert verified_receipt["status"] == "PASS_INDEPENDENT_CACHED_SCORE_NUMERIC_CHECK"
        assert verified_receipt["source_and_product_pins_unchanged"] is True
        assert verified_receipt["features"] == s["features"]
        assert verified_receipt["run_summary"]["sha256"] == digest(run / "RUN_SUMMARY.json")
    claims = []

    def claim(label, value, location, source="RUN_SUMMARY.json"):
        claims.append({"claim": label, "reported_value": str(value),
                       "source_file": str(run / source), "source_location": location,
                       "verification_status": "VERIFIED" if verified else "NEEDS REVIEW",
                       "source_readback_status": "MATCH",
                       "discrepancy": "None", "required_action":
                       "Independent numerical verification separately bound" if verified else "Await independent numerical verification",
                       "evidence_status": "post-outcome descriptive implementation; not a new confirmatory gate"})
        return value

    claim("Retained arrays", s["sample_n"], "sample_n")
    claim("Complete program rows", s["features"], "features")
    claim("Full/reduced ranks", f'{s["rank_full"]}/{s["rank_reduced"]}', "rank_full,rank_reduced")
    claim("OLS residual degrees of freedom", s["residual_df"], "residual_df")
    claim("Synthetic fixtures passed", s["synthetic_tests"]["fixtures"], "synthetic_tests.fixtures")
    lines = ["# GSE16710：原评分与无交互、普通线性及秩基线的实际比较", "",
             "2026-10-14。源文件数值读回核验通过；" +
             ("独立数值核验已绑定，见下方回执。" if verified else "**独立数值核验尚待完成，本页为计算后草稿。**"), "",
             "这次实际执行补齐锚数据集评分层面的基线比较，不重新生成评分、不重跑原主模型、不改变任何门控。全部程序均保留，没有新增剔除、插补或评分修正。"
             "同一输入上的系数一致是代数兼容检查，不是独立重复验证。研究总目标仍在推进，不能据此宣布整个 SAP 或项目完成。", "",
             "## 输入、执行及验证范围", "",
             f'- 原评分：Reactome/PROGENy/DoRothEA 三个缓存 RDS；共 {s["features"]:,} 行、{s["sample_n"]} 个保留阵列。',
             "- 四组为 Sham21=6、Sham60=4、Spi21=8、Spi60=8；沿用既有 QC 排除 GSM418682。没有本次新删除。原 27 行账本与 YAML 均保留。",
             "- 原交互对比为 `(Spi60−Spi21)−(Sham60−Sham21)`。完整/受限设计秩为 4/3，普通 OLS 残差自由度为 22。",
             f'- {s["synthetic_tests"]["fixtures"]} 项合成夹具通过，实际进程退出码 0。源文件在执行前后 SHA256 相同；本报告生成时再次读回核对。',
             "- r01 在预检矩阵附带属性比较处停止；确认元素、样本和列名完全相同后修正比较实现。r02/r03 为代码审查前的预检，均无真实基线效应计算；只有 r04 执行了本次真实基线。旧尝试保留。",
             "- 原始 CEL、旧失败错误正文、其他数据集模型均不在本次读取/运行范围。旧来源清单仅作版本账本，不代表本次重新验证原始 CEL。", "",
             "## 普通 OLS 与保留的 moderated 评分模型", "",
             "表中 q 是各自方法、各自完整家族内 BH 结果，不跨方法取并集。普通 OLS 的 P 与逐点 95% CI 以独立、同方差高斯误差为工作假设；本次不能证明这些假设成立。", "",
             "| 家族 | 完整程序数 | OLS/保留模型 q≤0.05 数 | P 值排序 Spearman | 前 min(20,n) 交集 |",
             "|---|---:|---:|---:|---:|"]
    for fam, c in collections.items():
        m = c["moderated_vs_OLS"]
        lines.append(f'| {fam} | {c["n"]} | {c["OLS_nominal_q05_n"]}/{c["retained_moderated_q05_n"]} | '
                     f'{m["rho"]:.6f} | {m["top_overlap"]}/{m["top_k"]} |')
        for key in ("n", "OLS_nominal_q05_n", "retained_moderated_q05_n", "q05_overlap_n",
                    "native_effect_compared_n", "ordinary_t_compared_n", "native_effect_max_abs_difference",
                    "ordinary_t_max_abs_difference", "f_t_squared_max_abs_difference"):
            claim(f"{fam} {key}", c[key], f"collections.{fam}.{key}")
        for key in ("rho", "paired_n", "top_k", "top_left", "top_right", "top_overlap"):
            claim(f"{fam} moderated_vs_OLS {key}", m[key], f"collections.{fam}.moderated_vs_OLS.{key}")
    lines.extend(["", "全部原量纲点估计均在预定的 `1e-8*(1+abs(reference))` 容差内一致；机器精度差异见运行摘要，不作为新的科学精度或独立复制主张。"
                  "每行的原量纲效应、SE、名义逐点 95% CI、P、家族内 BH q、四组均值及保留 moderated 结果均在完整结果表，不以只列显著条目替代效应与不确定性报告。"
                  "PROGENy 只有 14 个程序，前 min(20,n) 等于全家族，因此其 14/14 交集没有筛选稳定性证据。"
                  "高排序相关也不证明误差校准、同方差、动物独立性或生物学真实性；零发现不能解释为真实效应为零。", "",
                  "## 无交互与秩比较", "",
                  "无交互模型保留 injury/time 主效应，只将交互固定为零；它不是生物学无变化、健康或无损伤。受限模型的四组拟合均值和 SSE 已逐程序输出。"
                  "完整模型 SSE 的降低是样本内嵌套拟合的结构性结果，不是预测优势。F(1,22)=普通 t² 检查通过，它们是同一检验的两个形式，不能当作两份独立 null 证据。固定零没有人为 P/q/CI。", "",
                  "每个程序在同一 26 个阵列内作平均并列秩转换，再比较四组秩均值；这改变了 estimand，不能解释为原评分相同效应。没有输出秩 P/q/CI，也没有借用原评分 Freedman–Lane 结果为秩检验校准。", "",
                  "| 家族 | 原量纲/秩效应同向（非零配对） | 秩或原效应为零 | 缺失配对 | 绝对效应排序 Spearman |",
                  "|---|---:|---:|---:|---:|"])
    for fam, c in collections.items():
        lines.append(f'| {fam} | {c["rank_effect_direction_agree_n"]}/{c["rank_effect_paired_n"]} | '
                     f'{c["rank_effect_zero_n"]} | {c["rank_effect_missing_n"]} | {c["absolute_native_rank_effect_rho"]:.6f} |')
        for key in ("rank_effect_direction_agree_n", "rank_effect_paired_n", "rank_effect_zero_n",
                    "rank_effect_missing_n", "absolute_native_rank_effect_rho", "tied_score_feature_n"):
            claim(f"{fam} {key}", c[key], f"collections.{fam}.{key}")
    cm = collections["Reactome"]["camera_vs_OLS"]
    lines.extend(["", "没有原评分并列值；Reactome 的 3 个零秩对比来自组平均秩差为零，不是删去或插补的程序。方向比较的分母排除了精确零对比，完整结果仍保留这些行。", "",
                  "## CAMERA 不是同一个 null", "",
                  f'Reactome 的 {cm["paired_n"]} 个键全部匹配，无缺失/额外键。CAMERA 与普通评分 OLS 的 P 值排序 Spearman={cm["rho"]:.6f}，'
                  f'前 {cm["top_k"]} 交集为 {cm["top_overlap"]}/{cm["top_k"]}，方向同向为 {cm["direction_agree_n"]}/{cm["direction_paired_n"]}。'
                  "这是 competitive probe-level 与 sample-score 模型的对照，不能把相关程度当作同一估计量的互相验证，也不补造 CAMERA 系数或 CI。", "",
                  "## 方差与影响诊断", "",
                  "| 家族 | 四组样本方差最大/最小比中位数 | 同比最大值 | 最大 Cook 距离 |",
                  "|---|---:|---:|---:|"])
    for key, value in cm.items():
        if key != "null":
            claim(f"Reactome CAMERA vs OLS {key}", value, f"collections.Reactome.camera_vs_OLS.{key}")
    for fam, c in collections.items():
        v, k = c["variance_ratio_quantiles"], c["max_cook_quantiles"]
        lines.append(f'| {fam} | {v["50%"]:.3f} | {v["100%"]:.3f} | {k["100%"]:.3f} |')
        for metric, field, percentile in (("variance ratio median", "variance_ratio_quantiles", "50%"),
                                           ("variance ratio maximum", "variance_ratio_quantiles", "100%"),
                                           ("maximum Cook distance", "max_cook_quantiles", "100%")):
            claim(f"{fam} {metric}", c[field][percentile], f"collections.{fam}.{field}.{percentile}")
    lines.extend(["", "观察到组间样本方差不齐和影响程度差异，不能把名义 t 区间/P/BH 值升级为已经证明校准的推断。"
                  "小组样本方差本身也不稳定；这里没有以任意阈值判定错误、排除阵列或自动更换模型。完整表保留各组方差、最大学生化残差、Cook 距离及对应阵列 ID。", "",
                  "## 内容验证与仍待核实的元数据", "",
                  "Validation status：**NEEDS REVIEW（历史单位/批次边界）**；本次矩阵结构和代码数值检查通过。没有已确认的关键内容错误，没有新增删除。", "",
                  "核验前后原 metadata 均为 27 行；沿用既定 QC 的评分分析集仍为 26 阵列。本次核验删除原记录 0 行，评分行保留 1,771/1,771。", "",
                  "| 域 | 结果 | 类别/严重度 | 处理 |", "|---|---|---|---|",
                  "| 样本/特征/设计 | 27 行源账本，原 QC 后 26 阵列，1,771 完整评分行；无键重复、错配、非有限评分或本次无效推断行 | informational | 保留所有原记录；结果表不筛选 |",
                  "| biological_unit_id | 记录为 GSM 阵列 ID，不能独立认证动物身份或排除伪重复 | needs review / warning | 不宣称已认证的独立动物 |",
                  "| batch、biosample、SRA | 账本缺失，不能查证未记录批次或生物学独立性 | needs review / warning | 不猜测或自动批次校正 |",
                  "| QC 时间顺序 | outcome_blind 标签是保留文件的记载，不是独立历史认证 | needs review / warning | 保留来源限定 |",
                  "| 统计极端值 | 方差和影响诊断已报告，不自动归为 data error | possible outlier / needs review | 不删除、不截尾、不改变主模型 |", "",
                  "## Robustness 与 claim 边界", "",
                  "| 比较 | 实际结果 | 解释 | 对主结论的影响 |", "|---|---|---|---|",
                  "| 无交互参照 | 完整/受限拟合和代数检查完成 | 约束交互为零，不是独立生物学阴性对照 | 不救回 S1，不证明无效 |",
                  "| 普通/保留 moderated 评分模型 | 原量纲系数一致，排序较近，各自完整家族 q≤0.05 均为零 | 同一输入的方差处理比较，名义假设仍有限 | 无新增 admitted program；原门控不变 |",
                  "| 秩/原量纲 | 大部分但非全部非零效应同向，尺度与 estimand 不同 | 仅描述性；没有经过校准的秩交互推断 | 不能用作新的显著性支持 |",
                  "| CAMERA/评分 OLS | 匹配键下排序、前列和方向不完全相同 | competitive 与评分层面 null 不同 | 不互相替代或择优解释 |", "",
                  "| Proposed claim | Evidence class | Support | Allowed wording | Risk / additional evidence |",
                  "|---|---|---|---|---|",
                  "| 本次锚评分基线已实际计算 | computational / descriptive | SUPPORTED within scope | 同一缓存评分的三类比较已完成并保留全量结果 | 必须绑定独立数值回执；不能泛指全部 SAP |",
                  "| 程序已被独立重复验证，或复杂模型更优 | replication / prediction | UNSUPPORTED | 系数一致仅为代数兼容；无预测校准结论 | 独立样本与预设外部验证/校准尚需 |",
                  "| 没有 q≤0.05 就证明没有生物学效应 | biological absence | UNSUPPORTED | 本模型/家族门槛下未出现新统计发现 | 需要适当精度/等效性设计；不能由不显著推出 |",
                  "| 这些计算证明临床疗效或机制 | causal / mechanistic / clinical | UNSUPPORTED | 只描述本组阵列上的计算比较 | 需要独立因果、干预和临床证据 |", "",
                  "可保留计算范围内的描述。不得用独立复制、预测优越、临床有效或机制已证实措辞。本页没有此类 blocking overclaim，但模型假设和历史 provenance 不因此变为已验证。", "",
                  "## 当前完成对象与下一步", "",
                  "已完成的是 GSE16710 所有既有程序的评分层面 no-interaction/OLS/描述性 rank baseline。尚未完成其他数据集对应的源范围基线、RNA-seq 标签 null、正式秩校准、完整多规格与独立外部验证。"
                  "此处不补做或擅自批准 GSE277064 的额外 voom/DESeq2 分析，也不改变候选集合。完整研究目标继续保持 active。", "",
                  "现有 v5 Word 本轮未改动。待后续控制节点有完整结果，再将这些已核验的基线结果合并到后续版本；不把旧 controls review 覆盖成整体已完成。", "",
                  "## 可复核文件", "",
                  f'- [全量逐程序效应、CI、P/q 与诊断]({(run / "BASELINE_RESULTS.tsv").as_posix()})',
                  f'- [运行摘要]({(run / "RUN_SUMMARY.json").as_posix()})',
                  f'- [输入与 SHA256]({(run / "INPUT_MANIFEST.tsv").as_posix()})',
                  f'- [执行前注册对象]({(run / "REGISTRY.json").as_posix()})',
                  f'- [固定实现定义]({(root / "protocol/GSE16710_ANCHOR_BASELINE_COMPARISON_20261014.md").as_posix()})',
                  f'- [独立方法审查]({(root / "results/diagnostics/GSE16710_anchor_baseline_method_review_20261014/INDEPENDENT_METHOD_REVIEW.md").as_posix()})',
                  f'- [r04 静态代码审查]({(root / "results/diagnostics/GSE16710_anchor_baseline_method_review_20261014/INDEPENDENT_CODE_REVIEW_R04.md").as_posix()})'])
    if verified:
        lines.append(f'- [独立数值核验]({verification.as_posix()})；SHA256 `{digest(verification)}`。')
    lines.append(f'- [逐项数值来源与核验清单]({(out / "NUMERIC_CLAIM_INVENTORY.tsv").as_posix()})。')
    lines.extend(["", "完整评分导出和中间文件位于 R 盘；C 盘仅保留必要源代码、定义及小型报告。没有云端下载或原始数据迁移。", "",
                  "生成器从实际 JSON/TSV 读取上述数值；原量纲每程序效应与不确定性见全量表。逐项数值来源清单单独保存，不把计算后探索性细化写成原先预设的确认性结果。", ""])
    target.write_text("\n".join(lines), encoding="utf-8")
    claims_path = out / "NUMERIC_CLAIM_INVENTORY.tsv"
    with claims_path.open("w",encoding="utf-8",newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(claims[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(claims)
    receipt = {"status": "SOURCE_VALUES_VERIFIED_INDEPENDENT_REVIEW_BOUND" if verified else "DRAFT_AWAITING_INDEPENDENT_NUMERIC_REVIEW",
               "generated_utc": datetime.now(timezone.utc).isoformat(), "source_summary_sha256": digest(run / "RUN_SUMMARY.json"),
               "source_manifest_sha256": digest(run / "INPUT_MANIFEST.tsv"), "all_source_pins_read_back": True,
               "claim_rows":len(claims), "sources": len(before),
               "report": {"path": str(target), "sha256": digest(target), "bytes": target.stat().st_size},
               "claim_inventory": {"path":str(claims_path), "sha256":digest(claims_path)},
               "verification": {"path":str(verification), "sha256":digest(verification)} if verified else None,
               "scientific_provenance_status":"NEEDS_REVIEW; numerical acceptance does not authenticate animals/batches/history"}
    (out / "REPORT_BUILD_RECEIPT.json").write_text(json.dumps(receipt,indent=2),encoding="utf-8")
    print(json.dumps({"status":receipt["status"],"claim_rows":len(claims),"report":str(target)},ensure_ascii=False))


if __name__ == "__main__":
    main()
