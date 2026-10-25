# 纯计算研究计划与论文大纲

**项目：** 芍药甘草汤与脑卒中后痉挛  
**版本：** Route B2 / Protocol draft 1.0  
**日期：** 2026-08-17  
**依据：** `纯计算路线穿透式审查与重构_2026-08-17.md`  
**当前状态：** 方案已完成方法学重构；进入 protocol lock、数据账本与可复现分析阶段。

---

## 1. 固定研究定位

### 1.1 研究类型

本研究是：

> **a prospectively specified, phenotype-anchored, perturbation-grounded and exposure-constrained computational prioritization study**

不是疗效研究、临床验证、因果机制研究、target-engagement study 或纯网络药理学研究。

### 1.2 推荐题目

主标题：

> **Perturbation-grounded prioritization of Shaoyao–Gancao-responsive neurospinal programs relevant to post-stroke spasticity**

方法学分支通过 benchmark 后才可考虑：

> **Phenotype-anchored and perturbation-calibrated prioritization of intervention-responsive programs under sparse disease-specific evidence**

副标题或摘要中将 SGD/JGT–PSS 说明为 case study；未通过 benchmark 时不得把流程包装为新方法。

### 1.3 核心研究问题

> Which molecular programs associated with progressive experimental spasticity recur during chronic post-stroke spinal remodeling, and which of these programs show directionally coherent modulation by a chemically characterized Shaoyao–Gancao-related whole-formula perturbome and by an independent paeoniflorin–liquiritin intervention in injured spinal cord?

### 1.4 证据表述边界

使用：

- *associated with*；
- *recurred across models*；
- *directionally coherent modulation*；
- *perturbation-grounded prioritization*；
- *candidate program/regulator*；
- *preclinical computational hypothesis*。

禁止：

- *SGD alleviates/treats post-stroke spasticity*；
- *elucidated/confirmed/validated mechanism*；
- *causal target*；
- *brain–muscle synchronization/coupling*；
- *spatial shunting*；
- *virtual knockout validation*；
- *MD-confirmed binding*；
- *synergistic action*，除非未来获得真实组合效应数据。

---

## 2. 研究域与分析优先级

### Primary domain：神经性高反射相关 neurospinal program

证据链：

1. GSE16710：进行性实验性痉挛；
2. GSE128623：慢性卒中后颈髓时间过程；
3. GSE277064：卒中后颈髓 V2a/astrocyte/microglia 细胞类型与时相；
4. GSE289929：1:1 hot-water CEM JGT 的 PC12 真实剂量扰动；
5. GSE205906：受损脊髓中的 PF–liquiritin 真实干预；
6. GSE85871：同平台多成分经验零分布。

### Secondary domain：被动肌肉/ECM 重塑

证据链：

1. GSE31243：跨病种慢性痉挛/挛缩肌肉背景；
2. PRJNA630495：卒中后急性骨骼肌/NMJ/ECM 背景；
3. GSE227494：1:1 水提全方的 C2C12 扰动；
4. GSE227495：地塞米松肌萎缩中口服 whole-formula、两剂量、两块肌肉真实扰动。

此域若不能形成独立复现，只进入补充材料或未来研究，不影响 primary domain 的完成。神经域和肌肉域不得计算一个联合“耦合分数”。

### Boundary domain：随意运动恢复

GSE270766/GSE168934 只用于说明候选是否同时随真实运动恢复干预改变，以及区分“抗痉挛相关”与“一般恢复相关”。不作为 PSS 验证。

---

## 3. Work packages 与审查门

## WP0：Protocol lock、证据本体和项目治理

### 任务

1. 冻结 primary/secondary questions、estimands、datasets、contrasts、gene-set families 和 stop rules；
2. 建立 PSS 三域 ontology：neural hyperreflexia、voluntary motor control、passive tissue remodeling；
3. 建立逐样本 provenance ledger 和 evidence-source dependency matrix；
4. 建立 formula identity/PK/direct-assay ledger；
5. 在 OSF 或等价平台公开归档带时间戳的分析方案、实际锁定时间和分析时间线；由于结果查看已经开始，不得追溯性声称 preregistration；
6. 固定代码环境、容器、随机种子和目录结构。

### Gate 0

- 所有核心 accession、GSM/BioSample/SRA、对照、时间、剂量和统计单位可追溯；
- GSE298414 的 PC12 绝对剂量冲突明确标记；
- 研究标题、摘要目标和术语不再承诺机制或疗效。

### 产出

- `protocol.md`；
- `statistical_analysis_plan.md`；
- `datasets.tsv`；
- `formula_identity.tsv`；
- `evidence_dependency.tsv`；
- 可复现环境锁文件。

---

## WP1：数据获取、逐样本 QC 与设计恢复

### 任务

1. 优先下载官方 processed counts/matrices 和完整 sample metadata；
2. 对核心数据恢复原研究设计、独立动物/培养重复、批次和配对关系；
3. 完成 library size、mapping、PCA、sample correlation、outlier、sex/组织/时间一致性检查；
4. processed-count pilot 通过后，只重处理核心条件的 FASTQ；
5. 生成 dataset card，记录可做和不可做的 contrast。

### 关键样本规则

#### JGT 1:1、hot-water、CEM 主扰动

- **PC12 / GSE289929：** vehicle GSM8801731–GSM8801733；water JGT/JGW03 GSM8801734–GSM8801742；官方逐 GSM 剂量为 100/20/4 μg/mL；
- **C2C12 / GSE227494：** vehicle GSM7101313–GSM7101315；1:1 water JGT GSM7101316–GSM7101324；剂量按官方 metadata；
- **HepG2 / GSE299063 batch A：** matched vehicle GSM9032224/GSM9032226/GSM9032228；JGW03 GSM9032099–GSM9032101、GSM9032105–GSM9032107、GSM9032111–GSM9032113；HepG2 仅作代谢/组织不相关对照。

不得跨 accession 借用对照。513-profile 补充表不含 GSE227494 的 6 个 control 和 GSE295069 的 3 个 control；六个 GEO accession 实际合计 522 个 GSM。样本账本以 GSM/BioSample/SRA 为准，而不是论文总数或 accession 数。

论文、补充表和旧 GEO 对 vehicle 的 DMSO 终浓度仍有 0.02%/0.1% 表述不一致；主分析只使用实验内、同 accession、同批次且由提交者标记的 vehicle，不自行跨批次统一或借用 non-treated。

#### GSE298414

- 新 PC12 大系列不包含 canonical JGW03 CEM 1:1；
- 论文方法写 100/20/4 μg/mL，GEO/补充表写 500/100/20 μg/mL；
- 未获作者澄清前，只允许 high/mid/low ordinal sensitivity，不作绝对剂量、效力或临床暴露解释。

### Gate 1

- 设计和分析单位无法恢复的数据从推断分析移除；
- pooled library、单切片或细胞级伪重复数据降为描述性 reference；
- 核心样本未出现无法解释的 metadata/表达冲突。

---

## WP2：痉挛表型锚与慢性卒中脊髓复现

### 2.1 GSE16710 主分析

模型：

```text
expression ~ injury + time + injury:time
```

主 contrast：

```text
(Spi60 − Spi21) − (Sham60 − Sham21)
```

敏感性：

- Spi21 vs Sham21；
- Spi60 vs Sham60；
- leave-21d-out / leave-60d-out；
- 旧芯片 probe annotation 与当前一对一映射敏感性；
- 样本 influence diagnostics。

结果单位：gene effect、Reactome/SynGO program effect、PROGENy/DoRothEA activity；不以 DEG cutoff Venn 为主要结果。

### 2.2 GSE128623

主比较：

- chronic：D28/D42 vs control；
- phase control：D4/D7 vs control；
- D14 作为过渡时相。

目的：识别不只属于急性泛炎症、而在慢性卒中后颈髓持续的程序。

### 2.3 GSE277064

- V2a interneurons、astrocytes、microglia 分开建模；
- control/D3/D7/D14/D28 分开估计；
- D14 运动训练 2/4/8 h 属另一问题，不能与自然时间过程混合；
- n=3 的结果重视 effect size、CI、rank stability，不用夸张 gene-level discovery；
- 候选模块在 D14/D28 的细胞类型定位为 held-out 验证。

### Gate 2：疾病程序分级

- **Grade S1：** GSE16710 表型锚方向稳定；
- **Grade S2：** S1＋GSE128623 慢性同向，且不是仅急性炎症；
- **Grade S3：** S2＋GSE277064 至少一个预设细胞类型/慢性时相同向。

只有 S2–S3 进入 whole-formula 主检验。允许没有候选通过。

---

## WP3：全方扰动、制备因素与组织背景

### 3.1 主检验

每个 cell line 独立分析：

```text
expression/module_score ~ vehicle + ordered_dose + batch
```

主输出：

- high dose vs matched vehicle；
- log-dose 或 ordered-dose trend；
- 三重复内方向稳定性；
- 低 RIN 排除敏感性；
- 对 Gate 2 锁定模块的 disease-direction reversal/compatibility。

“reversal”必须写成 transcriptomic direction reversal，不等同疗效。

### 3.2 配方与提取方式的内部分析

在元数据完整的 C2C12/HepG2 中：

- ratios：2:1、1:1、1:2；
- solvent：water vs 70% ethanol；
- extraction：CEM vs IEM；
- single herbs：3:0/0:3；
- cell line 分层，不直接合并表达矩阵。

使用 hierarchical/factorial model 估计 preparation-dependent response，但同一资源的不同条件只算内部重复/稳健性。

### 3.3 非加和响应

如比较 CEM 与两个单味响应的加权线性期望，表述为：

> *departure from a linear transcriptomic expectation* 或 *non-additive transcriptomic response*

不得表述为 synergy，因为化学提取率、细胞暴露和剂量非线性没有被分离，也没有 Bliss/Loewe/HSA 真实组合效应。

### Gate 3：全方响应分级

- **J1：** PC12 1:1 water CEM 高剂量方向与疾病程序相反或相容；
- **J2：** J1＋剂量顺序稳定；
- **J3：** J2＋C2C12/HepG2 背景显示可解释的组织特异性，而不是三细胞普遍应激；
- **J4：** J2–J3 在比例/提取方式敏感性中稳健。

仅 J2 以上进入独立受损组织校准。

---

## WP4：独立真实扰动与经验特异性

### 4.1 GSE205906

预设 contrasts：

- disease: SNI vs sham；
- intervention: SNI＋PF–LQ vs SNI。

只检验 WP2–WP3 锁定的模块，避免同一数据发现再验证。报告每个模块的方向、效应量、CI/方向概率、rank stability；样本 n=3 时不把 nominal P 当确证。

GEO 个别样本存在人 GRCh38 的错误 genome 标注风险；最终结果不直接接受该处理矩阵，须从 raw reads 重新比对至大鼠 mRatBN7.2，并在 dataset card 中记录核查过程。

### 4.2 GSE85871

预先指定：

- paeoniflorin；
- liquiritin；
- glycyrrhizic acid。

以其余同平台成分建立 empirical null，并对 signature magnitude 匹配。报告：

- empirical percentile；
- empirical P；
- bootstrap rank；
- top-5%/top-10% stability；
- 随机方向/随机同规模基因集 null。

该数据只回答 specificity，不回答组织相关性或 PSS 疗效。

### 4.3 肌肉药理桥

GSE227495：

- gastrocnemius 与 soleus 分开并在 animal-level 建模；
- disease vs normal；
- JGT03 1:1 的 100/400 mg/kg vs DEX vehicle；
- ratio 2:1/1:1/1:2 与 dose interaction；
- resveratrol 作为原研究 positive comparator，而非 PSS 阳性对照。

可选 reference 仅用于注释：GSE167186/Human Muscle Ageing Cell Atlas 用于人肌肉/FAP 细胞定位，GSE243077/NeMO 用于健康人脊髓层板/细胞类型定位；GSE309478 可作跨情境 FAP/ECM 真实扰动参考。上述资源均不计为 PSS 或 SGD 疾病验证。

### Gate 4：最终证据级别

- **E1：** S2/S3＋J2；
- **E2：** E1＋GSE205906 独立受损脊髓方向支持；
- **E3：** E2＋GSE85871 经验特异性；
- **E4：** E3＋人体/动物真实暴露与直接 assay 证据。

主文候选程序原则上只报告 E2–E4；E1 作为探索性。0 个 E2 候选时不强行选靶点，论文转为负结果/证据边界资源或停止结果型投稿。

---

## WP5：暴露、制剂身份和靶点证据

### 5.1 Formula identity ledger

每个干预必须记录：

- botanical species、药用部位、炮制；
- ratio、CEM/IEM、water/ethanol；
- extraction yield、dose、vehicle、treatment time；
- HPLC/LC-MS marker 和批次；
- JGT、SGD、Shakuyakukanzoto、total glucosides、PF–LQ combination 的不可互换性。

### 5.2 Exposure ledger

人体口服血浆优先记录：

- albiflorin；
- paeoniflorin；
- glycycoumarin；
- isoliquiritigenin；
- glycyrrhetic acid；
- glycyrrhetic acid 3-O-monoglucuronide。

为每个实体保存 exact name、InChIKey、parent/metabolite、total/unbound Cmax、time、matrix、dose/formulation 和文献来源。未测脑/脊髓/肌肉暴露者明确标注 unknown。

### 5.3 Target ledger

靶点来源分级：

1. 直接 biochemical/biophysical binding；
2. 直接 functional assay；
3. curated high-confidence target；
4. expression/perturbation association；
5. database/ML prediction。

只有 1–2 级且与 E2–E3 程序一致的 pair 才能称 direct target evidence。预测项单独保留，不与 direct assay 合并计分。

### Gate 5

- compound 有实测体内暴露；
- target 有直接 assay；
- concentration plausibility 不明显违背 unbound exposure；
- formula identity 可追踪；
- 最多 0–3 compounds、0–2 targets。

无 pair 通过时，结构模块取消。

---

## WP6：门控式可选分析

### 6.1 空间转录组

仅在 E2 候选提出明确的皮层 lesion-border/white-matter hypothesis 时启动 PRJNA952594/PXD041388：

- animal/slide 为 n；
- predefined regions；
- cell2location＋RCTD 双方法；
- tissue histology/held-out markers；
- mixed model 或 spatial block bootstrap；
- 同一实验的空间和蛋白只算一个来源。

否则删除空间模块。

### 6.2 Virtual perturbation

默认不做。只有 TF 候选、多个动物、轨迹、scATAC/motif prior 和 context-matched 真实扰动均存在时，才允许 CellOracle/scTenifoldKnk 的补充 benchmark；必须优于 no-change、mean、linear/additive baseline。失败则完全删除。

### 6.3 Docking/MD

仅在 Gate 5 通过后，对最多 1 个主 pair 和 1 个备选 pair作结构解释：

- experimental holo PDB first；
- native-ligand redocking；
- known active/inactive/decoy；
- second docking engine/cross-docking；
- ≥3 independent MD replicas、收敛与 uncertainty；
- 只写 structural plausibility。

当前项目预算与算力不为此预留，除非上游结果触发。

---

## WP7：稳健性、负对照与复现

### 必做

- leave-one-study-out；
- leave-one-species-out；
- leave-one-timepoint-out；
- edgeR QL vs DESeq2/limma-voom；
- Reactome/SynGO vs PROGENy/DoRothEA；
- one-to-one ortholog vs pathway-only transfer；
- low-RIN exclusion；
- influence diagnostics；
- matched random gene sets；
- degree/expression/gene-set-size matched nulls；
- sample-label permutations；
- JGT hot-water CEM primary vs all-preparation secondary；
- GSE298414 excluded vs ordinal-only sensitivity；
- candidate selection frequency、Kendall tau、top-k Jaccard、specification curve。

### 证据独立性审计

每个结论必须标记：

- 是否共享 GSM/BioSample；
- 是否来自同一实验、动物或测序批次；
- 是否使用相同 gene-set/interaction database；
- 是否是同一模型输出的下游加工；
- 是否真正构成外部验证。

---

## 4. 统计分析计划摘要

### 4.1 主要方法

| 数据类型 | 主方法 | 敏感性 | 推断单位 |
|---|---|---|---|
| GSE16710 microarray | limma full factorial | annotation/remapping sensitivity | animal |
| bulk RNA-seq | edgeR quasi-likelihood | DESeq2/limma-voom | animal/independent culture |
| RiboTag | cell-type-stratified edgeR QL | shrinkage/rank stability | animal/library |
| sc/snRNA（条件性） | pseudobulk edgeR/muscat/dreamlet | mixed model | animal/donor |
| spatial（条件性） | region-level mixed model | spatial block bootstrap | animal/slide |
| pathway | camera/roast as primary | fgseaMultilevel | predefined family |
| regulator activity | decoupleR/PROGENy/DoRothEA | alternate regulon sources | sample/contrast |
| cross-study | direction/effect-size evidence matrix | multilevel meta-analysis only if commensurate | study |
| empirical specificity | matched empirical null/permutation | bootstrap | compound signature |

### 4.2 报告标准

- effect size、95% CI、exact n、BH FDR；
- 对小 n 结果报告不确定性，不以 p<0.05 二元化；
- 连续剂量不无依据二分；
- 每个多重检验 family 预先定义；
- 不把 cell/spot/frame 作为 n；
- 不把 SHAP、feature importance、GRN edge 或 ligand–receptor score 写成因果。

### 4.3 跨物种

- Ensembl 一对一直系同源映射为基因层主方案；
- pathway 层为主要跨物种比较；
- 记录映射率、缺失与一对多敏感性；
- 不把跨物种一致性写成人类 causal transport。

---

## 5. 可复现工程结构

```text
project/
├── README.md
├── protocol/
│   ├── protocol.md
│   ├── statistical_analysis_plan.md
│   └── deviations.md
├── metadata/
│   ├── datasets.tsv
│   ├── samples.tsv
│   ├── formula_identity.tsv
│   ├── exposure.tsv
│   ├── targets.tsv
│   └── evidence_dependency.tsv
├── config/
│   ├── contrasts.yaml
│   ├── gene_sets.yaml
│   └── thresholds.yaml
├── workflow/
│   ├── Snakefile or main.nf
│   └── envs/
├── R/
├── python/
├── data/
│   ├── raw/
│   ├── interim/
│   └── derived/
├── results/
│   ├── qc/
│   ├── effects/
│   ├── validation/
│   └── figures/
├── manuscript/
└── containers/
```

原始数据只读；所有中间与最终结果可由 workflow 重建。数据下载保存 checksum、URL、日期和 license。代码环境使用 renv＋conda-lock/uv lock，容器固定 OS 和主要软件版本。

---

## 6. 论文 Results 逻辑与图表计划

所有 Results 标题为预设结构，不预填“显著”“证实”等结果。

### Figure 1. Phenotype ontology, evidence graph and study design

- PSS 三域；
- 数据集—表型—组织—干预—证据等级；
- evidence-source dependency；
- claim ceiling。

### Figure 2. Progressive experimental-spasticity programs and chronic post-stroke spinal recurrence

- GSE16710 interaction effects；
- GSE128623 time-course effects；
- cross-model pathway direction and uncertainty；
- acute-stage negative control。

### Figure 3. Cell-type and temporal localization in post-stroke cervical spinal cord

- GSE277064 V2a/astrocyte/microglia；
- D3–D28 trajectory；
- locked module effects and CIs；
- no cell-as-replicate statistics。

### Figure 4. Whole-formula perturbation across dose, preparation and cellular context

- canonical PC12/C2C12 1:1 water CEM；
- ordered-dose response；
- ratio/solvent/extraction sensitivity；
- HepG2 negative/context control；
- metadata conflict panel and exclusion decision。

### Figure 5. Independent injured-spinal-cord perturbation and empirical compound specificity

- GSE205906 disease and intervention contrasts；
- candidate module direction；
- GSE85871 empirical null/percentile；
- positive/negative control diagnostics。

### Figure 6. Robustness, evidence grades and falsifiable hypotheses

- LOSO/multiverse specification curve；
- evidence grade E1–E4；
- dependency-aware evidence matrix；
- 0–3 final hypotheses or explicit zero-result outcome。

### Optional/Supplementary muscle domain

- GSE31243/PRJNA630495 background；
- GSE227494/GSE227495 whole-formula response；
- gastrocnemius vs soleus and dose/ratio effects；
- strictly separate from neurospinal main inference。

### Tables

1. Dataset and sample-level provenance；
2. Formula identity, preparation and exposure evidence；
3. Prespecified estimands, contrasts and statistical models；
4. Candidate program evidence matrix and claim level；
5. Critical exclusions, metadata conflicts and protocol deviations。

---

## 7. 英文论文大纲

## Abstract

结构：Background → Gap → Objective → Design/Data → Main analyses → Results（完成分析后填写）→ Evidence-bounded conclusion。

摘要禁止使用“confirmed mechanism”“therapeutic efficacy”“validated target”。

## Introduction

### Paragraph 1: Clinical and mechanistic heterogeneity of post-stroke spasticity

- neural hyperreflexia、voluntary motor control、passive tissue remodeling；
- MAS/一般运动恢复不等于痉挛机制；
- direct PSS molecular data scarcity。

### Paragraph 2: Why Shaoyao–Gancao requires intervention-identity and exposure constraints

- whole decoction、total glucosides、JGT/Kampo、single compounds 不可互换；
- 血浆暴露与靶组织暴露区别；
- 传统网络药理学的证据循环问题。

### Paragraph 3: Opportunity created by real perturbation datasets

- 2026 513-profile JGT resource；
- injured spinal PF–LQ；
- post-stroke cervical-spinal time/cell resources；
- need for dependency-aware evidence gating。

### Paragraph 4: Objective and claim boundary

- 明确 association/prioritization；
- 主要 neurospinal、次要 muscle domain；
- preregistration, holdout and negative controls。

## Methods

1. Study design, preregistration and evidence ontology  
2. Systematic dataset search and eligibility  
3. Sample-level provenance and dependency audit  
4. Formula identity, chemistry and exposure curation  
5. Preprocessing and quality control  
6. Experimental-spasticity anchor analysis  
7. Post-stroke spinal time-course and cell-type analysis  
8. Whole-formula perturbome analysis  
9. Injured-spinal-cord intervention analysis  
10. Compound-specific empirical null analysis  
11. Secondary muscle-domain analysis  
12. Cross-species mapping and program-level integration  
13. Robustness, multiverse and negative controls  
14. Evidence grading and stopping rules  
15. Code, data and ethics/governance。

## Results

1. Dataset audit altered the admissible evidence architecture  
2. Programs associated with experimental spasticity progression  
3. Recurrence and cellular localization during chronic post-stroke spinal remodeling  
4. Context- and preparation-dependent whole-formula perturbation  
5. Concordance in an independent injured-spinal-cord intervention  
6. Empirical specificity against unrelated herbal-compound signatures  
7. Robustness, evidence grades and retained/terminated hypotheses  
8. Secondary muscle-domain findings（仅通过 gate 时）。

## Discussion

1. Principal findings at their correct evidence level；
2. Biological interpretation by PSS domain；
3. What real perturbation adds beyond network pharmacology；
4. Formula identity, preparation and exposure implications；
5. Generalizability to sparse-phenotype complex interventions；
6. Limitations：no direct PSS＋SGD cohort、cross-disease transport、cell lines、small n、metadata conflicts、no tissue PK/target engagement；
7. Falsifiable predictions and future standardized clinical PK/PD trial。

## Conclusion

仅概括可复现程序和真实扰动的一致性，不宣称疗效或确证机制。

---

## 8. 目标期刊分支

### Branch A：生物学/系统药理主文

触发条件：至少一个 E2–E3 程序、跨数据稳健、真实全方和受损脊髓方向收敛。

候选：npj Systems Biology and Applications；强结果可尝试 Communications Biology 或 iScience。

### Branch B：方法/benchmark 主文

触发条件：通用框架在预注册 benchmark 中明显优于传统 network-pharmacology、exposure-agnostic 和简单 transcript-reversal baselines，并开放软件。

候选：PLOS Computational Biology、Bioinformatics、Briefings in Bioinformatics。

### Branch C：负结果/证据地图

触发条件：无 E2 候选，但数据审计、metadata、证据依赖和可重复资源具有独立价值。

改写为 evidence map/resource/protocol paper；不降低阈值制造阳性结果。

---

## 9. 里程碑与完成判据

| 里程碑 | 完成产物 | Go/No-go |
|---|---|---|
| M0 Protocol lock | protocol、SAP、术语/claim 表 | 关键 metadata 未解决则不分析 |
| M1 Provenance/QC | 样本账本、dataset cards、QC report | 伪重复/pooled 数据降级 |
| M2 Disease programs | GSE16710＋GSE128623＋GSE277064 effects | 无 S2 则停止药物机制叙事 |
| M3 Formula perturbome | canonical JGT effects＋prep sensitivity | 无 J2 则不称 formula-responsive program |
| M4 Independent calibration | GSE205906＋GSE85871＋muscle secondary | 无 E2 则不选 target |
| M5 Robustness | LOSO、multiverse、nulls、dependency audit | 不稳定候选删除 |
| M6 Target/exposure gate | exposure/direct assay ledger | 无 direct pair 则取消结构模块 |
| M7 Manuscript | 英文初稿、图表、补充、代码与数据声明 | 独立方法学审稿通过后投稿 |

---

## 10. 当前正式执行决策

1. 原路线 B 终止；Route B2 生效；
2. GNG12、PLOD3、6NPR、伪 sham、RRHO 同步、LogBB 空间分流和默认 MD 全部退出主链；
3. GSE16710、GSE128623、GSE277064、GSE289929、GSE205906 为 primary evidence chain；
4. GSE227495 为重要外周真实药理桥；
5. spatial/virtual KO/docking/MD 均为结果触发式可选模块；
6. 下一步从 WP0–WP1 开始，先形成可审计项目骨架、样本账本和 processed-count pilot。
