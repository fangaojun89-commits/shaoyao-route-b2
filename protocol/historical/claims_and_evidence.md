# Route B2 claims and evidence boundary

**Status:** protocol-level claim lock, 2026-08-17  
**Applies to:** all analyses, figures, abstracts, manuscripts and presentations produced under `route_b2`.

## 1. Non-negotiable scientific scope

Route B2 is a **preclinical computational prioritization study**. It asks whether molecular programs associated with experimental spasticity recur in independent post-stroke spinal-remodelling data and whether prespecified programs show directionally coherent responses in real, but cross-context, JGT and paeoniflorin-liquiritin perturbation datasets.

There is no dataset in the current route that jointly observes a standardized clinical SGD intervention, a directly quantified human PSS phenotype, target-tissue pharmacokinetics, molecular target engagement and a clinical outcome. Consequently, there is no estimand that can be interpreted as the treatment effect of SGD on human PSS and no complete causal mechanism can be identified.

The formula and exposure ledgers are authoritative for naming:

- `../metadata/formula_identity.tsv`
- `../metadata/exposure.tsv`

JGT, Chinese SGD, Japanese Shakuyakukanzoto/TJ-68, total glucosides and the paeoniflorin-liquiritin combination are **non-interchangeable intervention identities**. A synonym must never be substituted merely because the preparations share two source herbs or selected constituents.

## 2. PSS phenotype ontology: three domains must remain separate

The current domain model follows the 2026 American Heart Association scientific statement: <https://pmc.ncbi.nlm.nih.gov/articles/PMC13063836/>.

| Domain | Operational meaning | Preferred future measurement | Claim restriction in Route B2 |
|---|---|---|---|
| Neural hyperreflexia / involuntary muscle activity | Velocity-dependent stretch-reflex hyperexcitability, clonus and other involuntary activation attributable to altered descending and spinal control | Randomized-velocity stretch-evoked EMG and torque, tonic stretch-reflex threshold (TSRT), H-reflex recruitment/rate-dependent depression as complementary physiology | Experimental-spasticity and neurospinal-program language is allowed only when anchored to an actual spasticity phenotype; general stroke, FMA, NIHSS, infarct size or MAS alone does not create a PSS molecular signature |
| Impaired voluntary motor control | Weakness, loss of selective motor control and abnormal flexor/extensor synergies | FMA/ARAT, strength and task-based kinematics/EMG, analysed separately from stretch responses | Improvement or molecular change in this domain cannot be renamed antispastic efficacy |
| Passive tissue remodelling | Muscle shortening, altered fascicle/tendon properties, extracellular-matrix remodelling, fibrosis and non-neural mechanical stiffness | Slow-stretch passive torque, range of motion, ultrasound/shear-wave elastography and prespecified structural measures | Muscle/ECM findings form a secondary domain; they cannot be pooled with neural findings into an unvalidated “brain-muscle coupling” score |

MAS is a clinically useful ordinal resistance measure but mixes neural and non-neural contributions. It must not serve as the sole mechanistic gold standard.

## 3. Evidence classes

| Code | Evidence class | Examples in Route B2 | Maximum claim |
|---|---|---|---|
| E1 | Phenotype-linked observational association | Progression contrast in a spinal-cord-injury dataset with an experimental spasticity phenotype | “associated with experimental spasticity progression” |
| E2 | Cross-model recurrence/localisation | Reproducible program direction in independent chronic post-stroke spinal datasets or cell-type-resolved spinal data | “recurred during post-stroke spinal remodelling” or “localised to…”; not PSS-specificity |
| E3 | Real cross-context perturbation | Study-specific JGT cell perturbation, oral JGT in a non-stroke muscle model, or PF-LQ in SNI spinal cord | “directionally modulated in the stated model and preparation”; not clinical efficacy |
| E4 | Human systemic exposure | Six analytes detected in plasma after exact TJ-68 dosing in healthy adults | “detected in plasma after oral TJ-68”; not target-tissue exposure or SGD equivalence |
| E5 | Direct target evidence | Prespecified biochemical/biophysical binding or direct functional assay for the exact chemical entity and target | Direct-target wording only for the tested entity-target pair and assay context |
| E6 | Same-product clinical causal evidence | Randomized, adequately controlled trial of the locked SGD product with objective domain-specific PSS outcomes | Not available in Route B2 |

Multiple outputs derived from the same expression matrix, STRING/GO annotations or shared experimental resource are dependent evidence. Agreement among enrichment, network centrality and virtual perturbation does not upgrade evidence class.

## 4. Wording that may be used

The following wording is permissible only after the corresponding prospectively specified, internally locked analysis and robustness gate has passed. This wording does not imply public preregistration.

| Intended statement | Permissible English wording | Required evidence |
|---|---|---|
| Experimental anchor | “The program was associated with progression of experimental spasticity.” | E1, with the animal as the analysis unit and the prespecified injury-by-time contrast |
| Stroke recurrence | “The program recurred with a concordant direction during chronic post-stroke spinal remodelling.” | E2 in independent post-stroke datasets; report heterogeneity and uncertainty |
| Cell-type context | “The candidate program was enriched/localised in the prespecified spinal cell type and time window.” | Cell-type-resolved data with biological-replicate-aware inference |
| JGT response | “A chemically characterised 1:1 hot-water co-extracted JGT preparation directionally modulated the candidate program in differentiated PC12 cells.” | E3 from the exact JGT record; stable across replicates and prespecified dose analysis |
| PF-LQ bridge | “The program showed directionally coherent modulation by a paeoniflorin-liquiritin intervention in SNI-injured rat spinal cord.” | E3 from GSE205906; the primary article reports 3.98 mg/mL liquiritin, whereas the three GEO PLM samples are labelled 3.99, 3.98 and 3.100 mg/mL, respectively. This internally inconsistent metadata permits directional group analysis only until resolved. |
| Muscle bridge | “Orally administered JGT altered the prespecified module in dexamethasone-atrophy muscle.” | E3 from the exact GSE227495 preparation and tissue; identify it as a non-stroke model |
| Human exposure | “Albiflorin, paeoniflorin, glycycoumarin, isoliquiritigenin, glycyrrhetic acid and glycyrrhetic acid 3-O-monoglucuronide were detected in plasma after single oral administration of TJ-68 to healthy Japanese adults.” | E4; cite the exact commercial product, study population and primary source |
| Overall inference | “The analysis prioritises perturbation-grounded candidate programs for future PSS pharmacology and clinical PK/PD studies.” | Passed evidence gates with explicit cross-disease and preparation limitations |

Preferred verbs and modifiers are: `associated with`, `recurred across models`, `directionally coherent`, `perturbation-grounded prioritisation`, `candidate`, `preclinical`, `cross-context`, and `requires prospective validation`.

## 5. Wording that must not be used

Do not write, in the title, abstract, results, discussion, graphical abstract or press-style summary:

- “SGD treats/alleviates/reverses post-stroke spasticity.”
- “We elucidated/confirmed/validated the mechanism of SGD in PSS.”
- “The identified genes constitute a PSS-specific molecular signature.”
- “JGT, SGD and Shakuyakukanzoto are the same intervention.”
- “Total glucosides represent the whole decoction.”
- “PF-LQ validates the whole-formula mechanism or proves herbal synergy.”
- “Plasma detection proves brain, brainstem, spinal-cord or skeletal-muscle exposure.”
- “The compound reaches an effective concentration” unless the relevant total/unbound concentration, tissue, timing and potency have all been directly verified.
- “The candidate is a direct target” on the basis of expression overlap, network proximity, target prediction, docking or MD alone.
- “Virtual knockout validates necessity” or “MD confirms binding.”
- “Brain-muscle axis coupling”, “brain-muscle synchronisation” or “spatial shunting” without paired, time-aligned and biologically replicated data capable of identifying that relationship.
- “Synergistic action” without a real combination-effect design and an explicit additivity reference model.

Negative or null results do not permit post hoc relaxation of the identity or claim rules.

## 6. Exposure-specific rules

1. The human PK record concerns Tsumura TJ-68, lot H16612, in healthy Japanese adults; it is not automatically transferable to a Chinese SGD batch, Korean JGT extract, stroke population or PSS population.
2. The six human PK entities are recorded as analytes, not as proven active ingredients. Detection is evidence of systemic plasma exposure only.
3. Glycyrrhetic acid and glycyrrhetic acid 3-O-monoglucuronide are glycyrrhizic-acid-derived metabolites; trace glycyrrhetic acid in the formulation does not justify treating all circulating glycyrrhetic acid as an administered parent compound.
4. This ledger deliberately does not transcribe numeric PK values. Any concentration, AUC, half-life, dose-proportionality or potency comparison must be re-extracted from the primary table with unit, summary measure, denominator and uncertainty verified.
5. Unknown brain, brainstem, spinal-cord, skeletal-muscle and unbound exposure must remain `unknown`, not be imputed from OB, logBB, CNS-MPO or target-prediction scores.
6. A compound-target structural module is not admissible unless a relevant measured exposure and level-1/2 direct target record exist.

## 7. Claim-evidence decision table

| Claim | Current direct support | Decision |
|---|---|---|
| A molecular program tracks experimental spasticity progression | Potentially identifiable in the designated spasticity dataset | Test under the prospectively locked interaction model; retain as E1 only if robust |
| The same program occurs during chronic post-stroke spinal remodelling | Potentially identifiable across independent stroke spinal datasets | Test as cross-model recurrence; do not call it PSS-specific |
| Exact JGT preparations perturb the program in PC12/C2C12 or mouse muscle | Direct perturbation data exist in the stated contexts | Report preparation-, dose-, cell- and model-specific effects only |
| PF-LQ directionally perturbs the program in injured spinal cord | Direct SNI intervention transcriptome exists | Use as held-out cross-context calibration, not formula validation |
| Six entities enter human plasma after oral formula administration | Direct evidence exists for TJ-68 in healthy adults | Report for exact product/population only; no tissue or efficacy inference |
| SGD has a treatment effect on human PSS | No route-B2 counterfactual or same-product clinical outcome data | Prohibited |
| A compound directly binds and mediates an SGD effect through a selected target | No complete exposure-binding-engagement-phenotype chain | Prohibited unless future direct evidence is added |
| The two herbs or PF-LQ act synergistically | No valid combination-effect experiment with an additivity model in Route B2 | Prohibited |

## 8. Future clinical target product profile

The computational study should hand forward a falsifiable target product profile, not a claim that the product is already effective.

### 8.1 Product identity and quality

- Lock the exact national/pharmacopoeial identity of each crude drug, botanical species, plant part and processing method. If Paeoniae Radix Alba and prepared Glycyrrhizae Radix et Rhizoma are selected, state this prospectively rather than back-translating from JGT or TJ-68.
- Lock crude-drug ratio, combined versus separate extraction, solvent, temperature/time, concentration, drying method, extraction yield, excipients, final dosage form and storage.
- Authenticate both herbs, retain voucher specimens, archive batch certificates and specify contaminant/adulterant, pesticide, heavy-metal, microbial and mycotoxin limits.
- Establish a chemical fingerprint and release specifications for prespecified marker compounds; document batch-to-batch acceptance criteria before enrolment.
- Perform a bridging comparison against the JGT and TJ-68 materials used as computational/PK references. Shared herb names alone are not an equivalence demonstration.

### 8.2 Intended population and treatment context

- Adults with imaging-confirmed stroke and a prespecified post-stroke interval.
- Require an objective velocity-dependent stretch phenotype for entry into a neural-hyperreflexia trial; do not recruit solely by MAS.
- Record lesion location/CST burden, baseline weakness and synergy, passive stiffness/contracture, concomitant antispastic medication, botulinum-toxin timing and rehabilitation dose.
- Define whether the product is adjunctive to standardised rehabilitation and use an indistinguishable placebo where feasible.

### 8.3 Primary clinical estimand and outcomes

- Primary estimand: between-group difference in change at the prespecified time point for a velocity-dependent stretch-evoked neural response, adjusted for baseline and the repeated-measure structure.
- Candidate primary measurements: randomized-velocity stretch EMG/torque decomposition or TSRT.
- Neural secondary measures: H-reflex recruitment and rate-dependent depression, clonus and prespecified electrophysiology.
- Passive-tissue secondary measures: slow-stretch passive torque, range of motion, ultrasound/shear-wave elastography and muscle architecture.
- Voluntary-control outcomes: FMA/ARAT, strength and task performance, analysed as a separate domain rather than evidence of antispasticity.
- MAS/MTS may be included as clinical descriptors but must not be the sole mechanistic endpoint.

### 8.4 PK, pharmacodynamics and safety

- Run an analytical pilot on the final clinical batch before selecting PK analytes. The six TJ-68 plasma entities in `exposure.tsv` are candidates for assay development, not guaranteed markers of the final SGD product.
- Prespecify sampling times, total versus unbound plasma measurements, lower limit of quantification, metabolite identity and exposure-response analyses. Do not claim CNS or muscle penetration unless those matrices are directly measured under an ethically justified design.
- Select no more than 1-2 target-engagement hypotheses and 3-5 measurable biomarkers after the computational gates pass. A biomarker must map to a PSS domain and have an assay suitable for the planned population.
- Monitor potassium, blood pressure, oedema/body weight, renal and hepatic indices and concomitant drugs relevant to licorice-associated pseudoaldosteronism; define stopping and dose-modification rules.

### 8.5 Prospective claim ladder

| Future evidence obtained | Highest defensible claim |
|---|---|
| Locked product plus quantified plasma analytes | Systemic exposure to the exact product |
| Randomized difference in an objective domain-specific PSS endpoint | Clinical effect of the exact product on that endpoint and population |
| Exposure-response plus prespecified pharmacodynamic biomarker change | PK/PD association consistent with a candidate pathway |
| Direct binding/functional evidence plus in-vivo target engagement | Direct target evidence for the tested entity-target pair |
| Target perturbation or other valid causal mediation showing loss/rescue of the same phenotype | Mechanistic necessity/mediation, limited to the tested model and assumptions |

Only the corresponding evidence rung may be claimed; higher rungs cannot be inferred from lower ones.

## 9. Core sources

- AHA PSS scientific statement and three-domain framework: <https://pmc.ncbi.nlm.nih.gov/articles/PMC13063836/>
- Multidimensional JGT perturbome and formula/preparation metadata: <https://www.nature.com/articles/s41597-026-06759-6>
- JGT oral intervention in dexamethasone-atrophy muscle: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227495>
- PF-LQ intervention in SNI spinal cord: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205906>
- Human TJ-68 PK: <https://doi.org/10.1002/jps.24596>
- Total-glucosides clinical/sEMG report: <https://pubmed.ncbi.nlm.nih.gov/24515912/>
- Specific SGD analytical preparation used as an identity exemplar: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9055985/>
- PF-LQ primary article and the conflicting sequenced-sample metadata: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9489943/> and GEO samples GSM6235298–GSM6235300.
