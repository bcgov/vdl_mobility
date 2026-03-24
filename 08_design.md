# Occupational Mobility Research Design
*Working Notes — Research Proposal*

---

## Research Question

Do skill distance or hierarchical distance better rationalize occupational mobility among continuously employed workers in Canada, using linked 2016–2021 census data at the 5-digit NOC level?

---

## Methodology

- Sinkhorn optimal transport with entropic regularization
- Empirical joint distribution P(i,j) from linked census, with 2016 and 2021 occupational distributions imposed as hard margin constraints
- Goodness of fit: KL divergence between observed and predicted transport plan under each cost matrix
- ε set separately for each cost matrix via conditional median matching, reported across grid [0.5×, 2×] baseline at every test

---

## Cost Matrices

### C_hier — Tree Metric on NOC Structure (scaled by dividing by 4)

| Condition | Raw Cost | Scaled (/4) |
|---|---|---|
| 5-digit match | 0 | 0 |
| 4-digit match | 1 | 0.25 |
| 3-digit match | 2 | 0.5 |
| 2-digit match | 3 | 0.75 |
| 1-digit match, \|ΔTEER\| = 1–5 | 4–8 | 1.0–2.0 |
| Cross-sector (otherwise) | 9 | 2.25 |

*Note: 2nd digit of NOC 2021 is TEER, so 3-digit match implies TEER match. Max |ΔTEER| = 5, giving max within-sector cost of 8 (scaled: 2.0) — cross-sector penalty of 9 (scaled: 2.25) is only marginally above maximum within-sector cost.*

### C_skill — Continuous Euclidean Distance (scaled by dividing by 9)

- Euclidean distance in 10-dimensional PCA subspace of 161-dimensional O*NET space
- O*NET dimensions: skills, abilities, knowledge, work activities
- Scaled by dividing by 8th percentile value (= 9) — matching the 8th percentile of C_hier non-9 values after scaling

*Scaling rationale: conditional median matching — median of C_hier non-9 values = 4; 8th percentile of full C_skill distribution = 9. Dividing by these values anchors both distributions at 1 for their lower meaningful mass.*

*O*NET crosswalk: many O*NET codes mapping to one NOC is acceptable — NOC is an aggregate by design and the centroid skill vector is legitimate. One O*NET code mapping to many NOCs is the problematic case — affected NOCs receive identical skill vectors, making C_skill uninformative for those pairs.*

---

## Two-Stage Identification Structure

### Stage 1: Does any distance-based model work?

Education change used as a classifier sorting workers into groups where distance-based models are theoretically appropriate or not — not a causal claim.

**Clean group (all three conditions required):**

1. No change in highest attainment between 2016 and 2021 — origin NOC not stale
2. CoW = employee in both 2016 and 2021 — self-employed NOC reporting too noisy for distance-based model assumptions
3. 2016 NOC (translated to NOC 2021) TEER within ±1 of credential-implied TEER — origin NOC not mismeasured

*TEER 0 excluded from condition 3 — experience-based pathways to management make credential matching uninformative for this group. Condition 3 requires NOC 2016 → NOC 2021 translation to be applied before the filter; contingent on what Stats Can provides on data access.*

| Group | Expected Fit | Reason |
|---|---|---|
| Satisfies all clean group conditions | Good — main clean group | All assumptions satisfied |
| Change in highest attainment | Poor | Origin NOC stale |
| TEER-incommensurate in 2016 | Poor | Origin NOC mismeasured |
| Self-employed either period | Poor | NOC reporting noisy |
| Multiple violations | Poor — saturation expected | No prior on interaction; once origin NOC is noisy, additional violations likely redundant |

### Stage 2: Skill or hierarchy — conditional on framework being valid

Within the clean group, horse race between C_skill and C_hier with heterogeneity analyses by destination NOC characteristics, all with directional priors.

| Split | Directional Prior |
|---|---|
| Destination TEER | C_skill dominates at TEER 0 and 5; C_hier dominates in middle — inverted-U of occupational specificity |
| Destination gating (KL divergence of education distribution vs overall, quartiled) | C_hier dominates Q4 (most gated); C_skill dominates Q1 (least gated) |

**Destination TEER predictions:**

- **TEER 0 (management):** workers arrive from diverse sectors because general managerial capital is portable — C_hier's cross-sector penalty too harsh, C_skill should fit better
- **TEER 5:** workers arrive from diverse sectors because skill barriers are low — same prediction, opposite mechanism
- **TEER 2–3 (middle):** arrivals constrained by credential and sector — C_hier should fit better

*Destination gating measure: KL divergence of destination NOC's education distribution relative to the overall education distribution across all NOCs — data-driven identification of credential-locked occupations without relying on external regulatory classifications.*

Theoretically motivated cost structure retained as fixed throughout; the pattern of failures is the finding.

---

## Robustness Check

**Exclude 1:many O*NET to NOC 2021 mappings**

Main result retains all NOCs in the clean group. This check excludes any transition where either the origin or destination NOC has a 1:many O*NET mapping — i.e. where one O*NET code is assigned to multiple NOCs, giving those NOCs identical skill vectors and making C_skill uninformative for those pairs. This biases toward C_hier in the main result, since C_hier retains its discriminating power while C_skill does not. Excluding these pairs on both origin and destination is conservative but internally consistent:

- If C_hier wins in the main result, the robustness check removes a bias that was working in C_hier's favour — C_hier winning on the restricted sample is a stronger result
- If C_skill wins in the main result, C_skill is winning despite artificial compression of some of its distances — also a strong result

ε sensitivity absorbed into reporting protocol: every KL divergence test reported across grid [0.5×, 2×] baseline. Invariance of ranking across this range is a stability condition baked into every test, not a separate check.

---

## Sample

- 1/16 long form census, linked 2016–2021 via Statistics Canada Virtual Data Lab (custom linkage)
- Conditioned on employment (NOC reported) in both waves — definitional restriction, not a limitation. Estimand is occupational mobility among continuously employed workers.
- May 2021 reference week (mid-pandemic): pandemic displacement treated as a stress test of cost geometry rather than a threat to the design
- Highest educational attainment with CIP at 4-digit level expected to be available
- NOC vintage and internal concordance to be confirmed on data access; TEER assignments follow Stats Can's internal NOC 2021 coding

---

## Proposed Paper Structure

1. Does distance rationalize mobility? (Stage 1 — clean vs dirty group)
2. Skill or hierarchy? (Stage 2 — main horse race result within clean group)
3. Where does each metric dominate and why? (Stage 2 heterogeneity by destination TEER and destination gating)
4. Where does the framework break down and why? (Stage 1 failures — interpretable against theory)

---

## Notes on Gender

Gender differences explored descriptively in appendix only — no directional prior on which metric should dominate by gender, and occupational segregation argument is largely absorbed by the margin constraints in the Sinkhorn formulation. Not a main analysis.
