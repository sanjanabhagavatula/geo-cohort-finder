# Drug Response Meta-Analysis

## Purpose

Analyze publicly available omics evidence related to patient response to a cancer therapy.

Given a therapy and cancer type, identify relevant public datasets, extract and validate study metadata, determine available omics modalities and response phenotypes, analyze responder versus non-responder molecular differences when appropriate data are available, and synthesize evidence across studies.

The skill should preserve uncertainty and distinguish between unavailable information, negative findings, and evidence of non-response.

---

## Input

The user may provide a natural-language research question.

Example:

```text
Analyze pembrolizumab response in HNSCC.
```

The agent should extract relevant criteria from the request.

Possible criteria include:

- therapy or drug
- cancer or disease type
- omics modality
- treatment-response phenotype
- sample type
- treatment timepoint
- patient population
- other user-specified biological constraints

Not all criteria are required.

If a criterion is not specified, it should remain unrestricted unless it can be safely expanded from an explicitly provided concept.

For example:

```text
Input:
"Analyze pembrolizumab response in HNSCC."

Interpretation:

drug = pembrolizumab
disease = HNSCC
omics = unspecified
objective = treatment response
```

Do not invent missing biological or clinical constraints.

---

## Workflow

```text
 │
 │ "Analyze pembrolizumab response in HNSCC"
 ▼
LLM / Agent
 │
 │ reads SKILL.md
 │
 │ decides:
 │   - drug = pembrolizumab
 │   - disease = HNSCC
 │   - find public datasets
 │   - classify modalities
 │   - analyze response
 ▼
drug_response_meta.py
 │
 ├── dataset discovery
 ├── metadata extraction
 ├── data download
 ├── modality-specific processing
 ├── differential analysis
 ├── pathway analysis
 └── cross-study synthesis
 │
 ▼
output/
 ├── report.md
 ├── result.json
 ├── tables/
 ├── figures/
 └── reproducibility/
```

---

## Agent Responsibilities

The agent is responsible for interpreting the user's biological question and determining the appropriate analysis objective.

The agent should:

1. Identify the requested therapy.
2. Identify the requested cancer or disease context.
3. Identify any requested omics modality.
4. Identify whether treatment response, resistance, survival, or another phenotype is being investigated.
5. Preserve unspecified criteria as unspecified.
6. Invoke `drug_response_meta.py` with the appropriate parameters.
7. Interpret the structured results returned by the workflow.
8. Clearly communicate supporting evidence, conflicting evidence, missing information, and limitations.

The agent must not infer that a dataset contains a particular modality, treatment, response phenotype, or clinical variable unless this is established from retrieved metadata.

---

## Main Workflow

Run:

```bash
python drug_response_meta.py \
    --drug "<therapy>" \
    --disease "<cancer>"
```

Optional arguments may include:

```bash
--omics "<modality>"
--phenotype "<phenotype>"
--timepoint "<timepoint>"
--sample-type "<sample_type>"
```

Example:

```bash
python drug_response_meta.py \
    --drug "pembrolizumab" \
    --disease "HNSCC"
```

The workflow should perform the following stages.

### 1. Dataset Discovery

Search supported public repositories for studies relevant to the requested therapy and disease.

Search terms may be expanded to include established synonyms and related terminology when appropriate.

Candidate datasets should be recorded before determining suitability.

Dataset discovery should favor human patient-derived studies when the user's question concerns patient treatment response.

Do not assume that a dataset is suitable merely because its title or description mentions the requested therapy or cancer.

### 2. Metadata Extraction

For each candidate dataset, retrieve available study-level and sample-level metadata.

Where available, identify:

- dataset accession
- repository
- study title
- publication
- organism
- cancer type
- treatment
- omics modality
- number of patients
- number of samples
- sample type
- treatment timepoint
- response phenotype
- response definition
- response labels
- survival information
- availability of processed data
- availability of raw data

Metadata that cannot be established should be recorded as `unknown` or `null`.

Missing metadata must never be interpreted as a negative result.

### 3. Data Download

Download data only for datasets that contain sufficient evidence to address the user's research question.

For treatment-response analysis, preference should be given to datasets with:

- confirmed exposure to the requested therapy
- identifiable patient or sample groups
- usable molecular data
- confirmed response or outcome annotations

When practical, prefer analysis-ready processed data over reprocessing raw sequencing or mass-spectrometry data.

Record dataset identifiers, source files, download locations, and relevant version information for reproducibility.

### 4. Modality-Specific Processing

Determine the omics modality from retrieved metadata and route the dataset to an appropriate analysis workflow.

Possible modalities include:

- bulk RNA-seq
- single-cell RNA-seq
- microarray expression
- genomics
- proteomics
- phosphoproteomics
- epigenomics
- spatial omics
- multi-omics

Do not apply an analysis workflow designed for one modality to another modality.

If the modality is unsupported, retain the dataset in the discovery report but do not perform an inappropriate downstream analysis.

### 5. Differential Analysis

When appropriate response groups are available, compare molecular features between response groups.

Possible comparisons include:

```text
Responder vs Non-responder
CR/PR vs SD/PD
Sensitive vs Resistant
```

Response definitions must be derived from the source study.

Do not silently convert clinical outcome categories into binary response groups.

Any harmonization must be documented in the output and reproducibility record.

The analysis should record:

- original response labels
- derived analysis groups, if applicable
- sample counts per group
- statistical method
- thresholds
- significant features
- effect sizes
- uncertainty

### 6. Pathway Analysis

Where appropriate, perform pathway or gene-set analysis using the results of the molecular comparison.

The output should distinguish between:

- statistically supported pathway findings
- descriptive trends
- unsupported biological interpretations

Pathway findings should remain linked to the dataset and comparison from which they were derived.

### 7. Cross-Study Synthesis

If multiple suitable studies are available, compare findings across studies.

Evaluate:

- concordant molecular signals
- discordant findings
- recurring pathways
- differences in patient populations
- differences in treatment regimens
- differences in response definitions
- differences in omics modalities
- differences in sampling timepoints

Do not treat heterogeneous studies as directly comparable without documenting relevant differences.

Cross-study synthesis should identify where independent evidence converges and where evidence remains inconsistent or insufficient.

---

## Response Phenotype Handling

Treatment response may be represented differently across studies.

Examples include:

```text
CR / PR / SD / PD
Responder / Non-responder
Sensitive / Resistant
Progression-free survival
Overall survival
Duration of response
```

These outcomes are not automatically interchangeable.

The workflow should preserve the original phenotype definition and document any transformations used for downstream analysis.

For example:

```json
{
  "response_definition": "RECIST",
  "original_labels": ["CR", "PR", "SD", "PD"],
  "analysis_groups": {
    "responder": ["CR", "PR"],
    "nonresponder": ["SD", "PD"]
  }
}
```

If response metadata are absent, report response as unknown.

Do not interpret missing response information as non-response.

---

## Structured Output

The workflow should generate a machine-readable `result.json`.

Example structure:

```json
{
  "query": {
    "drug": "pembrolizumab",
    "disease": "HNSCC",
    "objective": "treatment response",
    "omics": null
  },

  "datasets": [
    {
      "accession": "GSEXXXXX",
      "repository": "GEO",
      "modality": "RNA-seq",
      "sample_count": 50,
      "therapy_confirmed": true,
      "response_available": true,
      "response_definition": "RECIST",
      "suitability": "suitable"
    }
  ],

  "analyses": {
    "differential_analysis": {},
    "pathway_analysis": {},
    "cross_study_synthesis": {}
  },

  "limitations": []
}
```

Use `null` or `unknown` for information that cannot be established.

Do not fabricate values to complete the schema.

---

## Human-Readable Report

Generate `output/report.md` summarizing the analysis.

The report should include:

### Research Question

State the interpreted therapy, disease, phenotype, and any additional constraints.

### Dataset Discovery

List candidate datasets and explain why they were identified.

### Dataset Suitability

For each dataset, report whether the required disease, treatment, modality, and response information were confirmed.

### Patient Response Evidence

Summarize response information available from the included datasets without generalizing beyond the studied populations.

### Molecular Differences

Summarize differential molecular features identified between relevant response groups.

### Pathways

Summarize supported pathway-level findings.

### Cross-Study Evidence

Describe findings that replicate or conflict across datasets.

### Limitations and Missing Evidence

Explicitly report unavailable metadata, small cohorts, heterogeneous response definitions, unsupported modalities, inaccessible data, or other limitations.

---

## Output Files

```text
output/
├── report.md
├── result.json
├── tables/
├── figures/
└── reproducibility/
```

### `report.md`

Human-readable synthesis of the analysis.

### `result.json`

Machine-readable representation of the query, datasets, analyses, and limitations.

### `tables/`

Structured results such as:

- dataset summaries
- sample metadata
- differential analysis results
- pathway enrichment results
- cross-study comparisons

### `figures/`

Generated visualizations where appropriate.

### `reproducibility/`

Record sufficient information to reproduce the workflow, including where possible:

- original query
- interpreted parameters
- dataset accessions
- search queries
- downloaded files
- analysis parameters
- software/tool versions
- commands executed

---

## Evidence and Safety Rules

This skill is intended for research and educational use.

It must not claim that an individual patient will or will not respond to a therapy.

Results should be described in the context of the populations and datasets from which they were derived.

Association does not establish causation.

A molecular feature associated with response should not automatically be described as a predictive biomarker.

A feature identified in one dataset should not automatically be described as validated.

Missing evidence is not evidence of absence.

Conflicting evidence should be reported rather than resolved by unsupported inference.

All dataset accessions and source information should be retained so users can inspect the underlying evidence.

---

## MVP Behavior

If the full analysis workflow is unavailable, the skill should still perform:

```text
research question
        ↓
dataset discovery
        ↓
metadata extraction
        ↓
dataset suitability assessment
        ↓
structured dataset report
```

Failure or unavailability of downstream differential or pathway analysis should not prevent reporting suitable public datasets.

---

## Future Extensions

Potential extensions include:

- additional public omics repositories
- automated response-label harmonization
- integration with existing ClawBio RNA-seq skills
- integration with existing ClawBio single-cell skills
- integration with proteomics analysis skills
- survival analysis
- biomarker replication across independent cohorts
- meta-analysis
- automatic identification of associated publications
- clinical-trial integration