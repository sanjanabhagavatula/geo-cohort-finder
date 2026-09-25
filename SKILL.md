---
name: geo-cohort-finder
description: >-
  Find analyzable human patient cohorts in NCBI GEO from a natural-language
  research question. Expands drug and disease synonyms, searches GEO, parses
  sample-level metadata, and returns a dataframe of candidate studies with
  patient counts, sample type, treatment exposure and treatment-response
  labels. Discovery and triage only; hands accessions to downstream skills.
license: MIT
metadata:
  version: "0.1.0"
  author: TODO-team-names
  domain: transcriptomics
  tags:
    - GEO
    - dataset discovery
    - cohort
    - treatment response
    - metadata
    - oncology
  inputs:
    - name: query
      type: string
      format:
        - text
      description: Natural-language research question, or explicit --drug/--disease flags
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Human-readable discovery and suitability report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable query, cohorts, and limitations
    - name: cohorts
      type: file
      format:
        - tsv
      description: One row per GEO series
    - name: samples
      type: file
      format:
        - tsv
      description: One row per GEO sample, with verbatim characteristics
  dependencies:
    python: ">=3.10"
    packages:
      - pandas>=2.0
      - requests>=2.31
  demo_data:
    - path: demo_data/esummary_hnscc.json
      description: Cached GEO esummary records, offline
    - path: demo_data/samples_hnscc.json
      description: Cached sample-level metadata, offline
  endpoints:
    cli: python skills/geo-cohort-finder/geo_cohort_finder.py --query "{query}" --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🔎"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: requests
    trigger_keywords:
      - find GEO datasets
      - public datasets for
      - responder vs non-responder
      - treatment response dataset
      - GEO cohort
      - GSE search
      - patient cohorts for
      - is there public data on
---

# 🔎 geo-cohort-finder

Find analyzable human patient cohorts in NCBI GEO, and say plainly which ones
are not analyzable and why.

## Trigger

Use this skill when a user asks whether public data exists for a biological
question — a drug, a cancer type, an omics modality, a treatment-response
comparison — or asks which GEO studies could support such an analysis.

Do **not** use this skill to download expression data or to run differential
expression, pathway or cross-study analysis. See **Integration with Bio
Orchestrator**.

## Why This Exists

GEO's own search covers series-level text only: title, summary and overall
design. The facts that decide whether a study is analyzable live in per-sample
`characteristics_ch1` fields, which that index does not reach. So the GEO
website cannot answer the question researchers actually ask — *are these
patients, how many of them, were they treated, and is there an outcome?*

The gap is large. A representative query for head and neck squamous carcinoma
(`Homo sapiens`, GSE entry type, expression profiling by high-throughput
sequencing) returns 784 series. Of 304 sampled, 121 had at least 20 samples,
and 40 of those 121 matched cell-line or perturbation keywords in their title
or summary alone — keyword matching undercounts, so the true proportion of
non-patient studies is higher. Sample count is also not patient count: the
largest hit by samples was a single-cell study of roughly 18 patients.

A researcher who filters on sample count alone therefore ranks single-cell
studies first, counts a six-patient study as large, and reads cell-line
experiments as patient cohorts.

## Core Capabilities

1. Parse a natural-language question into a typed query.
2. Expand drug and disease terms from an editable synonym table.
3. Search GEO via NCBI E-utilities and retrieve series-level summaries.
4. Fetch sample-level metadata for surviving candidates.
5. Classify sample type, patient count, treatment exposure and response labels
   against explicit rules.
6. Emit a suitability verdict per series, with the reason and the verbatim
   metadata string behind it.

## Scope

**One skill, one task.** This skill discovers and characterizes GEO series. It
does not download expression data, process raw reads, run differential
expression, perform pathway analysis, or synthesize findings across studies.
If the user wants those, route to `article-data-fetcher`,
`nfcore-rnaseq-wrapper`, `rnaseq-de` or `pathway-enricher`.

## Input Formats

A natural-language research question:

```text
Analyze pembrolizumab response in HNSCC.
List GEO datasets on HNSC with transcriptomic data.
```

Or explicit flags. Criteria are optional and independent.

The agent extracts criteria from the request. Possible criteria include
therapy or drug, cancer or disease type, omics modality,
treatment-response phenotype, sample type, treatment timepoint, patient
population, and other user-specified biological constraints.

**If a criterion is not specified, it remains unrestricted unless it can be
safely expanded from an explicitly provided concept. Do not invent missing
biological or clinical constraints.**

The parsed query is printed back before the search runs:

```text
drug:          pembrolizumab → [Keytruda, MK-3475, anti-PD-1, pembro]
cancer_type:   HNSCC → [head and neck squamous, SCCHN, oral cavity carcinoma, ...]
omics:         unspecified
objective:     treatment response
organism:      Homo sapiens
min_patients:  20
sample_type:   patient tumor
```

## Workflow

```text
"Analyze pembrolizumab response in HNSCC"
        ▼
   LLM / Agent  ── slot-filling only ──▶ typed query
        ▼
   geo_cohort_finder.py
        │
        ├── 1. synonym expansion        (rules table)
        ├── 2. broad search             (esearch, one call)
        ├── 3. series summaries         (esummary, batched)
        ├── 4. deep fetch               (sample characteristics, survivors only)
        └── 5. suitability scoring      (rules table)
        ▼
   output/  report.md · result.json · tables/ · reproducibility/
        ▼
   ── chains to ──▶ article-data-fetcher → nfcore-rnaseq-wrapper → rnaseq-de
```

### Stage boundaries

Stages have very different costs and failure modes, so each writes its output
to disk and later stages read from it. A run can resume from any stage.

| Stage | Method | Cost | Output |
|---|---|---|---|
| Broad search | `esearch` on `gds`, `GSE[ETYP]` | one call | `tables/candidates.tsv` |
| Series summaries | `esummary` on `gds`, batched | ~1 call per 100 | `tables/candidates.tsv` |
| Deep fetch | per-series sample metadata | seconds each | `tables/samples.tsv` |
| Scoring | local rules | instant | `tables/cohorts.tsv` |

Every network response is cached on disk keyed by accession. Re-runs do not
re-fetch.

## CLI Reference

```bash
# Standard usage
python geo_cohort_finder.py --query "pembrolizumab response in HNSCC" --output output/

# Explicit criteria
python geo_cohort_finder.py --drug pembrolizumab --disease HNSCC --omics rna-seq

# Discovery without a drug
python geo_cohort_finder.py --disease HNSC --omics transcriptomic --min-patients 30

# Demo mode (cached records, no network)
python geo_cohort_finder.py --demo

# Resume from a completed stage
python geo_cohort_finder.py --resume --stage characterize

# Via ClawBio runner
python clawbio.py run geo-cohort-finder --demo
```

## Demo

```bash
python clawbio.py run geo-cohort-finder --demo
```

Runs the full pipeline against cached GEO records for an HNSCC transcriptomics
query. No network access, no credentials, no user files. Produces the complete
output tree so the report format can be inspected offline.

## Algorithm / Methodology

### Search construction

The broad search targets the `gds` database, restricted to `GSE[ETYP]` and
`Homo sapiens[Organism]`, with expanded disease and drug terms joined by `OR`
and modality mapped to a `[DataSet Type]` filter.

Drug terms are **not** used to confirm treatment. A drug named in a series
summary is a candidate signal only. **Do not assume that a dataset is suitable
merely because its title or description mentions the requested therapy or
cancer.** Confirmation comes from sample-level metadata in the deep fetch.

### Domain decisions

These are the scientific rules this skill encodes. They are versioned here,
not inferred at runtime. Changing behavior means editing these tables.

#### Configurable defaults

Printed in the report header on every run.

```text
MIN_PATIENTS_TOTAL      20
MIN_PATIENTS_PER_ARM    10
MAX_SERIES_DEEP_FETCH   200
EXCLUDE_SAMPLE_TYPES    cell line, PDX, organoid, xenograft
ORGANISM                Homo sapiens
```

#### Sample-type classification

Applied to sample-level metadata first, falling back to series text.

| Pattern in `source_name` / `characteristics` | Class |
|---|---|
| `tumor biopsy`, `patient`, `primary tumor`, `FFPE`, `resection` | `patient_tumor` |
| `cell line`, named lines (`CAL27`, `FaDu`, `SCC-25`, `HN-SCC-*`) | `cell_line` |
| `PDX`, `xenograft`, `organoid` | `model` |
| `siRNA`, `shRNA`, `knockout`, `overexpress`, `treated ... for Nh` | `perturbation` |
| no match | `NEEDS_REVIEW` |

#### Response label mapping

Matched case-insensitively against sample characteristics. The verbatim source
string and its field name are retained for every mapped label.

| Source string | Mapped |
|---|---|
| `CR`, `complete response` | `CR` |
| `PR`, `partial response` | `PR` |
| `SD`, `stable disease` | `SD` |
| `PD`, `progressive disease`, `progression` | `PD` |
| `R`, `responder`, `response: yes`, `benefit` | `RESPONDER_UNSPECIFIED` |
| `NR`, `non-responder`, `nonresponder`, `response: no` | `NONRESPONDER_UNSPECIFIED` |
| `sensitive` / `resistant` | `SENSITIVE` / `RESISTANT` |
| anything else | `NEEDS_REVIEW` |

**An unmatched string is never guessed.** `NEEDS_REVIEW` means unparsed, not
absent.

#### Response group derivation

RECIST categories are reported raw, and **both** conventional binarizations are
reported as separate columns. The skill does not choose between them.

```text
strict:  responder = CR + PR           nonresponder = SD + PD
dcb:     responder = CR + PR + SD*     nonresponder = PD
         (* durable stable disease, where duration is recorded)
```

Stable disease is the pivot: a patient stable for 14 months falls in opposite
groups under the two rules, and published studies use both. **Do not silently
convert clinical outcome categories into binary response groups.** Any
harmonization is documented in the output and the reproducibility record.

Where a study defines its own response groups, those are preserved as the
primary labels and the derived groups are marked as derived.

#### Suitability verdicts

| Verdict | Criteria |
|---|---|
| `SUITABLE` | therapy confirmed in sample metadata AND response labels present AND `n_patients >= MIN_PATIENTS_TOTAL` AND both arms `>= MIN_PATIENTS_PER_ARM` AND `sample_type == patient_tumor` |
| `SUITABLE_NO_DRUG_FILTER` | as above, drug not requested |
| `INSUFFICIENT_N` | all criteria met except patient counts |
| `NO_RESPONSE_LABELS` | patient cohort, therapy confirmed, no outcome annotation |
| `SURVIVAL_ONLY` | survival recorded, no response call |
| `WRONG_SAMPLE_TYPE` | `cell_line`, `model` or `perturbation` only |
| `THERAPY_UNCONFIRMED` | drug in series text only, absent from sample metadata |
| `NEEDS_REVIEW` | metadata present but unparsed |

**Metadata that cannot be established is recorded as `unknown` or `null`.
Missing metadata is never interpreted as a negative result.** A
`NO_RESPONSE_LABELS` verdict means the labels were not found in the deposited
metadata — they may exist in the paper or its supplements.

### Patient counts

`n_samples` is the deposited sample count. `n_patients` is derived from
patient or subject identifiers in sample characteristics where present, and is
otherwise `unknown` — never assumed equal to `n_samples`.

Paired designs (pre- and on-treatment biopsies) and single-cell studies both
inflate sample count relative to patient count, in opposite magnitudes.
Filtering on `n_samples` is therefore not a proxy for cohort size.

## Example Queries

```text
Find GEO datasets on HNSC with transcriptomic data.
Which public cohorts have pembrolizumab responders and non-responders?
Is there public methylation data on cisplatin resistance in head and neck cancer?
Show me patient RNA-seq cohorts of at least 50 people in bladder cancer.
```

## Example Output

````text
# GEO Cohort Finder Report

## Research Question
drug: pembrolizumab (+4 synonyms) · disease: HNSCC (+6 synonyms)
omics: unspecified · objective: treatment response
defaults: MIN_PATIENTS_TOTAL=20, MIN_PATIENTS_PER_ARM=10

## Search
784 series matched the broad query. 121 passed the sample-count prefilter.
121 deep-fetched.

## Verdicts
SUITABLE                 4
INSUFFICIENT_N          11
NO_RESPONSE_LABELS      38
SURVIVAL_ONLY            6
WRONG_SAMPLE_TYPE       47
THERAPY_UNCONFIRMED      9
NEEDS_REVIEW             6

## Suitable cohorts
| accession | n_pat | n_samp | platform | response field | CR/PR/SD/PD | strict R/NR |
|-----------|-------|--------|----------|----------------|-------------|-------------|
| GSEnnnnnn |    36 |     72 | Illumina | response       | 3/8/11/14   | 11 / 25     |

## Limitations
- 6 series carry unparsed response strings (NEEDS_REVIEW), listed in tables/samples.tsv.
- n_patients could not be established for 19 series; reported as unknown.
- Verdicts reflect deposited metadata only, not the source publications.
````

## Output Structure

```text
output/
├── report.md
├── result.json
├── tables/
│   ├── candidates.tsv     one row per series, search stage
│   ├── cohorts.tsv        one row per series, scored
│   └── samples.tsv        one row per sample, verbatim characteristics
└── reproducibility/
    ├── commands.sh
    ├── query.json         parsed query and expanded synonyms
    ├── esearch_terms.txt  exact search string sent to NCBI
    └── environment.yml
```

### `tables/cohorts.tsv`

```text
gse_accession, title, pubmed_id, n_samples, n_patients,
organism, platform, omics_type, processing_level,
drug_matched, drug_evidence_field, drug_evidence_string,
sample_type, biopsy_timing,
has_response_labels, response_field, response_counts_raw,
n_responder_strict, n_nonresponder_strict,
n_responder_dcb, n_nonresponder_dcb,
has_survival, verdict, verdict_reason, ftp_link
```

### `tables/samples.tsv`

```text
gse_accession, gsm_accession, patient_id, source_name,
raw_characteristics, response_raw, response_mapped, timepoint, sample_type
```

`raw_characteristics` is the unmodified deposited string. Every derived column
can be audited against it.

### `result.json`

```json
{
  "query": {
    "drug": "pembrolizumab",
    "drug_synonyms": ["Keytruda", "MK-3475", "anti-PD-1"],
    "disease": "HNSCC",
    "objective": "treatment response",
    "omics": null,
    "defaults": {"min_patients_total": 20, "min_patients_per_arm": 10}
  },
  "search": {
    "esearch_term": "...",
    "n_matched": 784,
    "n_prefiltered": 121,
    "n_deep_fetched": 121
  },
  "cohorts": [
    {
      "accession": "GSEXXXXXX",
      "repository": "GEO",
      "modality": "RNA-seq",
      "n_samples": 72,
      "n_patients": 36,
      "sample_type": "patient_tumor",
      "therapy_confirmed": true,
      "therapy_evidence": {"field": "treatment", "value": "pembrolizumab 200mg q3w"},
      "response_available": true,
      "response_definition": "RECIST",
      "original_labels": ["CR", "PR", "SD", "PD"],
      "analysis_groups": {
        "strict": {"responder": ["CR", "PR"], "nonresponder": ["SD", "PD"]},
        "dcb": {"responder": ["CR", "PR", "SD"], "nonresponder": ["PD"]}
      },
      "verdict": "SUITABLE",
      "verdict_reason": "therapy confirmed; 36 patients; 11/25 strict"
    }
  ],
  "limitations": []
}
```

**Use `null` or `unknown` for information that cannot be established. Do not
fabricate values to complete the schema.**

## Human-Readable Report

`output/report.md` contains: the interpreted research question and the
defaults in force; the search terms and match counts; a verdict summary; the
suitable cohorts table; and a limitations section that explicitly reports
unavailable metadata, small cohorts, heterogeneous response definitions,
unsupported modalities and inaccessible data.

## Dependencies

- Python >= 3.10, `pandas`, `requests`
- Network access to `eutils.ncbi.nlm.nih.gov` (see `docs/data-handling.md`)
- Optional `NCBI_API_KEY` environment variable

No credentials are required. An NCBI API key raises the E-utilities rate limit
from 3 to 10 requests per second and is strongly recommended for interactive
use.

## Gotchas

- **GEO's search index is series-level.** Response labels, treatment and
  patient identifiers live in per-sample `characteristics_ch1` and are
  invisible to `esearch`. The deep fetch is not optional.
- **`n_samples` is not `n_patients`.** Paired timepoint designs double-count;
  single-cell series report thousands of cells for a few dozen patients.
- **One GSE can bundle several cohorts** — different treatments, platforms, or
  a cell-line arm beside a patient arm. Heterogeneous characteristics within a
  series are flagged rather than collapsed to one verdict.
- **SuperSeries and SubSeries duplicate each other.** Deduplicate on the
  relations field or counts inflate.
- **A drug named in a series summary was often not administered** in that
  study — it may be background, motivation or a cited comparison.
- **Biopsy timing is the silent confounder.** On-treatment samples differ from
  pre-treatment samples for reasons unrelated to response. Timing is reported
  wherever recorded, and `unknown` otherwise.
- **Free-text response strings are irregular**: `response: PR`,
  `best overall response: partial response`, `RECIST: 1`, `benefit: yes`, `R`.
  The mapping table is the product; extend it rather than loosening matching.
- **Processing level varies.** Supplementary deposits may be raw counts, FPKM,
  TPM or already-normalized values. A count-based downstream tool cannot
  consume an FPKM-only deposit.
- **Rate limits bite during demos.** Cache everything; use an API key.

## Safety

This skill is intended for research and educational use. It is not a medical
device and does not provide clinical decision support.

It must not claim that an individual patient will or will not respond to a
therapy. Results should be described in the context of the populations and
datasets from which they were derived.

Association does not establish causation. A molecular feature associated with
response should not automatically be described as a predictive biomarker. A
feature identified in one dataset should not automatically be described as
validated.

**Missing evidence is not evidence of absence.** Conflicting evidence should be
reported rather than resolved by unsupported inference.

All dataset accessions and source metadata strings are retained so users can
inspect the underlying evidence.

## Agent Boundary

The agent's role is **slot-filling and interpretation only**:

1. Extract criteria from the user's question.
2. Preserve unspecified criteria as unspecified.
3. Invoke `geo_cohort_finder.py` with those parameters.
4. Interpret the structured results.
5. Communicate supporting evidence, missing information and limitations.

**The agent must not infer that a dataset contains a particular modality,
treatment, response phenotype or clinical variable unless this is established
from retrieved metadata.** Classification, mapping and scoring are performed by
the rules tables above, not by the model. Identical input produces identical
output.

This skill reports which datasets exist and whether they are analyzable. It
does not assess whether a study is well powered for a given effect size, and a
`SUITABLE` verdict is not a claim that the cohort is adequate for genome-wide
discovery. Public treatment-response cohorts are typically 20–80 patients.

## Integration with Bio Orchestrator

`geo-cohort-finder` is a discovery front door. It emits accessions; other
skills consume them.

```text
geo-cohort-finder
  ├─▶ article-data-fetcher      retrieve deposited files for an accession
  ├─▶ nfcore-rnaseq-wrapper     FASTQ/BAM → counts
  ├─▶ rnaseq-de                 responder vs non-responder contrast
  ├─▶ de-summary                rank and summarize results
  ├─▶ pathway-enricher          gene-set enrichment
  ├─▶ scrna-orchestrator        single-cell series
  ├─▶ methylation-clock         methylation series
  └─▶ clinical-trial-finder     trials for the same drug and indication
```

Routing rule: questions about *which data exists* come here; questions about
*what the data shows* route onward.

## Maintenance

The synonym, sample-type, response-mapping and suitability tables are the
maintained surface of this skill. Extend them rather than relaxing matching
logic. Every table change should be accompanied by a gold-set run.

A gold set of hand-curated series with known correct verdicts lives in
`tests/gold_set.tsv`, including deliberate negatives: cell-line-only series,
combination-therapy cohorts, on-treatment-only designs and underpowered
studies. Verdict accuracy against the gold set is the skill's benchmark.

## Future Extensions

- Additional repositories: ArrayExpress, cBioPortal, dbGaP, EGA
- Automated response-label harmonization across studies
- Survival analysis and biomarker replication across independent cohorts
- Cross-study meta-analysis of concordant molecular signals
- Automatic identification of associated publications and supplements
- Clinical-trial integration

Cross-study synthesis is deliberately out of scope for this skill. Batch
effects are confounded with study, so pooled analysis across GEO series
recovers laboratory of origin as readily as biology; synthesis belongs in a
dedicated skill with an explicit meta-analytic model.

## Citations

- Barrett T, et al. NCBI GEO: archive for functional genomics data sets.
  *Nucleic Acids Res.* 2013;41(D1):D991–D995.
- Sayers EW, et al. Database resources of the National Center for
  Biotechnology Information. *Nucleic Acids Res.* 2024.
- Eisenhauer EA, et al. New response evaluation criteria in solid tumours:
  revised RECIST guideline (version 1.1). *Eur J Cancer.* 2009;45(2):228–247.
