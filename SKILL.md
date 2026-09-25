---
name: geo-cohort-finder
description: >-
  Find analyzable human patient cohorts in NCBI GEO from a natural-language
  research question. Records the question and its parsed query, expands drug
  and disease synonyms, searches GEO, parses sample-level metadata from GEO
  series matrix files, cross-checks the linked publication, and returns a
  dataframe of candidate studies with patient counts, sample type, treatment
  exposure and treatment-response labels. Discovery and triage only; hands
  accessions to downstream skills.
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
    - name: question
      type: string
      format:
        - text
      description: The user's natural-language research question, recorded verbatim
      required: true
    - name: query
      type: flags
      format:
        - cli
      description: Typed query converted from the question by the agent (--drug, --disease, --omics, --objective, ...)
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
      - requests>=2.31
  demo_data:
    - path: demo_data/cache/
      description: Cached network responses for the HNSCC transcriptomics demo query, offline
  endpoints:
    cli: python skills/geo-cohort-finder/geo_cohort_finder.py --question "{question}" --disease "{disease}" --output {output_dir}
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

1. Record the natural-language question verbatim alongside the typed query the
   agent converted it into.
2. Expand drug and disease terms from an editable synonym table **before** the
   search runs, and record every term used and where it came from.
3. Search GEO via NCBI E-utilities and retrieve series-level summaries.
4. Fetch sample-level metadata for surviving candidates from GEO series matrix
   files.
5. Cross-check each series against its linked publication where one exists.
6. Classify sample type, patient count, treatment exposure and response labels
   against explicit rules.
7. Emit a suitability verdict per series, with the reason and the verbatim
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

### Question → query conversion (recorded)

The agent converts the natural-language question into a typed query and passes
both to the script: the question verbatim via `--question`, the criteria via
flags. The script never re-interprets the question; it records it. Both are
written to `reproducibility/query.json` and printed in the report header, so a
reader can check the conversion.

```text
question:      "Analyze pembrolizumab response in HNSCC."
drug:          pembrolizumab
disease:       HNSCC
omics:         unspecified
objective:     response          (response | discovery; see Suitability verdicts)
organism:      Homo sapiens
sample_type:   patient_tumor
```

### Search term expansion (recorded)

Before the search runs, each drug and disease term is expanded to its synonyms.
Expansions come from the editable table `rules/synonyms.tsv`; the agent may add
further terms with `--extra-terms`. Every term is recorded with its source, and
the expansion is printed before any network call:

```text
disease: HNSCC
  HNSCC                                     user
  head and neck cancer                      rules/synonyms.tsv
  head and neck squamous cell carcinoma     rules/synonyms.tsv
  oral cancer                               rules/synonyms.tsv
  SCCHN                                     rules/synonyms.tsv
drug: pembrolizumab
  pembrolizumab                             user
  Keytruda                                  rules/synonyms.tsv
  MK-3475                                   rules/synonyms.tsv
```

The same list goes to `reproducibility/query.json` and to the exact search
string in `reproducibility/esearch_terms.txt`. A term not in the table and not
supplied by the agent is not searched.

## Workflow

```text
"Analyze pembrolizumab response in HNSCC"
        ▼
   LLM / Agent  ── question → typed query ──▶ --question + flags
        ▼
   geo_cohort_finder.py
        │
        ├── 0. record                   (question + typed query → query.json)
        ├── 1. expand                   (rules/synonyms.tsv, recorded)
        ├── 2. search                   (esearch, one call)
        ├── 3. summarize                (esummary, batched)
        ├── 4. characterize             (series matrix, first MAX_SERIES only)
        ├── 5. publications             (PubMed / Europe PMC cross-check)
        └── 6. score                    (rules tables)
        ▼
   output/  report.md · result.json · tables/ · reproducibility/
        ▼
   ── chains to ──▶ article-data-fetcher → nfcore-rnaseq-wrapper → rnaseq-de
```

### Stage boundaries

Stages have very different costs and failure modes. Every network response is
cached, so `--resume` re-runs all stages from the recorded query and the cache,
repeating only requests that never completed.

| Stage (`--stage`) | Method | Cost | Output |
|---|---|---|---|
| `record` | local | instant | `reproducibility/query.json` |
| `expand` | `rules/synonyms.tsv` | instant | `reproducibility/query.json`, `esearch_terms.txt` |
| `search` | `esearch` on `gds`, `GSE[ETYP]` | one call | `tables/candidates.tsv` |
| `summarize` | `esummary` on `gds`, batched | ~1 call per 100 | `tables/candidates.tsv` |
| `characterize` | series matrix file per series | seconds each | `tables/samples.tsv` |
| `publications` | PubMed + Europe PMC per linked PMID | ~3 calls per series | `tables/publications.tsv` |
| `score` | local rules | instant | `tables/cohorts.tsv` |

**Series cap.** Only the first `MAX_SERIES` candidates (default **10**) that
pass the prefilter are characterized and cross-checked. Series with a linked
PubMed ID are taken first, then the rest; within each group, in the order
`esearch` returns them (newest first). The ordering, the cap and the number of
candidates not characterized are all stated in the report. Raise the cap with
`--max-series` once the pipeline is validated.

**Prefilter.** Before the cap is applied, candidates are dropped only for
`n_samples < MIN_PATIENTS_TOTAL` (a patient contributes at least one sample, so
this cannot drop a qualifying cohort) and for being a SuperSeries whose
SubSeries are already in the candidate list.

**Rate limit.** No NCBI API key is used. All E-utilities calls are throttled to
at most 3 requests per second (≥ 0.34 s apart) and send the `tool` parameter
(and `email` if `--email` or `NCBI_EMAIL` is set), as NCBI requests.

Every network response is cached on disk keyed by accession. Re-runs do not
re-fetch.

## CLI Reference

```bash
# Standard usage: question recorded verbatim, typed query as flags
python geo_cohort_finder.py \
    --question "Analyze pembrolizumab response in HNSCC." \
    --drug pembrolizumab --disease HNSCC --objective response \
    --output output/

# Discovery without a drug
python geo_cohort_finder.py \
    --question "List GEO datasets on HNSC with transcriptomic data." \
    --disease HNSC --omics transcriptomic --objective discovery --min-patients 30

# Add search terms beyond the synonym table (recorded as source: agent)
python geo_cohort_finder.py --question "..." --disease HNSCC \
    --extra-terms "oropharyngeal carcinoma;laryngeal carcinoma" \
    --drug pembrolizumab --extra-drug-terms "anti-PD-1"

# Characterize more than the default 10 series
python geo_cohort_finder.py --question "..." --disease HNSCC --max-series 50

# Demo mode (cached records, no network)
python geo_cohort_finder.py --demo

# Resume: reuse the recorded query and every cached response in output/
python geo_cohort_finder.py --resume --output output/

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
`Homo sapiens[Organism]`. The expanded disease terms are joined by `OR`, as are
the expanded drug terms; the disease group and drug group are joined by `AND`.
Multi-word terms are quoted. Modality is mapped to a `[DataSet Type]` filter.
The exact string sent is saved to `reproducibility/esearch_terms.txt`.

Drug terms are **not** used to confirm treatment. A drug named in a series
summary is a candidate signal only. **Do not assume that a dataset is suitable
merely because its title or description mentions the requested therapy or
cancer.** Confirmation comes from sample-level metadata in the `characterize`
stage.

### Domain decisions

These are the scientific rules this skill encodes. They are versioned here,
not inferred at runtime. Changing behavior means editing these tables.

#### Configurable defaults

Printed in the report header on every run.

```text
MIN_PATIENTS_TOTAL      20
MIN_PATIENTS_PER_ARM    10
MAX_SERIES              10
NCBI_MAX_REQ_PER_SEC    3
EXCLUDE_SAMPLE_TYPES    cell line, PDX, organoid, xenograft
ORGANISM                Homo sapiens
```

#### Sample-type classification

Applied to sample-level metadata first, falling back to series text. Rows are
tested **in precedence order**; the first matching row wins. Model and cell-line
patterns come before patient patterns so that, for example, "patient-derived
xenograft" is classed as `model`, not `patient_tumor`.

| Precedence | Pattern in `source_name` / `characteristics` | Class |
|---|---|---|
| 1 | `PDX`, `patient-derived xenograft`, `xenograft`, `organoid` | `model` |
| 2 | a `cell line:` characteristic whose value is a line identifier (`HN-SCC-151`, `JHU-06`); the phrase `cell line` in a value; a name in `rules/cell_lines.tsv` | `cell_line` |
| 3 | `siRNA`, `shRNA`, `knockout`, `overexpress`, `treated ... for Nh` | `perturbation` |
| 4 | `peripheral blood`, `blood`, `PBMC`, `plasma`, `serum`, `saliva`, `adjacent normal`, `normal tissue/mucosa`, `lymph node` | `patient_other` |
| 5 | `tumor biopsy`, `patient`, `primary tumor`, `FFPE`, `resection`, `tumor`, `carcinoma`, `cancer`, `neoplasm` | `patient_tumor` |
| — | no match | `NEEDS_REVIEW` |

Patterns are matched against characteristic **values**, source name and title,
never against characteristic keys: submitters misuse keys (`cell line: T cells`
on patient tumor samples). A `cell line:` key counts only when its value looks
like a line identifier, not a cell type. `patient_tumor` and `patient_other`
both count as patient data; the split is reported in `sample_type`. Samples
whose stated organism differs from `--organism` are excluded before
classification.

Named cell lines (`CAL27`, `FaDu`, `SCC-25`, ...) live in `rules/cell_lines.tsv`,
not in this document, so other cancer types can be added without editing code.

#### Response label mapping

Response labels are matched **only within response fields**, and **only as a
whole value**. Substring matching across all characteristics is not used: it
would read `PD` inside `anti-PD-1` or `PD-L1` as progressive disease, and
`progression` inside `progression-free survival` as a response call.

1. **Find response fields.** A characteristic is a response field only if its
   key (the text before `:` in `characteristics_ch1`) matches, case-insensitively,
   one of: `response`, `best response`, `best overall response`, `BOR`,
   `RECIST`, `clinical response`, `treatment response`, `responder`,
   `response status`, `clinical benefit`. The list lives in
   `rules/response_fields.tsv`. Keys containing `survival`, `PFS`, `OS`, `time`,
   `days`, `months` or `PD-1`/`PD-L1` are never response fields.
2. **Match the whole value.** The value (the text after `:`), trimmed and
   case-folded, must equal a source string below exactly. `PR` matches
   `PR`; it does not match `PR (confirmed)` or `PRE`.
3. **Otherwise `NEEDS_REVIEW`.** A response field whose value does not match
   exactly is recorded verbatim with `response_mapped = NEEDS_REVIEW`. Extend
   the table; do not loosen the match.

The verbatim value and its field name are retained for every sample.

| Source string (whole value) | Mapped |
|---|---|
| `CR`, `complete response` | `CR` |
| `PR`, `partial response` | `PR` |
| `SD`, `stable disease` | `SD` |
| `PD`, `progressive disease` | `PD` |
| `R`, `responder`, `yes` (in a response field) | `RESPONDER_UNSPECIFIED` |
| `NR`, `non-responder`, `nonresponder`, `no` (in a response field) | `NONRESPONDER_UNSPECIFIED` |
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

Response duration is not parsed in this version, so the `dcb` counts are
`unknown` whenever any SD is present, and equal `strict` when there is none.
Study-defined labels (`RESPONDER_UNSPECIFIED`/`SENSITIVE` vs
`NONRESPONDER_UNSPECIFIED`/`RESISTANT`) fill the strict arms directly.
`NOT_AVAILABLE` values (`NA`, `not evaluable`, ...) are excluded from counts.

Stable disease is the pivot: a patient stable for 14 months falls in opposite
groups under the two rules, and published studies use both. **Do not silently
convert clinical outcome categories into binary response groups.** Any
harmonization is documented in the output and the reproducibility record.

Where a study defines its own response groups, those are preserved as the
primary labels and the derived groups are marked as derived.

#### Suitability verdicts

What counts as suitable depends on what was asked. The agent sets `--objective`
when converting the question:

| Objective | Set when the question... | Example |
|---|---|---|
| `response` | asks about response, resistance, responders, or benefit from a therapy | "pembrolizumab response in HNSCC" |
| `discovery` | asks only which data exist for a disease, modality or therapy | "GEO datasets on HNSC with transcriptomic data" |

Criteria are applied **in the precedence order below; the first verdict whose
condition holds is assigned.** Every failed criterion, not just the first, is
listed in `verdict_reason`.

| Order | Verdict | Condition | Applies to |
|---|---|---|---|
| 0 | `WRONG_ORGANISM` | no sample is from `--organism` (e.g. a mouse SubSeries of a human SuperSeries) | both |
| 1 | `WRONG_SAMPLE_TYPE` | every sample is `cell_line`, `model` or `perturbation` | both |
| 2 | `THERAPY_UNCONFIRMED` | a drug was requested and no sample-level field names it (series text only) | both, if drug given |
| 3 | `SURVIVAL_ONLY` | no response field, but a survival field is present | `response` |
| 4 | `NO_RESPONSE_LABELS` | no response field and no survival field | `response` |
| 5 | `NEEDS_REVIEW` | any sample type is `NEEDS_REVIEW`; or (`response`) any response value is `NEEDS_REVIEW`; or `n_patients` is `unknown` | both |
| 6 | `INSUFFICIENT_N` | `n_patients < MIN_PATIENTS_TOTAL`; or (`response`) an arm `< MIN_PATIENTS_PER_ARM` | both |
| 7 | `SUITABLE` | none of the above | both |

So for a `discovery` query, `SUITABLE` means: patient tumor samples, drug
confirmed if one was requested, the requested modality, and at least
`MIN_PATIENTS_TOTAL` patients. Response labels are not required; whether they
exist is still reported in `has_response_labels`. For a `response` query,
`SUITABLE` additionally requires parsed response labels with both arms of the
`strict` grouping at least `MIN_PATIENTS_PER_ARM`.

**An unknown patient count is `NEEDS_REVIEW`, never a pass.** `n_samples` is not
substituted for `n_patients` to reach a verdict.

**Metadata that cannot be established is recorded as `unknown` or `null`.
Missing metadata is never interpreted as a negative result.** A
`NO_RESPONSE_LABELS` verdict means the labels were not found in the deposited
metadata — they may exist in the paper or its supplements.

### Patient counts

`n_samples` is the deposited sample count. `n_patients` is derived in this
order, and `n_patients_basis` records which rule applied:

| Basis | Rule |
|---|---|
| `patient_id` | every patient sample has a patient/subject/donor ID; count distinct IDs |
| `assumed_one_per_sample` | **no** sample has a patient ID, **and** there is no timepoint field, **and** the series is not single-cell; count patient samples |
| `unknown` | anything else: some samples lack IDs, repeated sampling, or single-cell |

The assumption is stated in `verdict_reason` and the report whenever it is
used. Most GEO bulk studies deposit one sample per patient without an ID
field; the assumption fails for undeclared paired designs, which is why it is
never silent.

Paired designs (pre- and on-treatment biopsies) and single-cell studies both
inflate sample count relative to patient count, in opposite magnitudes.
Filtering on `n_samples` is therefore not a proxy for cohort size.

### Sample metadata source: GEO series matrix

Sample-level metadata comes from the series matrix file(s):

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnnnnn/GSExxxxxx/matrix/GSExxxxxx_series_matrix.txt.gz
```

where `GSEnnnnnn` is the accession with its last three digits replaced by
`nnn`. Only the header block (lines starting with `!`) is parsed; the
expression table is not read. From it the skill takes, per sample:
`!Sample_geo_accession`, `!Sample_title`, `!Sample_source_name_ch1`,
`!Sample_organism_ch1`, every `!Sample_characteristics_ch1` row (split into
`key: value`), `!Sample_platform_id`, and `!Series_pubmed_id`.

- **Multiple platforms.** A series run on more than one platform has one matrix
  file per platform (`GSExxxxxx-GPLyyyy_series_matrix.txt.gz`). All are fetched
  from the `matrix/` directory listing and their samples combined, with
  `platform` kept per sample.
- **SuperSeries** are detected from GEO's standard summary sentence ("This
  SuperSeries is composed of the SubSeries listed below"). They get verdict
  `SUPERSERIES` with the SubSeries listed in `verdict_reason`, are not counted
  against `MAX_SERIES`, and their SubSeries are characterized next, counted
  against `MAX_SERIES`.
- **Missing matrix** (rare, older or unusual series): fall back to the GEO text
  view `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSExxxxxx&targ=gsm&form=text&view=brief`
  and record `metadata_source = geo_text` instead of `series_matrix`.

### Publication cross-check

GEO metadata is often incomplete: response labels, treatment details and
patient counts frequently appear only in the paper. Where a series links to a
publication, the skill retrieves it and checks it against the deposited
metadata.

1. **Find the paper.** PubMed IDs come from `esummary` and
   `!Series_pubmed_id`. A series with none gets `pub_status = no_linked_pmid`
   and is otherwise unaffected.
2. **Retrieve.** For each PMID: PubMed record (title, journal, year, DOI,
   abstract) via `efetch`; PMCID via `elink`; if the article is open access,
   full text from Europe PMC
   (`https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML`).
   `pub_text_level` records what was obtained: `full_text`, `abstract_only`
   or `none`. Paywalled full text is never scraped.
3. **Check** (deterministic text search over what was retrieved, using the same
   expanded term lists):
   - `pub_mentions_accession`: the GSE accession appears in the text.
   - `pub_drug_mentioned`: any expanded drug term appears.
   - `pub_disease_mentioned`: any expanded disease term appears.
   - `pub_response_terms`: response vocabulary found (`RECIST`, `responder`,
     `objective response`, `CR`/`PR`/`SD`/`PD` as whole words).
   - `pub_patient_count_sentences`: sentences matching
     `\b\d+\s+(patients|subjects|participants|individuals)\b`, recorded verbatim.
     No number is parsed out of them automatically.
4. **Flag disagreements** in `pub_flags`, for example:
   - `THERAPY_IN_PAPER_ONLY`: verdict is `THERAPY_UNCONFIRMED` but the paper
     names the drug.
   - `RESPONSE_IN_PAPER_ONLY`: verdict is `NO_RESPONSE_LABELS` but the paper uses
     response vocabulary; labels may be in its supplements.
   - `ACCESSION_NOT_IN_PAPER`: full text retrieved but the GSE is not
     mentioned; the PubMed link may be to a related rather than source paper.

**The publication check annotates; it never changes a verdict.** Verdicts
reflect deposited GEO metadata so they can be reproduced exactly. Flags tell the
user where the paper may resolve a gap. The agent may read the retrieved text to
explain a flag, and must label anything it reports from that reading as its own
interpretation of the paper, with the sentence quoted.

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
question: "Analyze pembrolizumab response in HNSCC."
drug: pembrolizumab (+2 synonyms) · disease: HNSCC (+4 synonyms)
omics: unspecified · objective: response
defaults: MIN_PATIENTS_TOTAL=20, MIN_PATIENTS_PER_ARM=10, MAX_SERIES=10

## Search terms
disease: HNSCC [user] · head and neck cancer · head and neck squamous cell
carcinoma · oral cancer · SCCHN [rules/synonyms.tsv]
drug: pembrolizumab [user] · Keytruda · MK-3475 [rules/synonyms.tsv]

## Search
N series matched the query. M passed the prefilter.
10 characterized (first 10 in esearch order); M − 10 not characterized.

## Verdicts
SUITABLE                 1
INSUFFICIENT_N           1
NO_RESPONSE_LABELS       3
WRONG_SAMPLE_TYPE        3
THERAPY_UNCONFIRMED      1
NEEDS_REVIEW             1

## Suitable cohorts
| accession | paper | n_pat | n_samp | platform | response field | CR/PR/SD/PD | strict R/NR |
|-----------|-------|-------|--------|----------|----------------|-------------|-------------|
| GSEnnnnnn | PMID nnnnnnnn | 36 | 72 | Illumina | response | 3/8/11/14 | 11 / 25 |

## Publication flags
| accession | verdict | pub_text_level | flag |
|-----------|---------|----------------|------|
| GSEnnnnnn | NO_RESPONSE_LABELS | full_text | RESPONSE_IN_PAPER_ONLY |
| GSEnnnnnn | THERAPY_UNCONFIRMED | abstract_only | THERAPY_IN_PAPER_ONLY |

## Limitations
- Only 10 of M candidate series were characterized.
- 1 series carries unparsed response strings (NEEDS_REVIEW), listed in tables/samples.tsv.
- Verdicts reflect deposited GEO metadata; publication flags mark where the paper may add information.
````

## Output Structure

**The main output is one table, `geo_cohorts.tsv`: one row per GEO dataset,
no per-patient fields**, sorted with `SUITABLE` first. Columns:

```text
gse_accession, geo_link, title, paper_url, doi,
n_samples, n_patients, n_patients_basis, organism, platform, data_type,
sample_types, treatment, timepoints, sample_groups,
response_counts, n_responders, n_nonresponders, survival_data,
metadata_fields, verdict, verdict_reason, publication_flags
```

`treatment` and `sample_groups` give sample counts per value of each
characteristic (e.g. `treatment: pembrolizumab=20, placebo=12`,
`hpv: positive=12, negative=28`), skipping ID-like fields and fields with more
than 10 distinct values. `metadata_fields` lists every characteristic recorded
for the dataset (age, sex, HPV status, ...). Everything else below is
supporting detail.

```text
output/
├── geo_cohorts.tsv        MAIN TABLE, one row per dataset
├── report.html            self-contained HTML report: question, typed query, term
│                          expansion with sources, GEO query string, search counts,
│                          and the main table (sortable, filterable; paper link or N/A)
├── report.md
├── result.json
├── tables/                supporting detail
│   ├── candidates.tsv     one row per series, search stage
│   ├── cohorts.tsv        one row per series, scored
│   ├── samples.tsv        one row per sample, verbatim characteristics
│   └── publications.tsv   one row per linked PMID, cross-check results
├── cache/                 raw network responses (matrix headers, PubMed, Europe PMC)
└── reproducibility/
    ├── commands.sh
    ├── query.json         question verbatim, typed query, expanded terms with sources
    ├── esearch_terms.txt  exact search string sent to NCBI
    ├── environment.yml
    └── checksums.sha256
```

### `tables/cohorts.tsv`

```text
gse_accession, title, objective, metadata_source,
pubmed_id, pubmed_url, doi, pmcid,
n_samples, n_patients,
organism, platform, omics_type, processing_level,
drug_matched, drug_evidence_field, drug_evidence_string,
sample_type, biopsy_timing,
has_response_labels, response_field, response_counts_raw,
n_responder_strict, n_nonresponder_strict,
n_responder_dcb, n_nonresponder_dcb,
has_survival, verdict, verdict_reason,
pub_status, pub_text_level, pub_flags, ftp_link
```

`pubmed_url` is `https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/`; `doi` and
`pmcid` come from the PubMed record. All three are `null` when there is no
linked publication.

### `tables/publications.tsv`

```text
gse_accession, pubmed_id, pmcid, doi, title, journal, year,
pub_text_level, pub_mentions_accession, pub_drug_mentioned,
pub_disease_mentioned, pub_response_terms, pub_patient_count_sentences,
pub_flags
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
  "question": "Analyze pembrolizumab response in HNSCC.",
  "query": {
    "drug": "pembrolizumab",
    "disease": "HNSCC",
    "objective": "response",
    "omics": null,
    "defaults": {"min_patients_total": 20, "min_patients_per_arm": 10, "max_series": 10}
  },
  "search_terms": {
    "disease": [
      {"term": "HNSCC", "source": "user"},
      {"term": "head and neck cancer", "source": "rules/synonyms.tsv"},
      {"term": "head and neck squamous cell carcinoma", "source": "rules/synonyms.tsv"},
      {"term": "oral cancer", "source": "rules/synonyms.tsv"},
      {"term": "SCCHN", "source": "rules/synonyms.tsv"}
    ],
    "drug": [
      {"term": "pembrolizumab", "source": "user"},
      {"term": "Keytruda", "source": "rules/synonyms.tsv"},
      {"term": "MK-3475", "source": "rules/synonyms.tsv"}
    ]
  },
  "search": {
    "esearch_term": "...",
    "n_matched": null,
    "n_prefiltered": null,
    "n_characterized": 10,
    "n_not_characterized": null
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
      "verdict_reason": "therapy confirmed; 36 patients; 11/25 strict",
      "publication": {
        "pubmed_id": "nnnnnnnn",
        "pubmed_url": "https://pubmed.ncbi.nlm.nih.gov/nnnnnnnn/",
        "doi": "10.xxxx/xxxxx",
        "pub_text_level": "full_text",
        "pub_mentions_accession": true,
        "pub_flags": []
      }
    }
  ],
  "limitations": []
}
```

**Use `null` or `unknown` for information that cannot be established. Do not
fabricate values to complete the schema.**

## Human-Readable Report

`output/report.md` contains: the question verbatim and the typed query it was
converted to; the defaults in force; every search term with its source; match
counts and how many series were not characterized because of `MAX_SERIES`; a
verdict summary; the suitable cohorts table with paper links; the publication
flags; and a limitations section that explicitly reports
unavailable metadata, small cohorts, heterogeneous response definitions,
unsupported modalities and inaccessible data.

## Dependencies

- Python >= 3.10, `requests` (everything else is standard library);
  `pytest` for the tests
- Network access to `eutils.ncbi.nlm.nih.gov`, `ftp.ncbi.nlm.nih.gov` (over
  HTTPS) and `www.ebi.ac.uk` (Europe PMC)
- Optional contact email via `--email` or `NCBI_EMAIL`, sent to NCBI as asked
  in its usage policy

No credentials and no NCBI API key are used. E-utilities calls are throttled to
3 requests per second, the limit without a key. At the default
`MAX_SERIES = 10` a full run makes on the order of 50 network requests.

## Gotchas

- **GEO's search index is series-level.** Response labels, treatment and
  patient identifiers live in per-sample `characteristics_ch1` and are
  invisible to `esearch`. The `characterize` stage is not optional.
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
- **`PD` is also a drug target.** `anti-PD-1`, `PD-L1 status` and
  `PD-1 inhibitor` appear in treatment fields of exactly the cohorts this skill
  looks for. Substring matching would call every such sample progressive
  disease; hence response matching is restricted to response fields and whole
  values.
- **"Patient" appears in model names.** "Patient-derived xenograft" and
  "patient-derived organoid" contain `patient`; sample-type rules are applied in
  precedence order with models first.
- **The linked PubMed ID is not always the source paper.** Some series link a
  later reanalysis or a companion paper. `ACCESSION_NOT_IN_PAPER` flags this
  when full text is available.
- **Processing level varies.** Supplementary deposits may be raw counts, FPKM,
  TPM or already-normalized values. A count-based downstream tool cannot
  consume an FPKM-only deposit.
- **Rate limits bite during demos.** Without an API key NCBI allows 3 requests
  per second; exceeding it returns HTTP 429. Cache everything and throttle.

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

1. Convert the user's question into a typed query: drug, disease, omics,
   sample type, and `objective` (`response` or `discovery`).
2. Preserve unspecified criteria as unspecified.
3. Optionally propose extra search terms via `--extra-terms`; these are
   recorded with source `agent`.
4. Invoke `geo_cohort_finder.py` with the question verbatim (`--question`) and
   the typed query as flags, so the conversion is recorded.
5. Interpret the structured results, including publication flags.
6. Communicate supporting evidence, missing information and limitations.

**The agent must not infer that a dataset contains a particular modality,
treatment, response phenotype or clinical variable unless this is established
from retrieved metadata.** Classification, mapping and scoring are performed by
the rules tables above, not by the model. Given the same cached network
responses, identical input produces identical output; GEO itself changes over
time, so the run date is recorded.

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
maintained surface of this skill. The editable ones live in `rules/`:
`synonyms.tsv`, `cell_lines.tsv`, `response_fields.tsv` and
`response_labels.tsv`. Extend them rather than relaxing matching
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
- Parsing response labels and patient counts from publication supplements
- Raising `MAX_SERIES` beyond 10 once verdicts are validated against the gold set
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
