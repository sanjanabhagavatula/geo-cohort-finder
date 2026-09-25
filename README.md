# geo-cohort-finder

Find analyzable human patient cohorts in NCBI GEO — and say plainly which ones
are **not** analyzable, and why.

A [ClawBio](https://github.com/ClawBio/ClawBio) skill.

## The problem

GEO's search index covers series-level text only: title, summary, overall
design. The facts that decide whether a study is analyzable — is this a
patient or a cell line, how many people, were they treated, is there an
outcome — live in per-sample `characteristics_ch1` fields that the index
cannot reach.

So the GEO website cannot answer the question researchers actually ask.

A representative HNSCC transcriptomics query returns **784 series**. Of 304
sampled, 121 had ≥20 samples — and 40 of those matched cell-line keywords in
their title alone. Sample count is not patient count: the largest hit by
samples was a single-cell study of ~18 patients.

## What it does

Natural-language question → typed query → GEO search → sample-level metadata
→ a scored dataframe of candidate cohorts.

```bash
python geo_cohort_finder.py --demo
python geo_cohort_finder.py --drug pembrolizumab --disease HNSCC --omics rna-seq
python geo_cohort_finder.py --disease HNSC --omics transcriptomic --min-patients 30
```

Output: `report.md`, `result.json`, `tables/cohorts.tsv`, `tables/samples.tsv`,
`reproducibility/`.

## What makes it more than a search wrapper

**Synonym expansion.** Searching `HNSCC` misses studies written as "SCCHN" or
"oral cavity carcinoma"; searching `pembrolizumab` misses "Keytruda" and
"MK-3475". Expansion runs from an editable table.

**Columns the search index does not have.** Sample type, patient count vs
sample count, treatment exposure, response labels — all parsed from
sample-level metadata.

**Both response binarizations, side by side.** RECIST collapses two common
ways, and stable disease is the pivot. On `GSE159067` (CR 5, PR 6, SD 27,
PD 64) the strict rule gives 11 responders and the durable-benefit rule gives
38 — from identical data. The skill reports both and endorses neither.

**Evidence tiers for treatment.** Sample metadata often records only the class
("immunotherapy") and never the agent. That earns `SUITABLE_CLASS_ONLY`,
never a silent upgrade to a confirmed drug match.

## What it does not do

Discovery and triage only. It does not download expression data, process raw
reads, run differential expression, or synthesize across studies — those are
`article-data-fetcher`, `nfcore-rnaseq-wrapper`, `rnaseq-de` and
`pathway-enricher`.

A `SUITABLE` verdict means labels exist and counts clear the thresholds. It is
**not** a claim that the cohort is adequately powered. Public treatment-response
cohorts are typically 20–80 patients.

Missing metadata is never read as a negative result. `NEEDS_REVIEW` means
unparsed, not absent.

## Layout

```
skills/geo-cohort-finder/
├── SKILL.md              spec: rules tables, verdicts, gotchas
├── geo_cohort_finder.py
├── demo_data/            cached GEO records; --demo runs offline
└── tests/gold_set.tsv    hand-curated verdicts, the benchmark
```

## Requirements

Python ≥3.10, `requests`. An `NCBI_API_KEY` environment variable raises the
E-utilities rate limit from 3 to 10 requests/second.
