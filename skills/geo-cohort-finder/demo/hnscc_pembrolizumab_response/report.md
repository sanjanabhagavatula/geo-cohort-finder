# GEO Cohort Finder Report

Run: 2026-09-25T17:55:21+00:00 · skill v0.1.0

## Research Question
question: "Which HNSCC cohorts have pembrolizumab responders and non-responders?"
drug: pembrolizumab · disease: HNSCC · omics: unspecified · objective: response
defaults: MIN_PATIENTS_TOTAL=20, MIN_PATIENTS_PER_ARM=10, MAX_SERIES=10

## Search terms
- disease: HNSCC [user] · HNSC [rules/synonyms.tsv] · SCCHN [rules/synonyms.tsv] · head and neck cancer [rules/synonyms.tsv] · head and neck squamous cell carcinoma [rules/synonyms.tsv] · squamous cell carcinoma of the head and neck [rules/synonyms.tsv] · oral cancer [rules/synonyms.tsv]
- drug: pembrolizumab [user] · Keytruda [rules/synonyms.tsv] · MK-3475 [rules/synonyms.tsv] · lambrolizumab [rules/synonyms.tsv] · immune checkpoint inhibitor [rules/synonyms.tsv] · immune checkpoint blockade [rules/synonyms.tsv] · ICI [rules/synonyms.tsv] · anti-PD-1 [rules/synonyms.tsv] · anti-PD-L1 [rules/synonyms.tsv] · anti-CTLA-4 [rules/synonyms.tsv] · nivolumab [rules/synonyms.tsv] · atezolizumab [rules/synonyms.tsv] · durvalumab [rules/synonyms.tsv] · ipilimumab [rules/synonyms.tsv]

## Search
45 series matched the query (45 returned). 31 passed the prefilter (n_samples >= 20).
10 characterized (series with a linked paper first, then esearch order); 21 not characterized because of MAX_SERIES.

## Verdicts
- THERAPY_UNCONFIRMED: 4
- NO_RESPONSE_LABELS: 3
- SURVIVAL_ONLY: 2
- WRONG_SAMPLE_TYPE: 1
- SUPERSERIES: 1

## All characterized series

| accession | verdict | paper | n_pat | n_samp | sample types | reason |
|---|---|---|---|---|---|---|
| [GSE310856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310856) | THERAPY_UNCONFIRMED | [PMID 41749511](https://pubmed.ncbi.nlm.nih.gov/41749511/) | 11 | 33 | patient_other:33 | neither the drug nor its class named in any sample-level field; no response or survival field; 11 patients < 20 |
| [GSE301741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741) | THERAPY_UNCONFIRMED | [PMID 41923630](https://pubmed.ncbi.nlm.nih.gov/41923630/) | unknown | 58 | NEEDS_REVIEW:58 | neither the drug nor its class named in any sample-level field; no response or survival field; 58 samples with unclassified sample type; n_patients unknown (no patient IDs and repeated sampling or single-cell design) |
| [GSE307471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307471) | WRONG_SAMPLE_TYPE | [PMID 41482523](https://pubmed.ncbi.nlm.nih.gov/41482523/) | unknown | 33 | cell_line:33 | all samples are cell_line/model/perturbation; 9 samples from another organism excluded; neither the drug nor its class named in any sample-level field; no response or survival field |
| [GSE286827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286827) | NO_RESPONSE_LABELS | [PMID 41045934](https://pubmed.ncbi.nlm.nih.gov/41045934/) | 29 | 87 | NEEDS_REVIEW:87 | no response or survival field; 87 samples with unclassified sample type |
| [GSE291246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291246) | NO_RESPONSE_LABELS | [PMID 40314973](https://pubmed.ncbi.nlm.nih.gov/40314973/) | unknown | 35 | patient_tumor:29;patient_other:4;NEEDS_REVIEW:2 | no response or survival field; 2 samples with unclassified sample type; n_patients unknown (no patient IDs and repeated sampling or single-cell design) |
| [GSE288199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288199) | THERAPY_UNCONFIRMED | [PMID 40086437](https://pubmed.ncbi.nlm.nih.gov/40086437/) | 180 | 180 | patient_tumor:180 | neither the drug nor its class named in any sample-level field; no response or survival field |
| [GSE281729](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281729) | NO_RESPONSE_LABELS | [PMID 39585339](https://pubmed.ncbi.nlm.nih.gov/39585339/) | 67 | 67 | patient_tumor:67 | no response or survival field |
| [GSE277573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277573) | SURVIVAL_ONLY | [PMID 39558036](https://pubmed.ncbi.nlm.nih.gov/39558036/) | 100 | 100 | patient_other:100 | survival field present, no response field |
| [GSE259280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259280) | SURVIVAL_ONLY | [PMID 39049036](https://pubmed.ncbi.nlm.nih.gov/39049036/) | 34 | 375 | NEEDS_REVIEW:191;patient_tumor:184 | survival field present, no response field; 191 samples with unclassified sample type |
| [GSE234138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234138) | SUPERSERIES | [PMID 38383563](https://pubmed.ncbi.nlm.nih.gov/38383563/) | unknown | 114 |  | SuperSeries; see SubSeries GSE233980,GSE234136 |
| [GSE233980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233980) | THERAPY_UNCONFIRMED | [PMID 38383563](https://pubmed.ncbi.nlm.nih.gov/38383563/) | 23 | 23 | NEEDS_REVIEW:23 | neither the drug nor its class named in any sample-level field; no response or survival field; 23 samples with unclassified sample type |

## Publication flags

| accession | PMID | pub_text_level | flags |
|---|---|---|---|
| GSE310856 | 41749511 | full_text | THERAPY_IN_PAPER_ONLY |
| GSE301741 | 41923630 | full_text | THERAPY_IN_PAPER_ONLY |
| GSE286827 | 41045934 | full_text | RESPONSE_IN_PAPER_ONLY |
| GSE291246 | 40314973 | full_text | RESPONSE_IN_PAPER_ONLY |
| GSE288199 | 40086437 | full_text | THERAPY_IN_PAPER_ONLY |
| GSE281729 | 39585339 | full_text | RESPONSE_IN_PAPER_ONLY |
| GSE233980 | 38383563 | full_text | THERAPY_IN_PAPER_ONLY |

## Limitations
- Only 10 of 31 prefiltered series were characterized (MAX_SERIES=10).
- n_patients could not be established for 3 series; reported as unknown.
- n_patients for 4 series assumes one sample per patient (no patient IDs deposited); see n_patients_basis.
- Verdicts reflect deposited GEO metadata only; publication flags mark where the paper may add information.

_Research use only. Not a medical device; not clinical decision support._
