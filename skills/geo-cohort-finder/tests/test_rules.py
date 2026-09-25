"""Offline tests for the rules and scoring in geo_cohort_finder.py (no network)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geo_cohort_finder as g  # noqa: E402


@pytest.fixture(scope="module")
def rules():
    return g.Rules()


def sample(gsm="GSM1", source="HNSCC tumor", title="", **chars):
    return {"gsm_accession": gsm, "title": title, "source_name": source, "organism": "Homo sapiens",
            "platform": "GPL1", "characteristics": {k.replace("_", " "): v for k, v in chars.items()}}


def query(objective="response", drug="pembrolizumab", min_total=2, min_arm=1):
    return {"drug": drug, "disease": "HNSCC", "omics": None, "objective": objective,
            "defaults": {"min_patients_total": min_total, "min_patients_per_arm": min_arm, "max_series": 10}}


# --- response matching ------------------------------------------------------

def test_anti_pd1_treatment_is_not_progressive_disease(rules):
    s = g.annotate_sample(sample(treatment="anti-PD-1"), rules, None)
    assert s["response_field"] == ""
    assert s["response_mapped"] == ""


def test_pdl1_field_is_never_a_response_field(rules):
    assert not g.is_response_field("PD-L1 response", rules)
    assert not g.is_response_field("progression-free survival", rules)
    assert not g.is_response_field("response time (days)", rules)


def test_response_fields_recognized(rules):
    for key in ["response", "Best Overall Response", "RECIST", "bor", "response (RECIST 1.1)"]:
        assert g.is_response_field(key, rules), key


def test_whole_value_match_only(rules):
    assert g.map_response("PR", rules) == "PR"
    assert g.map_response("  partial response ", rules) == "PR"
    assert g.map_response("PR (confirmed)", rules) == "NEEDS_REVIEW"
    assert g.map_response("PRE", rules) == "NEEDS_REVIEW"
    assert g.map_response("RECIST: 1", rules) == "NEEDS_REVIEW"
    assert g.map_response("NA", rules) == "NOT_AVAILABLE"


# --- sample type precedence ---------------------------------------------------

@pytest.mark.parametrize("source,expected", [
    ("patient-derived xenograft", "model"),
    ("patient-derived organoid", "model"),
    ("FaDu", "cell_line"),
    ("HNSCC cell line", "cell_line"),
    ("tumor biopsy, patient 3", "patient_tumor"),
    ("peripheral blood", "patient_other"),
    ("adjacent normal mucosa", "patient_other"),
    ("Head and neck mucosal squamous cell carcinoma", "patient_tumor"),  # GSE268014
    ("H&N squamous cel carcinoma", "patient_tumor"),                     # GSE333537
    ("unlabelled", "NEEDS_REVIEW"),
])
def test_sample_type_precedence(rules, source, expected):
    assert g.classify_sample_type(sample(source=source), rules) == expected


def test_cell_line_key_with_line_id(rules):
    assert g.classify_sample_type(sample(source="x", cell_line="UM-SCC-1"), rules) == "cell_line"
    assert g.classify_sample_type(sample(source="HN-SCC-151", cell_line="HN-SCC-151"), rules) == "cell_line"
    assert g.classify_sample_type(sample(source="JHU-06", cell_line="JHU-06"), rules) == "cell_line"


def test_cell_line_key_misused_for_cell_type(rules):
    # GSE296954 / GSE296867: patient samples with "cell line: T cells"
    assert g.classify_sample_type(
        sample(source="tumor tissue", tissue="tumor tissue", cell_line="T cells"), rules) == "patient_tumor"
    assert g.classify_sample_type(
        sample(source="peripheral blood", tissue="peripheral blood", cell_line="T cells"), rules) == "patient_other"


# --- expansion ------------------------------------------------------------------

def test_hnscc_expansion_recorded_with_sources(rules):
    terms = rules.expand("disease", "HNSCC")
    assert terms[0] == {"term": "HNSCC", "source": "user"}
    got = {t["term"] for t in terms}
    assert {"head and neck cancer", "head and neck squamous cell carcinoma", "oral cancer"} <= got
    assert all(t["source"] == "rules/synonyms.tsv" for t in terms[1:])


def test_unknown_term_searches_user_term_only(rules):
    assert rules.expand("disease", "made-up cancer") == [{"term": "made-up cancer", "source": "user"}]


# --- verdicts ---------------------------------------------------------------------

def _annotated(rules, specs, drug="pembrolizumab"):
    drug_re = g._term_regex(["pembrolizumab", "Keytruda"]) if drug else None
    return [g.annotate_sample(s, rules, drug_re) for s in specs]


def test_suitable_response_cohort(rules):
    smp = _annotated(rules, [
        sample("GSM1", patient_id="P1", treatment="pembrolizumab", response="PR"),
        sample("GSM2", patient_id="P2", treatment="pembrolizumab", response="PD"),
    ])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["verdict"] == "SUITABLE"
    assert (row["n_responder_strict"], row["n_nonresponder_strict"]) == (1, 1)


def test_no_patient_ids_assumes_one_per_sample(rules):
    smp = _annotated(rules, [sample("GSM1", treatment="pembrolizumab", response="PR"),
                             sample("GSM2", treatment="pembrolizumab", response="PD")])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["n_patients"] == 2
    assert row["n_patients_basis"] == "assumed_one_per_sample"
    assert row["verdict"] == "SUITABLE"
    assert "assumed" in row["verdict_reason"]


def test_no_patient_ids_with_timepoints_is_unknown(rules):
    smp = _annotated(rules, [sample("GSM1", treatment="pembrolizumab", response="PR", timepoint="pre"),
                             sample("GSM2", treatment="pembrolizumab", response="PR", timepoint="on")])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["n_patients"] == "unknown"
    assert row["verdict"] == "NEEDS_REVIEW"


def test_no_patient_ids_single_cell_is_unknown(rules):
    smp = _annotated(rules, [sample("GSM1"), sample("GSM2")], drug=None)
    row = g.summarize_series("GSE1", {"summary": "single-cell RNA-seq of HNSCC"}, smp,
                             query(objective="discovery", drug=None))
    assert row["n_patients"] == "unknown"


def test_partial_patient_ids_is_unknown(rules):
    smp = _annotated(rules, [sample("GSM1", patient_id="P1"), sample("GSM2")], drug=None)
    row = g.summarize_series("GSE1", {}, smp, query(objective="discovery", drug=None))
    assert row["n_patients"] == "unknown"


def test_wrong_organism(rules):
    # GSE333739: mouse MOC1 tumours inside a human SuperSeries
    mouse = [dict(sample(f"GSM{i}", patient_id=f"M{i}"), organism="Mus musculus") for i in range(3)]
    smp = [g.annotate_sample(s, rules, None, "Homo sapiens") for s in mouse]
    row = g.summarize_series("GSE1", {}, smp, query(objective="discovery", drug=None))
    assert row["verdict"] == "WRONG_ORGANISM"


def test_wrong_sample_type_takes_precedence(rules):
    smp = _annotated(rules, [sample("GSM1", source="FaDu cell line", response="PR")])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["verdict"] == "WRONG_SAMPLE_TYPE"
    assert "THERAPY_UNCONFIRMED" not in row["verdict"]
    assert "drug not named" in row["verdict_reason"]  # all failures still listed


def test_discovery_does_not_require_response(rules):
    smp = _annotated(rules, [sample("GSM1", patient_id="P1"), sample("GSM2", patient_id="P2")], drug=None)
    row = g.summarize_series("GSE1", {}, smp, query(objective="discovery", drug=None))
    assert row["verdict"] == "SUITABLE"


def test_response_objective_without_labels(rules):
    smp = _annotated(rules, [sample("GSM1", patient_id="P1", treatment="pembrolizumab",
                                    overall_survival_months="12")])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["verdict"] == "SURVIVAL_ONLY"


def test_dcb_unknown_when_stable_disease_present(rules):
    smp = _annotated(rules, [
        sample("GSM1", patient_id="P1", treatment="pembrolizumab", response="SD"),
        sample("GSM2", patient_id="P2", treatment="pembrolizumab", response="PR"),
    ])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["n_responder_dcb"] == "unknown"


# --- main table -------------------------------------------------------------------

def test_main_table_has_no_patient_level_columns(tmp_path, rules):
    smp = _annotated(rules, [
        sample("GSM1", patient_id="P1", treatment="pembrolizumab", hpv="positive", response="PR"),
        sample("GSM2", patient_id="P2", treatment="pembrolizumab", hpv="negative", response="PD"),
        sample("GSM3", patient_id="P3", treatment="placebo", hpv="negative", response="PD"),
    ])
    row = g.summarize_series("GSE1", {}, smp, query())
    assert row["treatment"] == "treatment: pembrolizumab=2, placebo=1"
    assert "hpv: negative=2, positive=1" in row["sample_groups"]
    assert "patient id" not in row["metadata_fields"]

    g.write_main_table(tmp_path / "geo_cohorts.tsv", [row])
    header = (tmp_path / "geo_cohorts.tsv").read_text().splitlines()[0].split("\t")
    assert not any("patient_id" in h for h in header)
    assert {"gse_accession", "paper_url", "n_samples", "organism", "treatment"} <= set(header)


def test_html_report(tmp_path, rules):
    smp = _annotated(rules, [sample("GSM1", patient_id="P1"), sample("GSM2", patient_id="P2")], drug=None)
    with_paper = g.summarize_series("GSE1", {"title": "Study <one> & co"}, smp, query("discovery", None))
    with_paper.update(pubmed_id="12345", doi="10.1/x", ftp_link="https://geo/GSE1")
    no_paper = g.summarize_series("GSE2", {"title": "Study two"}, smp, query("discovery", None))
    no_paper.update(pubmed_id=None, ftp_link="https://geo/GSE2")
    record = {"question": "HNSCC transcriptomics?", "run_date": "2026-09-25",
              "query": {**query("discovery", None), "organism": "Homo sapiens", "sample_type": "patient_tumor"},
              "search_terms": {"disease": rules.expand("disease", "HNSCC")}, "gds_types": []}
    search = {"esearch_term": '"HNSCC"[All Fields]', "n_matched": 5, "n_returned": 5,
              "n_prefiltered": 2, "n_characterized": 2, "n_not_characterized": 0}
    g.write_html_report(tmp_path, record, search, [with_paper, no_paper], [])
    page = (tmp_path / "report.html").read_text()
    assert "https://pubmed.ncbi.nlm.nih.gov/12345/" in page
    assert ">N/A<" in page
    assert "Study &lt;one&gt; &amp; co" in page          # escaped
    assert "head and neck squamous cell carcinoma" in page  # term expansion shown
    assert "rules/synonyms.tsv" in page


# --- parsing ---------------------------------------------------------------------

MATRIX = (
    '!Series_title\t"A study"\n'
    '!Series_pubmed_id\t"12345"\n'
    '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
    '!Sample_source_name_ch1\t"tumor"\t"tumor"\n'
    '!Sample_characteristics_ch1\t"patient id: P1"\t"patient id: P2"\n'
    '!Sample_characteristics_ch1\t"response: CR"\t"response: PD"\n'
)


def test_parse_matrix_header():
    series, samples = g.parse_matrix_header(MATRIX)
    assert series["!Series_pubmed_id"] == ["12345"]
    assert [s["gsm_accession"] for s in samples] == ["GSM1", "GSM2"]
    assert samples[1]["characteristics"] == {"patient id": "P2", "response": "PD"}


def test_matrix_dir_url():
    assert g.matrix_dir_url("GSE100866").endswith("/GSE100nnn/GSE100866/matrix/")
    assert g.matrix_dir_url("GSE999").endswith("/GSEnnn/GSE999/matrix/")


def test_esearch_term_quotes_and_groups(rules):
    rec = {"search_terms": {"disease": rules.expand("disease", "HNSCC")},
           "query": {"organism": "Homo sapiens"}, "gds_types": []}
    term = g.build_esearch_term(rec)
    assert '"head and neck cancer"[All Fields]' in term
    assert '"gse"[Entry Type]' in term and '"Homo sapiens"[Organism]' in term
