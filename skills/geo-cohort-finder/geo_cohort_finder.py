#!/usr/bin/env python3
"""geo-cohort-finder: find analyzable human patient cohorts in NCBI GEO.

Discovery and triage only. Emits accessions for downstream ClawBio skills.
All scientific decisions live in the RULES tables below and in SKILL.md --
never in model weights. Identical input produces identical output.
"""

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# --------------------------------------------------------------------------
# Configurable defaults -- printed in every report header.
# --------------------------------------------------------------------------
DEFAULTS = {
    "MIN_PATIENTS_TOTAL": 20,
    "MIN_PATIENTS_PER_ARM": 10,
    "MAX_SERIES_DEEP_FETCH": 40,
    "ORGANISM": "Homo sapiens",
}

# --------------------------------------------------------------------------
# RULES: synonym expansion. Extend these rather than loosening matching.
# --------------------------------------------------------------------------
DRUG_SYNONYMS = {
    "pembrolizumab": ["Keytruda", "MK-3475", "anti-PD-1", "pembro"],
    "nivolumab": ["Opdivo", "BMS-936558", "anti-PD-1"],
    "atezolizumab": ["Tecentriq", "MPDL3280A", "anti-PD-L1"],
    "cisplatin": ["cis-platinum", "CDDP", "platinum"],
    "cetuximab": ["Erbitux", "anti-EGFR"],
}

DISEASE_SYNONYMS = {
    "hnscc": ["head and neck squamous", "HNSCC", "SCCHN", "HNSC",
              "oral cavity carcinoma", "oropharyngeal carcinoma",
              "laryngeal carcinoma"],
    "hnsc": ["head and neck squamous", "HNSCC", "SCCHN", "HNSC"],
    "melanoma": ["melanoma", "cutaneous melanoma", "SKCM"],
    "bladder": ["bladder cancer", "urothelial carcinoma", "BLCA"],
    "lung": ["lung cancer", "NSCLC", "non-small cell lung", "LUAD", "LUSC"],
}

# Drug-class fallback. Sample metadata frequently records only the class
# ("immunotherapy", "anti-PD-1") and never the agent. Class evidence is
# reported as class_only -- weaker than a named drug, and never silently
# upgraded to one.
DRUG_CLASSES = {
    "pembrolizumab": ["immunotherapy", "immune checkpoint", "checkpoint inhibitor",
                      "anti-pd-1", "pd-1", "pd1", "ICI"],
    "nivolumab": ["immunotherapy", "immune checkpoint", "checkpoint inhibitor",
                  "anti-pd-1", "pd-1", "pd1", "ICI"],
    "atezolizumab": ["immunotherapy", "immune checkpoint", "anti-pd-l1", "pd-l1"],
    "cisplatin": ["chemotherapy", "platinum", "chemoradiation"],
    "cetuximab": ["anti-egfr", "egfr inhibitor", "targeted therapy"],
}

OMICS_FILTERS = {
    "rna-seq": "expression profiling by high throughput sequencing",
    "transcriptomic": "expression profiling by high throughput sequencing",
    "array": "expression profiling by array",
    "methylation": "methylation profiling by genome tiling array",
}

# --------------------------------------------------------------------------
# RULES: sample-type classification.
# --------------------------------------------------------------------------
SAMPLE_TYPE_RULES = [
    ("cell_line", r"cell line|cell-line|CAL27|CAL-27|FaDu|SCC-?\d|HN-?SCC-?\d|"
                  r"UM-?SCC|Detroit ?562|HeLa|passage \d"),
    ("model", r"\bPDX\b|xenograft|organoid|\bmouse\b|murine"),
    ("perturbation", r"siRNA|shRNA|knock ?down|knock ?out|CRISPR|overexpress|"
                     r"transfect|vector|treated with .{0,30}for \d+ ?h"),
    ("patient_tumor", r"tumor|tumour|biopsy|patient|primary|FFPE|resection|"
                      r"specimen|surgical|\bTx\b"),
]

# --------------------------------------------------------------------------
# RULES: response labels. Field name is matched loosely because real GEO
# fields look like "best response on immunotherapy (recist)", not "response".
# --------------------------------------------------------------------------
RESPONSE_FIELD_PAT = re.compile(
    r"resp|recist|benefit|outcome|sensitiv|resistan|irrecist", re.I)
SURVIVAL_FIELD_PAT = re.compile(r"\bos\b|\bpfs\b|surviv|progression.free", re.I)
TIMEPOINT_FIELD_PAT = re.compile(r"time ?point|\bstate\b|treatment status|visit", re.I)
PATIENT_ID_PAT = re.compile(r"patient|subject|donor|case ?id|participant", re.I)
# Fields whose names match above but which hold clinical values, not IDs.
PATIENT_ID_EXCLUDE = re.compile(
    r"diagnos|histolog|status|outcome|age|sex|gender|stage|grade|type|site", re.I)

RESPONSE_VALUE_RULES = [
    ("CR", r"^\s*(CR|complete response|complete remission)\s*$"),
    ("PR", r"^\s*(PR|partial response|partial remission)\s*$"),
    ("SD", r"^\s*(SD|stable disease|stable)\s*$"),
    ("PD", r"^\s*(PD|progressive disease|progression|progressive)\s*$"),
    # Merged categories used by real studies (Riaz et al. deposit "PRCR").
    ("PR_OR_CR", r"^\s*(PRCR|PR/CR|CR/PR|CRPR|PR or CR)\s*$"),
    ("SD_OR_PD", r"^\s*(SDPD|SD/PD|PD/SD)\s*$"),
    ("RESPONDER_UNSPECIFIED", r"^\s*(R|responder|response|yes|benefit|DCB)\s*$"),
    ("NONRESPONDER_UNSPECIFIED",
     r"^\s*(NR|non-?responder|no response|no|NDB)\s*$"),
    ("SENSITIVE", r"^\s*sensitive\s*$"),
    ("RESISTANT", r"^\s*resistant\s*$"),
    # Explicitly recorded as unknown by the submitter. Distinct from
    # NEEDS_REVIEW, which means this skill could not parse the string.
    ("UNKNOWN_RECORDED", r"^\s*(UNK|unknown|NA|N/?A|not evaluable|NE|not assessed)\s*$"),
]

TIMEPOINT_VALUE_RULES = [
    ("pre_treatment", r"pre-?tx|pre-?treat|baseline|\bB1\b|screening|naive"),
    ("on_treatment", r"on-?tx|on-?treat|during|cycle \d|week \d"),
    ("post_treatment", r"post-?tx|post-?treat|after|progression|relapse"),
]

RECIST = {"CR", "PR", "SD", "PD"}

# RULES: patient identifiers encoded in sample titles. Deposited titles carry
# the patient far more often than any characteristics field does, and unlike a
# count quoted in an abstract this yields a per-sample mapping -- so repeated
# timepoints from one patient collapse correctly.
PATIENT_TITLE_RULES = [
    r"\bpatient (?:number |no\.? |#)?(\d{1,4})\b",
    r"\b(?:subject|donor|case)[ _#-]?(\w{1,6})\b",
    r"^(Pt[ _-]?\d{1,4})\b",
    r"^(P\d{1,4})\b",
    r"^([A-Z]{1,4}[-_]?\d{1,4})[ _-]",
]


def patient_from_title(title):
    """Extract a patient identifier from a sample title, or '' if none."""
    for pat in PATIENT_TITLE_RULES:
        m = re.search(pat, title, re.I)
        if m:
            return re.sub(r"[ _-]", "", m.group(1)).upper()
    return ""


# --------------------------------------------------------------------------
# Query parsing
# --------------------------------------------------------------------------
def expand(term, table):
    """Expand a term via a synonym table; always includes the term itself."""
    if not term:
        return []
    key = term.strip().lower()
    out = [term]
    for k, syns in table.items():
        if k in key or key in k:
            out.extend(syns)
    return sorted(set(out), key=str.lower)


def parse_query(args):
    q = {
        "drug": args.drug,
        "disease": args.disease,
        "omics": args.omics,
        "objective": "treatment response" if args.drug else "dataset discovery",
        "organism": DEFAULTS["ORGANISM"],
        "min_patients": args.min_patients,
    }
    q["drug_synonyms"] = expand(args.drug, DRUG_SYNONYMS)
    q["disease_synonyms"] = expand(args.disease, DISEASE_SYNONYMS)
    q["drug_class_terms"] = DRUG_CLASSES.get(
        (args.drug or "").strip().lower(), []) if args.drug else []
    return q


def build_esearch_term(q):
    parts = []
    if q["disease_synonyms"]:
        parts.append("(" + " OR ".join(
            f'"{s}"[All Fields]' for s in q["disease_synonyms"]) + ")")
    parts.append(f'"{q["organism"]}"[Organism]')
    parts.append("GSE[Entry Type]")
    if q["omics"]:
        f = OMICS_FILTERS.get(q["omics"].lower())
        if f:
            parts.append(f'"{f}"[DataSet Type]')
    return " AND ".join(parts)


# --------------------------------------------------------------------------
# Network, with on-disk cache
# --------------------------------------------------------------------------
class Fetcher:
    def __init__(self, cache_dir, offline=False):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.api_key = os.environ.get("NCBI_API_KEY")
        self.calls = 0

    def _cached(self, key, fn):
        path = self.cache / key
        if path.exists():
            return path.read_text()
        if self.offline:
            raise SystemExit(
                f"offline mode: no cached record for {key}. "
                f"Run without --demo, or add the file to demo_data/.")
        text = fn()
        path.write_text(text)
        # 3 req/s without a key, 10 with one.
        time.sleep(0.12 if self.api_key else 0.36)
        self.calls += 1
        return text

    def _get(self, url):
        req = urllib.request.Request(url, headers={"User-Agent": "geo-cohort-finder/0.1"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.read().decode("utf-8", "replace")

    def esearch(self, term, retmax=400):
        key = "esearch_" + hashlib.md5(term.encode()).hexdigest()[:12] + ".json"
        params = {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"}
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
        return json.loads(self._cached(key, lambda: self._get(url)))

    def esummary(self, uids):
        key = "esummary_" + hashlib.md5(",".join(uids).encode()).hexdigest()[:12] + ".json"
        params = {"db": "gds", "id": ",".join(uids), "retmode": "json"}
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params)
        return json.loads(self._cached(key, lambda: self._get(url)))

    def pubmed(self, pmids):
        """Abstracts only. Full text is frequently paywalled; abstracts are not."""
        key = "pubmed_" + "_".join(pmids) + ".txt"
        params = {"db": "pubmed", "id": ",".join(pmids),
                  "rettype": "abstract", "retmode": "text"}
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params)
        return self._cached(key, lambda: self._get(url))

    def pmids_for(self, gse):
        """Resolve a GSE accession to its linked PubMed IDs via esummary."""
        res = self.esearch(f"{gse}[Accession] AND GSE[Entry Type]", retmax=1)
        uids = res["esearchresult"]["idlist"]
        if not uids:
            return []
        d = self.esummary(uids)["result"]
        return [str(x) for x in d[uids[0]].get("pubmedids", [])]

    def samples(self, gse):
        """Sample-level SOFT records. This is the only place response labels live."""
        key = f"{gse}_gsm.txt"
        params = {"acc": gse, "targ": "gsm", "form": "text", "view": "brief"}
        url = f"{GEO_ACC}?" + urllib.parse.urlencode(params)
        return self._cached(key, lambda: self._get(url))


# --------------------------------------------------------------------------
# Parsing sample metadata
# --------------------------------------------------------------------------
def parse_soft(text):
    """Split a SOFT dump into per-sample dicts of {field: value}."""
    samples, cur = [], None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gsm": line.split("=")[-1].strip(), "chars": [], "title": "",
                   "source": ""}
        elif cur is None:
            continue
        elif line.startswith("!Sample_characteristics_ch1"):
            v = line.split("=", 1)[-1].strip()
            if ":" in v:
                f, val = v.split(":", 1)
                cur["chars"].append((f.strip(), val.strip()))
            else:
                cur["chars"].append(("", v))
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[-1].strip()
        elif line.startswith("!Sample_source_name_ch1"):
            cur["source"] = line.split("=", 1)[-1].strip()
    if cur:
        samples.append(cur)
    return samples


def classify_sample_type(sample):
    blob = " ".join([sample["title"], sample["source"]] +
                    [f"{f}: {v}" for f, v in sample["chars"]])
    for label, pat in SAMPLE_TYPE_RULES:
        if re.search(pat, blob, re.I):
            return label
    return "NEEDS_REVIEW"


def map_response(value):
    for label, pat in RESPONSE_VALUE_RULES:
        if re.match(pat, value, re.I):
            return label
    return "NEEDS_REVIEW"


def map_timepoint(value):
    for label, pat in TIMEPOINT_VALUE_RULES:
        if re.search(pat, value, re.I):
            return label
    return "unknown"


def characterize(gse, soft_text, query):
    """Derive one cohort record from sample-level metadata."""
    samples = parse_soft(soft_text)
    rows, resp_field, tp_field, pid_field = [], None, None, None

    for f, _ in [c for s in samples for c in s["chars"]]:
        if resp_field is None and RESPONSE_FIELD_PAT.search(f):
            resp_field = f
        if tp_field is None and TIMEPOINT_FIELD_PAT.search(f):
            tp_field = f
        if (pid_field is None and PATIENT_ID_PAT.search(f)
                and not PATIENT_ID_EXCLUDE.search(f)):
            pid_field = f

    if pid_field:
        vals = [dict(s["chars"]).get(pid_field, "") for s in samples]
        nonblank = [v for v in vals if v]
        distinct = len(set(nonblank))
        too_uniform = len(samples) > 2 and distinct < max(2, 0.2 * len(samples))
        too_long = nonblank and (sum(len(v) for v in nonblank) / len(nonblank)) > 20
        if too_uniform or too_long:
            pid_field = None

    has_survival = any(SURVIVAL_FIELD_PAT.search(f)
                       for s in samples for f, _ in s["chars"])

    drug_hits = 0
    class_hits = 0
    drug_terms = [d.lower() for d in query.get("drug_synonyms", [])]
    class_terms = [c.lower() for c in query.get("drug_class_terms", [])]
    patients = set()

    for s in samples:
        d = dict(s["chars"])
        raw = "; ".join(f"{f}: {v}" for f, v in s["chars"])
        resp_raw = d.get(resp_field, "") if resp_field else ""
        pid = (d.get(pid_field, "") if pid_field else "") or \
            patient_from_title(s["title"])
        row = {
            "gse_accession": gse,
            "gsm_accession": s["gsm"],
            "patient_id": pid,
            "patient_id_source": ("characteristics" if pid_field and
                                  d.get(pid_field) else
                                  "title" if pid else "none"),
            "source_name": s["source"],
            "raw_characteristics": raw,
            "response_raw": resp_raw,
            "response_mapped": map_response(resp_raw) if resp_raw else "",
            "timepoint": map_timepoint(d.get(tp_field, "")) if tp_field else "unknown",
            "sample_type": classify_sample_type(s),
        }
        if row["patient_id"]:
            patients.add(row["patient_id"])
        blob = (raw + " " + s["title"] + " " + s["source"]).lower()
        if drug_terms and any(t in blob for t in drug_terms):
            drug_hits += 1
        if class_terms and any(t in blob for t in class_terms):
            class_hits += 1
        rows.append(row)

    counts = {}
    for r in rows:
        if r["response_mapped"] and r["response_mapped"] != "NEEDS_REVIEW":
            counts[r["response_mapped"]] = counts.get(r["response_mapped"], 0) + 1
    needs_review = sum(1 for r in rows if r["response_mapped"] == "NEEDS_REVIEW")
    recorded_unknown = counts.pop("UNKNOWN_RECORDED", 0)

    types = {}
    for r in rows:
        types[r["sample_type"]] = types.get(r["sample_type"], 0) + 1
    dominant_type = max(types, key=types.get) if types else "NEEDS_REVIEW"

    n_samples = len(rows)
    n_patients = len(patients) if patients else None

    strict_r = (counts.get("CR", 0) + counts.get("PR", 0)
                + counts.get("PR_OR_CR", 0)
                + counts.get("RESPONDER_UNSPECIFIED", 0)
                + counts.get("SENSITIVE", 0))
    strict_nr = (counts.get("SD", 0) + counts.get("PD", 0)
                 + counts.get("SD_OR_PD", 0)
                 + counts.get("NONRESPONDER_UNSPECIFIED", 0)
                 + counts.get("RESISTANT", 0))
    dcb_r = (counts.get("CR", 0) + counts.get("PR", 0)
             + counts.get("PR_OR_CR", 0) + counts.get("SD", 0))
    dcb_nr = counts.get("PD", 0)

    cohort = {
        "accession": gse,
        "n_samples": n_samples,
        "n_patients": n_patients,
        "n_patients_source": ("characteristics" if pid_field else
                              "sample titles" if n_patients else "unavailable"),
        "sample_type": dominant_type,
        "sample_type_breakdown": types,
        "mixed_cohort": len([t for t, c in types.items() if c >= 3]) > 1,
        "response_field": resp_field,
        "response_counts_raw": counts,
        "response_needs_review": needs_review,
        "response_recorded_unknown": recorded_unknown,
        "n_responder_strict": strict_r,
        "n_nonresponder_strict": strict_nr,
        "n_responder_dcb": dcb_r,
        "n_nonresponder_dcb": dcb_nr,
        "has_survival": has_survival,
        "therapy_confirmed": drug_hits > 0 if drug_terms else None,
        "therapy_evidence_level": ("drug_named" if drug_hits else
                                   "class_only" if class_hits else "none"),
        "therapy_evidence_n": drug_hits or class_hits,
        "timepoint_field": tp_field,
        "pmids": [],
        "pub_drug_evidence": "none",
        "pub_quote": "",
        "pub_n_patients_mentions": "",
    }
    cohort["verdict"], cohort["verdict_reason"] = score(cohort, query)
    return cohort, rows


PUB_N_PATIENTS_PAT = re.compile(
    r"(\d{2,4})\s+patients?\b", re.I)


def scan_publication(text, query):
    """Scan abstract text for drug evidence and a reported patient count.

    Publication evidence is a separate, weaker tier than sample-level
    metadata: it describes the study, not any individual deposited sample.
    It is never merged into the sample-level counts.
    """
    out = {"pub_drug_evidence": "none", "pub_quote": "",
           "pub_n_patients_mentions": ""}
    if not text:
        return out
    flat = " ".join(text.split())
    for term in query.get("drug_synonyms", []):
        m = re.search(r"[^.]*\b" + re.escape(term) + r"\b[^.]*\.", flat, re.I)
        if m:
            out["pub_drug_evidence"] = "drug_named"
            out["pub_quote"] = m.group(0).strip()[:300]
            break
    if out["pub_drug_evidence"] == "none":
        for term in query.get("drug_class_terms", []) + ["PD-1/PD-L1", "PD-L1"]:
            m = re.search(r"[^.]*" + re.escape(term) + r"[^.]*\.", flat, re.I)
            if m:
                out["pub_drug_evidence"] = "class_only"
                out["pub_quote"] = m.group(0).strip()[:300]
                break
    # Report every "N patients" mention rather than picking one. An abstract
    # routinely covers several cohorts (e.g. "102 and 82 patients with HNSCC
    # or NSCLC"), so any single extracted number is as likely to belong to a
    # different cohort as to this series. These are unverified mentions, not
    # a patient count, and are never used in scoring.
    counts = [x for x in PUB_N_PATIENTS_PAT.findall(flat)]
    if counts:
        out["pub_n_patients_mentions"] = ",".join(dict.fromkeys(counts))
    return out


def score(c, query):
    """RULES: suitability verdicts. See SKILL.md for the criteria table."""
    eff_n = c["n_patients"] if c["n_patients"] else c["n_samples"]
    n_label = "patients" if c["n_patients"] else "samples (patient count unknown)"

    if c["sample_type"] in ("cell_line", "model", "perturbation"):
        return "WRONG_SAMPLE_TYPE", f"dominant sample type is {c['sample_type']}"
    if c["sample_type"] == "NEEDS_REVIEW":
        return "NEEDS_REVIEW", "sample type could not be classified from metadata"

    has_resp = bool(c["response_counts_raw"])
    n_labeled = sum(c["response_counts_raw"].values())
    n_unparsed = c.get("response_needs_review", 0)
    if has_resp and n_unparsed > 0.1 * max(n_labeled, 1):
        return ("NEEDS_REVIEW",
                f"{n_unparsed} of {n_labeled + n_unparsed} response strings "
                f"are unparsed; arm counts would be incomplete. Extend the "
                f"response mapping table before using this cohort.")
    drug_asked = query.get("drug") is not None

    if drug_asked and c["therapy_evidence_level"] == "none":
        if c.get("pub_drug_evidence", "none") != "none":
            return ("THERAPY_PUBLICATION_ONLY",
                    f"treatment described in the linked publication "
                    f"({c['pub_drug_evidence']}) but absent from sample-level "
                    f"metadata; samples cannot be assigned to arms")
        return ("THERAPY_UNCONFIRMED",
                "neither the drug nor its class appears in sample-level "
                "metadata; series text alone is not confirmation")
    if not has_resp:
        if c["has_survival"]:
            return "SURVIVAL_ONLY", "survival recorded but no response call found"
        return ("NO_RESPONSE_LABELS",
                "no response annotation found in deposited sample metadata")
    if eff_n < DEFAULTS["MIN_PATIENTS_TOTAL"]:
        return "INSUFFICIENT_N", f"{eff_n} {n_label} < {DEFAULTS['MIN_PATIENTS_TOTAL']}"
    if min(c["n_responder_strict"], c["n_nonresponder_strict"]) < \
            DEFAULTS["MIN_PATIENTS_PER_ARM"]:
        return ("INSUFFICIENT_N",
                f"arms {c['n_responder_strict']}/{c['n_nonresponder_strict']} "
                f"below {DEFAULTS['MIN_PATIENTS_PER_ARM']} per arm (strict rule)")
    if not drug_asked:
        verdict = "SUITABLE_NO_DRUG_FILTER"
    elif c["therapy_evidence_level"] == "drug_named":
        verdict = "SUITABLE"
    else:
        verdict = "SUITABLE_CLASS_ONLY"
    return verdict, (f"{eff_n} {n_label}; strict {c['n_responder_strict']}/"
                     f"{c['n_nonresponder_strict']}; dcb {c['n_responder_dcb']}/"
                     f"{c['n_nonresponder_dcb']}")


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------
def write_tsv(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                for c in cols) + "\n")


def render_report(q, term, stats, cohorts):
    L = []
    L.append("# GEO Cohort Finder Report\n")
    L.append("## Research Question\n")
    L.append(f"- **drug**: {q['drug'] or 'unspecified'}"
             + (f" (+{len(q['drug_synonyms']) - 1} synonyms)" if q["drug"] else ""))
    L.append(f"- **disease**: {q['disease'] or 'unspecified'}"
             + (f" (+{len(q['disease_synonyms']) - 1} synonyms)" if q["disease"] else ""))
    L.append(f"- **omics**: {q['omics'] or 'unspecified'}")
    L.append(f"- **objective**: {q['objective']}")
    L.append(f"- **defaults**: " + ", ".join(f"{k}={v}" for k, v in DEFAULTS.items()))
    L.append("\n## Search\n")
    L.append(f"```\n{term}\n```\n")
    L.append(f"- {stats['n_matched']} series matched the broad query")
    L.append(f"- {stats['n_deep_fetched']} deep-fetched for sample-level metadata\n")

    L.append("## Verdicts\n")
    tally = {}
    for c in cohorts:
        tally[c["verdict"]] = tally.get(c["verdict"], 0) + 1
    L.append("| verdict | n |\n|---|---|")
    for v, n in sorted(tally.items(), key=lambda x: -x[1]):
        L.append(f"| `{v}` | {n} |")

    good = [c for c in cohorts if c["verdict"].startswith("SUITABLE")]
    L.append("\n## Suitable cohorts\n")
    if not good:
        L.append("_None met the criteria. See limitations._")
    else:
        L.append("| accession | n_samples | n_patients | sample type | response field "
                 "| raw counts | strict R/NR | DCB R/NR |\n"
                 "|---|---|---|---|---|---|---|---|")
        for c in good:
            raw = ", ".join(f"{k}:{v}" for k, v in sorted(c["response_counts_raw"].items()))
            L.append(f"| {c['accession']} | {c['n_samples']} | "
                     f"{c['n_patients'] or 'unknown'} | {c['sample_type']} | "
                     f"`{c['response_field']}` | {raw} | "
                     f"{c['n_responder_strict']}/{c['n_nonresponder_strict']} | "
                     f"{c['n_responder_dcb']}/{c['n_nonresponder_dcb']} |")

    L.append("\n## All screened series\n")
    L.append("| accession | n_samples | sample type | verdict | reason |\n|---|---|---|---|---|")
    for c in sorted(cohorts, key=lambda x: x["verdict"]):
        L.append(f"| {c['accession']} | {c['n_samples']} | {c['sample_type']} | "
                 f"`{c['verdict']}` | {c['verdict_reason']} |")

    L.append("\n## Limitations\n")
    nr = sum(c["response_needs_review"] for c in cohorts)
    unk = sum(1 for c in cohorts if c["n_patients"] is None)
    mixed = [c["accession"] for c in cohorts if c["mixed_cohort"]]
    L.append(f"- `n_patients` could not be established for {unk} of {len(cohorts)} "
             "series; reported as unknown, never assumed equal to `n_samples`.")
    L.append(f"- {nr} sample-level response strings were unparsed (`NEEDS_REVIEW`). "
             "Unparsed is not absent.")
    if mixed:
        L.append(f"- Mixed cohorts detected ({', '.join(mixed)}): one series bundles "
                 "several sample types. A single verdict may not apply to all samples.")
    L.append("- Verdicts reflect deposited metadata only. Response labels may exist "
             "in the source publication or its supplements and not in GEO.")
    L.append("- Response group counts depend on the binarization rule. Both are "
             "reported; this skill does not choose between them.")
    L.append("- A `SUITABLE` verdict means labels exist and counts clear the "
             "thresholds. It is not a claim that the cohort is adequately powered.")
    L.append("\n## Next steps\n")
    L.append("Chain to `article-data-fetcher` to retrieve deposited files, then "
             "`rnaseq-de` for the responder vs non-responder contrast.")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Find analyzable patient cohorts in GEO.")
    p.add_argument("--query", help="natural-language question (agent fills slots)")
    p.add_argument("--drug")
    p.add_argument("--disease")
    p.add_argument("--omics")
    p.add_argument("--min-patients", type=int, default=DEFAULTS["MIN_PATIENTS_TOTAL"])
    p.add_argument("--max-fetch", type=int, default=DEFAULTS["MAX_SERIES_DEEP_FETCH"])
    p.add_argument("--accessions", help="comma-separated GSEs; skips the search stage")
    p.add_argument("--output", default="output")
    p.add_argument("--with-publications", action="store_true",
                   help="fetch linked PubMed abstracts for supporting evidence")
    p.add_argument("--demo", action="store_true", help="cached records, no network")
    args = p.parse_args()

    here = Path(__file__).parent
    if args.demo:
        args.drug = args.drug or "pembrolizumab"
        args.disease = args.disease or "HNSCC"
        args.omics = args.omics or "rna-seq"
        args.accessions = args.accessions or "GSE159067,GSE200996,GSE288199"
        args.with_publications = True
        cache = here / "demo_data"
    else:
        cache = Path(args.output) / ".cache"

    DEFAULTS["MIN_PATIENTS_TOTAL"] = args.min_patients
    DEFAULTS["MAX_SERIES_DEEP_FETCH"] = args.max_fetch

    q = parse_query(args)
    term = build_esearch_term(q)
    fetcher = Fetcher(cache, offline=args.demo)

    print("Parsed query:")
    for k in ("drug", "disease", "omics", "objective", "organism"):
        print(f"  {k:12s} {q[k]}")
    if q["drug_synonyms"]:
        print(f"  {'expanded':12s} {', '.join(q['drug_synonyms'])}")
    if q["disease_synonyms"]:
        print(f"  {'expanded':12s} {', '.join(q['disease_synonyms'])}")
    print()

    # Stage 1-2: discovery
    if args.accessions:
        accs = [a.strip() for a in args.accessions.split(",")]
        n_matched = len(accs)
        print(f"Using {n_matched} supplied accessions (search stage skipped).")
    else:
        print(f"esearch: {term}")
        res = fetcher.esearch(term)["esearchresult"]
        n_matched = int(res["count"])
        uids = res["idlist"]
        print(f"  {n_matched} series matched; retrieving summaries for {len(uids)}")
        accs, summaries = [], {}
        for i in range(0, len(uids), 100):
            d = fetcher.esummary(uids[i:i + 100])["result"]
            for u in d["uids"]:
                summaries[d[u]["accession"]] = d[u]
        ranked = sorted(summaries.values(),
                        key=lambda r: -int(r.get("n_samples", 0)))
        accs = [r["accession"] for r in ranked][:args.max_fetch]
        print(f"  deep-fetching top {len(accs)} by sample count")

    # Stage 3-4: deep fetch and scoring
    cohorts, sample_rows = [], []
    for i, gse in enumerate(accs, 1):
        try:
            soft = fetcher.samples(gse)
            c, rows = characterize(gse, soft, q)
        except Exception as e:                      # noqa: BLE001
            print(f"  [{i}/{len(accs)}] {gse}: FETCH FAILED ({e})")
            continue
        if args.with_publications:
            try:
                pmids = fetcher.pmids_for(gse)
                if pmids:
                    c["pmids"] = pmids
                    c.update(scan_publication(fetcher.pubmed(pmids), q))
            except Exception as e:                  # noqa: BLE001
                c["pub_quote"] = f"lookup failed: {e}"
            c["verdict"], c["verdict_reason"] = score(c, q)
        cohorts.append(c)
        sample_rows.extend(rows)
        print(f"  [{i}/{len(accs)}] {gse:<12} n={c['n_samples']:<5} "
              f"{c['sample_type']:<14} {c['verdict']}")

    # Output
    out = Path(args.output)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "reproducibility").mkdir(parents=True, exist_ok=True)

    write_tsv(out / "tables" / "cohorts.tsv", cohorts, [
        "accession", "n_samples", "n_patients", "n_patients_source", "sample_type", "mixed_cohort",
        "therapy_confirmed", "therapy_evidence_level", "therapy_evidence_n",
        "response_field",
        "response_counts_raw", "response_needs_review", "response_recorded_unknown",
        "n_responder_strict", "n_nonresponder_strict",
        "n_responder_dcb", "n_nonresponder_dcb",
        "has_survival", "timepoint_field", "pmids", "pub_drug_evidence",
        "pub_n_patients_mentions", "pub_quote", "verdict", "verdict_reason"])
    write_tsv(out / "tables" / "samples.tsv", sample_rows, [
        "gse_accession", "gsm_accession", "patient_id", "source_name",
        "patient_id_source", "raw_characteristics", "response_raw", "response_mapped",
        "timepoint", "sample_type"])

    stats = {"n_matched": n_matched, "n_deep_fetched": len(cohorts)}
    (out / "report.md").write_text(render_report(q, term, stats, cohorts))
    (out / "result.json").write_text(json.dumps(
        {"query": q, "search": {"esearch_term": term, **stats},
         "cohorts": cohorts,
         "limitations": [
             "n_patients is unknown where no patient identifier is deposited.",
             "NEEDS_REVIEW means unparsed, not absent.",
             "Verdicts reflect deposited metadata only, not source publications.",
             "Both response binarizations are reported; neither is endorsed.",
         ]}, indent=2, default=str))
    (out / "reproducibility" / "commands.sh").write_text(
        "#!/usr/bin/env bash\n# exact replay\n"
        + " ".join(f'"{a}"' if " " in a else a for a in sys.argv) + "\n")
    (out / "reproducibility" / "query.json").write_text(json.dumps(q, indent=2))
    (out / "reproducibility" / "esearch_terms.txt").write_text(term + "\n")

    tally = {}
    for c in cohorts:
        tally[c["verdict"]] = tally.get(c["verdict"], 0) + 1
    print("\nVerdicts: " + ", ".join(f"{v}={n}" for v, n in sorted(tally.items())))
    print(f"Wrote {out}/report.md, result.json, tables/, reproducibility/")


if __name__ == "__main__":
    main()
