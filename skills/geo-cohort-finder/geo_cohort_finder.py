#!/usr/bin/env python3
"""geo-cohort-finder: find analyzable human patient cohorts in NCBI GEO.

Stages (see SKILL.md): record -> expand -> search -> summarize ->
characterize -> publications -> score. Every network response is cached on
disk; classification and verdicts come from the tables in rules/, not from a
model.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import os
import platform
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import requests

SKILL_DIR = Path(__file__).resolve().parent
RULES_DIR = SKILL_DIR / "rules"
DEMO_CACHE = SKILL_DIR / "demo_data" / "cache"
VERSION = "0.1.0"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
TOOL_NAME = "geo-cohort-finder"

DEFAULTS = {
    "min_patients_total": 20,
    "min_patients_per_arm": 10,
    "max_series": 10,
    "ncbi_max_req_per_sec": 3,
    "organism": "Homo sapiens",
}

SUPERSERIES_SENTENCE = "this superseries is composed of the subseries"

DEMO_ARGS = {
    "question": "List GEO datasets on HNSCC with transcriptomic data.",
    "disease": "HNSCC",
    "omics": "transcriptomic",
    "objective": "discovery",
}


# --------------------------------------------------------------------------
# Rules tables
# --------------------------------------------------------------------------

def _read_tsv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _norm(text: str) -> str:
    """Case-fold, trim, collapse whitespace and underscores."""
    return re.sub(r"[\s_]+", " ", text.strip().lower())


def _term_regex(terms: list[str]) -> Optional[re.Pattern]:
    """Whole-word, case-insensitive regex matching any of the terms."""
    if not terms:
        return None
    alts = "|".join(re.escape(t) for t in sorted(terms, key=len, reverse=True))
    return re.compile(rf"(?<![A-Za-z0-9])(?:{alts})(?![A-Za-z0-9])", re.IGNORECASE)


class Rules:
    def __init__(self, rules_dir: Path = RULES_DIR):
        self.synonyms = _read_tsv(rules_dir / "synonyms.tsv")
        self.cell_lines = [r["name"] for r in _read_tsv(rules_dir / "cell_lines.tsv")]
        self.response_fields = [_norm(r["field"]) for r in _read_tsv(rules_dir / "response_fields.tsv")]
        self.response_labels = {_norm(r["source"]): r["mapped"]
                                for r in _read_tsv(rules_dir / "response_labels.tsv")}
        self.omics = {_norm(r["omics"]): [t.strip() for t in r["gds_types"].split(";")]
                      for r in _read_tsv(rules_dir / "omics.tsv")}
        # Sample metadata often records only the therapeutic class and never
        # the agent, so a class match is tracked as its own weaker tier.
        self.drug_classes = {_norm(r["drug"]): [c.strip() for c in r["class_terms"].split(";") if c.strip()]
                             for r in _read_tsv(rules_dir / "drug_classes.tsv")}
        self.cell_line_re = _term_regex(self.cell_lines)

    def expand(self, kind: str, term: str) -> list[dict]:
        """Return [{term, source}] for a user term; user term first."""
        out = [{"term": term, "source": "user"}]
        key = _norm(term)
        for row in self.synonyms:
            if row["kind"] != kind:
                continue
            syns = [s.strip() for s in row["synonyms"].split(";") if s.strip()]
            if key == _norm(row["canonical"]) or key in {_norm(s) for s in syns}:
                for s in syns:
                    if _norm(s) not in {_norm(o["term"]) for o in out}:
                        out.append({"term": s, "source": "rules/synonyms.tsv"})
        return out


# Sample-type rules, tested in precedence order; first match wins. Patterns
# are matched against characteristic values, source name and title, never
# against characteristic keys (submitters misuse keys, e.g. "cell line: T cells").
SAMPLE_TYPE_RULES = [
    ("model", re.compile(r"\b(pdx|patient[- ]derived xenografts?|xenografts?|organoids?)\b", re.I)),
    ("cell_line", re.compile(r"\bcell[ _-]?lines?\b", re.I)),
    ("perturbation", re.compile(
        r"\b(sirna|shrna|knock[- ]?out|overexpress\w*)\b|\btreated\b.{0,60}\bfor\s+\d+\s*h(ours?|rs?)?\b", re.I)),
    ("patient_other", re.compile(
        r"\b(peripheral blood|whole blood|blood|pbmcs?|plasma|serum|saliva|"
        r"adjacent normal|normal adjacent|normal (tissue|mucosa|epithelium)|lymph nodes?)\b", re.I)),
    ("patient_tumor", re.compile(
        r"\b(tumou?r biops(y|ies)|patients?|primary tumou?rs?|ffpe|resections?|tumou?rs?|biops(y|ies)|"
        r"carcinomas?|cancers?|neoplasms?|tumou?r tissue)\b", re.I)),
]
PATIENT_TYPES = {"patient_tumor", "patient_other"}
CELL_LINE_KEY = re.compile(r"^cell[ _-]?line$", re.I)
# A "cell line:" value counts only if it looks like a line identifier
# (HN-SCC-151, JHU-06, HSC-2, MOC1), not a cell type ("T cells").
LINE_ID_VALUE = re.compile(r"^(?=.*\d)[A-Za-z0-9][A-Za-z0-9:\-/. ]{1,20}$")
GENERIC_CELL_VALUE = re.compile(r"\b(cells?|primary|patient|tissue|blood|none|n/?a)\b", re.I)
SINGLE_CELL = re.compile(r"single[- ]cell|scrna|snrna|\b10x\b|chromium", re.I)
PATIENT_ID_KEY = re.compile(
    r"^(patient|subject|donor|individual|case)([ _-]?(id|identifier|number|no\.?|#))?$", re.I)
# Hard excludes: quantities that are not a response call, whatever else the
# field name says. "pd-1"/"pd-l1" are deliberately NOT here -- Hugo et al.
# (GSE78220) name the field "anti-pd-1 response", a genuine RECIST call, and
# excluding the marker name discarded it. A field naming the marker without a
# response token fails the token test anyway.
RESPONSE_KEY_EXCLUDE = re.compile(
    r"surviv|\bpfs\b|\bos\b|\bdfs\b|\btime\b|\bdays?\b|\bmonths?\b|\bdate\b|duration", re.I)
SURVIVAL_KEY = re.compile(
    r"surviv|\bos\b|\bpfs\b|\bdfs\b|vital[ _]status|\bdeath\b|\bdeceased\b|follow[ -]?up", re.I)
TIMEPOINT_KEY = re.compile(
    r"time ?point|treatment[ _]status|biopsy[ _]timing|\bvisit\b|\bcycle\b|pre[/_-]?post", re.I)
TREATMENT_KEY = re.compile(r"treat|therap|\bdrug|\bagent|regimen|\barm\b|intervention", re.I)
ID_LIKE_KEY = re.compile(r"\b(id|identifier|barcode|sample name|title|accession)\b", re.I)
PATIENT_COUNT_SENTENCE = re.compile(
    r"[^.]*\b\d+\s+(patients|subjects|participants|individuals)\b[^.]*\.", re.I)
PUB_RESPONSE_TERMS = ["RECIST", "responder", "responders", "non-responder", "non-responders",
                      "objective response"]
# Bare RECIST abbreviations: case-sensitive, and not part of "PD-1"/"PD-L1".
PUB_RECIST_ABBREV = re.compile(r"(?<![A-Za-z0-9-])(CR|PR|SD|PD)(?![A-Za-z0-9-])")


# --------------------------------------------------------------------------
# Network: throttled, cached, optionally offline
# --------------------------------------------------------------------------

class CacheMiss(RuntimeError):
    pass


class Fetcher:
    """GET with a global throttle (<= 3 req/s, no API key) and a disk cache."""

    def __init__(self, cache_dir: Path, email: Optional[str], offline: bool = False,
                 max_per_sec: float = 3.0):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.email = email
        self.offline = offline
        self.min_interval = 1.0 / max_per_sec + 0.01
        self._last = 0.0
        self.session = requests.Session()
        self.session.headers["User-Agent"] = f"{TOOL_NAME}/{VERSION}"
        self.n_requests = 0

    def _key(self, url: str, params: Optional[dict]) -> Path:
        blob = url + "?" + json.dumps(params or {}, sort_keys=True)
        return self.cache_dir / hashlib.sha1(blob.encode()).hexdigest()

    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()

    def _params(self, url: str, params: Optional[dict]) -> Optional[dict]:
        if url.startswith(EUTILS):
            params = dict(params or {})
            params["tool"] = TOOL_NAME
            if self.email:
                params["email"] = self.email
        return params

    def get_text(self, url: str, params: Optional[dict] = None,
                 gz_header_only: bool = False) -> Optional[str]:
        """Return response text, or None for 404. Cached, including misses.

        gz_header_only: stream a gzipped series matrix and keep only the
        metadata header, stopping at !series_matrix_table_begin.
        """
        path = self._key(url, params)
        if path.with_suffix(".miss").exists():
            return None
        if path.exists():
            return path.read_text(encoding="utf-8")
        if self.offline:
            raise CacheMiss(f"not in cache (offline mode): {url} {params or ''}")

        for attempt in range(4):
            self._throttle()
            self.n_requests += 1
            try:
                resp = self.session.get(url, params=self._params(url, params),
                                        timeout=60, stream=gz_header_only)
            except requests.RequestException as exc:
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
                continue
            if resp.status_code == 404:
                path.with_suffix(".miss").touch()
                return None
            if resp.status_code in (429, 500, 502, 503, 504):
                resp.close()
                time.sleep(2 ** (attempt + 1))
                continue
            resp.raise_for_status()
            if gz_header_only:
                text = _read_matrix_header(resp)
            else:
                text = resp.text
            path.write_text(text, encoding="utf-8")
            return text
        raise RuntimeError(f"giving up after retries: {url}")


def _read_matrix_header(resp: requests.Response) -> str:
    lines = []
    with gzip.GzipFile(fileobj=resp.raw) as gz:
        for line in io.TextIOWrapper(gz, encoding="utf-8", errors="replace"):
            if line.startswith("!series_matrix_table_begin"):
                break
            lines.append(line)
    resp.close()
    return "".join(lines)


# --------------------------------------------------------------------------
# Stage 0/1: record and expand
# --------------------------------------------------------------------------

def build_query(args: argparse.Namespace, rules: Rules) -> dict:
    query = {
        "drug": args.drug,
        "disease": args.disease,
        "omics": args.omics,
        "objective": args.objective,
        "organism": args.organism,
        "sample_type": "patient_tumor",
        "defaults": {
            "min_patients_total": args.min_patients,
            "min_patients_per_arm": args.min_per_arm,
            "max_series": args.max_series,
            "ncbi_max_req_per_sec": DEFAULTS["ncbi_max_req_per_sec"],
        },
    }
    terms: dict[str, list[dict]] = {}
    if args.disease:
        terms["disease"] = rules.expand("disease", args.disease)
        for t in _split_terms(args.extra_terms):
            terms["disease"].append({"term": t, "source": "agent"})
    if args.drug:
        terms["drug"] = rules.expand("drug", args.drug)
        query["drug_class_terms"] = rules.drug_classes.get(_norm(args.drug), [])
        for t in _split_terms(args.extra_drug_terms):
            terms["drug"].append({"term": t, "source": "agent"})

    gds_types: list[str] = []
    if args.omics:
        gds_types = rules.omics.get(_norm(args.omics), [])
        if not gds_types:
            raise SystemExit(f"--omics '{args.omics}' not in rules/omics.tsv; "
                             f"known: {', '.join(sorted(rules.omics))}")
    return {
        "question": args.question,
        "query": query,
        "search_terms": terms,
        "gds_types": gds_types,
        "run_date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "skill_version": VERSION,
    }


def _split_terms(raw: Optional[str]) -> list[str]:
    return [t.strip() for t in (raw or "").split(";") if t.strip()]


def build_esearch_term(record: dict) -> str:
    def group(terms: list[dict]) -> str:
        return "(" + " OR ".join(f'"{t["term"]}"[All Fields]' for t in terms) + ")"

    parts = []
    for kind in ("disease", "drug"):
        if record["search_terms"].get(kind):
            parts.append(group(record["search_terms"][kind]))
    parts.append('"gse"[Entry Type]')
    parts.append(f'"{record["query"]["organism"]}"[Organism]')
    if record["gds_types"]:
        parts.append("(" + " OR ".join(f'"{t}"[DataSet Type]' for t in record["gds_types"]) + ")")
    return " AND ".join(parts)


def print_expansion(record: dict) -> None:
    print(f'question: "{record["question"]}"')
    for k in ("drug", "disease", "omics", "objective", "organism"):
        print(f"  {k + ':':<11}{record['query'][k] if record['query'][k] else 'unspecified'}")
    for kind, terms in record["search_terms"].items():
        print(f"{kind} terms:")
        for t in terms:
            print(f"  {t['term']:<45}{t['source']}")
    if len(record["search_terms"].get("disease", [])) == 1 and record["query"]["disease"]:
        print("  (no synonyms in rules/synonyms.tsv for this disease; searching the user term only)")


# --------------------------------------------------------------------------
# Stage 2/3: search and summarize
# --------------------------------------------------------------------------

def esearch(fetcher: Fetcher, term: str) -> tuple[int, list[str]]:
    text = fetcher.get_text(f"{EUTILS}/esearch.fcgi",
                            {"db": "gds", "term": term, "retmax": 10000, "retmode": "json"})
    res = json.loads(text)["esearchresult"]
    return int(res["count"]), res.get("idlist", [])


def esummary(fetcher: Fetcher, uids: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(uids), 100):
        batch = uids[i:i + 100]
        text = fetcher.get_text(f"{EUTILS}/esummary.fcgi",
                                {"db": "gds", "id": ",".join(batch), "retmode": "json"})
        result = json.loads(text)["result"]
        for uid in batch:  # keep esearch order
            doc = result.get(uid)
            if not doc:
                continue
            acc = doc.get("accession", "")
            if not acc.startswith("GSE"):
                continue
            out.append({
                "gse_accession": acc,
                "title": doc.get("title", ""),
                "summary": doc.get("summary", ""),
                "gdstype": doc.get("gdstype", ""),
                "n_samples": int(doc.get("n_samples") or 0),
                "pubmed_ids": [str(p) for p in doc.get("pubmedids", [])],
                "taxon": doc.get("taxon", ""),
                "platform": "GPL" + str(doc.get("gpl", "")).replace(";", ";GPL") if doc.get("gpl") else "",
                "is_superseries": SUPERSERIES_SENTENCE in doc.get("summary", "").lower(),
            })
    return out


# --------------------------------------------------------------------------
# Stage 4: characterize (series matrix header)
# --------------------------------------------------------------------------

def matrix_dir_url(acc: str) -> str:
    digits = acc[3:]
    stub = "GSE" + (digits[:-3] if len(digits) > 3 else "") + "nnn"
    return f"{GEO_FTP}/{stub}/{acc}/matrix/"


def fetch_series_metadata(fetcher: Fetcher, acc: str) -> dict:
    """Return {'series': {...}, 'samples': [...], 'metadata_source': str}."""
    listing = fetcher.get_text(matrix_dir_url(acc))
    files = sorted(set(re.findall(r'href="([^"]*_series_matrix\.txt\.gz)"', listing or "")))
    if files:
        series: dict = {}
        samples: list[dict] = []
        for fname in files:
            header = fetcher.get_text(matrix_dir_url(acc) + fname, gz_header_only=True)
            s, smp = parse_matrix_header(header or "")
            for k, v in s.items():
                series.setdefault(k, v)
            samples.extend(smp)
        return {"series": series, "samples": samples, "metadata_source": "series_matrix"}

    text = fetcher.get_text(GEO_ACC, {"acc": acc, "targ": "gsm", "form": "text", "view": "brief"})
    if text:
        return {"series": {}, "samples": parse_geo_text_samples(text), "metadata_source": "geo_text"}
    return {"series": {}, "samples": [], "metadata_source": "unavailable"}


def _split_characteristic(cell: str) -> tuple[str, str]:
    if ":" in cell:
        k, v = cell.split(":", 1)
        return k.strip(), v.strip()
    return "characteristic", cell.strip()


def _add_char(chars: dict, key: str, value: str) -> None:
    k, i = key, 2
    while k in chars:
        k, i = f"{key} ({i})", i + 1
    chars[k] = value


def parse_matrix_header(text: str) -> tuple[dict, list[dict]]:
    series: dict = {}
    rows: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if not line.startswith("!"):
            continue
        parts = line.rstrip("\n").split("\t")
        key = parts[0]
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key.startswith("!Series_"):
            series.setdefault(key, []).extend(v for v in vals if v)
        elif key.startswith("!Sample_"):
            rows.append((key, vals))

    accs = next((v for k, v in rows if k == "!Sample_geo_accession"), [])
    samples = [{"gsm_accession": a, "title": "", "source_name": "", "organism": "",
                "platform": "", "characteristics": {}} for a in accs]
    simple = {"!Sample_title": "title", "!Sample_source_name_ch1": "source_name",
              "!Sample_organism_ch1": "organism", "!Sample_platform_id": "platform"}
    for key, vals in rows:
        for i, v in enumerate(vals[:len(samples)]):
            if key in simple:
                samples[i][simple[key]] = v
            elif key == "!Sample_characteristics_ch1" and v:
                _add_char(samples[i]["characteristics"], *_split_characteristic(v))
    return series, samples


def parse_geo_text_samples(text: str) -> list[dict]:
    samples: list[dict] = []
    cur: Optional[dict] = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm_accession": line.split("=", 1)[1].strip(), "title": "", "source_name": "",
                   "organism": "", "platform": "", "characteristics": {}}
            samples.append(cur)
            continue
        if cur is None or " = " not in line:
            continue
        key, val = line.split(" = ", 1)
        if key == "!Sample_title":
            cur["title"] = val
        elif key == "!Sample_source_name_ch1":
            cur["source_name"] = val
        elif key == "!Sample_organism_ch1":
            cur["organism"] = val
        elif key == "!Sample_platform_id":
            cur["platform"] = val
        elif key == "!Sample_characteristics_ch1":
            _add_char(cur["characteristics"], *_split_characteristic(val))
    return samples


def fetch_subseries(fetcher: Fetcher, acc: str) -> list[str]:
    text = fetcher.get_text(GEO_ACC, {"acc": acc, "targ": "self", "form": "text", "view": "brief"})
    return re.findall(r"!Series_relation = SuperSeries of:\s*(GSE\d+)", text or "")


# --------------------------------------------------------------------------
# Classification (rules only)
# --------------------------------------------------------------------------

def _is_line_id(value: str, rules: Rules) -> bool:
    v = value.strip()
    if rules.cell_line_re and rules.cell_line_re.search(v):
        return True
    return bool(LINE_ID_VALUE.match(v)) and not GENERIC_CELL_VALUE.search(v)


def classify_sample_type(sample: dict, rules: Rules) -> str:
    chars = sample["characteristics"]
    text = " | ".join([sample.get("source_name", ""), sample.get("title", "")] + list(chars.values()))
    for cls, pattern in SAMPLE_TYPE_RULES:
        if cls == "cell_line":
            if any(CELL_LINE_KEY.match(k) and _is_line_id(v, rules) for k, v in chars.items()) or \
                    pattern.search(text) or (rules.cell_line_re and rules.cell_line_re.search(text)):
                return cls
        elif pattern.search(text):
            return cls
    return "NEEDS_REVIEW"


def is_response_field(key: str, rules: Rules) -> bool:
    k = _norm(key)
    if RESPONSE_KEY_EXCLUDE.search(k):
        return False
    return any(k == f or k.startswith(f + " ") or k.startswith(f + "(")
               or k.endswith(" " + f) for f in rules.response_fields)


def map_response(value: str, rules: Rules) -> str:
    """Whole-value match only; anything unmatched is NEEDS_REVIEW."""
    return rules.response_labels.get(_norm(value), "NEEDS_REVIEW")


# Patient identifiers encoded in sample titles. Submitters record the patient
# in the title far more often than in a characteristics field, and unlike a
# count quoted in an abstract this gives a per-sample mapping -- so repeated
# timepoints from one patient collapse instead of being counted twice.
PATIENT_TITLE_RULES = [
    re.compile(r"\bpatient (?:number |no\.? |#)?(\d{1,4})\b", re.I),
    re.compile(r"\b(?:subject|donor|case)[ _#-]?(\w{1,6})\b", re.I),
    re.compile(r"^(Pt[ _-]?\d{1,4})(?!\d)", re.I),
    re.compile(r"^(P\d{1,4})(?!\d)"),
]


def patient_id(sample: dict) -> Optional[str]:
    for k, v in sample["characteristics"].items():
        if PATIENT_ID_KEY.match(k.strip()) and v.strip():
            return v.strip()
    for pat in PATIENT_TITLE_RULES:
        m = pat.search(sample.get("title", "") or "")
        if m:
            return re.sub(r"[ _-]", "", m.group(1)).upper()
    return None


def annotate_sample(sample: dict, rules: Rules, drug_re: Optional[re.Pattern],
                    organism: str = DEFAULTS["organism"],
                    drug_class_re: Optional[re.Pattern] = None) -> dict:
    chars = sample["characteristics"]
    org = sample.get("organism", "")
    resp_fields = [(k, v) for k, v in chars.items() if is_response_field(k, rules)]
    resp_key, resp_raw = resp_fields[0] if resp_fields else ("", "")
    drug_field, drug_string = "", ""
    if drug_re:
        for k, v in list(chars.items()) + [("source_name", sample.get("source_name", ""))]:
            if drug_re.search(f"{k}: {v}"):
                drug_field, drug_string = k, v
                break
    class_field, class_string = "", ""
    if drug_class_re:
        for k, v in list(chars.items()) + [("source_name", sample.get("source_name", ""))]:
            if drug_class_re.search(f"{k}: {v}"):
                class_field, class_string = k, v
                break
    timepoints = [v for k, v in chars.items() if TIMEPOINT_KEY.search(k)]
    return {
        **sample,
        "sample_type": classify_sample_type(sample, rules),
        # Unknown organism is not a mismatch; only a stated different organism is.
        "organism_match": (not org) or _norm(organism) in _norm(org),
        "patient_id": patient_id(sample),
        "response_field": resp_key,
        "response_raw": resp_raw,
        "response_mapped": map_response(resp_raw, rules) if resp_fields else "",
        "drug_field": drug_field,
        "drug_string": drug_string,
        "drug_class_field": class_field,
        "drug_class_string": class_string,
        "has_survival": any(SURVIVAL_KEY.search(k) for k in chars),
        "timepoint": "; ".join(timepoints),
    }


def _value_counts(samples: list[dict], key_filter, max_levels: int = 10, max_len: int = 400) -> str:
    """'key: a=3, b=2; key2: ...' for characteristic keys passing key_filter.

    Keys whose values are all distinct (IDs, ages) or have more than
    max_levels levels are skipped, so the result reads as groups.
    """
    by_key: dict[str, Counter] = {}
    for s in samples:
        for k, v in s["characteristics"].items():
            if key_filter(k) and not PATIENT_ID_KEY.match(k.strip()) and not ID_LIKE_KEY.search(k):
                by_key.setdefault(k, Counter())[v] += 1
    parts = []
    for k, c in by_key.items():
        if 1 <= len(c) <= max_levels and (len(c) < len(samples) or len(samples) <= 2):
            parts.append(f"{k}: " + ", ".join(f"{v}={n}" for v, n in c.most_common()))
    out = "; ".join(parts)
    return out if len(out) <= max_len else out[:max_len - 3] + "..."


def _metadata_fields(samples: list[dict]) -> str:
    keys = {k for s in samples for k in s["characteristics"]
            if not PATIENT_ID_KEY.match(k.strip()) and not ID_LIKE_KEY.search(k)}
    return "; ".join(sorted(keys, key=str.lower))


RECIST = {"CR", "PR", "SD", "PD"}
STRICT_R = {"CR", "PR", "PR_OR_CR", "RESPONDER_UNSPECIFIED", "SENSITIVE"}
STRICT_NR = {"SD", "PD", "SD_OR_PD", "NONRESPONDER_UNSPECIFIED", "RESISTANT"}
EXCLUDED_TYPES = {"cell_line", "model", "perturbation"}


def summarize_series(acc: str, meta: dict, samples: list[dict], query: dict) -> dict:
    """Aggregate annotated samples to one series row and assign a verdict."""
    d = query["defaults"]
    objective = query["objective"]
    types = Counter(s["sample_type"] for s in samples)
    right_org = [s for s in samples if s.get("organism_match", True)]
    usable = [s for s in right_org if s["sample_type"] not in EXCLUDED_TYPES]

    # Patient count: from patient IDs on every usable sample; otherwise assume
    # one sample per patient only when nothing suggests repeated sampling
    # (no timepoint field, no patient IDs at all, not single-cell).
    ids = [s["patient_id"] for s in usable]
    series_text = " ".join([meta.get("title", ""), meta.get("summary", ""), meta.get("gdstype", "")] +
                           [v for s in usable for v in s["characteristics"].values()])
    n_patients: Optional[int] = None
    basis = "unknown"
    if usable and all(ids):
        from_title = any(not any(PATIENT_ID_KEY.match(k.strip()) for k in s["characteristics"])
                         for s in usable)
        n_patients = len(set(ids))
        basis = "patient_id_from_title" if from_title else "patient_id"
    elif usable and not any(ids) and not any(s["timepoint"] for s in usable) \
            and not SINGLE_CELL.search(series_text):
        n_patients, basis = len(usable), "assumed_one_per_sample"

    # Response labels, counted per patient where patient IDs exist.
    resp_samples = [s for s in usable if s["response_field"]]
    has_response = bool(resp_samples)
    unparsed = [s for s in resp_samples if s["response_mapped"] == "NEEDS_REVIEW"]
    per_unit: dict[str, set] = {}
    for s in resp_samples:
        if s["response_mapped"] in ("NOT_AVAILABLE", "NEEDS_REVIEW"):
            continue
        unit = s["patient_id"] or s["gsm_accession"]
        per_unit.setdefault(unit, set()).add(s["response_mapped"])
    conflicts = [u for u, labs in per_unit.items() if len(labs) > 1]
    labels = Counter(next(iter(l)) for u, l in per_unit.items() if len(l) == 1)
    n_r = sum(labels[k] for k in STRICT_R)
    n_nr = sum(labels[k] for k in STRICT_NR)
    # dcb puts durable SD with responders; duration is not parsed, so dcb is
    # only defined when no SD is present (then it equals strict).
    dcb_known = labels["SD"] == 0 and labels["SD_OR_PD"] == 0
    drug_hits = [s for s in samples if s["drug_field"]]
    class_hits = [s for s in samples if s.get("drug_class_field")]

    row = {
        "gse_accession": acc,
        "title": meta.get("title", ""),
        "objective": objective,
        "metadata_source": meta.get("metadata_source", ""),
        "n_samples": len(samples) or meta.get("n_samples"),
        "n_patients": n_patients if n_patients is not None else "unknown",
        "n_patients_basis": basis,
        "organism": ";".join(sorted({s["organism"] for s in samples if s["organism"]})),
        "platform": ";".join(sorted({s["platform"] for s in samples if s["platform"]})),
        "omics_type": meta.get("gdstype", ""),
        "drug_matched": bool(drug_hits) if query["drug"] else "",
        "drug_evidence_field": drug_hits[0]["drug_field"] if drug_hits else "",
        "drug_evidence_string": drug_hits[0]["drug_string"] if drug_hits else "",
        "drug_evidence_level": ("drug_named" if drug_hits else
                                "class_only" if class_hits else "none") if query["drug"] else "",
        "drug_class_evidence_string": class_hits[0]["drug_class_string"] if class_hits else "",
        "sample_type": ";".join(f"{k}:{v}" for k, v in types.most_common()),
        "heterogeneous_sample_types": len(types) > 1,
        "biopsy_timing": ";".join(sorted({s["timepoint"] for s in samples if s["timepoint"]})) or "unknown",
        "has_response_labels": has_response,
        "response_field": ";".join(sorted({s["response_field"] for s in resp_samples})),
        "response_counts_raw": ";".join(f"{k}:{v}" for k, v in sorted(labels.items())),
        "n_responder_strict": n_r if has_response else "",
        "n_nonresponder_strict": n_nr if has_response else "",
        "n_responder_dcb": (n_r if dcb_known else "unknown") if has_response else "",
        "n_nonresponder_dcb": (n_nr if dcb_known else "unknown") if has_response else "",
        "has_survival": any(s["has_survival"] for s in usable),
        "treatment": _value_counts(usable, lambda k: bool(TREATMENT_KEY.search(k))) or "not recorded",
        "sample_groups": _value_counts(usable, lambda k: not TREATMENT_KEY.search(k)),
        "metadata_fields": _metadata_fields(usable),
    }

    # Verdict: precedence order, first match wins; all failures listed.
    failures: list[tuple[str, str]] = []
    if samples and not right_org:
        failures.append(("WRONG_ORGANISM", f"no samples from {query.get('organism') or DEFAULTS['organism']}"))
    elif samples and not usable:
        failures.append(("WRONG_SAMPLE_TYPE", "all samples are cell_line/model/perturbation"))
    if len(right_org) < len(samples) and right_org:
        failures.append(("NEEDS_REVIEW", f"{len(samples) - len(right_org)} samples from another organism excluded"))
    if query["drug"] and not drug_hits and not class_hits:
        failures.append(("THERAPY_UNCONFIRMED",
                         "neither the drug nor its class named in any sample-level field"))
    if objective == "response" and not has_response:
        if row["has_survival"]:
            failures.append(("SURVIVAL_ONLY", "survival field present, no response field"))
        else:
            failures.append(("NO_RESPONSE_LABELS", "no response or survival field"))
    if not samples:
        failures.append(("NEEDS_REVIEW", "sample metadata unavailable"))
    if types.get("NEEDS_REVIEW"):
        failures.append(("NEEDS_REVIEW", f"{types['NEEDS_REVIEW']} samples with unclassified sample type"))
    if objective == "response" and unparsed:
        failures.append(("NEEDS_REVIEW", f"{len(unparsed)} unparsed response values"))
    if objective == "response" and conflicts:
        failures.append(("NEEDS_REVIEW", f"{len(conflicts)} patients with conflicting response labels"))
    if usable and n_patients is None:
        why = "some samples lack a patient ID" if any(ids) else \
            "no patient IDs and repeated sampling or single-cell design"
        failures.append(("NEEDS_REVIEW", f"n_patients unknown ({why})"))
    if n_patients is not None and n_patients < d["min_patients_total"]:
        failures.append(("INSUFFICIENT_N", f"{n_patients} patients < {d['min_patients_total']}"))
    if objective == "response" and has_response and min(n_r, n_nr) < d["min_patients_per_arm"]:
        failures.append(("INSUFFICIENT_N", f"strict arms {n_r}/{n_nr} < {d['min_patients_per_arm']}"))

    order = ["WRONG_ORGANISM", "WRONG_SAMPLE_TYPE", "THERAPY_UNCONFIRMED", "SURVIVAL_ONLY",
             "NO_RESPONSE_LABELS", "NEEDS_REVIEW", "INSUFFICIENT_N"]
    if failures:
        row["verdict"] = min((f[0] for f in failures), key=order.index)
        row["verdict_reason"] = "; ".join(f[1] for f in failures)
    else:
        class_only = bool(query["drug"]) and not drug_hits and bool(class_hits)
        row["verdict"] = "SUITABLE_CLASS_ONLY" if class_only else "SUITABLE"
        bits = [f"{n_patients} patients" +
                (" (assumed one sample per patient)" if basis == "assumed_one_per_sample" else "")]
        if query["drug"]:
            bits.insert(0, "therapy class confirmed, agent unestablished"
                        if class_only else "therapy confirmed")
        if has_response:
            bits.append(f"{n_r}/{n_nr} strict")
        row["verdict_reason"] = "; ".join(bits)
    return row


# --------------------------------------------------------------------------
# Stage 5: publications
# --------------------------------------------------------------------------

def fetch_pubmed(fetcher: Fetcher, pmid: str) -> dict:
    text = fetcher.get_text(f"{EUTILS}/efetch.fcgi", {"db": "pubmed", "id": pmid, "retmode": "xml"})
    rec = {"pubmed_id": pmid, "title": "", "journal": "", "year": "", "doi": "", "pmcid": "", "abstract": ""}
    if not text:
        return rec
    root = ET.fromstring(text)
    art = root.find(".//PubmedArticle")
    if art is None:
        return rec
    rec["title"] = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
    rec["journal"] = art.findtext(".//Journal/Title", "")
    rec["year"] = art.findtext(".//Journal/JournalIssue/PubDate/Year", "")
    rec["abstract"] = " ".join("".join(a.itertext()) for a in art.findall(".//Abstract/AbstractText"))
    for aid in art.findall(".//PubmedData/ArticleIdList/ArticleId"):
        if aid.get("IdType") == "doi":
            rec["doi"] = (aid.text or "").strip()
        elif aid.get("IdType") == "pmc":
            rec["pmcid"] = (aid.text or "").strip()
    return rec


def fetch_fulltext(fetcher: Fetcher, pmcid: str) -> Optional[str]:
    """Open-access full text from NCBI PMC (efetch db=pmc), or None.

    NCBI PMC is used rather than Europe PMC because www.ebi.ac.uk rejects the
    TLS handshake from FIPS-mode OpenSSL on some HPC systems (e.g. BU SCC).
    Articles outside the open-access subset come back without a <body>; those
    return None and the caller falls back to the abstract.
    """
    try:
        xml = fetcher.get_text(f"{EUTILS}/efetch.fcgi",
                               {"db": "pmc", "id": pmcid.removeprefix("PMC"), "retmode": "xml"})
    except (requests.RequestException, RuntimeError) as exc:
        # Full text is optional; keep the PubMed record and fall back to the abstract.
        print(f"  full text unavailable for {pmcid}: {type(exc).__name__}")
        return None
    if not xml:
        return None
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return None
    # Search the whole article (data-availability statements sit outside <body>),
    # but only when a body exists, i.e. the full text was actually returned.
    return " ".join(root.itertext()) if root.find(".//body") is not None else None


def check_publication(acc: str, pmid: str, fetcher: Fetcher, record: dict, verdict: str) -> dict:
    rec = fetch_pubmed(fetcher, pmid)
    full = fetch_fulltext(fetcher, rec["pmcid"]) if rec["pmcid"] else None
    if full:
        level, text = "full_text", full
    elif rec["abstract"]:
        level, text = "abstract_only", rec["title"] + " " + rec["abstract"]
    else:
        level, text = "none", rec["title"]

    def mentioned(kind: str) -> Optional[bool]:
        terms = [t["term"] for t in record["search_terms"].get(kind, [])]
        rx = _term_regex(terms)
        return bool(rx.search(text)) if rx else None

    resp_terms: list[str] = []
    if level != "none":
        resp_terms = sorted({m.group(0).lower() for m in _term_regex(PUB_RESPONSE_TERMS).finditer(text)})
        # Bare CR/PR/SD/PD only count as a set of at least three (PR alone is
        # often progesterone receptor).
        abbrevs = {m.group(1) for m in PUB_RECIST_ABBREV.finditer(text)}
        if len(abbrevs) >= 3:
            resp_terms += sorted(abbrevs)
    mentions_acc = bool(re.search(rf"\b{acc}\b", text)) if level == "full_text" else None
    drug_hit = mentioned("drug")

    flags = []
    if verdict == "THERAPY_UNCONFIRMED" and drug_hit:
        flags.append("THERAPY_IN_PAPER_ONLY")
    if verdict == "NO_RESPONSE_LABELS" and resp_terms:
        flags.append("RESPONSE_IN_PAPER_ONLY")
    if level == "full_text" and mentions_acc is False:
        flags.append("ACCESSION_NOT_IN_PAPER")

    return {
        "gse_accession": acc,
        "pubmed_id": pmid,
        "pmcid": rec["pmcid"],
        "doi": rec["doi"],
        "title": rec["title"],
        "journal": rec["journal"],
        "year": rec["year"],
        "pub_text_level": level,
        "pub_mentions_accession": "" if mentions_acc is None else mentions_acc,
        "pub_drug_mentioned": "" if drug_hit is None else drug_hit,
        "pub_disease_mentioned": "" if mentioned("disease") is None else mentioned("disease"),
        "pub_response_terms": ";".join(resp_terms),
        "pub_patient_count_sentences": " || ".join(
            m.group(0).strip() for m in PATIENT_COUNT_SENTENCE.finditer(text))[:2000] if level != "none" else "",
        "pub_flags": ";".join(flags),
    }


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

# Main output: one row per GEO dataset, no per-patient fields.
MAIN_COLS = [
    ("gse_accession", "gse_accession"), ("geo_link", "ftp_link"), ("title", "title"),
    ("paper_url", "pubmed_url"), ("doi", "doi"),
    ("n_samples", "n_samples"), ("n_patients", "n_patients"), ("n_patients_basis", "n_patients_basis"),
    ("organism", "organism"), ("platform", "platform"), ("data_type", "omics_type"),
    ("sample_types", "sample_type"), ("treatment", "treatment"), ("timepoints", "biopsy_timing"),
    ("sample_groups", "sample_groups"), ("therapy_evidence", "drug_evidence_level"),
    ("response_counts", "response_counts_raw"),
    ("n_responders_strict", "n_responder_strict"), ("n_nonresponders_strict", "n_nonresponder_strict"),
    ("n_responders_dcb", "n_responder_dcb"), ("n_nonresponders_dcb", "n_nonresponder_dcb"),
    ("survival_data", "has_survival"), ("metadata_fields", "metadata_fields"),
    ("verdict", "verdict"), ("verdict_reason", "verdict_reason"), ("publication_flags", "pub_flags"),
]
VERDICT_RANK = ["SUITABLE", "SUITABLE_CLASS_ONLY", "INSUFFICIENT_N", "NEEDS_REVIEW", "SURVIVAL_ONLY", "NO_RESPONSE_LABELS",
                "THERAPY_UNCONFIRMED", "WRONG_SAMPLE_TYPE", "WRONG_ORGANISM", "SUPERSERIES"]


def _sorted_cohorts(cohorts: list[dict]) -> list[dict]:
    rank = {v: i for i, v in enumerate(VERDICT_RANK)}
    return sorted(cohorts, key=lambda c: rank.get(c.get("verdict", ""), len(rank)))


def write_main_table(path: Path, cohorts: list[dict]) -> None:
    rows = []
    for c in _sorted_cohorts(cohorts):
        row = {name: c.get(src) for name, src in MAIN_COLS}
        row["paper_url"] = row["paper_url"] or "N/A"
        row["doi"] = row["doi"] or "N/A"
        rows.append(row)
    write_tsv(path, rows, [name for name, _ in MAIN_COLS])


VERDICT_HELP = {
    "SUITABLE": "Meets every criterion for the question asked.",
    "SUITABLE_CLASS_ONLY": "Meets every criterion, but sample metadata names only the "
                           "drug class (e.g. 'immunotherapy'), never the agent. The cohort "
                           "received some drug of that class; which one is unestablished.",
    "INSUFFICIENT_N": "Right kind of data, too few patients (or too few per response arm).",
    "NEEDS_REVIEW": "Metadata present but could not be parsed or counted; check by hand.",
    "SURVIVAL_ONLY": "Survival recorded but no response call.",
    "NO_RESPONSE_LABELS": "No response or survival annotation in the deposited metadata.",
    "THERAPY_UNCONFIRMED": "Drug named only in series text, not in sample metadata.",
    "WRONG_SAMPLE_TYPE": "Only cell lines, models or perturbation experiments.",
    "WRONG_ORGANISM": "No samples from the requested organism.",
    "SUPERSERIES": "Container series; its SubSeries are listed separately.",
}

HTML_CSS = """
:root{--bg:#f7f7f5;--card:#fff;--ink:#1d1d1b;--muted:#6b6b66;--line:#e3e2dd;--accent:#1f5f8b;
--ok:#1e7a46;--okbg:#e3f3e9;--warn:#8a5a00;--warnbg:#fbf0d9;--bad:#9b2c2c;--badbg:#f8e3e3;--nbg:#ecebe7}
@media (prefers-color-scheme:dark){:root{--bg:#161615;--card:#1f1f1d;--ink:#ecebe7;--muted:#a3a29c;
--line:#34332f;--accent:#7cb7e0;--ok:#7fd3a0;--okbg:#1d3326;--warn:#f0c46a;--warnbg:#3a2f16;
--bad:#f09a9a;--badbg:#3b1f1f;--nbg:#2a2a27}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1400px;margin:0 auto;padding:28px 20px 60px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:0 0 12px}
.sub{color:var(--muted);margin:0 0 24px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin:0 0 18px}
.q{font-size:17px;font-style:italic;margin:0 0 14px}
.kv{display:grid;grid-template-columns:max-content 1fr;gap:4px 18px;margin:0}
.kv dt{color:var(--muted)}.kv dd{margin:0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:18px}
.stats{display:flex;flex-wrap:wrap;gap:10px}
.stat{background:var(--nbg);border-radius:8px;padding:8px 14px;min-width:120px}
.stat b{display:block;font-size:20px}.stat span{color:var(--muted);font-size:12px}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12.5px}
pre{background:var(--nbg);border-radius:8px;padding:10px 12px;white-space:pre-wrap;word-break:break-word;margin:0}
table{border-collapse:collapse;width:100%}
th,td{text-align:left;vertical-align:top;padding:7px 10px;border-bottom:1px solid var(--line)}
th{font-size:12px;text-transform:uppercase;letter-spacing:.03em;color:var(--muted);font-weight:600}
.terms td:last-child{color:var(--muted)}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--card)}
.main th{position:sticky;top:0;background:var(--card);cursor:pointer;white-space:nowrap;z-index:1}
.main th:hover{color:var(--ink)}.main td{font-size:13px}
.main td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.main td.wide{min-width:220px;max-width:340px}.main td.title{min-width:240px;max-width:360px}
.main tr:hover td{background:var(--nbg)}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
.na{color:var(--muted)}.note{color:var(--muted);font-size:12px}
.badge{display:inline-block;padding:2px 8px;border-radius:999px;font-size:11.5px;font-weight:600;white-space:nowrap}
.b-ok{background:var(--okbg);color:var(--ok)}.b-warn{background:var(--warnbg);color:var(--warn)}
.b-bad{background:var(--badbg);color:var(--bad)}.b-n{background:var(--nbg);color:var(--muted)}
.toolbar{display:flex;gap:12px;align-items:center;margin:0 0 10px;flex-wrap:wrap}
.toolbar input{flex:1;min-width:200px;max-width:360px;padding:7px 10px;border:1px solid var(--line);
border-radius:8px;background:var(--card);color:var(--ink);font:inherit}
ul{margin:0;padding-left:20px}.legend td:first-child{white-space:nowrap}
footer{color:var(--muted);font-size:12px;margin-top:24px}
"""

HTML_JS = """
(function(){
  var t=document.getElementById('cohorts');if(!t)return;
  var tb=t.tBodies[0],f=document.getElementById('filter'),cnt=document.getElementById('count');
  function upd(){var q=f.value.toLowerCase(),n=0;
    for(var r of tb.rows){var s=r.textContent.toLowerCase().indexOf(q)>=0;r.hidden=!s;if(s)n++;}
    cnt.textContent=n+' of '+tb.rows.length+' datasets';}
  f.addEventListener('input',upd);upd();
  Array.prototype.forEach.call(t.tHead.rows[0].cells,function(th,i){
    th.addEventListener('click',function(){
      var asc=th.dataset.dir!=='asc';th.dataset.dir=asc?'asc':'desc';
      var rows=Array.prototype.slice.call(tb.rows);
      rows.sort(function(a,b){var x=a.cells[i].dataset.v||a.cells[i].textContent,
        y=b.cells[i].dataset.v||b.cells[i].textContent,nx=parseFloat(x),ny=parseFloat(y);
        var c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);return asc?c:-c;});
      rows.forEach(function(r){tb.appendChild(r);});});});
})();
"""


def _badge(verdict: str) -> str:
    cls = {"SUITABLE": "b-ok", "SUITABLE_CLASS_ONLY": "b-ok", "INSUFFICIENT_N": "b-warn", "NEEDS_REVIEW": "b-warn",
           "SURVIVAL_ONLY": "b-warn", "NO_RESPONSE_LABELS": "b-warn",
           "THERAPY_UNCONFIRMED": "b-bad", "WRONG_SAMPLE_TYPE": "b-bad",
           "WRONG_ORGANISM": "b-bad"}.get(verdict, "b-n")
    return f'<span class="badge {cls}" title="{html.escape(VERDICT_HELP.get(verdict, ""))}">' \
           f'{html.escape(verdict)}</span>'


def _paper_cell(c: dict) -> str:
    pmids = [p for p in str(c.get("pubmed_id") or "").split(";") if p]
    if not pmids:
        return '<span class="na">N/A</span>'
    links = [f'<a href="https://pubmed.ncbi.nlm.nih.gov/{html.escape(p)}/" target="_blank" '
             f'rel="noopener">PMID {html.escape(p)}</a>' for p in pmids]
    if c.get("doi"):
        links.append(f'<a class="note" href="https://doi.org/{html.escape(c["doi"])}" target="_blank" '
                     f'rel="noopener">DOI</a>')
    return "<br>".join(links)


def _cell(value, cls: str = "") -> str:
    v = "" if value is None else str(value)
    body = html.escape(v) if v else '<span class="na">—</span>'
    return f'<td class="{cls}">{body}</td>' if cls else f"<td>{body}</td>"


def write_html_report(out: Path, record: dict, search: dict, cohorts: list[dict], pubs: list[dict]) -> None:
    q, d = record["query"], record["query"]["defaults"]
    e = html.escape
    verdicts = Counter(c["verdict"] for c in cohorts)

    query_rows = "".join(
        f"<dt>{e(label)}</dt><dd>{e(str(q[key])) if q.get(key) else '<span class=na>unspecified</span>'}</dd>"
        for label, key in [("Disease", "disease"), ("Drug", "drug"), ("Data type", "omics"),
                           ("Objective", "objective"), ("Organism", "organism"), ("Sample type", "sample_type")])
    defaults = (f"<dt>Min. patients</dt><dd>{d['min_patients_total']}</dd>"
                f"<dt>Min. per arm</dt><dd>{d['min_patients_per_arm']}</dd>"
                f"<dt>Max. series</dt><dd>{d['max_series']}</dd>")

    term_blocks = []
    for kind, terms in record["search_terms"].items():
        rows = "".join(f"<tr><td>{e(t['term'])}</td><td>{e(t['source'])}</td></tr>" for t in terms)
        term_blocks.append(f'<div><h2>{e(kind.capitalize())} terms ({len(terms)})</h2>'
                           f'<table class="terms"><thead><tr><th>Term</th><th>Source</th></tr></thead>'
                           f'<tbody>{rows}</tbody></table></div>')
    if record.get("gds_types"):
        rows = "".join(f"<tr><td>{e(t)}</td><td>rules/omics.tsv</td></tr>" for t in record["gds_types"])
        term_blocks.append(f'<div><h2>Data type filter</h2><table class="terms"><thead><tr><th>GEO DataSet Type'
                           f'</th><th>Source</th></tr></thead><tbody>{rows}</tbody></table></div>')

    stats = "".join(f'<div class="stat"><b>{v}</b><span>{e(k)}</span></div>' for k, v in [
        ("series matched", search["n_matched"]),
        (f"passed prefilter (≥{d['min_patients_total']} samples)", search["n_prefiltered"]),
        ("characterized", search["n_characterized"]),
        ("not characterized (cap)", search["n_not_characterized"])])
    verdict_stats = "".join(f'<div class="stat"><b>{n}</b><span>{_badge(v)}</span></div>'
                            for v, n in sorted(verdicts.items(), key=lambda x: VERDICT_RANK.index(x[0])
                                               if x[0] in VERDICT_RANK else 99))

    head = ["GEO", "Title", "Paper", "Samples", "Patients", "Organism", "Data type", "Sample types",
            "Treatment", "Therapy evidence", "Timepoints", "Sample groups", "Response",
            "Survival", "Verdict", "Reason"]
    body_rows = []
    for c in _sorted_cohorts(cohorts):
        acc = e(c["gse_accession"])
        n_pat = c.get("n_patients")
        pat_note = ('<br><span class="note">assumed 1 per sample</span>'
                    if c.get("n_patients_basis") == "assumed_one_per_sample" else "")
        resp = c.get("response_counts_raw") or ""
        if c.get("n_responder_strict") not in ("", None):
            resp = (f"{resp}<br><span class=\"note\">strict R {c['n_responder_strict']} / "
                    f"NR {c['n_nonresponder_strict']} · durable-benefit "
                    f"{c.get('n_responder_dcb')} / {c.get('n_nonresponder_dcb')}</span>")
        lvl = c.get("drug_evidence_level") or ""
        lvl_html = {
            "drug_named": '<span class="badge b-ok">drug named</span>',
            "class_only": '<span class="badge b-warn" title="Sample metadata names only the '
                          'therapeutic class, never the agent.">class only</span>',
            "none": '<span class="badge b-bad">not in metadata</span>',
        }.get(lvl, '<span class="na">&mdash;</span>')
        body_rows.append(
            "<tr>"
            f'<td class="mono"><a href="{e(c.get("ftp_link", ""))}" target="_blank" rel="noopener">{acc}</a></td>'
            + _cell(c.get("title"), "title")
            + f"<td>{_paper_cell(c)}</td>"
            + f'<td class="num">{e(str(c.get("n_samples") or ""))}</td>'
            + f'<td class="num" data-v="{e(str(n_pat if isinstance(n_pat, int) else -1))}">'
              f'{e(str(n_pat))}{pat_note}</td>'
            + _cell(c.get("organism"))
            + _cell(c.get("omics_type"))
            + _cell(c.get("sample_type"))
            + _cell(c.get("treatment"), "wide")
            + f"<td>{lvl_html}</td>"
            + _cell(c.get("biopsy_timing") if c.get("biopsy_timing") != "unknown" else "")
            + _cell(c.get("sample_groups"), "wide")
            + (f"<td>{resp}</td>" if resp else _cell(""))
            + _cell("yes" if c.get("has_survival") is True else ("no" if c.get("has_survival") is False else ""))
            + f'<td data-v="{VERDICT_RANK.index(c["verdict"]) if c["verdict"] in VERDICT_RANK else 99}">'
              f'{_badge(c["verdict"])}</td>'
            + _cell(c.get("verdict_reason"), "wide")
            + "</tr>")

    flagged = [p for p in pubs if p.get("pub_flags")]
    flag_html = ("<table><thead><tr><th>GEO</th><th>PMID</th><th>Text checked</th><th>Flags</th></tr></thead><tbody>"
                 + "".join(f"<tr><td class=mono>{e(p['gse_accession'])}</td><td>{e(p['pubmed_id'])}</td>"
                           f"<td>{e(p['pub_text_level'])}</td><td>{e(p['pub_flags'])}</td></tr>" for p in flagged)
                 + "</tbody></table>") if flagged else '<p class="na">No publication flags.</p>'

    lim = []
    if search["n_not_characterized"]:
        lim.append(f"Only {search['n_characterized']} of {search['n_prefiltered']} prefiltered series were "
                   f"characterized (MAX_SERIES={d['max_series']}); series with a linked paper first.")
    assumed = sum(1 for c in cohorts if c.get("n_patients_basis") == "assumed_one_per_sample")
    if assumed:
        lim.append(f"Patient counts for {assumed} datasets assume one sample per patient (no patient IDs "
                   "deposited).")
    failed = sum(1 for c in cohorts if c.get("pub_status") == "check_failed")
    if failed:
        lim.append(f"Publication check failed for {failed} datasets; paper links are still shown.")
    lim.append("Verdicts reflect deposited GEO metadata only; the linked paper may contain more.")
    legend = "".join(f"<tr><td>{_badge(v)}</td><td>{e(h)}</td></tr>" for v, h in VERDICT_HELP.items())

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GEO Cohort Finder — {e(q.get('disease') or q.get('drug') or '')}</title><style>{HTML_CSS}</style></head>
<body><div class="wrap">
<h1>GEO Cohort Finder report</h1>
<p class="sub">Run {e(record['run_date'])} · geo-cohort-finder v{VERSION} · {fetcher_note(search)}</p>

<section class="card"><h2>Question</h2>
<p class="q">“{e(record['question'] or '')}”</p>
<div class="grid"><dl class="kv">{query_rows}</dl><dl class="kv">{defaults}</dl></div></section>

<section class="card"><h2>Search term expansion</h2>
<p class="note" style="margin-top:-6px">Every term searched, and where it came from. Terms are expanded before
the search runs; nothing outside this list is searched.</p>
<div class="grid">{''.join(term_blocks)}</div>
<h2 style="margin-top:18px">GEO query sent</h2><pre>{e(search['esearch_term'])}</pre></section>

<section class="card"><h2>Search results</h2><div class="stats">{stats}</div>
<h2 style="margin-top:16px">Verdicts</h2><div class="stats">{verdict_stats}</div></section>

<section><h2>Datasets</h2>
<div class="toolbar"><input id="filter" type="search" placeholder="Filter rows…"><span id="count" class="note"></span>
<span class="note">Click a column header to sort.</span></div>
<div class="tablewrap"><table class="main" id="cohorts"><thead><tr>{''.join(f'<th>{h}</th>' for h in head)}</tr></thead>
<tbody>{''.join(body_rows)}</tbody></table></div></section>

<div class="grid" style="margin-top:18px">
<section class="card"><h2>Publication flags</h2>{flag_html}</section>
<section class="card"><h2>Limitations</h2><ul>{''.join(f'<li>{e(x)}</li>' for x in lim)}</ul></section>
</div>
<section class="card"><h2>Verdict legend</h2><table class="legend"><tbody>{legend}</tbody></table></section>

<footer>Research use only. Not a medical device; not clinical decision support. Full data:
<code>geo_cohorts.tsv</code>, <code>result.json</code>, <code>tables/</code>; reproduce with
<code>reproducibility/commands.sh</code>.</footer>
</div><script>{HTML_JS}</script></body></html>
"""
    (out / "report.html").write_text(page, encoding="utf-8")


def fetcher_note(search: dict) -> str:
    return f"{search['n_characterized']} datasets characterized"


COHORT_COLS = [
    "gse_accession", "title", "objective", "metadata_source",
    "pubmed_id", "pubmed_url", "doi", "pmcid",
    "n_samples", "n_patients", "n_patients_basis", "organism", "platform", "omics_type", "processing_level",
    "drug_matched", "drug_evidence_level", "drug_evidence_field", "drug_evidence_string",
    "drug_class_evidence_string",
    "sample_type", "heterogeneous_sample_types", "biopsy_timing",
    "has_response_labels", "response_field", "response_counts_raw",
    "n_responder_strict", "n_nonresponder_strict", "n_responder_dcb", "n_nonresponder_dcb",
    "has_survival", "treatment", "sample_groups", "metadata_fields", "verdict", "verdict_reason",
    "pub_status", "pub_text_level", "pub_flags", "ftp_link",
]
SAMPLE_COLS = ["gse_accession", "gsm_accession", "patient_id", "organism", "source_name", "raw_characteristics",
               "response_field", "response_raw", "response_mapped", "timepoint", "sample_type"]
PUB_COLS = ["gse_accession", "pubmed_id", "pmcid", "doi", "title", "journal", "year",
            "pub_text_level", "pub_mentions_accession", "pub_drug_mentioned", "pub_disease_mentioned",
            "pub_response_terms", "pub_patient_count_sentences", "pub_flags"]
CANDIDATE_COLS = ["gse_accession", "title", "gdstype", "n_samples", "pubmed_ids", "platform",
                  "is_superseries", "prefilter", "characterized"]


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def write_report(out: Path, record: dict, search: dict, cohorts: list[dict], pubs: list[dict]) -> None:
    q, d = record["query"], record["query"]["defaults"]
    lines = ["# GEO Cohort Finder Report", "",
             f"Run: {record['run_date']} · skill v{VERSION}", "",
             "## Research Question",
             f'question: "{record["question"]}"',
             f"drug: {q['drug'] or 'unspecified'} · disease: {q['disease'] or 'unspecified'} · "
             f"omics: {q['omics'] or 'unspecified'} · objective: {q['objective']}",
             f"defaults: MIN_PATIENTS_TOTAL={d['min_patients_total']}, "
             f"MIN_PATIENTS_PER_ARM={d['min_patients_per_arm']}, MAX_SERIES={d['max_series']}", "",
             "## Search terms"]
    for kind, terms in record["search_terms"].items():
        lines.append(f"- {kind}: " + " · ".join(f"{t['term']} [{t['source']}]" for t in terms))
    lines += ["", "## Search",
              f"{search['n_matched']} series matched the query "
              f"({search['n_returned']} returned). {search['n_prefiltered']} passed the prefilter "
              f"(n_samples >= {d['min_patients_total']}).",
              f"{search['n_characterized']} characterized (series with a linked paper first, "
              f"then esearch order); "
              f"{search['n_not_characterized']} not characterized because of MAX_SERIES.", "",
              "## Verdicts"]
    for v, n in Counter(c["verdict"] for c in cohorts).most_common():
        lines.append(f"- {v}: {n}")

    lines += ["", "## All characterized series", "",
              "| accession | verdict | paper | n_pat | n_samp | sample types | reason |",
              "|---|---|---|---|---|---|---|"]
    for c in cohorts:
        paper = f"[PMID {c['pubmed_id']}]({c['pubmed_url']})" if c.get("pubmed_id") else "none"
        lines.append(f"| [{c['gse_accession']}]({c['ftp_link']}) | {c['verdict']} | {paper} | "
                     f"{c['n_patients']} | {c['n_samples']} | {c['sample_type']} | {c['verdict_reason']} |")

    flagged = [p for p in pubs if p["pub_flags"]]
    lines += ["", "## Publication flags", ""]
    if flagged:
        lines += ["| accession | PMID | pub_text_level | flags |", "|---|---|---|---|"]
        lines += [f"| {p['gse_accession']} | {p['pubmed_id']} | {p['pub_text_level']} | {p['pub_flags']} |"
                  for p in flagged]
    else:
        lines.append("None.")

    lim = [f"Only {search['n_characterized']} of {search['n_prefiltered']} prefiltered series "
           f"were characterized (MAX_SERIES={d['max_series']})."] if search["n_not_characterized"] else []
    unk = sum(1 for c in cohorts if c["n_patients"] == "unknown" and c["verdict"] != "SUPERSERIES")
    if unk:
        lim.append(f"n_patients could not be established for {unk} series; reported as unknown.")
    assumed = sum(1 for c in cohorts if c.get("n_patients_basis") == "assumed_one_per_sample")
    if assumed:
        lim.append(f"n_patients for {assumed} series assumes one sample per patient "
                   "(no patient IDs deposited); see n_patients_basis.")
    nr = sum(1 for c in cohorts if c["verdict"] == "NEEDS_REVIEW")
    if nr:
        lim.append(f"{nr} series need manual review; see verdict_reason and tables/samples.tsv.")
    lim.append("Verdicts reflect deposited GEO metadata only; publication flags mark where the paper "
               "may add information.")
    lines += ["", "## Limitations"] + [f"- {x}" for x in lim]
    lines += ["", "_Research use only. Not a medical device; not clinical decision support._", ""]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def write_reproducibility(out: Path, record: dict, esearch_term: str, argv: list[str]) -> None:
    rep = out / "reproducibility"
    rep.mkdir(parents=True, exist_ok=True)
    (rep / "query.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    (rep / "esearch_terms.txt").write_text(esearch_term + "\n", encoding="utf-8")
    cmd = " ".join(_shell_quote(a) for a in [sys.executable, str(Path(__file__).resolve())] + argv)
    (rep / "commands.sh").write_text(f"#!/usr/bin/env bash\n# run {record['run_date']}\n{cmd}\n",
                                     encoding="utf-8")
    (rep / "environment.yml").write_text(
        "name: geo-cohort-finder\n"
        f"python: \"{platform.python_version()}\"\n"
        f"platform: \"{platform.platform()}\"\n"
        "packages:\n"
        f"  requests: \"{requests.__version__}\"\n", encoding="utf-8")
    sums = []
    for p in sorted([out / "geo_cohorts.tsv", out / "report.html", out / "report.md", out / "result.json"] +
                    list((out / "tables").glob("*.tsv"))):
        if p.exists():
            sums.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}")
    (rep / "checksums.sha256").write_text("\n".join(sums) + "\n", encoding="utf-8")


def _shell_quote(s: str) -> str:
    return s if re.fullmatch(r"[\w./=:,-]+", s) else "'" + s.replace("'", "'\\''") + "'"


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

def run(args: argparse.Namespace, argv: list[str]) -> int:
    rules = Rules()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    # Stage 0/1: record + expand
    if args.resume:
        record = json.loads((out / "reproducibility" / "query.json").read_text())
    else:
        record = build_query(args, rules)
    term = build_esearch_term(record)
    print_expansion(record)
    print(f"\nesearch term:\n  {term}\n")
    (out / "reproducibility").mkdir(exist_ok=True)
    (out / "reproducibility" / "query.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    cache_dir = Path(args.cache_dir) if args.cache_dir else out / "cache"
    fetcher = Fetcher(cache_dir, email=args.email, offline=args.offline,
                      max_per_sec=DEFAULTS["ncbi_max_req_per_sec"])
    query = record["query"]
    d = query["defaults"]

    # Stage 2/3: search + summarize. Explicit accessions skip the search and
    # are looked up directly -- used by the gold-set harness and for targeted
    # re-runs, where the question is whether one known series scores correctly.
    if getattr(args, "accessions", None):
        accs = [a.strip() for a in args.accessions.split(",") if a.strip()]
        uids = []
        for a in accs:
            _, found = esearch(fetcher, f"{a}[Accession] AND gse[Entry Type]")
            uids.extend(found[:1])
        n_matched = len(uids)
        print(f"accessions: {len(accs)} supplied; search stage skipped")
        candidates = esummary(fetcher, uids) if uids else []
    else:
        n_matched, uids = esearch(fetcher, term)
        print(f"search: {n_matched} series matched; summarizing {len(uids)}")
        candidates = esummary(fetcher, uids)
    for c in candidates:
        c["prefilter"] = "pass" if c["n_samples"] >= d["min_patients_total"] else "n_samples_below_min"
        c["characterized"] = False
    passing = [c for c in candidates if c["prefilter"] == "pass"]
    by_acc = {c["gse_accession"]: c for c in candidates}

    # Stage 4: characterize, first MAX_SERIES only; SuperSeries expand to SubSeries.
    drug_re = _term_regex([t["term"] for t in record["search_terms"].get("drug", [])])
    drug_class_re = _term_regex(query.get("drug_class_terms", []))
    # Order: series with a linked PubMed ID first, then the rest; esearch order within each.
    queue = [c["gse_accession"] for c in passing if c["pubmed_ids"]] + \
            [c["gse_accession"] for c in passing if not c["pubmed_ids"]]
    seen: set[str] = set()
    cohorts: list[dict] = []
    all_samples: list[dict] = []
    n_char = 0
    while queue and n_char < d["max_series"]:
        acc = queue.pop(0)
        if acc in seen:
            continue
        seen.add(acc)
        meta = by_acc.get(acc, {"gse_accession": acc})
        if meta.get("is_superseries"):
            subs = fetch_subseries(fetcher, acc)
            cohorts.append({"gse_accession": acc, "title": meta.get("title", ""),
                            "objective": query["objective"], "verdict": "SUPERSERIES",
                            "verdict_reason": "SuperSeries; see SubSeries " + ",".join(subs),
                            "n_samples": meta.get("n_samples"), "n_patients": "unknown",
                            "sample_type": "", "pubmed_ids": meta.get("pubmed_ids", [])})
            queue = [s for s in subs if s not in seen] + queue
            continue
        print(f"characterize [{n_char + 1}/{d['max_series']}] {acc}")
        try:
            fetched = fetch_series_metadata(fetcher, acc)
        except (requests.RequestException, RuntimeError, OSError, EOFError) as exc:
            print(f"  metadata fetch failed: {exc}")
            fetched = {"series": {}, "samples": [], "metadata_source": f"failed: {type(exc).__name__}"}
        samples = [annotate_sample(s, rules, drug_re, query["organism"], drug_class_re)
                   for s in fetched["samples"]]
        series = fetched["series"]
        meta = {**meta, "metadata_source": fetched["metadata_source"]}
        meta.setdefault("title", (series.get("!Series_title") or [""])[0])
        meta.setdefault("gdstype", ";".join(series.get("!Series_type", [])))
        pmids = list(dict.fromkeys(meta.get("pubmed_ids", []) + series.get("!Series_pubmed_id", [])))
        row = summarize_series(acc, meta, samples, query)
        row["pubmed_ids"] = pmids
        cohorts.append(row)
        for s in samples:
            all_samples.append({**s, "gse_accession": acc, "raw_characteristics":
                                " | ".join(f"{k}: {v}" for k, v in s["characteristics"].items())})
        if acc in by_acc:
            by_acc[acc]["characterized"] = True
        n_char += 1

    # Stage 5: publications (annotate only; verdicts unchanged)
    pubs: list[dict] = []
    for c in cohorts:
        c["ftp_link"] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={c['gse_accession']}"
        pmids = c.pop("pubmed_ids", [])
        if not pmids:
            c.update({"pub_status": "no_linked_pmid", "pubmed_id": None, "pubmed_url": None,
                      "doi": None, "pmcid": None, "pub_text_level": "none", "pub_flags": ""})
            continue
        checks = []
        for p in pmids:
            try:
                checks.append(check_publication(c["gse_accession"], p, fetcher, record, c["verdict"]))
            except (requests.RequestException, RuntimeError, ET.ParseError) as exc:
                print(f"  publication check failed for PMID {p}: {exc}")
        if not checks:
            c.update({"pub_status": "check_failed", "pubmed_id": ";".join(pmids),
                      "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/",
                      "doi": None, "pmcid": None, "pub_text_level": "none", "pub_flags": ""})
            continue
        pubs.extend(checks)
        first = checks[0]
        c.update({"pub_status": "checked", "pubmed_id": ";".join(pmids),
                  "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/",
                  "doi": first["doi"] or None, "pmcid": first["pmcid"] or None,
                  "pub_text_level": first["pub_text_level"],
                  "pub_flags": ";".join(sorted({f for p in checks for f in p["pub_flags"].split(";") if f}))})

    # Stage 6: score is done per series above; write outputs.
    search = {"esearch_term": term, "n_matched": n_matched, "n_returned": len(uids),
              "n_prefiltered": len(passing), "n_characterized": n_char,
              "n_not_characterized": max(len(passing) - len(seen & {c["gse_accession"] for c in passing}), 0)}
    write_main_table(out / "geo_cohorts.tsv", cohorts)
    tables = out / "tables"
    write_tsv(tables / "candidates.tsv",
              [{**c, "pubmed_ids": ";".join(c["pubmed_ids"])} for c in candidates], CANDIDATE_COLS)
    write_tsv(tables / "cohorts.tsv", cohorts, COHORT_COLS)
    write_tsv(tables / "samples.tsv", all_samples, SAMPLE_COLS)
    write_tsv(tables / "publications.tsv", pubs, PUB_COLS)
    result = {"question": record["question"], "query": query, "search_terms": record["search_terms"],
              "search": search, "run_date": record["run_date"], "cohorts": cohorts,
              "publications": pubs, "network_requests": fetcher.n_requests}
    (out / "result.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    write_report(out, record, search, cohorts, pubs)
    write_html_report(out, record, search, cohorts, pubs)
    write_reproducibility(out, record, term, argv)

    print(f"\n{n_char} series characterized; {fetcher.n_requests} network requests.")
    for v, n in Counter(c["verdict"] for c in cohorts).most_common():
        print(f"  {v:<22}{n}")
    print(f"main table: {out / 'geo_cohorts.tsv'}")
    print(f"report: {out / 'report.html'} (also report.md)")
    return 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Find analyzable patient cohorts in NCBI GEO.")
    p.add_argument("--question", help="the user's natural-language question, recorded verbatim")
    p.add_argument("--accessions", help="comma-separated GSE accessions; skips the search stage")
    p.add_argument("--disease")
    p.add_argument("--drug")
    p.add_argument("--omics", help="see rules/omics.tsv")
    p.add_argument("--objective", choices=["response", "discovery"], default="discovery")
    p.add_argument("--organism", default=DEFAULTS["organism"])
    p.add_argument("--extra-terms", help="extra disease terms, ';'-separated (source: agent)")
    p.add_argument("--extra-drug-terms", help="extra drug terms, ';'-separated (source: agent)")
    p.add_argument("--min-patients", type=int, default=DEFAULTS["min_patients_total"])
    p.add_argument("--min-per-arm", type=int, default=DEFAULTS["min_patients_per_arm"])
    p.add_argument("--max-series", type=int, default=DEFAULTS["max_series"])
    p.add_argument("--email", default=os.environ.get("NCBI_EMAIL"),
                   help="contact email sent to NCBI (or set NCBI_EMAIL)")
    p.add_argument("--output", default="output")
    p.add_argument("--cache-dir", help="default: <output>/cache")
    p.add_argument("--offline", action="store_true", help="use cache only; fail on a cache miss")
    p.add_argument("--resume", action="store_true", help="reuse <output>/reproducibility/query.json and cache")
    p.add_argument("--demo", action="store_true", help="run the bundled HNSCC query from demo_data/cache, offline")
    args = p.parse_args(argv)

    if args.demo:
        for k, v in DEMO_ARGS.items():
            setattr(args, k, v)
        args.cache_dir = str(DEMO_CACHE)
        args.offline = True
        if args.output == "output":
            args.output = "output_demo"
    if not args.resume:
        if not args.question:
            p.error("--question is required (the user's question, verbatim)")
        if not (args.disease or args.drug):
            p.error("give at least one of --disease or --drug")
    return args


def main(argv: Optional[list[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    try:
        return run(args, argv)
    except CacheMiss as exc:
        print(f"error: {exc}\n(demo cache not populated yet; see demo_data/README.md)", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
