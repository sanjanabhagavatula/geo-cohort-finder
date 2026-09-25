# Demo output: HNSCC pembrolizumab response

Example output of a live run on `main` @ `69a6e40` (2026-09-25), kept so the
report format can be inspected without running anything. Start with
`report.html` (open in a browser) or `geo_cohorts.tsv` (the main table).

```text
question:  "Which HNSCC cohorts have pembrolizumab responders and non-responders?"
flags:     --disease HNSCC --drug pembrolizumab --objective response
defaults:  MIN_PATIENTS_TOTAL=20, MIN_PATIENTS_PER_ARM=10, MAX_SERIES=10
result:    45 series matched; 10 characterized; 0 SUITABLE
```

The raw network cache is not included. Re-run with
`reproducibility/commands.sh`; results will differ as GEO changes.

## Known issues visible in this output

- **Drug expansion pulls in the whole ICI class.** "pembrolizumab" also
  searched nivolumab, atezolizumab, durvalumab, ipilimumab and anti-CTLA-4,
  because pembrolizumab appears in the `immune checkpoint inhibitor` row of
  `rules/synonyms.tsv` and `expand()` merges every row containing the term.
  GSE281729 (nivolumab) was therefore searched for a pembrolizumab question.
- **`bestresponse` (no space) is not recognised as a response field.**
  GSE259280 has `bestresponse: PD=202, PR=76, SD=67, CR=30` but is reported
  `SURVIVAL_ONLY`.
- **The 10-series cap excludes known cohorts.** GSE200996 matched the search
  but was not characterized. GSE159067 does not match the search terms.
- **GSE288199** has no treatment or response fields in its sample
  characteristics; response may be encoded in sample titles.
