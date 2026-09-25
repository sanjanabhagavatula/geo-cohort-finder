# Demo cache

`--demo` runs the bundled query offline from `demo_data/cache/`:

```text
question:  "List GEO datasets on HNSCC with transcriptomic data."
disease:   HNSCC   omics: transcriptomic   objective: discovery
```

The cache is populated once from a live run and committed:

```bash
python geo_cohort_finder.py \
    --question "List GEO datasets on HNSCC with transcriptomic data." \
    --disease HNSCC --omics transcriptomic --objective discovery \
    --cache-dir demo_data/cache --output /tmp/geo_demo_seed
```

Re-seed whenever the demo query, `rules/synonyms.tsv` or `rules/omics.tsv`
changes, since those change the search string and therefore the cache keys.
