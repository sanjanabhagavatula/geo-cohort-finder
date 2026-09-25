#!/usr/bin/env bash
# run 2026-09-25T17:55:21+00:00 on main @ 69a6e40, Python 3.13.8
# Absolute cluster paths replaced with repo-relative ones; run from the repo root.
python3 skills/geo-cohort-finder/geo_cohort_finder.py --question 'Which HNSCC cohorts have pembrolizumab responders and non-responders?' --disease HNSCC --drug pembrolizumab --objective response --output output/hnscc_pembrolizumab_response
