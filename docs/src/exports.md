# Export files

## Exports

| File | Content |
|---|---|
| `assignments.csv` | one row per input record (original ID), assignment, lineage, flags, selected process IDs |
| `bold_matches.csv` | every (ASV, reference) match: identity, coverage, length, reference taxonomy, `identified_by`, reference sequence hash, LCA decision |
| `reference_usage.csv` | the three counts per reference record |
| `taxonomic_attribution.csv` | statistics per identifier |
| `analysis_provenance.json` | analysis ID, date (UTC), input hashes, BOLD source and retrieval dates (per sequence: queried or reused from the checkpoint), parameters, LCA and software versions, selected process IDs, assignments, limitations |
| `used_process_ids.tsv` | process IDs in LCA sets, with the citation columns |

Cells starting with `= + - @` get a leading apostrophe (spreadsheet formula
injection).
