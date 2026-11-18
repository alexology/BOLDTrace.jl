# BOLD data and the LCA

## BOLD data sources

Checked on 2026-11-04 against the published OpenAPI documents:

- **Identification**: the official Barcode ID API
  (`https://id.boldsystems.org/openapi.json`): `POST /submission` with a
  FASTA file, `GET /submission/status/{sub_id}`, `GET
  /submission/results/{sub_id}`. Limits: 3 submissions/minute and 40
  requests/minute; BOLDassignment.jl spaces the requests, retries and keeps a
  checkpoint per sequence. No undocumented endpoint is used. Modes: rapid
  (≥ 94 %, 25 hits, ≤ 1000 sequences per submission), genus and species
  (≥ 90 %, 50 hits), exhaustive (≥ 75 %, 100 hits).
- **Reference metadata**: the portal API
  (`https://portal.boldsystems.org/openapi.json`, `/api/query` and
  `/api/documents/{query_id}`, BCDM records): reference sequence (`nuc`),
  `identified_by`, `identification_method`, projects/datasets
  (`bold_recordset_code_arr`), voucher type, museum ID, storing institution,
  upload date. The portal returns **one record per specimen and marker**;
  BOLDTrace uses the record of the marker of the match (e.g. `COI-5P`).
- **Not available** from the APIs: a version or snapshot of the reference
  libraries, the alignment software and its version, bibliographic
  information. The provenance report records the retrieval dates and states
  these limitations: the date of the run, for every sequence whether its hits
  were queried now or reused from the checkpoint (with the date of the
  checkpoint file, i.e. of the original download), and the date of the
  portal metadata. Reused hits are also listed among the limitations.

Per-sequence coverage is the engine's `qcovhsp` (percent of the query covered
by the alignment); alignment length is its `length`.

## LCA (`BOLDTrace.LCA`, version 1.1.1)

An independent module (`lca(matches, LCAParams(...))`, no BOLD dependency).
For each ASV:

1. matches below `min_identity` (default 90 %), `min_coverage` (80 %) or
   `min_align_length` (off) are not admissible;
2. the best admissible match is the one with the highest identity;
3. admissible matches within `max_delta` (1 percentage point) of the best
   enter the LCA set, closest first, at most `max_matches` (100); ties are
   broken by identifier, and truncation within ties is flagged;
4. ranks are visited from phylum to species; at each rank the most frequent
   name among the selected references identified at that rank is assigned if
   at least `min_agreement` percent of them carry it. `min_agreement = 100`
   (default) is the strict LCA. A lower value, e.g. 95, tolerates a few
   discordant references, typically misidentified records: with 99 records
   of one species and one of another, the species is assigned with 99 %
   agreement. Discordant references are reported as *dissenting* (rank in
   `dissent_rank`), leave the deeper ranks, and give no support or credit to
   the assignment. Ties never resolve; `min_agreement` must be above 50;
5. references not identified at a rank (missing name, open nomenclature such
   as `Genus sp.`, `cf.`, `aff.`, interim codes) do not vote at that rank or
   deeper but still support the ranks above; if no selected reference is
   identified at a rank, the assignment stops there;
6. `rank_cap = true` (default) limits the rank by the identity of the best
   match: species ≥ 97, genus ≥ 95, family ≥ 90, order ≥ 85, class ≥ 75,
   phylum ≥ 50 % (as BOLDigger3/BOLDassignment), each threshold configurable
   (`rank_thresholds`). The lineage stops above the first rank whose
   threshold is not reached, so it never skips a rank. The uncapped LCA rank
   is reported (`lca_rank`).

Flags: `discordant:<rank>`, `majority:<rank>` (assigned with dissenting
references), `incomplete_taxonomy:<rank>`,
`inconsistent_hierarchy:<rank>` (one name under different parents),
`identical_reference_sequences` and
`identical_reference_sequences_discordant:<rank>` (identical reference
sequences, by hash of the portal sequence, with different names),
`truncated_at_max_matches`, `reference_classification_conflict` (a
reference record classified differently by hits of the same analysis, e.g. a
hit reused from an older checkpoint and a new one: each hit keeps its own
classification in the LCA, in `bold_matches.csv` and in the stored events;
`reference_usage.csv` lists all the classifications, and the analysis
reports a warning). Status: `assigned`, `no_matches`,
`no_admissible_match`, `unresolved`, and `result_missing` (flag
`bold_result_missing`) when a finished BOLD submission did not return the
result of a sequence: its matches are unknown, which is not the same as no
match. The analysis lists it among its errors and in the provenance
(`hits_retrieval.by_sequence`); the sequence is not checkpointed, so running
the analysis again queries it. An analysis with missing results is **not
recorded** in the usage statistics (`record_analysis!` returns `reason =
"incomplete"`), so the complete run is recorded when the missing sequences
are retrieved.

The rank is the depth of agreement among the selected references, **not** a
statistical confidence. An assignment is *accepted* when it reaches
`accept_rank` (default `:phylum`, i.e. any assignment).
