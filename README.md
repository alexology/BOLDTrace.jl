# BOLDTrace.jl

Identification of amplicon sequence variants (ASVs) through
[BOLD Systems v5](https://id.boldsystems.org), lowest common ancestor (LCA)
assignment, traceability of the BOLD reference records used and
documentation of the people who identified the reference specimens. It
includes a web interface built on [Genie.jl](https://genieframework.com).

This package is part of an experimental effort to use AI to port and develop
bioinformatics and quantitative ecology tools in Julia, with the goal of
creating a unified environment for high-performance scientific computing.
High performance is a goal, not a measured guarantee.

```
FASTA → BOLD identification engine → matches → LCA → assignment → usage counters → export
```

## Install

Not registered. BOLDTrace depends on
[BOLDassignment.jl](https://github.com/alexology/BOLDassignment.jl), which is
not registered either:

```julia
using Pkg
Pkg.add(url = "https://github.com/alexology/BOLDassignment.jl")
Pkg.develop(path = "path/to/BOLDTrace.jl")
```

With Julia ≥ 1.11, `[sources]` in `Project.toml` points to a sibling checkout
`../BOLDassignment.jl`.

## Web application

```
julia --threads=2 examples/run_app.jl
```

opens `http://127.0.0.1:8080` with seven sections:

- **Upload**: multi-FASTA, validation, number of sequences and of unique
  sequences, errors and warnings, file and sequence-set hashes;
- **Analysis**: BOLD library and search mode, every LCA parameter (identity,
  coverage, distance from the best match, number of matches, alignment
  length, minimum agreement, rank cap and its threshold for each rank),
  progress and errors;
- **Assignments**: ASV ID, taxon, rank, best identity, coverage, matches in
  the LCA, agreement and dissenting references; click a row to see every
  reference, included or excluded, with the reason. Downloads: assignments,
  all hits, and the hits of the selected ASV;
- **Reference usage**: process ID, taxon, identifier and the three usage
  counts (download);
- **Taxonomic attribution**: statistics per taxonomic identifier (download);
- **Statistics**: personal (project/user) and global counters;
- **Export**: the six files listed below.

`serve(; database, checkpoint, host, port, allowed_hosts, ...)` runs the
server; `serve(...; async = true)` returns it for `shutdown`. Analyses run
one at a time (the BOLD rate limits apply to the whole client); at most
`max_pending` (20) may wait, and queued or running analyses are never
dropped from memory. `shutdown(server; timeout = Inf)` cancels the queued
analyses, lets the running one finish and be recorded, then closes the
database. `serve(...; defaults = LCAParams(...), accept_rank = :genus)` sets
the parameters proposed by the page; each analysis can still change them.
The default host accepts only local connections and there is no authentication: on a
shared server use a reverse proxy with HTTPS and `allowed_hosts`. Request
bodies above `max_body` (64 MiB) are refused before being read; responses
carry a restrictive Content Security Policy.

## Use from Julia

```julia
using BOLDTrace

input = read_asvs("asvs.fasta")                    # validation, duplicates, hashes
a = run_analysis(input; project = "lake2026", user = "AL",
                 params = LCAParams(min_identity = 90, min_coverage = 80, max_delta = 1),
                 db = :animal_public, mode = :rapid, checkpoint = "bold_checkpoint")

assignment_table(a)          # one row per input record
match_table(a)               # every (ASV, reference) match with the LCA decision
reference_usage(a)           # retrieved / LCA / contribution counts per record
taxonomic_attribution(a)     # per identified_by value

store = open_store("boldtrace.sqlite")
record_analysis!(store, a)   # not recorded twice
usage_stats(store; project = "lake2026")   # personal
usage_stats(store)                         # global
export_analysis(a, "results/")
```

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
  information. The provenance report records the retrieval date and states
  these limitations.

Per-sequence coverage is the engine's `qcovhsp` (percent of the query covered
by the alignment); alignment length is its `length`.

## LCA (`BOLDTrace.LCA`, version 1.1.0)

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
   (`rank_thresholds`). The uncapped LCA rank is reported (`lca_rank`).

Flags: `discordant:<rank>`, `majority:<rank>` (assigned with dissenting
references), `incomplete_taxonomy:<rank>`,
`inconsistent_hierarchy:<rank>` (one name under different parents),
`identical_reference_sequences` and
`identical_reference_sequences_discordant:<rank>` (identical reference
sequences, by hash of the portal sequence, with different names),
`truncated_at_max_matches`. Status: `assigned`, `no_matches`,
`no_admissible_match`, `unresolved`.

The rank is the depth of agreement among the selected references, **not** a
statistical confidence. An assignment is *accepted* when it reaches
`accept_rank` (default `:phylum`, i.e. any assignment).

## Reference usage and attribution

For each reference record (`process_id` + marker) and the unique ASVs of an
analysis:

- **retrieved usage**: ASVs whose BOLD matches include the record;
- **LCA participation**: ASVs whose LCA set includes it;
- **assignment contribution**: ASVs with an accepted assignment whose LCA set
  includes it and which it supports (not dissenting).

They describe participation; they do not show that a record was
indispensable for an assignment.

Taxonomic attribution uses the verbatim BOLD `identified_by` field. Names
are not normalized or merged (similar names are not assumed to be the same
person; a value listing several people stays one entry); an empty `orcid`
column is reserved for a future normalization. No contribution is assigned
to the sequencing laboratory, the storing institution, the uploader or the
project manager. Columns: identifier, records, ASVs supported, LCA
participations, distinct taxa, and an **experimental** fractional credit
(each accepted ASV distributes 1 equally among the distinct identifiers of
the references supporting it). There are no rankings or merit indicators.

## Usage counters (SQLite)

`open_store(path)` creates a SQLite database with `analyses`,
`asv_results`, `bold_references`, `identification_events` and
`citation_links`. Events, not only counts, are stored: one event per
(unique ASV, reference record) of an analysis. An analysis with the same
fingerprint (project, user, sequence set, BOLD library and mode, LCA version
and parameters) is not recorded again, and event keys are unique.
`usage_stats` reports analyses, projects, events (retrieved, in LCA sets, in
accepted assignments), distinct ASVs, process IDs and reference sequences,
for a project/user/analysis (personal) or for the whole database (global).

Global statistics count only the identifications recorded by this
application, not all uses of BOLD worldwide.

The SQL is portable to PostgreSQL (TEXT/INTEGER/DOUBLE PRECISION,
`ON CONFLICT` upserts, ISO timestamps); a migration needs another
DBInterface driver and `$n` placeholders.

## Bibliographic links

`CitationLink(process_id, doi, relation, evidence_source, verification_status)`
with relations `original_barcode_publication`, `taxonomic_revision`,
`barcode_reuse_publication`, stored by `add_citation!`. Only the data
structure and the export of the process IDs used are implemented; there is
no literature search. A publication linked to a BOLD project does not
necessarily document every barcode of the project.

## Exports

| File | Content |
|---|---|
| `assignments.csv` | one row per input record (original ID), assignment, lineage, flags, selected process IDs |
| `bold_matches.csv` | every (ASV, reference) match: identity, coverage, length, reference taxonomy, `identified_by`, reference sequence hash, LCA decision |
| `reference_usage.csv` | the three counts per reference record |
| `taxonomic_attribution.csv` | statistics per identifier |
| `analysis_provenance.json` | analysis ID, date (UTC), input hashes, BOLD source and retrieval date, parameters, LCA and software versions, selected process IDs, assignments, limitations |
| `used_process_ids.tsv` | process IDs in LCA sets, with the citation columns |

Cells starting with `= + - @` get a leading apostrophe (spreadsheet formula
injection).

## Input

`read_asvs` reads FASTA with FASTX. Identifiers (first word of the header)
are kept. Sequences are normalized (no gaps or whitespace, upper case,
`U → T`) and hashed with SHA-256; identical sequences are analysed once and
mapped back to every input ID. Records with invalid characters, empty
sequences or an identifier reused for a different sequence are reported and
skipped without stopping the analysis; short sequences (< 100 bases, the
minimum overlap of the engine) and ambiguous ones are warnings. A file that
is not valid FASTA is rejected as a whole.

## Tests

```
julia --project=. --threads=4 test/runtests.jl
```

Simulated identification engine and portal (formats observed on BOLD),
the real BOLD response of BOLDassignment.jl, SQLite and the HTTP server; no
network access.

## License

MIT (see `LICENSE`). BOLD data retrieved with this package are subject to the
terms of use of BOLD Systems. The dnambientale.it logo in `src/assets` is the
one used by ASVAtlas.jl.
