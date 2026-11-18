# Attribution, counters and citations

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
Each event keeps a snapshot of the reference as used by its analysis
(sequence hash, taxon, rank, `identified_by`): `bold_references` holds the
latest metadata, but later analyses never change the counts of earlier ones.
Databases of schema version 1 are migrated when opened; their old events
take the snapshot from the reference catalogue at that moment.

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
