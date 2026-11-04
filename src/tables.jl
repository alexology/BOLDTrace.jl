# Result tables of one analysis. Every table is a vector of NamedTuples with
# the same fields, ready for CSV/JSON export.

const NOT_PROVIDED = "(not provided)"

# Columns of every table, also when it has no rows (stable export schemas).
assignment_columns(ranks) = (:asv_id, :sequence_hash, :status, :taxon, :rank, Symbol.(ranks)...,
                             :agreement, :lca_rank, :rank_capped, :best_identity, :best_coverage,
                             :n_matches, :n_admissible, :n_selected, :n_dissenting, :accepted,
                             :flags, :selected_process_ids, :error)
match_columns(ranks) = (:asv_ids, :sequence_hash, :process_id, :record_id, :marker, :bin_uri,
                        :reference_sequence_hash, :reference_taxon, :reference_rank,
                        (Symbol(:ref_, r) for r in ranks)..., :identified_by, :identity, :coverage,
                        :align_length, :included_in_lca, :decision, :dissent_rank,
                        :supports_assignment, :final_taxon, :final_rank, :accepted)
const REFERENCE_USAGE_COLUMNS = (:process_id, :record_id, :marker, :reference_sequence_hash, :taxon,
                                 :rank, :identified_by, :bin_uri, :recordsets, :retrieved_usage,
                                 :lca_participation, :assignment_contribution)
const ATTRIBUTION_COLUMNS = (:identifier, :orcid, :n_records, :n_asvs_supported,
                             :n_lca_participations, :n_taxa, :fractional_credit)

_join(v) = join(v, ";")
_lineage_value(res::LCAResult, r::Symbol) =
    (i = findfirst(p -> first(p) === r, res.lineage); i === nothing ? missing : res.lineage[i].second)

"""
    assignment_table(a::Analysis) -> Vector{NamedTuple}

One row per input record (original identifiers, input order), including
rejected records (`status = "input_error"`). Columns: `asv_id`,
`sequence_hash`, `status`, `taxon`, `rank`, one column per rank,
`lca_rank` (before the rank cap), `rank_capped`, `best_identity`,
`best_coverage`, `n_matches`, `n_admissible`, `n_selected`, `accepted`,
`flags`, `selected_process_ids`, `error`.
"""
function assignment_table(a::Analysis)
    rows = NamedTuple[]
    for rec in a.input.records
        if rec.error !== nothing
            push!(rows, _assignment_row(a, rec.id, missing, nothing, nothing, rec.error))
            continue
        end
        i = rec.unique_index
        push!(rows, _assignment_row(a, rec.id, rec.hash, i, a.results[i], missing))
    end
    return rows
end

function _assignment_row(a, id, hash, i, res, err)
    ranks = a.params.ranks
    lineage = NamedTuple{Tuple(ranks)}(Tuple(res === nothing ? missing : _lineage_value(res, r)
                                             for r in ranks))
    if res === nothing
        return merge((asv_id = id, sequence_hash = hash, status = "input_error", taxon = missing,
                      rank = missing), lineage,
                     (agreement = missing, lca_rank = missing, rank_capped = missing,
                      best_identity = missing, best_coverage = missing, n_matches = 0,
                      n_admissible = 0, n_selected = 0, n_dissenting = 0, accepted = false,
                      flags = "", selected_process_ids = "", error = err))
    end
    pids = selected_process_ids(a, i)
    return merge((asv_id = id, sequence_hash = hash, status = String(res.status),
                  taxon = something(res.taxon, missing),
                  rank = res.rank === nothing ? missing : String(res.rank)), lineage,
                 (agreement = isempty(res.agreement) ? missing : round(last(res.agreement); digits = 2),
                  lca_rank = res.lca_rank === nothing ? missing : String(res.lca_rank),
                  rank_capped = res.rank_capped, best_identity = res.best_identity,
                  best_coverage = res.best_coverage, n_matches = res.n_matches,
                  n_admissible = res.n_admissible, n_selected = res.n_selected,
                  n_dissenting = count(m -> m.included && !supports(a, i, m), a.matches[i]),
                  accepted = accepted(a, i), flags = _join(res.flags),
                  selected_process_ids = _join(pids), error = err))
end

"""
    selected_process_ids(a::Analysis, i) -> Vector{String}

Process IDs of the LCA set of unique sequence `i`, by decreasing identity
then process ID (independent of the order of the BOLD results).
"""
function selected_process_ids(a::Analysis, i::Integer)
    sel = sort([m for m in a.matches[i] if m.included];
               by = m -> (-m.identity, a.references[m.key].process_id, m.key))
    return unique!([a.references[m.key].process_id for m in sel])
end

"""
    match_table(a::Analysis) -> Vector{NamedTuple}

One row per (unique sequence, reference record): `asv_ids` (all input IDs
with that sequence), `sequence_hash`, `process_id`, `record_id`, `marker`,
`bin_uri`, `reference_sequence_hash`, `reference_taxon`,
`reference_rank`, one column per rank of the reference lineage,
`identified_by`, `identity`, `coverage`, `align_length`, `included_in_lca`,
`decision`, `final_taxon`, `final_rank`, `accepted`.
"""
function match_table(a::Analysis)
    rows = NamedTuple[]
    ranks = a.params.ranks
    for (i, u) in enumerate(a.input.unique), m in a.matches[i]
        ref = a.references[m.key]
        res = a.results[i]
        taxon, rank = reference_taxon(ref, ranks)
        lineage = NamedTuple{Tuple(Symbol(:ref_, r) for r in ranks)}(
            Tuple(get(ref.lineage, r, missing) for r in ranks))
        push!(rows, merge((asv_ids = _join(u.ids), sequence_hash = u.hash,
                           process_id = ref.process_id, record_id = ref.record_id,
                           marker = ref.marker, bin_uri = ref.bin_uri,
                           reference_sequence_hash = ref.sequence_hash,
                           reference_taxon = taxon, reference_rank = rank), lineage,
                          (identified_by = ref.identified_by, identity = m.identity,
                           coverage = m.coverage, align_length = m.align_length,
                           included_in_lca = m.included, decision = String(m.reason),
                           dissent_rank = m.dissent === nothing ? missing : String(m.dissent),
                           supports_assignment = supports(a, i, m),
                           final_taxon = something(res.taxon, missing),
                           final_rank = res.rank === nothing ? missing : String(res.rank),
                           accepted = accepted(a, i))))
    end
    return rows
end

"""
    reference_usage(a::Analysis) -> Vector{NamedTuple}

One row per reference record with three descriptive counts over the unique
sequences of the analysis:

- `retrieved_usage`: sequences whose BOLD matches include the record;
- `lca_participation`: sequences whose LCA set includes the record;
- `assignment_contribution`: sequences with an accepted assignment whose LCA
  set includes the record.

These counts describe participation; they do not show that a record was
necessary for an assignment. Sorted by decreasing participation.
"""
function reference_usage(a::Analysis)
    retrieved = Dict{String,Int}()
    lca_part = Dict{String,Int}()
    contrib = Dict{String,Int}()
    for i in eachindex(a.results), m in a.matches[i]
        retrieved[m.key] = get(retrieved, m.key, 0) + 1
        if m.included
            lca_part[m.key] = get(lca_part, m.key, 0) + 1
            accepted(a, i) && supports(a, i, m) && (contrib[m.key] = get(contrib, m.key, 0) + 1)
        end
    end
    rows = NamedTuple[]
    for (key, n) in retrieved
        ref = a.references[key]
        taxon, rank = reference_taxon(ref, a.params.ranks)
        push!(rows, (process_id = ref.process_id, record_id = ref.record_id, marker = ref.marker,
                     reference_sequence_hash = ref.sequence_hash, taxon = taxon, rank = rank,
                     identified_by = coalesce(ref.identified_by, NOT_PROVIDED),
                     bin_uri = ref.bin_uri, recordsets = ref.recordsets,
                     retrieved_usage = n, lca_participation = get(lca_part, key, 0),
                     assignment_contribution = get(contrib, key, 0)))
    end
    sort!(rows; by = r -> (-r.assignment_contribution, -r.lca_participation, -r.retrieved_usage,
                           r.process_id, coalesce(r.marker, "")))
    return rows
end

"""
    taxonomic_attribution(a::Analysis) -> Vector{NamedTuple}

Aggregate statistics per taxonomic identifier, i.e. the verbatim BOLD
`identified_by` value of the reference records. Names are **not**
normalized or merged: similar spellings stay separate, a value listing
several people stays one entry. Records without `identified_by` are grouped
under `"(not provided)"`. Contribution is never assigned to the sequencing
laboratory, the storing institution, the uploader or the project manager.

Columns: `identifier`, `orcid` (empty, reserved for a future normalization),
`n_records` (distinct records among the matches), `n_asvs_supported`
(unique sequences with an accepted assignment whose LCA set includes one of
the identifier's records), `n_lca_participations` (record × sequence pairs
in LCA sets), `n_taxa` (distinct taxa of the identifier's records) and
`fractional_credit`, an **experimental participation indicator**: each
accepted sequence distributes 1 equally among the distinct identifiers of
its LCA set. It is not a measure of scientific merit.
"""
function taxonomic_attribution(a::Analysis)
    records = Dict{String,Set{String}}()
    asvs = Dict{String,Set{Int}}()
    parts = Dict{String,Int}()
    taxa = Dict{String,Set{String}}()
    credit = Dict{String,Float64}()
    ident(ref) = coalesce(ref.identified_by, NOT_PROVIDED)
    for i in eachindex(a.results)
        ok = accepted(a, i)
        in_set = Set{String}()
        for m in a.matches[i]
            ref = a.references[m.key]
            who = ident(ref)
            push!(get!(records, who, Set{String}()), m.key)
            t, _ = reference_taxon(ref, a.params.ranks)
            ismissing(t) || push!(get!(taxa, who, Set{String}()), t)
            get!(asvs, who, Set{Int}())
            m.included || continue
            parts[who] = get(parts, who, 0) + 1
            supports(a, i, m) || continue      # dissenting references give no support
            push!(in_set, who)
            ok && push!(asvs[who], i)
        end
        if ok && !isempty(in_set)
            for who in in_set
                credit[who] = get(credit, who, 0.0) + 1 / length(in_set)
            end
        end
    end
    rows = NamedTuple[]
    for (who, recs) in records
        push!(rows, (identifier = who, orcid = missing, n_records = length(recs),
                     n_asvs_supported = length(asvs[who]),
                     n_lca_participations = get(parts, who, 0),
                     n_taxa = length(get(taxa, who, Set{String}())),
                     fractional_credit = round(get(credit, who, 0.0); digits = 6)))
    end
    sort!(rows; by = r -> (r.identifier == NOT_PROVIDED, r.identifier))
    return rows
end

"""
    analysis_summary(a::Analysis) -> NamedTuple

Counts of the analysis: unique sequences by status, accepted assignments,
assignments by rank, distinct references and process IDs.
"""
function analysis_summary(a::Analysis)
    status = Dict{String,Int}()
    by_rank = Dict{String,Int}()
    for r in a.results
        status[String(r.status)] = get(status, String(r.status), 0) + 1
        r.rank === nothing || (by_rank[String(r.rank)] = get(by_rank, String(r.rank), 0) + 1)
    end
    used = Set(m.key for ms in a.matches for m in ms if m.included)
    return (analysis_id = a.id, n_records = length(a.input.records),
            n_unique = length(a.input.unique),
            n_input_errors = count(r -> r.error !== nothing, a.input.records),
            status = status, by_rank = by_rank,
            n_accepted = count(i -> accepted(a, i), eachindex(a.results)),
            n_references = length(a.references), n_references_in_lca = length(used),
            n_process_ids = length(unique(r.process_id for r in values(a.references))),
            errors = copy(a.errors))
end
