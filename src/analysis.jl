# FASTA → BOLD identification → reference metadata → LCA.

const PORTAL_URL = "https://portal.boldsystems.org"

"""
    ReferenceRecord

A BOLD reference record used as a match. `key` is `"<process_id>|<marker>"`
(the portal has one record per specimen and marker). Fields from the
identification engine: `process_id`, `marker`, `bin_uri`, `taxid`, `lineage`.
Fields from the portal (`missing` when metadata were not retrieved):
`record_id`, `sequence_hash` (SHA-256 of the normalized `nuc`), `identified_by`
(verbatim), `identification`, `identification_method`, `identification_rank`,
`institution` (storing institution, never used for attribution),
`recordsets` (BOLD projects/datasets), `voucher_type`, `museum_id`,
`sequence_upload_date`, `species_reference`.
"""
struct ReferenceRecord
    key::String
    process_id::String
    marker::Union{Missing,String}
    bin_uri::Union{Missing,String}
    taxid::Union{Missing,String}
    lineage::Dict{Symbol,String}
    record_id::Union{Missing,String}
    sequence_hash::Union{Missing,String}
    identified_by::Union{Missing,String}
    identification::Union{Missing,String}
    identification_method::Union{Missing,String}
    identification_rank::Union{Missing,String}
    institution::Union{Missing,String}
    recordsets::Union{Missing,String}
    voucher_type::Union{Missing,String}
    museum_id::Union{Missing,String}
    sequence_upload_date::Union{Missing,String}
    species_reference::Union{Missing,String}
end

"""
    MatchRow

One match of a unique query sequence: the reference `key`, percent `identity`,
`coverage`, `align_length`, the LCA decision (`included`, `reason`) and
`dissent`, the rank at which a selected reference disagreed with an
assignment accepted by majority (`nothing` otherwise), and `lineage`, the
classification of the reference reported by BOLD with this hit (the one used
by the LCA). Hits of one analysis may classify the same reference
differently (e.g. one from an old checkpoint): see
[`reference_conflicts`](@ref).
"""
struct MatchRow
    key::String
    identity::Float64
    coverage::Union{Missing,Float64}
    align_length::Union{Missing,Int}
    included::Bool
    reason::Symbol
    dissent::Union{Nothing,Symbol}
    lineage::Dict{Symbol,String}
end

"""
    Analysis

Result of [`run_analysis`](@ref) or [`analyze_hits`](@ref): identifiers
(`id`, `fingerprint`), `created_at` (UTC), `project`, `user`, the `input`,
LCA `params`, `accept_rank`, data `source` (provenance of the BOLD data),
one `LCAResult` and the `MatchRow`s per unique sequence (`results`,
`matches`), the `references` by key, non-fatal `errors`, and the downloaded
data (`hits` from the identification engine, portal `records` or
`nothing`), kept so that [`reanalyze`](@ref) can recompute the LCA without
contacting BOLD again.
"""
struct Analysis
    id::String
    fingerprint::String
    created_at::DateTime
    project::String
    user::String
    input::InputSet
    params::LCAParams
    accept_rank::Symbol
    source::Dict{String,Any}
    results::Vector{LCAResult}
    matches::Vector{Vector{MatchRow}}
    references::Dict{String,ReferenceRecord}
    errors::Vector{String}
    hits::DataFrame
    records::Union{Nothing,DataFrame}
end

"""
    accepted(analysis, i) -> Bool

Whether the assignment of unique sequence `i` is accepted: assigned at
`analysis.accept_rank` or deeper.
"""
function accepted(a::Analysis, i::Integer)
    r = a.results[i]
    r.status === :assigned || return false
    ranks = a.params.ranks
    return findfirst(==(r.rank), ranks) >= findfirst(==(a.accept_rank), ranks)
end

"""
    supports(a::Analysis, i, m::MatchRow) -> Bool

Whether reference match `m` of unique sequence `i` supports its assignment:
in the LCA set and not dissenting at the assigned rank or above (a reference
dissenting only below the assigned rank, e.g. after the rank cap, supports
it).
"""
function supports(a::Analysis, i::Integer, m::MatchRow)
    m.included || return false
    m.dissent === nothing && return true
    r = a.results[i].rank
    r === nothing && return false
    ranks = a.params.ranks
    return findfirst(==(m.dissent), ranks) > findfirst(==(r), ranks)
end

"""
    params_dict(p::LCAParams) -> Dict

Parameters as a JSON-ready dictionary (sorted keys in the JSON output).
"""
params_dict(p::LCAParams) = Dict{String,Any}(
    "min_identity" => p.min_identity, "min_coverage" => p.min_coverage,
    "max_delta" => p.max_delta, "max_matches" => p.max_matches,
    "min_align_length" => p.min_align_length, "min_agreement" => p.min_agreement,
    "rank_cap" => p.rank_cap,
    "rank_thresholds" => Dict(String(k) => v for (k, v) in p.rank_thresholds),
    "ranks" => String.(p.ranks))

function _fingerprint(project, user, input, source, params, accept_rank)
    parts = [project, user, input.set_hash, string(source["database_code"]),
             JSON3.write(source["mode_parameters"]), LCA_VERSION,
             _canonical_json(params_dict(params)), String(accept_rank)]
    return bytes2hex(SHA.sha256(join(parts, '\n')))
end

_canonical_json(d::AbstractDict) =
    "{" * join((JSON3.write(string(k)) * ":" * _canonical_json(d[k]) for k in sort!(collect(keys(d)); by = string)), ",") * "}"
_canonical_json(v::AbstractVector) = "[" * join(_canonical_json.(v), ",") * "]"
_canonical_json(x) = JSON3.write(x)

_now_utc() = Dates.now(Dates.UTC)
_iso(t::DateTime) = Dates.format(t, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"

"""
    run_analysis(input::InputSet; project = "default", user = "", params = LCAParams(),
                 accept_rank = :phylum, db = :animal_public, mode = :rapid,
                 engine = IDEngine(), checkpoint = nothing, metadata = true,
                 portal_url = "https://portal.boldsystems.org",
                 progress = msg -> nothing) -> Analysis

Identify the unique sequences of `input` with the official BOLD v5
identification engine (`BOLDassignment.bold_identify`), retrieve the
metadata of the matched records from the BOLD portal API (`metadata = true`)
and compute the LCA of every sequence.

A failure of the portal does not stop the analysis: references keep the
fields of the identification engine, the failure is listed in `errors` and
in the provenance. `checkpoint` (a directory) stores the raw engine result of
every sequence and makes interrupted runs resumable; the provenance
(`source["hits_retrieval"]`) tells, for each sequence, whether its hits were
queried now or reused from the checkpoint and when they were retrieved
(checkpoint file date). `progress` receives
short status messages (used by the web interface).
"""
function run_analysis(input::InputSet; project::AbstractString = "default",
                      user::AbstractString = "", params::LCAParams = LCAParams(),
                      accept_rank::Symbol = :phylum, db = :animal_public, mode::Symbol = :rapid,
                      engine::IDEngine = IDEngine(), checkpoint = nothing, metadata::Bool = true,
                      portal_url::AbstractString = PORTAL_URL, progress = msg -> nothing)
    accept_rank in params.ranks || throw(ArgumentError("accept_rank :$accept_rank is not a rank"))
    code = BOLDassignment._database_code(db)
    mode_params = BOLDassignment._mode(mode)
    seqs = [u.sequence for u in input.unique]
    isempty(seqs) && throw(ArgumentError("no valid sequence to identify"))
    run_at = _now_utc()
    # Checkpoint files present before the call: a file still unchanged after
    # it was reused, and its date is the original retrieval from BOLD.
    files = checkpoint === nothing ? Dict{String,String}() :
        Dict(u.hash => joinpath(checkpoint, BOLDassignment._sequence_key(code, mode_params,
                                    normalize_sequence(u.sequence)) * ".jsonl") for u in input.unique)
    before = Dict(h => mtime(f) for (h, f) in files if isfile(f))
    progress("Submitting $(length(seqs)) unique sequence(s) to the BOLD identification engine")
    hits = bold_identify(seqs; db, mode, engine, checkpoint)
    hits_at = _now_utc()
    position = Dict(u.sequence => i for (i, u) in enumerate(input.unique))
    retrieval = _hits_retrieval(input, files, before, hits_at, _missing_results(input, hits, position))
    n_reused = retrieval["n_from_checkpoint"]
    n_reused > 0 && progress("Hits of $n_reused sequence(s) reused from the checkpoint " *
                             "(retrieved from BOLD between $(retrieval["checkpoint_oldest"]) and " *
                             "$(retrieval["checkpoint_newest"]))")
    errors = String[]
    records = nothing
    metadata_status = "not requested"
    metadata_at = nothing
    if metadata
        pids = unique!(collect(skipmissing(hits.process_id)))
        progress("Retrieving metadata of $(length(pids)) BOLD record(s) from the portal")
        try
            records = isempty(pids) ? DataFrame() : bold_records(pids; engine, base_url = portal_url)
            metadata_status = "retrieved"
            metadata_at = _iso(_now_utc())
        catch err
            err isa InterruptException && rethrow()
            msg = "BOLD portal metadata not retrieved: $(sprint(showerror, err))"
            push!(errors, first(msg, 1000))
            metadata_status = "failed"
        end
    end
    source = Dict{String,Any}(
        "service" => "BOLD Systems v5 identification engine (official Barcode ID API)",
        "engine_url" => engine.base_url,
        "database" => string(db), "database_code" => code,
        "mode" => string(mode),
        "mode_parameters" => Dict("mi" => mode_params.mi, "mo" => mode_params.mo,
                                  "maxh" => mode_params.maxh),
        "analysis_run_at" => _iso(run_at),
        "hits_retrieval" => retrieval,
        "checkpoint_configured" => checkpoint !== nothing,
        "checkpoint_used" => n_reused > 0,
        "portal_url" => metadata ? portal_url : nothing,
        "metadata_status" => metadata_status,
        "metadata_retrieved_at" => metadata_at,
    )
    progress("Computing LCA assignments")
    return analyze_hits(input, hits, records; project, user, params, accept_rank, source, errors)
end

# Reference keys whose hits carry more than one classification, with the
# distinct classifications (in order of first appearance).
function _conflicts(matches)
    seen = Dict{String,Vector{Dict{Symbol,String}}}()
    for ms in matches, m in ms
        v = get!(seen, m.key, Dict{Symbol,String}[])
        m.lineage in v || push!(v, m.lineage)
    end
    return Dict(k => v for (k, v) in seen if length(v) > 1)
end

_lineage_names(ls, ranks) = [coalesce(first(reference_taxon(l, ranks)), "unclassified") for l in ls]

"""
    reference_conflicts(a::Analysis) -> Dict{String,Vector{Dict{Symbol,String}}}

Reference records (by key) whose hits in `a` carry different
classifications, with the distinct classifications.
"""
reference_conflicts(a::Analysis) = _conflicts(a.matches)

"""
    is_incomplete(a::Analysis) -> Bool

Whether BOLD did not return the result of some sequence of `a` (status
`result_missing`). Incomplete analyses are not recorded in the usage
statistics; run the analysis again to query the missing sequences.
"""
is_incomplete(a::Analysis) = any(r -> r.status === :result_missing, a.results)

_RESULT_MISSING() = LCAResult(:result_missing, nothing, nothing, Pair{Symbol,String}[], Float64[],
                              nothing, false, missing, missing, 0, 0, 0, MatchDecision[],
                              Union{Nothing,Symbol}[], ["bold_result_missing"])

# Unique sequences whose result BOLD did not return (`result_received`
# false, BOLDassignment ≥ 0.1.2; older tables have no such column).
function _missing_results(input::InputSet, hits::AbstractDataFrame, position)
    absent = falses(length(input.unique))
    hasproperty(hits, :result_received) || return absent
    for (row, s) in enumerate(hits.sequence)
        i = get(position, normalize_sequence(s), 0)
        i == 0 || hits.result_received[row] || (absent[i] = true)
    end
    return absent
end

# Origin of the hits of every unique sequence: queried now, or reused from a
# checkpoint file written at an earlier date (file modification time, UTC).
function _hits_retrieval(input::InputSet, files, before, hits_at::DateTime, absent)
    by_sequence = Dict{String,Any}()
    saved = DateTime[]
    for (i, u) in enumerate(input.unique)
        t0 = get(before, u.hash, nothing)
        f = get(files, u.hash, nothing)
        if absent[i]
            by_sequence[u.hash] = Dict("source" => "missing_from_bold_result", "retrieved_at" => nothing)
        elseif t0 !== nothing && isfile(f) && mtime(f) == t0
            t = Dates.unix2datetime(t0)
            push!(saved, t)
            by_sequence[u.hash] = Dict("source" => "checkpoint", "retrieved_at" => _iso(t))
        else
            by_sequence[u.hash] = Dict("source" => "bold_identification_engine",
                                       "retrieved_at" => _iso(hits_at))
        end
    end
    n_cached = length(saved)
    n_absent = count(absent)
    n_new = length(input.unique) - n_cached - n_absent
    return Dict{String,Any}(
        "n_from_bold" => n_new, "n_from_checkpoint" => n_cached, "n_missing_from_bold" => n_absent,
        "retrieved_at" => n_new > 0 ? _iso(hits_at) : nothing,
        "checkpoint_oldest" => isempty(saved) ? nothing : _iso(minimum(saved)),
        "checkpoint_newest" => isempty(saved) ? nothing : _iso(maximum(saved)),
        "checkpoint_date_note" => "date of the checkpoint file (when the hits were saved), from the file system",
        "by_sequence" => by_sequence)
end

"""
    analyze_hits(input, hits, records = nothing; project = "default", user = "",
                 params = LCAParams(), accept_rank = :phylum, source = Dict(),
                 errors = String[]) -> Analysis

Offline part of [`run_analysis`](@ref): `hits` is a table returned by
`BOLDassignment.bold_identify` for the unique sequences of `input`, `records`
a table of portal records (`BOLDassignment.bold_records`) or `nothing`.
"""
function analyze_hits(input::InputSet, hits::AbstractDataFrame, records = nothing;
                      project::AbstractString = "default", user::AbstractString = "",
                      params::LCAParams = LCAParams(), accept_rank::Symbol = :phylum,
                      source::AbstractDict = Dict{String,Any}(), errors = String[])
    accept_rank in params.ranks || throw(ArgumentError("accept_rank :$accept_rank is not a rank"))
    source = Dict{String,Any}(source)
    get!(source, "database_code", "unknown")
    get!(source, "mode_parameters", Dict{String,Any}())
    portal = _portal_index(records)
    position = Dict(u.sequence => i for (i, u) in enumerate(input.unique))
    per_seq = [Int[] for _ in input.unique]
    absent = _missing_results(input, hits, position)
    for (row, s) in enumerate(hits.sequence)
        i = get(position, normalize_sequence(s), 0)
        i == 0 && continue
        ismissing(hits.similarity[row]) && continue
        push!(per_seq[i], row)
    end
    references = Dict{String,ReferenceRecord}()
    results = LCAResult[]
    matches = Vector{MatchRow}[]
    valid = BOLDassignment.is_valid_name
    for (i, rows) in enumerate(per_seq)
        seen = Set{String}()
        lca_matches = Match[]
        keys_i = String[]
        hashes = Union{Nothing,String}[]
        # Deterministic order, independent of the order of the BOLD results.
        sort!(rows; by = r -> (-hits.similarity[r], String(hits.process_id[r]),
                               coalesce(_col(hits, :marker, r), "")))
        for row in rows
            pid = String(hits.process_id[row])
            marker = _col(hits, :marker, row)
            key = string(pid, '|', coalesce(marker, ""))
            key in seen && continue
            push!(seen, key)
            lineage = Dict{Symbol,String}()
            for r in params.ranks
                hasproperty(hits, r) || continue
                v = hits[row, r]
                ismissing(v) || (lineage[r] = String(v))
            end
            ref = get!(references, key) do
                _reference(key, pid, marker, hits, row, lineage, get(portal, (pid, marker), nothing))
            end
            push!(lca_matches, Match(key, Float64(hits.similarity[row]),
                                     _col(hits, :query_coverage, row),
                                     _col(hits, :align_length, row), lineage))
            push!(keys_i, key)
            push!(hashes, coalesce(ref.sequence_hash, nothing))
        end
        res = absent[i] ? _RESULT_MISSING() :
              lca(lca_matches, params; valid_name = valid, sequence_hash = hashes)
        push!(results, res)
        push!(matches, [MatchRow(keys_i[k], lca_matches[k].identity, lca_matches[k].coverage,
                                 lca_matches[k].align_length, res.decisions[k].included,
                                 res.decisions[k].reason, res.dissent[k], lca_matches[k].lineage)
                         for k in eachindex(keys_i)])
    end
    conflicts = _conflicts(matches)
    if !isempty(conflicts)
        for (i, ms) in enumerate(matches)
            any(m -> haskey(conflicts, m.key), ms) || continue
            "reference_classification_conflict" in results[i].flags ||
                push!(results[i].flags, "reference_classification_conflict")
        end
        detail = join(sort!([string(references[k].process_id, " (", join(_lineage_names(ls, params.ranks), " vs "), ")")
                             for (k, ls) in conflicts]), "; ")
        msg = first("Warning: BOLD classified $(length(conflicts)) reference record(s) differently " *
                    "in hits of this analysis (e.g. hits reused from an older checkpoint): $detail. " *
                    "Each hit keeps its own classification in the LCA and in bold_matches.csv " *
                    "(flag reference_classification_conflict).", 1000)
        errors = collect(String, errors)
        msg in errors || push!(errors, msg)
    end
    if any(absent)
        ids = join((join(input.unique[i].ids, ";") for i in findall(absent)), ", ")
        msg = first("BOLD did not return the result of $(count(absent)) sequence(s) ($ids): " *
                    "status result_missing, not 'no match'; run the analysis again to query them", 1000)
        errors = collect(String, errors)
        msg in errors || push!(errors, msg)
    end
    created = _now_utc()
    fp = _fingerprint(String(project), String(user), input, source, params, accept_rank)
    return Analysis(string(uuid4()), fp, created, String(project), String(user), input, params,
                    accept_rank, source, results, matches, references, collect(String, errors),
                    DataFrame(hits; copycols = false),
                    records === nothing ? nothing : DataFrame(records; copycols = false))
end

"""
    reanalyze(a::Analysis; params = a.params, accept_rank = a.accept_rank,
              project = a.project, user = a.user) -> Analysis

New analysis computed from the data already downloaded for `a` (BOLD hits
and portal metadata): no request to BOLD. Use it to try other LCA
parameters. The provenance keeps the original data source and retrieval
date and records `reanalysis_of = a.id`. The BOLD library and search mode
cannot change (they decide which hits exist): a lower `min_identity` than
the minimum of the search mode (e.g. 94 % in rapid mode) finds no new hits.
"""
function reanalyze(a::Analysis; params::LCAParams = a.params, accept_rank::Symbol = a.accept_rank,
                   project::AbstractString = a.project, user::AbstractString = a.user)
    source = copy(a.source)
    source["reanalysis_of"] = a.id
    return analyze_hits(a.input, a.hits, a.records; project, user, params, accept_rank, source,
                        errors = copy(a.errors))
end

_col(df, name, row) = hasproperty(df, name) ? df[row, name] : missing

function _portal_index(records)
    index = Dict{Tuple{String,Union{Missing,String}},Any}()
    (records === nothing || nrow(records) == 0 || !hasproperty(records, :processid)) && return index
    for r in eachrow(records)
        pid = r.processid
        ismissing(pid) && continue
        marker = hasproperty(records, :marker_code) ? r.marker_code : missing
        index[(String(pid), ismissing(marker) ? missing : String(marker))] = r
    end
    return index
end

function _text_field(r, name)
    (r === nothing || !hasproperty(r, name)) && return missing
    v = r[name]
    (v === nothing || ismissing(v)) && return missing
    s = strip(string(v))
    # Arrays come from the portal as JSON text, e.g. ["UTOCP","DS-X"].
    if startswith(s, "[")
        try
            s = join((string(x) for x in JSON3.read(s)), "; ")
        catch
        end
    end
    return isempty(s) || s == "None" ? missing : String(s)
end

function _reference(key, pid, marker, hits, row, lineage, portal_row)
    nuc = _text_field(portal_row, :nuc)
    seq_hash = missing
    if !ismissing(nuc)
        seq_hash = try
            sequence_hash(nuc)
        catch
            missing
        end
    end
    return ReferenceRecord(key, pid, marker, _col(hits, :bin_uri, row), _col(hits, :taxid, row),
                           lineage, _text_field(portal_row, :record_id), seq_hash,
                           _text_field(portal_row, :identified_by),
                           _text_field(portal_row, :identification),
                           _text_field(portal_row, :identification_method),
                           _text_field(portal_row, :identification_rank),
                           _text_field(portal_row, :inst),
                           _text_field(portal_row, :bold_recordset_code_arr),
                           _text_field(portal_row, :voucher_type),
                           _text_field(portal_row, :museumid),
                           _text_field(portal_row, :sequence_upload_date),
                           _text_field(portal_row, :species_reference))
end

"""
    reference_taxon(ref::ReferenceRecord, ranks = DEFAULT_RANKS) -> (taxon, rank)
    reference_taxon(lineage::AbstractDict, ranks = DEFAULT_RANKS) -> (taxon, rank)

Deepest valid name of the reference lineage (`(missing, missing)` if none).
`ref.lineage` is the classification of the first hit of the record; the
classification of each hit is `MatchRow.lineage`.
"""
reference_taxon(ref::ReferenceRecord, ranks = DEFAULT_RANKS) = reference_taxon(ref.lineage, ranks)

function reference_taxon(lineage::AbstractDict, ranks = DEFAULT_RANKS)
    for r in reverse(collect(ranks))
        nm = get(lineage, r, nothing)
        nm !== nothing && BOLDassignment.is_valid_name(r, nm) && return (nm, String(r))
    end
    return (missing, missing)
end
