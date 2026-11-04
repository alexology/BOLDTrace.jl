# Bibliographic links (data structure only), provenance and export files.

"Relation types between a BOLD record and a publication."
const CITATION_RELATIONS = ("original_barcode_publication", "taxonomic_revision",
                            "barcode_reuse_publication")

"Verification states of a bibliographic link."
const CITATION_STATUSES = ("unverified", "verified", "rejected")

const USED_PROCESS_ID_COLUMNS = (:process_id, :record_id, :marker, :lca_participation,
                                 :assignment_contribution, :doi, :relation, :evidence_source,
                                 :verification_status)

"""
    CitationLink(process_id, doi, relation, evidence_source, verification_status = "unverified")

Link between a BOLD process ID and a publication, for a future independent
registry. `relation` is one of [`CITATION_RELATIONS`](@ref) and
`verification_status` one of [`CITATION_STATUSES`](@ref). A publication
linked to a BOLD project does not necessarily document every barcode of the
project: links are per process ID and carry their evidence. No literature
search is performed by this package.
"""
struct CitationLink
    process_id::String
    doi::String
    relation::String
    evidence_source::String
    verification_status::String
    function CitationLink(process_id, doi, relation, evidence_source,
                          verification_status = "unverified")
        isempty(strip(process_id)) && throw(ArgumentError("empty process ID"))
        occursin(r"^10\.\d{4,9}/\S+$", strip(doi)) || throw(ArgumentError("invalid DOI: $doi"))
        relation in CITATION_RELATIONS || throw(ArgumentError("unknown relation: $relation"))
        verification_status in CITATION_STATUSES ||
            throw(ArgumentError("unknown verification status: $verification_status"))
        new(strip(process_id), strip(doi), relation, evidence_source, verification_status)
    end
end

"""
    used_process_ids(a::Analysis; citations = CitationLink[]) -> Vector{NamedTuple}

One row per reference record that entered at least one LCA set: `process_id`,
`record_id`, `marker`, `lca_participation`, `assignment_contribution` and
the citation columns (`doi`, `relation`, `evidence_source`,
`verification_status`), empty unless `citations` holds links for that
process ID (one row per link).
"""
function used_process_ids(a::Analysis; citations::AbstractVector{CitationLink} = CitationLink[])
    by_pid = Dict{String,Vector{CitationLink}}()
    for c in citations
        push!(get!(by_pid, c.process_id, CitationLink[]), c)
    end
    rows = NamedTuple[]
    for u in reference_usage(a)
        u.lca_participation > 0 || continue
        base = (process_id = u.process_id, record_id = u.record_id, marker = u.marker,
                lca_participation = u.lca_participation,
                assignment_contribution = u.assignment_contribution)
        links = get(by_pid, u.process_id, CitationLink[])
        if isempty(links)
            push!(rows, merge(base, (doi = missing, relation = missing, evidence_source = missing,
                                     verification_status = missing)))
        else
            for c in links
                push!(rows, merge(base, (doi = c.doi, relation = c.relation,
                                         evidence_source = c.evidence_source,
                                         verification_status = c.verification_status)))
            end
        end
    end
    return rows
end

"""
    provenance(a::Analysis) -> Dict

Reproducibility report: analysis identifier and fingerprint, date (UTC),
project and user, input hashes, BOLD data source and retrieval dates (run
date, origin and date of the hits of each sequence: queried now or reused
from the checkpoint, date of the portal metadata), alignment software, LCA version and parameters, software versions, selected
process IDs per sequence, final assignments and known limitations.
"""
function provenance(a::Analysis)
    selected = Dict{String,Vector{String}}()
    assignments = Any[]
    for (i, u) in enumerate(a.input.unique)
        res = a.results[i]
        selected[u.hash] = selected_process_ids(a, i)
        push!(assignments, Dict("sequence_hash" => u.hash, "asv_ids" => u.ids,
                                "status" => String(res.status),
                                "taxon" => res.taxon,
                                "rank" => res.rank === nothing ? nothing : String(res.rank),
                                "accepted" => accepted(a, i), "flags" => res.flags))
    end
    limitations = String[
        "The BOLD API does not report a version or snapshot of the reference libraries; " *
        "the retrieval dates are recorded instead (hits per sequence, portal metadata) and results " *
        "may change as BOLD is updated.",
        "Alignments are computed by the BOLD identification engine on the BOLD servers; " *
        "the API does not report the alignment software or its version.",
        "Global statistics count only the analyses recorded by this application, not all uses of BOLD.",
        "The LCA rank is the depth of agreement among selected references, not a statistical confidence.",
    ]
    hr = get(a.source, "hits_retrieval", nothing)
    if hr isa AbstractDict && get(hr, "n_from_checkpoint", 0) > 0
        push!(limitations, "The BOLD hits of $(hr["n_from_checkpoint"]) sequence(s) were reused from " *
                           "the local checkpoint, retrieved from BOLD between $(hr["checkpoint_oldest"]) " *
                           "and $(hr["checkpoint_newest"]) (file dates); BOLD may have changed since. " *
                           "Per-sequence origin and date: bold_data_source.hits_retrieval.by_sequence.")
    end
    get(a.source, "metadata_status", "") == "retrieved" ||
        push!(limitations, "Reference metadata from the BOLD portal are missing " *
                           "(identified_by, reference sequences); attribution and identical-sequence checks are incomplete.")
    return Dict{String,Any}(
        "analysis_id" => a.id, "fingerprint" => a.fingerprint,
        "reanalysis_of" => get(a.source, "reanalysis_of", nothing),
        "created_at" => _iso(a.created_at), "project" => a.project, "user" => a.user,
        "input" => Dict("name" => a.input.name, "file_sha256" => a.input.file_sha256,
                        "sequence_set_hash" => a.input.set_hash,
                        "sequence_hash_method" => "SHA-256 of the normalized sequence (no gaps/whitespace, upper case, U→T)",
                        "n_records" => length(a.input.records),
                        "n_unique_sequences" => length(a.input.unique),
                        "n_input_errors" => count(r -> r.error !== nothing, a.input.records)),
        "bold_data_source" => a.source,
        "database_version" => "not provided by the BOLD API (see limitations)",
        "alignment_software" => Dict("name" => "BOLD identification engine (server side)",
                                     "version" => "not reported by the API"),
        "lca" => Dict("version" => LCA_VERSION, "parameters" => params_dict(a.params),
                      "accept_rank" => String(a.accept_rank)),
        "software" => Dict("BOLDTrace" => string(pkgversion(@__MODULE__)),
                           "BOLDassignment" => string(pkgversion(BOLDassignment)),
                           "julia" => string(VERSION)),
        "selected_process_ids" => selected,
        "assignments" => assignments,
        "errors" => a.errors,
        "limitations" => limitations,
    )
end

# --- CSV/TSV writing ---------------------------------------------------------

_cell(::Missing) = ""
_cell(::Nothing) = ""
_cell(x::Bool) = x ? "true" : "false"
_cell(x::AbstractFloat) = isinteger(x) && abs(x) < 1e15 ? string(Int(x)) : string(x)
_cell(x) = string(x)

function _csv_field(s::AbstractString)
    # Spreadsheet formula injection: prefix cells starting with = + - @.
    t = !isempty(s) && first(s) in ('=', '+', '-', '@') && tryparse(Float64, s) === nothing ?
        "'" * s : s
    return any(c -> c in (',', '"', '\n', '\r'), t) ? "\"" * replace(t, "\"" => "\"\"") * "\"" : t
end

function _tsv_field(s::AbstractString)
    t = replace(s, '\t' => ' ', '\n' => ' ', '\r' => ' ')
    return !isempty(t) && first(t) in ('=', '+', '-', '@') && tryparse(Float64, t) === nothing ?
           "'" * t : t
end

"""
    write_table(io_or_path, rows; delim = ',', columns = nothing)

Write a vector of NamedTuples as CSV (`delim = ','`, RFC 4180 quoting) or TSV
(`delim = '\\t'`, tabs and newlines replaced by spaces). Cells starting with
`= + - @` (not numbers) get a leading apostrophe against formula injection.
`columns`, when given, fixes the header and the order of the columns (also
with no rows); otherwise the fields of the first row are used.
"""
function write_table(io::IO, rows::AbstractVector; delim::Char = ',', columns = nothing)
    cols = columns !== nothing ? collect(Symbol, columns) :
           isempty(rows) ? Symbol[] : collect(keys(first(rows)))
    field = delim == ',' ? _csv_field : _tsv_field
    println(io, join((field(string(c)) for c in cols), delim))
    for r in rows
        println(io, join((field(_cell(getproperty(r, c))) for c in cols), delim))
    end
    return nothing
end

write_table(path::AbstractString, rows::AbstractVector; kwargs...) =
    open(io -> write_table(io, rows; kwargs...), path, "w")

"Names of the export files, in the order of the Export page."
const EXPORT_FILES = ("assignments.csv", "bold_matches.csv", "reference_usage.csv",
                      "taxonomic_attribution.csv", "analysis_provenance.json",
                      "used_process_ids.tsv")

"""
    export_file(a::Analysis, name; citations = CitationLink[], sequence_hash = nothing) -> String

Content of one export file (`name` in [`EXPORT_FILES`](@ref)). For
`bold_matches.csv`, `sequence_hash` keeps only the matches of one unique
sequence (all its input IDs), with the same columns.
"""
function export_file(a::Analysis, name::AbstractString;
                     citations::AbstractVector{CitationLink} = CitationLink[],
                     sequence_hash::Union{Nothing,AbstractString} = nothing)
    sequence_hash === nothing || name == "bold_matches.csv" ||
        throw(ArgumentError("sequence_hash applies only to bold_matches.csv"))
    io = IOBuffer()
    if name == "assignments.csv"
        write_table(io, assignment_table(a); columns = assignment_columns(a.params.ranks))
    elseif name == "bold_matches.csv"
        rows = match_table(a)
        if sequence_hash !== nothing
            any(u -> u.hash == sequence_hash, a.input.unique) ||
                throw(ArgumentError("unknown sequence hash"))
            rows = filter(r -> r.sequence_hash == sequence_hash, rows)
        end
        write_table(io, rows; columns = match_columns(a.params.ranks))
    elseif name == "reference_usage.csv"
        write_table(io, reference_usage(a); columns = REFERENCE_USAGE_COLUMNS)
    elseif name == "taxonomic_attribution.csv"
        write_table(io, taxonomic_attribution(a); columns = ATTRIBUTION_COLUMNS)
    elseif name == "analysis_provenance.json"
        JSON3.pretty(io, JSON3.write(provenance(a)))
        println(io)
    elseif name == "used_process_ids.tsv"
        write_table(io, used_process_ids(a; citations); delim = '\t',
                    columns = USED_PROCESS_ID_COLUMNS)
    else
        throw(ArgumentError("unknown export file: $name"))
    end
    return String(take!(io))
end

"""
    export_analysis(a::Analysis, dir; citations = CitationLink[]) -> Vector{String}

Write the six export files ([`EXPORT_FILES`](@ref)) to `dir` (created if
needed) and return their paths. Existing files are replaced.
"""
function export_analysis(a::Analysis, dir::AbstractString;
                         citations::AbstractVector{CitationLink} = CitationLink[])
    mkpath(dir)
    paths = String[]
    for name in EXPORT_FILES
        p = joinpath(dir, name)
        write(p, export_file(a, name; citations))
        push!(paths, p)
    end
    return paths
end
