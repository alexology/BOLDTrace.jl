"""
    BOLDTrace

Identification of ASVs through the BOLD Systems v5 identification engine,
LCA assignment, traceability of the BOLD reference records used, attribution
to the people who identified the reference specimens, usage counters and
reproducible exports, with a web interface built on Genie.jl.

Pipeline: [`read_asvs`](@ref) → [`run_analysis`](@ref) (BOLD engine, portal
metadata, [`LCA.lca`](@ref)) → tables ([`assignment_table`](@ref),
[`reference_usage`](@ref), [`taxonomic_attribution`](@ref)) →
[`record_analysis!`](@ref) → [`export_analysis`](@ref). [`serve`](@ref)
starts the web application.
"""
module BOLDTrace

using BOLDassignment
using BOLDassignment: IDEngine, bold_identify, bold_records, normalize_sequence
using DataFrames
using Dates
using DBInterface
using FASTX
using HTTP
using JSON3
using SHA
using SQLite
using UUIDs

include("lca.jl")
using .LCA

export LCA, LCAParams, LCAResult, Match, MatchDecision, lca, LCA_VERSION,
       InputSet, InputRecord, UniqueSequence, read_asvs, input_summary, sequence_hash,
       Analysis, ReferenceRecord, MatchRow, run_analysis, analyze_hits, reanalyze, accepted, supports, reference_taxon,
       assignment_table, match_table, reference_usage, taxonomic_attribution,
       analysis_summary,
       CitationLink, CITATION_RELATIONS, CITATION_STATUSES, used_process_ids, provenance,
       write_table, export_file, export_analysis, EXPORT_FILES,
       UsageStore, open_store, record_analysis!, usage_stats, list_analyses, add_citation!,
       citation_links, GLOBAL_STATS_NOTE,
       serve, shutdown

include("fasta.jl")
include("analysis.jl")
include("tables.jl")
include("export.jl")
include("store.jl")
include("server.jl")

end
