"""
    LCA

Lowest common ancestor of the reference matches of one query sequence.
Independent of BOLD: matches are plain values with identity, coverage,
alignment length and a lineage.
"""
module LCA

export LCAParams, Match, LCAResult, MatchDecision, lca, LCA_VERSION, DEFAULT_RANKS,
       DEFAULT_RANK_THRESHOLDS

"Version of the LCA algorithm, recorded in every analysis report."
const LCA_VERSION = "1.1.0"

"Ranks from the most general to the most specific."
const DEFAULT_RANKS = (:phylum, :class, :order, :family, :genus, :species)

"""
Minimum identity (percent) of the best admissible match for an assignment to
reach each rank when `rank_cap = true` (same values as BOLDigger3 and
BOLDassignment).
"""
const DEFAULT_RANK_THRESHOLDS = Dict{Symbol,Float64}(
    :phylum => 50.0, :class => 75.0, :order => 85.0,
    :family => 90.0, :genus => 95.0, :species => 97.0)

"""
    LCAParams(; min_identity = 90, min_coverage = 80, max_delta = 1.0,
              max_matches = 100, min_align_length = 0, min_agreement = 100,
              rank_cap = true, rank_thresholds = DEFAULT_RANK_THRESHOLDS,
              ranks = DEFAULT_RANKS)

Parameters of [`lca`](@ref):

- `min_identity`: minimum percent identity of an admissible match;
- `min_coverage`: minimum percent of the query covered by the match; a match
  with unknown coverage is excluded when `min_coverage > 0`;
- `max_delta`: maximum identity difference (percentage points) from the best
  admissible match for a match to enter the LCA set;
- `max_matches`: maximum size of the LCA set (closest matches first);
- `min_align_length`: minimum alignment length in bases (0 = not used; a
  match with unknown length is excluded when it is positive);
- `min_agreement`: minimum percent of the selected references identified at a
  rank that must carry the most frequent name for that rank to be assigned.
  `100` (default) is the strict LCA: one discordant reference stops the
  assignment. A lower value, e.g. `95`, tolerates a few discordant
  references (probable misidentifications in the database): they are
  flagged as dissenting and leave the deeper ranks. Must be above 50 (a
  strict majority); ties for the most frequent name always stop;
- `rank_cap`: when true, the assignment cannot go deeper than the deepest rank
  whose `rank_thresholds` value the best admissible identity reaches (e.g. no
  species below 97 %), even if the LCA set agrees further.
"""
struct LCAParams
    min_identity::Float64
    min_coverage::Float64
    max_delta::Float64
    max_matches::Int
    min_align_length::Int
    rank_cap::Bool
    rank_thresholds::Dict{Symbol,Float64}
    ranks::Vector{Symbol}
    min_agreement::Float64
end

function LCAParams(; min_identity::Real = 90, min_coverage::Real = 80, max_delta::Real = 1.0,
                   max_matches::Integer = 100, min_align_length::Integer = 0,
                   min_agreement::Real = 100, rank_cap::Bool = true,
                   rank_thresholds = DEFAULT_RANK_THRESHOLDS, ranks = DEFAULT_RANKS)
    (isfinite(min_agreement) && 50 < min_agreement <= 100) ||
        throw(ArgumentError("min_agreement must be within (50, 100]"))
    for (name, v) in (("min_identity", min_identity), ("min_coverage", min_coverage))
        (isfinite(v) && 0 <= v <= 100) || throw(ArgumentError("$name must be within [0, 100]"))
    end
    (isfinite(max_delta) && 0 <= max_delta <= 100) ||
        throw(ArgumentError("max_delta must be within [0, 100]"))
    max_matches >= 1 || throw(ArgumentError("max_matches must be at least 1"))
    min_align_length >= 0 || throw(ArgumentError("min_align_length must be non-negative"))
    rs = Symbol[Symbol(r) for r in ranks]
    (isempty(rs) || !allunique(rs)) && throw(ArgumentError("ranks must be non-empty and unique"))
    limits = Dict{Symbol,Float64}()
    for (k, v) in pairs(rank_thresholds)
        r = Symbol(k)
        r in rs || continue
        (v isa Real && isfinite(v) && 0 <= v <= 100) ||
            throw(ArgumentError("rank threshold for :$r must be within [0, 100]"))
        limits[r] = Float64(v)
    end
    rank_cap && any(r -> !haskey(limits, r), rs) &&
        throw(ArgumentError("rank_cap needs a threshold for every rank"))
    return LCAParams(Float64(min_identity), Float64(min_coverage), Float64(max_delta),
                     Int(max_matches), Int(min_align_length), rank_cap, limits, rs,
                     Float64(min_agreement))
end

"""
    Match(id, identity, coverage, align_length, lineage)

One reference match: an identifier (e.g. a BOLD process ID), percent
`identity`, percent query `coverage` and `align_length` (`missing` when
unknown), and `lineage`, a `Dict{Symbol,String}` from rank to name. A rank
absent from `lineage` is unknown for that reference.
"""
struct Match
    id::String
    identity::Float64
    coverage::Union{Missing,Float64}
    align_length::Union{Missing,Int}
    lineage::Dict{Symbol,String}
end

"""
    MatchDecision

Outcome of one match: `included` in the LCA set or not, and `reason`:
`:selected`, `:low_identity`, `:low_coverage`, `:unknown_coverage`,
`:short_alignment`, `:unknown_length`, `:outside_delta`, `:over_max_matches`.
"""
struct MatchDecision
    included::Bool
    reason::Symbol
end

"""
    LCAResult

- `status`: `:assigned`, `:no_matches` (nothing retrieved),
  `:no_admissible_match` (all matches filtered out) or `:unresolved` (the
  LCA set does not agree even at the first rank);
- `taxon`, `rank`: deepest assigned taxon and its rank (`nothing` unless
  assigned);
- `lineage`: assigned name for each rank down to `rank`;
- `agreement`: for each rank of `lineage`, percent of the selected references
  identified at that rank (and consistent with the ranks above) that carry
  the assigned name;
- `dissent`: one entry per input match: the rank at which a selected
  reference disagreed with an assignment accepted by majority
  (`min_agreement < 100`), otherwise `nothing`;
- `lca_rank`: deepest rank supported by the LCA set before the rank cap;
- `rank_capped`: whether the rank cap made the assignment shallower;
- `best_identity`, `best_coverage`: of the best admissible match;
- `n_matches`, `n_admissible`, `n_selected`;
- `decisions`: one [`MatchDecision`](@ref) per input match, same order;
- `flags`: notes such as `"incomplete_taxonomy:species"`,
  `"discordant:genus"`, `"majority:species"` (assigned with dissenting
  references), `"inconsistent_hierarchy:genus"`,
  `"identical_reference_sequences_discordant"`, `"truncated_at_max_matches"`.

The rank is the depth of agreement of the selected references, **not** a
statistical confidence.
"""
struct LCAResult
    status::Symbol
    taxon::Union{Nothing,String}
    rank::Union{Nothing,Symbol}
    lineage::Vector{Pair{Symbol,String}}
    agreement::Vector{Float64}
    lca_rank::Union{Nothing,Symbol}
    rank_capped::Bool
    best_identity::Union{Missing,Float64}
    best_coverage::Union{Missing,Float64}
    n_matches::Int
    n_admissible::Int
    n_selected::Int
    decisions::Vector{MatchDecision}
    dissent::Vector{Union{Nothing,Symbol}}
    flags::Vector{String}
end

_name(m::Match, r::Symbol) = get(m.lineage, r, nothing)

"""
    lca(matches, params = LCAParams(); valid_name = (rank, name) -> !isempty(strip(name)),
        sequence_hash = nothing) -> LCAResult

Lowest common ancestor of `matches` (a vector of [`Match`](@ref)):

1. matches below `min_identity`, `min_coverage` or `min_align_length` are not
   admissible;
2. the best admissible match is the one with the highest identity;
3. admissible matches within `max_delta` of the best form the candidate set,
   ordered by identity, then coverage (descending), then identifier, and cut
   to `max_matches`;
4. ranks are visited from the most general: at each rank the most frequent
   name among the selected references identified at that rank is assigned if
   it is carried by at least `min_agreement` percent of them (100 = all:
   **strict LCA**); otherwise, or on a tie, the assignment stops at the
   previous rank. References with another name are *dissenting*: they are
   recorded in `dissent` and leave the deeper ranks;
5. a reference not identified at a rank (name missing, or rejected by
   `valid_name`, e.g. `Genus sp.`) does not vote at that rank or deeper ones,
   but still supports the ranks above; if no selected reference is identified
   at a rank, the assignment stops there (incomplete taxonomies never force a
   deeper rank);
6. with `rank_cap`, the result is cut to the deepest rank allowed by the
   identity of the best admissible match.

`sequence_hash`, when given, is a vector with the hash of each reference
sequence (`nothing` when unknown): selected references with identical
sequences but different names at the assigned or next rank are flagged.
"""
function lca(matches::AbstractVector{Match}, params::LCAParams = LCAParams();
             valid_name = (rank, name) -> !isempty(strip(name)),
             sequence_hash::Union{Nothing,AbstractVector} = nothing)
    n = length(matches)
    sequence_hash === nothing || length(sequence_hash) == n ||
        throw(ArgumentError("sequence_hash must have one entry per match"))
    decisions = Vector{MatchDecision}(undef, n)
    dissent = Vector{Union{Nothing,Symbol}}(nothing, n)
    flags = String[]
    n == 0 && return LCAResult(:no_matches, nothing, nothing, Pair{Symbol,String}[], Float64[],
                               nothing, false, missing, missing, 0, 0, 0, decisions, dissent, flags)
    admissible = Int[]
    for (i, m) in enumerate(matches)
        (isfinite(m.identity) && 0 <= m.identity <= 100) ||
            throw(ArgumentError("match $(m.id): identity must be within [0, 100]"))
        ismissing(m.coverage) || (isfinite(m.coverage) && 0 <= m.coverage <= 100) ||
            throw(ArgumentError("match $(m.id): coverage must be missing or within [0, 100]"))
        ismissing(m.align_length) || m.align_length >= 0 ||
            throw(ArgumentError("match $(m.id): align_length must be missing or non-negative"))
        reason = if m.identity < params.min_identity
            :low_identity
        elseif params.min_coverage > 0 && ismissing(m.coverage)
            :unknown_coverage
        elseif !ismissing(m.coverage) && m.coverage < params.min_coverage
            :low_coverage
        elseif params.min_align_length > 0 && ismissing(m.align_length)
            :unknown_length
        elseif !ismissing(m.align_length) && m.align_length < params.min_align_length
            :short_alignment
        else
            :admissible
        end
        decisions[i] = MatchDecision(false, reason)
        reason === :admissible && push!(admissible, i)
    end
    if isempty(admissible)
        return LCAResult(:no_admissible_match, nothing, nothing, Pair{Symbol,String}[], Float64[],
                         nothing, false, missing, missing, n, 0, 0, decisions, dissent, flags)
    end
    order = sort(admissible; by = i -> (-matches[i].identity,
                                         -coalesce(matches[i].coverage, -1.0), matches[i].id))
    best = matches[order[1]]
    floor_identity = best.identity - params.max_delta
    candidates = Int[]
    for i in order
        if matches[i].identity >= floor_identity - 1e-9
            push!(candidates, i)
        else
            decisions[i] = MatchDecision(false, :outside_delta)
        end
    end
    selected = candidates[1:min(end, params.max_matches)]
    for i in candidates[length(selected)+1:end]
        decisions[i] = MatchDecision(false, :over_max_matches)
    end
    if length(candidates) > length(selected)
        cut = matches[candidates[length(selected)+1]].identity
        cut == matches[selected[end]].identity && push!(flags, "truncated_at_max_matches_within_ties")
        push!(flags, "truncated_at_max_matches")
    end
    for i in selected
        decisions[i] = MatchDecision(true, :selected)
    end

    lineage = Pair{Symbol,String}[]
    agreement = Float64[]
    # Still voting: identified at every rank so far and agreeing with the
    # names assigned so far.
    active = trues(length(selected))
    for r in params.ranks
        names = Dict{String,Int}()
        n_unidentified = 0
        for (k, i) in enumerate(selected)
            active[k] || continue
            nm = _name(matches[i], r)
            if nm === nothing || !valid_name(r, nm)
                active[k] = false
                n_unidentified += 1
                continue
            end
            names[strip(nm)] = get(names, strip(nm), 0) + 1
        end
        n_unidentified > 0 && push!(flags, "incomplete_taxonomy:$r")
        isempty(names) && break
        total = sum(values(names))
        top = maximum(values(names))
        if length(names) > 1
            push!(flags, "discordant:$r")
            # Ties never resolve; otherwise the majority must reach min_agreement.
            count(==(top), values(names)) == 1 || break
            100 * top / total >= params.min_agreement * (1 - 1e-12) || break
            push!(flags, "majority:$r")
        end
        name = first(k for (k, v) in names if v == top)
        for (k, i) in enumerate(selected)
            active[k] || continue
            if strip(_name(matches[i], r)) != name
                dissent[i] = r
                active[k] = false
            end
        end
        push!(lineage, r => name)
        push!(agreement, 100 * top / total)
    end
    _check_hierarchy!(flags, matches, selected, params.ranks, valid_name)
    sequence_hash === nothing ||
        _check_identical!(flags, matches, selected, sequence_hash, params.ranks, valid_name)

    lca_rank = isempty(lineage) ? nothing : first(last(lineage))
    capped = false
    if params.rank_cap && !isempty(lineage)
        allowed = [p for p in lineage if best.identity >= params.rank_thresholds[first(p)]]
        if length(allowed) < length(lineage)
            capped = true
            lineage = allowed
            agreement = agreement[1:length(allowed)]
        end
    end
    if isempty(lineage)
        return LCAResult(:unresolved, nothing, nothing, lineage, agreement, lca_rank, capped,
                         best.identity, best.coverage, n, length(admissible), length(selected),
                         decisions, dissent, flags)
    end
    return LCAResult(:assigned, last(lineage).second, last(lineage).first, lineage, agreement,
                     lca_rank, capped, best.identity, best.coverage, n, length(admissible),
                     length(selected), decisions, dissent, flags)
end

# The same name at a rank placed under different parents among the selected
# references (e.g. one genus in two families) signals inconsistent
# classifications in the reference database.
function _check_hierarchy!(flags, matches, selected, ranks, valid_name)
    for j in 2:length(ranks)
        r, parent = ranks[j], ranks[j-1]
        parents = Dict{String,Set{String}}()
        for i in selected
            nm, pn = _name(matches[i], r), _name(matches[i], parent)
            (nm === nothing || pn === nothing) && continue
            (valid_name(r, nm) && valid_name(parent, pn)) || continue
            push!(get!(parents, strip(nm), Set{String}()), strip(pn))
        end
        any(s -> length(s) > 1, values(parents)) && push!(flags, "inconsistent_hierarchy:$r")
    end
end

# Selected references with identical sequences carrying different names at
# some rank: identical barcodes cannot discriminate those taxa.
function _check_identical!(flags, matches, selected, hashes, ranks, valid_name)
    groups = Dict{String,Vector{Int}}()
    for i in selected
        h = hashes[i]
        (h === nothing || ismissing(h)) && continue
        push!(get!(groups, String(h), Int[]), i)
    end
    any(g -> length(g) > 1, values(groups)) && push!(flags, "identical_reference_sequences")
    # Most general rank at which some group of identical sequences disagrees.
    # (Separate loops: `break` in `for a in A, b in B` would leave both.)
    first_rank = 0
    for g in values(groups)
        length(g) > 1 || continue
        for (k, r) in enumerate(ranks)
            nms = Set{String}()
            for i in g
                nm = _name(matches[i], r)
                nm !== nothing && valid_name(r, nm) && push!(nms, strip(nm))
            end
            if length(nms) > 1
                first_rank = first_rank == 0 ? k : min(first_rank, k)
                break
            end
        end
    end
    first_rank > 0 && push!(flags, "identical_reference_sequences_discordant:$(ranks[first_rank])")
    return nothing
end

end # module LCA
