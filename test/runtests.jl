using BOLDTrace
using BOLDTrace.LCA
using BOLDassignment
using BOLDassignment: IDEngine
using DataFrames
using FASTX
using HTTP
using JSON3
using Test

include("fake.jl")

# Deterministic pseudo-random DNA (no dependency on the RNG of the Julia version).
function dna(seed, n = 150)
    x = UInt64(seed) * 0x9E3779B97F4A7C15 + 1
    io = IOBuffer()
    for _ in 1:n
        x = x * 6364136223846793005 + 1442695040888963407
        print(io, "ACGT"[(x >> 33) % 4 + 1])
    end
    return String(take!(io))
end

lineage(; phylum = "Arthropoda", class = "Insecta", order = "Coleoptera",
        family = "Scarabaeidae", genus = "Onthophagus", species = "Onthophagus taurus") =
    Dict{Symbol,String}(k => v for (k, v) in (:phylum => phylum, :class => class, :order => order,
                                              :family => family, :genus => genus,
                                              :species => species) if v !== nothing)
m(id, identity, lin; cov = 100.0, len = 300) = Match(id, identity, cov, len, lin)
valid = BOLDassignment.is_valid_name

@testset "BOLDTrace" begin

@testset "FASTA input" begin
    a, b, c = dna(1), dna(2), dna(3)
    wrapped = a[1:60] * "\n" * a[61:end]
    text = ">ASV1 sample=1\n$wrapped\n>ASV2\n$(lowercase(b))\n>ASV3\n$(a[1:20])-$(a[21:end])\n" *
           ">bad\nACGTXXACGT\n>ASV2\n$c\n>ASV4\nACGTACGT\n>ASV5\n$(replace(b, 'T' => 'U'))\n"
    input = read_asvs(Vector{UInt8}(codeunits(text)); name = "t.fasta")
    @test [r.id for r in input.records] == ["ASV1", "ASV2", "ASV3", "bad", "ASV2", "ASV4", "ASV5"]
    @test input.records[1].header == "ASV1 sample=1"
    @test length(input.unique) == 3
    @test input.unique[1].ids == ["ASV1", "ASV3"]          # gap removed: same sequence
    @test input.unique[2].ids == ["ASV2", "ASV5"]          # case and U→T normalized
    @test input.records[1].hash == sequence_hash(a) == sequence_hash(lowercase(a))
    @test input.records[1].hash == bytes2hex(BOLDTrace.SHA.sha256(a))
    @test occursin("invalid nucleotide", input.records[4].error)
    @test occursin("different sequence", input.records[5].error)
    @test any(occursin("8 bases, shorter than 100", w) for w in input.messages)
    s = input_summary(input)
    @test (s.n_records, s.n_valid, s.n_unique, s.n_duplicates, s.n_errors) == (7, 5, 3, 2, 2)
    # Set hash: independent of order, identifiers and wrapping.
    other = read_asvs(Vector{UInt8}(codeunits(">x\nACGTACGT\n>y\n$b\n>z\n$a\n")))
    @test other.set_hash == input.set_hash
    @test other.file_sha256 != input.file_sha256
    @test_throws ArgumentError read_asvs(UInt8[])
    @test_throws ArgumentError read_asvs(Vector{UInt8}(codeunits("no header\nACGT\n")))
    @test_throws ArgumentError read_asvs(Vector{UInt8}(codeunits(">a\nACGT\n>b\nACGT\n")); max_records = 1)
    path = joinpath(mktempdir(), "in.fasta")
    write(path, text)
    @test read_asvs(path).name == "in.fasta"
end

@testset "LCA cases" begin
    p = LCAParams(; min_identity = 90, min_coverage = 80, max_delta = 1.0)
    # Concordant at species.
    r = lca([m("a", 100, lineage()), m("b", 99.5, lineage())], p; valid_name = valid)
    @test (r.status, r.taxon, r.rank, r.n_selected) == (:assigned, "Onthophagus taurus", :species, 2)
    @test all(d -> d.included && d.reason == :selected, r.decisions)
    # Different species of one genus → genus.
    r = lca([m("a", 100, lineage()), m("b", 99.8, lineage(species = "Onthophagus illyricus"))], p;
            valid_name = valid)
    @test (r.taxon, r.rank) == ("Onthophagus", :genus)
    @test "discordant:species" in r.flags
    # Different genera of one family → family.
    r = lca([m("a", 99, lineage()), m("b", 98.5, lineage(genus = "Caccobius", species = "Caccobius schreberi"))],
            p; valid_name = valid)
    @test (r.taxon, r.rank) == ("Scarabaeidae", :family)
    # Discordant classification above: strict stop, no majority vote.
    r = lca([m("a", 99, lineage()), m("b", 99, lineage()),
             m("c", 99, lineage(family = "Geotrupidae", genus = "Onthophagus"))], p; valid_name = valid)
    @test (r.taxon, r.rank) == ("Coleoptera", :order)
    @test "discordant:family" in r.flags && "inconsistent_hierarchy:genus" in r.flags
    # Incomplete taxonomy: open nomenclature and missing ranks do not vote but
    # support higher ranks; if nothing is identified at a rank, stop there.
    r = lca([m("a", 100, lineage()), m("b", 100, lineage(species = "Onthophagus sp. 2")),
             m("c", 99.9, lineage(species = nothing))], p; valid_name = valid)
    @test (r.taxon, r.rank) == ("Onthophagus taurus", :species)
    @test "incomplete_taxonomy:species" in r.flags
    r = lca([m("a", 100, lineage(species = nothing)), m("b", 100, lineage(species = "Onthophagus cf. taurus"))],
            p; valid_name = valid)
    @test (r.taxon, r.rank) == ("Onthophagus", :genus)
    r = lca([m("a", 100, lineage(genus = nothing, species = nothing))], p; valid_name = valid)
    @test r.rank == :family
    # Identical reference sequences with different species.
    r = lca([m("a", 100, lineage()), m("b", 100, lineage(species = "Onthophagus illyricus"))], p;
            valid_name = valid, sequence_hash = ["h1", "h1"])
    @test r.rank == :genus
    @test "identical_reference_sequences" in r.flags
    @test "identical_reference_sequences_discordant:species" in r.flags
    # Singleton groups before the discordant pair (case found on real BOLD data:
    # O. illyricus TDAOE1457-22 has the sequence of an O. taurus record).
    for order in ([1, 2, 3, 4], [4, 3, 2, 1], [2, 4, 1, 3])
        ms = [m("a", 100, lineage()), m("b", 100, lineage()), m("c", 99.9, lineage()),
              m("d", 99.7, lineage(species = "Onthophagus illyricus"))]
        hs = ["h1", "h2", "h3", "h3"]
        r = lca(ms[order], p; valid_name = valid, sequence_hash = hs[order])
        @test "identical_reference_sequences_discordant:species" in r.flags
    end
    r = lca([m("a", 100, lineage()), m("b", 100, lineage(genus = "Caccobius", species = "Caccobius schreberi")),
             m("c", 100, lineage()), m("d", 100, lineage(species = "Onthophagus illyricus"))], p;
            valid_name = valid, sequence_hash = ["h1", "h1", "h2", "h2"])
    @test "identical_reference_sequences_discordant:genus" in r.flags     # most general rank
    @test !("identical_reference_sequences_discordant:species" in r.flags)
    r = lca([m("a", 100, lineage()), m("b", 100, lineage())], p; valid_name = valid,
            sequence_hash = ["h1", nothing])
    @test !any(startswith("identical"), r.flags)
    # No admissible match / no match.
    r = lca([m("a", 85, lineage()), m("b", 99, lineage(); cov = 50.0)], p)
    @test r.status == :no_admissible_match && r.taxon === nothing && r.n_admissible == 0
    @test [d.reason for d in r.decisions] == [:low_identity, :low_coverage]
    @test lca(Match[], p).status == :no_matches
    # Coverage or length unknown, minimum alignment length.
    r = lca([Match("a", 99, missing, 300, lineage())], p)
    @test r.decisions[1].reason == :unknown_coverage
    @test lca([Match("a", 99, missing, 300, lineage())], LCAParams(; min_coverage = 0)).status == :assigned
    q = LCAParams(; min_align_length = 200)
    @test lca([m("a", 99, lineage(); len = 150)], q).decisions[1].reason == :short_alignment
    @test lca([Match("a", 99, 100.0, missing, lineage())], q).decisions[1].reason == :unknown_length
    # Distance from the best and maximum number of matches.
    r = lca([m("a", 100, lineage()), m("b", 98.9, lineage(species = "Onthophagus illyricus"))], p;
            valid_name = valid)
    @test r.rank == :species && r.decisions[2].reason == :outside_delta
    r = lca([m("a", 100, lineage()), m("b", 98.9, lineage(species = "Onthophagus illyricus"))],
            LCAParams(; max_delta = 2); valid_name = valid)
    @test r.rank == :genus
    r = lca([m("c", 99.5, lineage()), m("a", 100, lineage()), m("b", 99.8, lineage())],
            LCAParams(; max_matches = 2))
    @test [d.reason for d in r.decisions] == [:over_max_matches, :selected, :selected]
    @test "truncated_at_max_matches" in r.flags
    # Rank cap: the LCA set agrees on the species but the best match is 96 %.
    r = lca([m("a", 96, lineage())], p; valid_name = valid)
    @test (r.rank, r.lca_rank, r.rank_capped) == (:genus, :species, true)
    @test lca([m("a", 96, lineage())], LCAParams(; rank_cap = false)).rank == :species
    r = lca([m("a", 91, lineage())], p)
    @test r.rank == :family && r.rank_capped
    # Parameter checks.
    @test_throws ArgumentError LCAParams(; min_identity = 101)
    @test_throws ArgumentError LCAParams(; max_matches = 0)
    @test_throws ArgumentError LCAParams(; max_delta = -1)
    @test_throws ArgumentError lca([m("a", 120, lineage())])
    # Deterministic order for ties: identifier breaks the tie.
    r = lca([m("z", 100, lineage()), m("y", 100, lineage())], LCAParams(; max_matches = 1))
    @test r.decisions == [MatchDecision(false, :over_max_matches), MatchDecision(true, :selected)]
    @test "truncated_at_max_matches_within_ties" in r.flags
end

@testset "Real BOLD response" begin
    # Identification engine response for CASBI15164-26 (Onthophagus taurus) and
    # a low-complexity sequence (public.tax-derep, rapid mode, 2026-10-23), the
    # fixture of BOLDassignment.jl.
    body = read(joinpath(@__DIR__, "data", "bold_results_taurus.jsonl"))
    seqs = [String(JSON3.read(l)[:sequence]) for l in eachline(IOBuffer(body))]
    input = read_asvs(Vector{UInt8}(codeunits(join((">q$i
$s
" for (i, s) in enumerate(seqs))))))
    hits = BOLDassignment._hits_table([BOLDassignment.normalize_sequence(s) for s in seqs], body)
    a = analyze_hits(input, hits)
    r = a.results[1]
    # Nine records at 100 %, one O. illyricus at 99.68 % within 1 point: genus.
    @test (r.status, r.taxon, r.rank) == (:assigned, "Onthophagus", :genus)
    @test "discordant:species" in r.flags && "incomplete_taxonomy:species" in r.flags
    @test r.n_matches == 14 && r.best_coverage == 100
    @test a.results[2].status == :no_matches
    # The two hits covering 91 % of the query are excluded with min_coverage = 95.
    b = analyze_hits(input, hits; params = LCAParams(; min_coverage = 95))
    @test count(m -> m.reason == :low_coverage, b.matches[1]) == 2
    # Without metadata the provenance states the limitation.
    @test any(occursin("metadata from the BOLD portal are missing", l) for l in provenance(a)["limitations"])
end

# Shared scenario: three ASVs (one duplicated), one without matches.
s1, s2, s3, s4 = dna(11), dna(12), dna(13), dna(14)
nuc(i) = dna(100 + i, 600)
catalog = Dict{String,Vector{Any}}(
    s1 => [("P1", "COI-5P", 0.0, 100, 150, tax("Onthophagus", "Onthophagus taurus")),
           ("P2", "COI-5P", 0.3, 100, 150, tax("Onthophagus", "Onthophagus taurus")),
           ("P3", "COI-5P", 0.5, 99, 149, tax("Onthophagus", "Onthophagus illyricus")),
           ("P9", "COI-5P", 3.0, 100, 150, tax("Onthophagus", "Onthophagus vacca"))],
    s2 => [("P4", "COI-5P", 0.0, 100, 150, tax("Caccobius", "Caccobius schreberi")),
           ("P1", "COI-5P", 4.0, 60, 90, tax("Onthophagus", "Onthophagus taurus"))],
    s3 => [("P5", "COI-5P", 3.5, 100, 150, tax("Aphodius", "Aphodius sp.")),
           ("P6", "COI-5P", 3.8, 100, 150, tax("Aphodius", "Aphodius fimetarius"))],
)
portal = Dict(
    "P1" => [doc("P1", "COI-3P", nuc(99), "Wrong Marker"), doc("P1", "COI-5P", nuc(1), "A. Rossi")],
    "P2" => [doc("P2", "COI-5P", nuc(1), "A. Rossi")],      # identical sequence to P1
    "P3" => [doc("P3", "COI-5P", nuc(3), "B. Bianchi")],
    "P4" => [doc("P4", "COI-5P", nuc(4), "Rossi, A.")],     # not merged with "A. Rossi"
    "P5" => [doc("P5", "COI-5P", nuc(5), nothing)],
    "P6" => [doc("P6", "COI-5P", nuc(6), "C. Verdi")],
    "P9" => [doc("P9", "COI-5P", nuc(9), "B. Bianchi")],
)
fasta = ">A1\n$s1\n>A2\n$s2\n>A3\n$s3\n>A1b\n$s1\n>A4\n$s4\n>bad\nACGTZ\n"

function scenario(; params = LCAParams(), accept_rank = :phylum, portal_fails = false,
                  project = "proj", user = "me")
    f = FakeBOLD(catalog, portal)
    f.portal_fails = portal_fails
    input = read_asvs(Vector{UInt8}(codeunits(fasta)); name = "asvs.fasta")
    a = run_analysis(input; project, user, params, accept_rank, engine = fake_engine(f),
                     portal_url = "https://portal.test")
    return a, f
end

@testset "Analysis with simulated BOLD services" begin
    a, f = scenario()
    @test count(startswith("POST"), f.requests) == 1     # unique sequences sent once
    @test length(a.input.unique) == 4
    st = Dict(u.ids[1] => a.results[i] for (i, u) in enumerate(a.input.unique))
    @test (st["A1"].taxon, st["A1"].rank) == ("Onthophagus", :genus)   # P3 within 1 point
    @test st["A1"].n_selected == 3 && st["A1"].n_matches == 4
    @test "identical_reference_sequences" in st["A1"].flags         # P1 and P2 share nuc
    @test (st["A2"].taxon, st["A2"].rank) == ("Caccobius schreberi", :species)
    @test (st["A3"].taxon, st["A3"].rank) == ("Aphodius", :genus)      # best 96.5 %: species capped; sp. ignored
    @test st["A3"].lca_rank == :species && st["A3"].rank_capped
    @test st["A4"].status == :no_matches
    # Marker-aware metadata: P1 takes the COI-5P document.
    ref = a.references["P1|COI-5P"]
    @test ref.identified_by == "A. Rossi" && ref.record_id == "P1.COI-5P"
    @test ref.sequence_hash == sequence_hash(nuc(1)) == a.references["P2|COI-5P"].sequence_hash
    @test ref.recordsets == "PROJ; DS-TEST" && ref.institution == "Some Institute"
    @test ismissing(a.references["P5|COI-5P"].identified_by)
    @test isempty(a.errors) && a.source["metadata_status"] == "retrieved"

    rows = assignment_table(a)
    @test [r.asv_id for r in rows] == ["A1", "A2", "A3", "A1b", "A4", "bad"]
    @test rows[1].taxon == rows[4].taxon == "Onthophagus" && rows[1].sequence_hash == rows[4].sequence_hash
    @test rows[1].selected_process_ids == "P1;P2;P3"
    @test rows[6].status == "input_error" && occursin("invalid nucleotide", rows[6].error)
    @test rows[5].status == "no_matches" && !rows[5].accepted
    @test ismissing(rows[1].species) && rows[1].genus == "Onthophagus"

    mt = match_table(a)
    @test length(mt) == 8
    p9 = only(r for r in mt if r.process_id == "P9")
    @test !p9.included_in_lca && p9.decision == "outside_delta" && p9.final_taxon == "Onthophagus"
    p1b = only(r for r in mt if r.process_id == "P1" && r.sequence_hash == sequence_hash(s2))
    @test p1b.decision == "low_coverage" && p1b.coverage == 60
    @test only(r for r in mt if r.process_id == "P1" && r.sequence_hash == sequence_hash(s1)).asv_ids == "A1;A1b"

    use = Dict(r.process_id => r for r in reference_usage(a))
    @test (use["P1"].retrieved_usage, use["P1"].lca_participation, use["P1"].assignment_contribution) == (2, 1, 1)
    @test (use["P9"].retrieved_usage, use["P9"].lca_participation) == (1, 0)
    @test use["P5"].identified_by == "(not provided)"

    attr = Dict(r.identifier => r for r in taxonomic_attribution(a))
    @test Set(keys(attr)) == Set(["A. Rossi", "Rossi, A.", "B. Bianchi", "C. Verdi", "(not provided)"])
    @test (attr["A. Rossi"].n_records, attr["A. Rossi"].n_asvs_supported,
           attr["A. Rossi"].n_lca_participations, attr["A. Rossi"].n_taxa) == (2, 1, 2, 1)
    @test attr["B. Bianchi"].n_records == 2 && attr["B. Bianchi"].n_lca_participations == 1
    @test attr["B. Bianchi"].n_taxa == 2
    # A1: identifiers {A. Rossi, B. Bianchi} → 1/2 each; A2: Rossi, A. → 1;
    # A3: {(not provided), C. Verdi} → 1/2 each.
    @test attr["A. Rossi"].fractional_credit == 0.5 && attr["Rossi, A."].fractional_credit == 1.0
    @test sum(r.fractional_credit for r in values(attr)) ≈ 3
    @test !any(r -> r.identifier in ("Some Institute", "Wrong Marker"), values(attr))
    @test ismissing(attr["C. Verdi"].orcid)

    @test !isdefined(BOLDTrace, :taxonomic_contribution_value)     # monetary value removed
    # Accepted assignments from the species only: one.
    a_sp, _ = scenario(; accept_rank = :species)
    @test [r.asv_id for r in assignment_table(a_sp) if r.accepted] == ["A2"]
    use_sp = Dict(r.process_id => r for r in reference_usage(a_sp))
    @test use_sp["P1"].lca_participation == 1 && use_sp["P1"].assignment_contribution == 0
    @test use_sp["P4"].assignment_contribution == 1

    s = analysis_summary(a)
    @test s.n_accepted == 3 && s.status["no_matches"] == 1 && s.n_input_errors == 1
    @test s.n_references == 7 && s.n_references_in_lca == 6

    # Portal failure: the analysis completes, the limitation is reported.
    b, _ = scenario(; portal_fails = true)
    @test b.source["metadata_status"] == "failed" && length(b.errors) == 1
    @test ismissing(b.references["P1|COI-5P"].identified_by)
    @test assignment_table(b)[1].taxon == "Onthophagus"
    @test any(occursin("metadata from the BOLD portal are missing", l) for l in provenance(b)["limitations"])
end

@testset "Minimum agreement" begin
    # 99 records of O. taurus and one O. illyricus, all at 100 %.
    many = [m("t$(lpad(k, 3, '0'))", 100, lineage()) for k in 1:99]
    push!(many, m("x", 100, lineage(species = "Onthophagus illyricus")))
    strict = lca(many, LCAParams(); valid_name = valid)
    @test strict.rank == :genus && all(isnothing, strict.dissent)
    r = lca(many, LCAParams(; min_agreement = 95); valid_name = valid)
    @test (r.taxon, r.rank) == ("Onthophagus taurus", :species)
    @test r.agreement == [100, 100, 100, 100, 100, 99]
    @test "majority:species" in r.flags && "discordant:species" in r.flags
    @test r.dissent[end] == :species && count(!isnothing, r.dissent) == 1
    # 90/10 is below 95 %; a tie never resolves.
    ten = [[m("t$k", 100, lineage()) for k in 1:9]; m("x", 100, lineage(species = "Onthophagus illyricus"))]
    @test lca(ten, LCAParams(; min_agreement = 95); valid_name = valid).rank == :genus
    @test lca(ten, LCAParams(; min_agreement = 90); valid_name = valid).rank == :species
    tie = [m("a", 100, lineage()), m("b", 100, lineage(species = "Onthophagus illyricus"))]
    @test lca(tie, LCAParams(; min_agreement = 51); valid_name = valid).rank == :genus
    # A dissenter at the family leaves the deeper ranks.
    fam = [[m("t$k", 100, lineage()) for k in 1:19];
           m("g", 100, lineage(family = "Geotrupidae", genus = "Geotrupes", species = "Geotrupes stercorarius"))]
    r = lca(fam, LCAParams(; min_agreement = 95); valid_name = valid)
    @test r.rank == :species && r.dissent[end] == :family && r.agreement[4] == 95
    # Unidentified references do not count in the agreement.
    inc = [m("a", 100, lineage()), m("b", 100, lineage(species = nothing)),
           m("c", 100, lineage(species = "Onthophagus sp.")),
           m("d", 100, lineage(species = "Onthophagus illyricus"))]
    @test lca(inc, LCAParams(; min_agreement = 60); valid_name = valid).rank == :genus   # 1 of 2
    @test_throws ArgumentError LCAParams(; min_agreement = 50)
    @test_throws ArgumentError LCAParams(; min_agreement = 101)

    # In an analysis: A1 has P1, P2 (O. taurus) and P3 (O. illyricus): 2 of 3.
    a, _ = scenario(; params = LCAParams(; min_agreement = 60))
    rows = Dict(r.asv_id => r for r in assignment_table(a))
    @test rows["A1"].taxon == "Onthophagus taurus" && rows["A1"].agreement ≈ 66.67
    @test rows["A1"].n_dissenting == 1 && occursin("majority:species", rows["A1"].flags)
    p3 = only(r for r in match_table(a) if r.process_id == "P3")
    @test p3.included_in_lca && p3.dissent_rank == "species" && !p3.supports_assignment
    use = Dict(r.process_id => r for r in reference_usage(a))
    @test use["P3"].lca_participation == 1 && use["P3"].assignment_contribution == 0
    attr = Dict(r.identifier => r for r in taxonomic_attribution(a))
    @test attr["B. Bianchi"].n_asvs_supported == 0 && attr["B. Bianchi"].fractional_credit == 0
    @test attr["A. Rossi"].fractional_credit == 1.0
    @test provenance(a)["lca"]["parameters"]["min_agreement"] == 60
    store = open_store()
    record_analysis!(store, a)
    @test usage_stats(store).n_contribution_events == 5 && usage_stats(store).n_lca_events == 6
    close(store)
    # A dissenter below the assigned rank (after the rank cap) still supports it.
    capped = [m("a", 96, lineage()), m("b", 96, lineage()), m("c", 96, lineage(species = "Onthophagus illyricus"))]
    rc = lca(capped, LCAParams(; min_agreement = 60); valid_name = valid)
    @test rc.rank == :genus && rc.rank_capped && rc.dissent[3] == :species
end

@testset "Provenance and exports" begin
    a, _ = scenario()
    p = provenance(a)
    @test p["analysis_id"] == a.id && p["fingerprint"] == a.fingerprint
    @test p["input"]["sequence_set_hash"] == a.input.set_hash
    @test p["bold_data_source"]["database_code"] == "public.tax-derep"
    @test p["bold_data_source"]["mode_parameters"]["mi"] == 0.94
    @test occursin("not provided", p["database_version"])
    @test p["lca"]["version"] == LCA_VERSION && p["lca"]["parameters"]["max_delta"] == 1.0
    @test p["selected_process_ids"][sequence_hash(s1)] == ["P1", "P2", "P3"]
    @test length(p["assignments"]) == 4
    dir = mktempdir()
    links = [CitationLink("P1", "10.1234/abcd.5", "original_barcode_publication", "manual")]
    paths = export_analysis(a, dir; citations = links)
    @test basename.(paths) == collect(EXPORT_FILES)
    js = JSON3.read(read(joinpath(dir, "analysis_provenance.json"), String))
    @test js[:analysis_id] == a.id
    lines = readlines(joinpath(dir, "assignments.csv"))
    @test startswith(lines[1], "asv_id,sequence_hash,status,taxon,rank,phylum")
    @test length(lines) == 7
    tsv = readlines(joinpath(dir, "used_process_ids.tsv"))
    @test split(tsv[1], '\t')[1:4] == ["process_id", "record_id", "marker", "lca_participation"]
    @test length(tsv) == 1 + 6
    p1 = split(only(filter(l -> startswith(l, "P1\t"), tsv)), '\t')
    @test "10.1234/abcd.5" in p1 && "original_barcode_publication" in p1
    @test length(readlines(joinpath(dir, "bold_matches.csv"))) == 9
    @test_throws ArgumentError export_file(a, "other.csv")
    @test_throws ArgumentError CitationLink("P1", "not-a-doi", "taxonomic_revision", "x")
    @test_throws ArgumentError CitationLink("P1", "10.1234/x", "cited", "x")
    # CSV quoting and formula injection.
    io = IOBuffer()
    BOLDTrace.write_table(io, [(a = "x,y", b = "=SUM(A1)", c = -1.5, d = missing, e = "say \"hi\"")])
    @test String(take!(io)) == "a,b,c,d,e\n\"x,y\",'=SUM(A1),-1.5,,\"say \"\"hi\"\"\"\n"
    io = IOBuffer()
    BOLDTrace.write_table(io, NamedTuple[]; columns = (:a, :b))
    @test String(take!(io)) == "a,b\n"
end

@testset "Usage database" begin
    path = joinpath(mktempdir(), "usage.sqlite")
    store = open_store(path)
    a, _ = scenario()
    r1 = record_analysis!(store, a)
    @test r1.recorded && r1.n_events == 8
    # The same analysis again (new run, same input/parameters/project): not recorded.
    a2, _ = scenario()
    @test a2.id != a.id && a2.fingerprint == a.fingerprint
    r2 = record_analysis!(store, a2)
    @test !r2.recorded && r2.analysis_id == a.id
    g = usage_stats(store)
    @test g.scope == "global" && g.n_analyses == 1 && g.n_events == 8
    @test g.n_lca_events == 6 && g.n_contribution_events == 6
    @test g.n_asvs == 3 && g.n_process_ids == 7 && g.n_process_ids_in_lca == 6
    @test g.n_reference_sequences == 6          # P1 and P2 share one sequence
    @test g.note == GLOBAL_STATS_NOTE
    # Another project with the same input is a separate analysis.
    b, _ = scenario(; project = "other", user = "you", params = LCAParams(; max_delta = 2))
    @test record_analysis!(store, b).recorded
    g = usage_stats(store)
    @test g.n_analyses == 2 && g.n_projects == 2 && g.n_events == 16 && g.n_asvs == 3
    mine = usage_stats(store; project = "proj")
    @test mine.scope == "personal" && mine.n_events == 8 && mine.n_projects == 1
    @test usage_stats(store; user = "you").n_analyses == 1
    @test usage_stats(store; analysis_id = a.id).n_events == 8
    @test usage_stats(store; project = "none").n_events == 0
    @test length(list_analyses(store)) == 2 && length(list_analyses(store; project = "other")) == 1
    add_citation!(store, CitationLink("P1", "10.1234/abcd.5", "taxonomic_revision", "paper"))
    add_citation!(store, CitationLink("P1", "10.1234/abcd.5", "taxonomic_revision", "paper", "verified"))
    add_citation!(store, CitationLink("P7", "10.5555/xyz", "barcode_reuse_publication", "x"))
    links = citation_links(store, ["P1"])
    @test length(links) == 1 && links[1].verification_status == "verified"
    @test length(citation_links(store)) == 2
    close(store)
    # Reopening keeps everything.
    store = open_store(path)
    @test usage_stats(store).n_events == 16
    close(store)
end

@testset "Web application" begin
    f = FakeBOLD(catalog, portal)
    dir = mktempdir()
    port = 18000 + rand(0:999)
    server = serve(; database = joinpath(dir, "web.sqlite"), port, async = true,
                   engine = () -> fake_engine(f), checkpoint = joinpath(dir, "ckpt"),
                   portal_url = "https://portal.test", max_body = 200_000)
    base = "http://127.0.0.1:$port"
    req(method, path; body = nothing, type = "application/json") =
        HTTP.request(method, base * path, body === nothing ? [] : ["Content-Type" => type],
                     body === nothing ? UInt8[] : body; status_exception = false, retry = false)
    jsonof(r) = JSON3.read(String(r.body))
    try
        r = req("GET", "/")
        @test r.status == 200 && occursin("BOLDTrace", String(r.body))
        @test occursin("default-src 'self'", HTTP.header(r, "Content-Security-Policy"))
        @test req("GET", "/assets/app.js").status == 200
        logo = req("GET", "/assets/dnambientale.png")
        @test logo.status == 200 && HTTP.header(logo, "Content-Type") == "image/png"
        @test logo.body == read(joinpath(pkgdir(BOLDTrace), "src", "assets", "dnambientale.png"))
        @test occursin("Powered by", String(req("GET", "/").body))
        o = jsonof(req("GET", "/api/options"))
        @test o[:defaults][:min_identity] == 90 && !haskey(o, :tcv_disclaimer)
        # Upload.
        @test req("POST", "/api/upload"; body = "x", type = "text/plain").status == 415
        @test req("POST", "/api/upload"; body = JSON3.write((fasta = "garbage",))).status == 400
        @test req("POST", "/api/upload"; body = JSON3.write((fasta = repeat("A", 300_000),))).status == 413
        u = jsonof(req("POST", "/api/upload"; body = JSON3.write((fasta = fasta, filename = "a.fasta"))))
        @test u[:n_records] == 6 && u[:n_unique] == 4 && u[:n_errors] == 1
        # Analysis.
        bad = (token = u[:token], params = (min_identity = 150,))
        @test req("POST", "/api/analyze"; body = JSON3.write(bad)).status == 400
        @test req("POST", "/api/analyze"; body = JSON3.write((token = "nope",))).status == 404
        j = jsonof(req("POST", "/api/analyze"; body = JSON3.write((token = u[:token], project = "web",
                                                                    params = (max_delta = 1,)))))
        id = j[:id]
        status = ""
        for _ in 1:200
            status = jsonof(req("GET", "/api/job?id=$id"))[:status]
            status in ("done", "failed") && break
            sleep(0.1)
        end
        job = jsonof(req("GET", "/api/job?id=$id"))
        @test status == "done"
        @test job[:summary][:n_accepted] == 3 && job[:recorded][:recorded] == true
        @test req("GET", "/api/job?id=unknown").status == 404
        as = jsonof(req("GET", "/api/assignments?id=$id"))
        @test length(as) == 6 && as[1][:taxon] == "Onthophagus"
        ms = jsonof(req("GET", "/api/matches?id=$id&hash=$(as[1][:sequence_hash])"))
        @test length(ms) == 4 && count(x -> x[:included_in_lca], ms) == 3
        @test req("GET", "/api/matches?id=$id&hash=zzz").status == 404
        @test length(jsonof(req("GET", "/api/references?id=$id"))) == 7
        @test length(jsonof(req("GET", "/api/attribution?id=$id"))) == 5
        @test req("GET", "/api/tcv?id=$id").status == 404
        # Hits of one ASV: same columns, only its rows, file named after its first ID.
        one = req("GET", "/api/export?id=$id&file=bold_matches.csv&hash=$(as[1][:sequence_hash])")
        @test one.status == 200
        @test occursin("bold_matches_A1.csv", HTTP.header(one, "Content-Disposition"))
        lines = split(strip(String(one.body)), '
')
        @test length(lines) == 5 && length(split(lines[1], ',')) == length(BOLDTrace.match_columns(BOLDTrace.LCA.DEFAULT_RANKS))
        @test req("GET", "/api/export?id=$id&file=bold_matches.csv&hash=zzz").status == 404
        @test req("GET", "/api/export?id=$id&file=assignments.csv&hash=$(as[1][:sequence_hash])").status == 400
        # Rank thresholds from the page: genus at 97 % moves A3 (best 96.5 %) to the family.
        @test req("POST", "/api/analyze"; body = JSON3.write((token = u[:token], params = (min_agreement = 50,)))).status == 400
        bad = (token = u[:token], params = (rank_thresholds = (tribe = 90,),))
        @test req("POST", "/api/analyze"; body = JSON3.write(bad)).status == 400
        bad = (token = u[:token], params = (rank_thresholds = (genus = 101,),))
        @test req("POST", "/api/analyze"; body = JSON3.write(bad)).status == 400
        jt = jsonof(req("POST", "/api/analyze"; body = JSON3.write((token = u[:token], project = "thr",
                                  params = (rank_thresholds = (genus = 97,),)))))
        for _ in 1:200
            jsonof(req("GET", "/api/job?id=$(jt[:id])"))[:status] in ("done", "failed") && break
            sleep(0.1)
        end
        a3 = only(r for r in jsonof(req("GET", "/api/assignments?id=$(jt[:id])")) if r[:asv_id] == "A3")
        @test a3[:rank] == "family" && a3[:rank_capped] == true
        for name in EXPORT_FILES
            e = req("GET", "/api/export?id=$id&file=$name")
            @test e.status == 200 && occursin(name, HTTP.header(e, "Content-Disposition"))
        end
        @test req("GET", "/api/export?id=$id&file=..%2Fsecret").status == 404
        stats = jsonof(req("GET", "/api/stats?project=web"))
        @test stats[:personal][:n_events] == 8 && stats[:global_][:n_events] == 16   # + "thr"
        # Running the same analysis again does not increase the counters.
        j2 = jsonof(req("POST", "/api/analyze"; body = JSON3.write((token = u[:token], project = "web",
                                                                     params = (max_delta = 1,)))))
        for _ in 1:200
            jsonof(req("GET", "/api/job?id=$(j2[:id])"))[:status] in ("done", "failed") && break
            sleep(0.1)
        end
        @test jsonof(req("GET", "/api/job?id=$(j2[:id])"))[:recorded][:recorded] == false
        @test jsonof(req("GET", "/api/stats"))[:global_][:n_events] == 16
        # Checkpoint reused: the second run sent no new submission.
        @test count(startswith("POST"), f.requests) == 1
        @test req("GET", "/missing").status == 404
    finally
        shutdown(server)
    end
end

@testset "Review regressions BT1-BT5" begin
    # BT5: coverage and length outside their domain are rejected.
    for cov in (NaN, Inf, -Inf, 101.0, -1.0)
        @test_throws ArgumentError lca([Match("a", 99, cov, 300, lineage())])
    end
    @test_throws ArgumentError lca([Match("a", 99, 100.0, -1, lineage())])
    @test lca([Match("a", 99, 0.0, 0, lineage())], LCAParams(; min_coverage = 0)).status == :assigned

    # BT4: same export schema with and without matches.
    nohit_fasta = ">N1\n$(dna(71))\n>N2\n$(dna(72))\n"
    f0 = FakeBOLD(Dict{String,Vector{Any}}(), Dict{String,Vector{Dict{String,Any}}}())
    empty_input = read_asvs(Vector{UInt8}(codeunits(nohit_fasta)); name = "nohit.fasta")
    e = run_analysis(empty_input; project = "empty", engine = fake_engine(f0),
                     portal_url = "https://portal.test")
    @test all(r -> r.status == :no_matches, e.results)
    a, _ = scenario()
    for name in EXPORT_FILES
        endswith(name, ".json") && continue
        h(x) = first(split(export_file(x, name), '\n'))
        @test h(e) == h(a)
    end
    @test length(split(first(split(export_file(e, "bold_matches.csv"), '\n')), ',')) == 26
    @test length(split(first(split(export_file(e, "reference_usage.csv"), '\n')), ',')) == 12
    @test length(split(first(split(export_file(e, "taxonomic_attribution.csv"), '\n')), ',')) == 7
    @test length(split(first(split(export_file(e, "used_process_ids.tsv"), '\n')), '\t')) == 9
    # Rows carry exactly the declared columns.
    ranks = a.params.ranks
    @test collect(keys(first(assignment_table(a)))) == collect(BOLDTrace.assignment_columns(ranks))
    @test collect(keys(first(match_table(a)))) == collect(BOLDTrace.match_columns(ranks))
    @test collect(keys(first(reference_usage(a)))) == collect(BOLDTrace.REFERENCE_USAGE_COLUMNS)
    @test collect(keys(first(taxonomic_attribution(a)))) == collect(BOLDTrace.ATTRIBUTION_COLUMNS)
    @test collect(keys(first(used_process_ids(a)))) == collect(BOLDTrace.USED_PROCESS_ID_COLUMNS)

    # BT3: analyses without matches count as analyses and projects.
    store = open_store()
    @test record_analysis!(store, e).recorded
    g = usage_stats(store)
    @test (g.n_analyses, g.n_projects, g.n_events) == (1, 1, 0)
    @test usage_stats(store; project = "empty").n_analyses == 1
    @test usage_stats(store; project = "other").n_analyses == 0
    close(store)

    # Server with an engine factory held by a gate: the first job stays running.
    function gated_server(dir; kwargs...)
        f = FakeBOLD(catalog, portal)
        gate = Channel{Nothing}(Inf)
        port = 19000 + rand(0:999)
        server = serve(; database = joinpath(dir, "gate.sqlite"), port, async = true,
                       engine = () -> (take!(gate); fake_engine(f)), checkpoint = nothing,
                       portal_url = "https://portal.test", kwargs...)
        return server, gate, "http://127.0.0.1:$port"
    end
    req(base, method, path; body = nothing) =
        HTTP.request(method, base * path, body === nothing ? [] : ["Content-Type" => "application/json"],
                     body === nothing ? UInt8[] : body; status_exception = false, retry = false,
                     readtimeout = 120)
    jsonof(r) = JSON3.read(String(r.body))
    function start_job(base, token, project)
        r = req(base, "POST", "/api/analyze"; body = JSON3.write((token = token, project = project)))
        return r.status == 200 ? jsonof(r)[:id] : r.status
    end
    status(base, id) = (r = req(base, "GET", "/api/job?id=$id"); r.status == 200 ? jsonof(r)[:status] : r.status)
    function wait_status(base, id, wanted; tries = 300)
        s = nothing
        for _ in 1:tries
            s = status(base, id)
            s in wanted && return s
            sleep(0.1)
        end
        return s
    end

    # BT2: running and queued jobs are not evicted; the queue is bounded.
    dir = mktempdir()
    server, gate, base = gated_server(dir; max_kept = 1, max_pending = 2)
    try
        token = jsonof(req(base, "POST", "/api/upload"; body = JSON3.write((fasta = fasta,))))[:token]
        id1 = start_job(base, token, "p1")
        @test wait_status(base, id1, ("running",)) == "running"
        id2 = start_job(base, token, "p2")
        @test status(base, id1) == "running"          # not evicted by the second job
        @test status(base, id2) == "queued"
        @test start_job(base, token, "p3") == 503     # max_pending reached
        put!(gate, nothing); put!(gate, nothing)
        @test wait_status(base, id1, ("done", "failed")) == "done"
        @test wait_status(base, id2, ("done", "failed")) == "done"
        # Finished jobs are evicted beyond max_kept.
        put!(gate, nothing)
        id3 = start_job(base, token, "p3")
        @test wait_status(base, id3, ("done", "failed")) == "done"
        @test status(base, id1) == 404
    finally
        put!(gate, nothing); put!(gate, nothing)
        @test shutdown(server)
    end

    # BT1: shutdown waits for the running job before closing the database;
    # queued jobs are cancelled.
    dir = mktempdir()
    server, gate, base = gated_server(dir)
    app = BOLDTrace._APP[]
    token = jsonof(req(base, "POST", "/api/upload"; body = JSON3.write((fasta = fasta,))))[:token]
    id1 = start_job(base, token, "q1")
    @test wait_status(base, id1, ("running",)) == "running"
    id2 = start_job(base, token, "q2")
    @test shutdown(server; timeout = 0.5) == false     # job still running
    @test !istaskdone(app.worker) && isopen(app.store.db)
    put!(gate, nothing)
    wait(app.worker)
    @test app.jobs[id1].status == :done && app.jobs[id1].recorded.recorded
    @test app.jobs[id2].status == :cancelled
    @test !isopen(app.store.db)
    store = open_store(joinpath(dir, "gate.sqlite"))
    @test usage_stats(store).n_analyses == 1 && usage_stats(store; project = "q1").n_events == 8
    close(store)
    @test BOLDTrace._APP[] === nothing || BOLDTrace._APP[].stopping[]
end

end
