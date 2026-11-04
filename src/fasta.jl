# Reading and validation of the multi-FASTA input.

"""
    InputRecord

One record of the input FASTA: original `id` (first word of the header),
full `header`, normalized `sequence` and its `hash` (SHA-256 of the normalized
sequence), `unique_index` (position in `InputSet.unique`), and `error` when the
record is not usable (`sequence`, `hash` and `unique_index` are then `nothing`).
"""
struct InputRecord
    id::String
    header::String
    sequence::Union{Nothing,String}
    hash::Union{Nothing,String}
    unique_index::Union{Nothing,Int}
    error::Union{Nothing,String}
end

"""
    UniqueSequence

A distinct normalized sequence of the input, its `hash` and the `ids` of every
input record carrying it, in input order.
"""
struct UniqueSequence
    hash::String
    sequence::String
    ids::Vector{String}
end

"""
    InputSet

Validated input: `records` in input order, `unique` sequences in order of first
appearance, `messages` (warnings that do not make a record unusable),
`file_sha256` of the raw bytes, `set_hash` (SHA-256 of the sorted unique
sequence hashes, independent of record order, IDs and line wrapping) and the
source `name`.
"""
struct InputSet
    name::String
    records::Vector{InputRecord}
    unique::Vector{UniqueSequence}
    messages::Vector{String}
    file_sha256::String
    set_hash::String
end

"""
    sequence_hash(seq) -> String

Stable hash of a nucleotide sequence: hexadecimal SHA-256 of the sequence
normalized with `BOLDassignment.normalize_sequence` (no whitespace or gaps,
upper case, `U` → `T`). Identical sequences always get the same hash,
whatever their case, gaps or line wrapping.
"""
sequence_hash(seq::AbstractString) = _hash_normalized(normalize_sequence(seq))
_hash_normalized(s::AbstractString) = bytes2hex(SHA.sha256(s))

"""
    read_asvs(path_or_io; name = basename(path), min_length = 100,
              max_records = 100_000) -> InputSet

Read and validate a (multi-)FASTA file with FASTX.

Every record keeps its original identifier. A record is rejected, with an
error message and without stopping the others, when its sequence is empty or
contains characters outside the IUPAC nucleotide alphabet, or when its
identifier repeats an earlier one with a different sequence. Warnings, which
keep the record, are given for sequences shorter than `min_length` (the BOLD
engine requires an overlap of at least 100 bases), repeated identifiers with
the same sequence and sequences that are more than 5 % ambiguous bases.
Identical sequences are recognized after normalization and grouped in
`unique`. A file that is not valid FASTA throws `ArgumentError`.
"""
function read_asvs(path::AbstractString; name::AbstractString = basename(path), kwargs...)
    isfile(path) || throw(ArgumentError("file not found: $path"))
    return read_asvs(read(path); name, kwargs...)
end

read_asvs(io::IO; name::AbstractString = "input.fasta", kwargs...) =
    read_asvs(read(io); name, kwargs...)

function read_asvs(bytes::AbstractVector{UInt8}; name::AbstractString = "input.fasta",
                   min_length::Integer = 100, max_records::Integer = 100_000)
    data = Vector{UInt8}(bytes)
    file_sha = bytes2hex(SHA.sha256(data))
    records = InputRecord[]
    uniques = UniqueSequence[]
    by_hash = Dict{String,Int}()
    seen_ids = Dict{String,Union{Nothing,String}}()  # id => hash of its first record
    messages = String[]
    isempty(strip(String(copy(data)))) && throw(ArgumentError("the FASTA file is empty"))
    reader = FASTA.Reader(IOBuffer(data))
    try
        n = 0
        for rec in reader
            n += 1
            n > max_records && throw(ArgumentError("more than $max_records records"))
            id = String(FASTA.identifier(rec))
            header = String(FASTA.description(rec))
            isempty(id) && (id = "record_$n")
            raw = String(FASTA.sequence(rec))
            seq = try
                normalize_sequence(raw)
            catch err
                err isa ArgumentError || rethrow()
                push!(records, InputRecord(id, header, nothing, nothing, nothing,
                                           "record $n ($id): $(err.msg)"))
                haskey(seen_ids, id) || (seen_ids[id] = nothing)
                continue
            end
            h = _hash_normalized(seq)
            if haskey(seen_ids, id)
                if seen_ids[id] == h
                    push!(messages, "record $n ($id): identifier repeated with the same sequence")
                else
                    push!(records, InputRecord(id, header, nothing, nothing, nothing,
                                               "record $n ($id): identifier already used for a different sequence"))
                    continue
                end
            else
                seen_ids[id] = h
            end
            length(seq) < min_length &&
                push!(messages, "record $n ($id): $(length(seq)) bases, shorter than $min_length")
            ambiguous = count(c -> !(c in "ACGT"), seq)
            ambiguous > 0.05 * length(seq) &&
                push!(messages, "record $n ($id): $ambiguous ambiguous bases")
            k = get(by_hash, h, 0)
            if k == 0
                push!(uniques, UniqueSequence(h, seq, String[]))
                k = by_hash[h] = length(uniques)
            end
            push!(uniques[k].ids, id)
            push!(records, InputRecord(id, header, seq, h, k, nothing))
        end
    catch err
        err isa ArgumentError && rethrow()
        throw(ArgumentError("not a valid FASTA file: $(sprint(showerror, err))"))
    finally
        close(reader)
    end
    isempty(records) && throw(ArgumentError("no FASTA record found"))
    set_hash = bytes2hex(SHA.sha256(join(sort!([u.hash for u in uniques]), '\n')))
    return InputSet(String(name), records, uniques, messages, file_sha, set_hash)
end

"""
    input_summary(input::InputSet) -> NamedTuple

Counts for the upload page: records, valid records, unique sequences,
duplicated records (valid records whose sequence occurs earlier), errors and
warnings.
"""
function input_summary(input::InputSet)
    valid = count(r -> r.error === nothing, input.records)
    errors = [r.error for r in input.records if r.error !== nothing]
    return (name = input.name, n_records = length(input.records), n_valid = valid,
            n_unique = length(input.unique), n_duplicates = valid - length(input.unique),
            n_errors = length(errors), errors = errors, warnings = copy(input.messages),
            file_sha256 = input.file_sha256, set_hash = input.set_hash)
end
