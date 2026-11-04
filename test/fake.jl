# Simulated BOLD services in the formats observed on 2026-10-23 (identification
# engine) and 2026-11-04 (portal documents, one per specimen and marker).

mutable struct FakeBOLD
    catalog::Dict{String,Vector{Any}}      # sequence => [(pid, marker, pdist, qcov, len, taxonomy)]
    portal::Dict{String,Vector{Dict{String,Any}}}   # pid => documents
    submissions::Dict{String,Vector{Pair{String,String}}}
    portal_fails::Bool
    requests::Vector{String}
    omit::Set{String}          # sequences left out of finished submissions
end
FakeBOLD(catalog, portal) = FakeBOLD(catalog, portal, Dict(), false, String[], Set{String}())

_bytes(s) = Vector{UInt8}(codeunits(s))

function fake_request(f::FakeBOLD, method, url, payload)
    push!(f.requests, "$method $url")
    if startswith(url, "https://portal.test")
        f.portal_fails && return (503, _bytes("unavailable"))
        if occursin("/api/query?", url)
            q = HTTP.unescapeuri(match(r"query=([^&]*)", url)[1])
            ids = [last(split(t, ':')) for t in split(q, ';')]
            return (200, _bytes(JSON3.write(Dict("query_id" => join(ids, ',')))))
        end
        m = match(r"/api/documents/([^?]*)\?length=(\d+)&start=(\d+)", url)
        ids = split(HTTP.unescapeuri(m[1]), ',')
        docs = reduce(vcat, [get(f.portal, i, Dict{String,Any}[]) for i in ids]; init = Dict{String,Any}[])
        n, start = parse(Int, m[2]), parse(Int, m[3])
        return (200, _bytes(JSON3.write(Dict("data" => docs[start+1:min(end, start + n)]))))
    end
    if method == "POST"
        records = Pair{String,String}[]
        for rec in FASTA.Reader(IOBuffer(payload))
            push!(records, String(FASTA.identifier(rec)) => String(FASTA.sequence(rec)))
        end
        sub = "sub$(length(f.submissions) + 1)"
        f.submissions[sub] = records
        return (200, _bytes("{\"sub_id\":\"$sub\"}"))
    end
    sub = last(split(url, '/'))
    haskey(f.submissions, sub) || return (404, UInt8[])
    occursin("/status/", url) && return (200, _bytes("{\"status\":\"complete\"}"))
    io = IOBuffer()
    for (id, seq) in f.submissions[sub]
        seq in f.omit && continue
        results = Dict{String,Any}()
        for (pid, marker, pdist, qcov, len, tax) in get(f.catalog, seq, [])
            results["$pid|$marker|BOLD:AAA0001|1|public|0"] =
                Dict("pdist" => pdist, "pident" => 100 - pdist, "qcovhsp" => qcov, "length" => len,
                     "taxonomy" => tax)
        end
        println(io, JSON3.write(Dict("seqid" => id, "sequence" => seq, "results" => results)))
    end
    return (200, take!(io))
end

fake_engine(f::FakeBOLD) =
    IDEngine(; base_url = "https://id.test", request = (m, u, p) -> fake_request(f, m, u, p),
             sleep = dt -> nothing, verbose = false, submit_interval = 0.0, request_interval = 0.0,
             status_wait = 0.0, retry_wait = 0.0, max_retries = 2,
             throttle = BOLDassignment._Throttle())

tax(genus, species; family = "Scarabaeidae") =
    Dict("phylum" => "Arthropoda", "class" => "Insecta", "order" => "Coleoptera",
         "family" => family, "genus" => genus, "species" => species)

doc(pid, marker, nuc, identified_by; species = "x") =
    Dict{String,Any}("processid" => pid, "marker_code" => marker, "record_id" => "$pid.$marker",
                     "nuc" => nuc, "identified_by" => identified_by, "inst" => "Some Institute",
                     "bold_recordset_code_arr" => ["PROJ", "DS-TEST"], "species" => species)
