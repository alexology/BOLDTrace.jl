# Web application built on Genie.jl. Page and static files come from a fixed
# list; the API answers JSON. Analyses run one at a time in a background task,
# because the BOLD rate limits apply to the whole client.

import Genie
import Sockets

const _ASSETS = joinpath(@__DIR__, "assets")
const _STATIC = Dict(
    "/" => ("index.html", "text/html; charset=utf-8"),
    "/assets/app.js" => ("app.js", "text/javascript; charset=utf-8"),
    "/assets/style.css" => ("style.css", "text/css; charset=utf-8"),
    "/assets/dnambientale.png" => ("dnambientale.png", "image/png"),
)
for (_, (file, _)) in _STATIC
    include_dependency(joinpath(_ASSETS, file))
end
const _STATIC_BODY = Dict(path => read(joinpath(_ASSETS, file)) for (path, (file, _)) in _STATIC)

const _SECURITY_HEADERS = [
    "Content-Security-Policy" => "default-src 'self'; img-src 'self' data:; style-src 'self'; " *
                                 "script-src 'self'; connect-src 'self'; base-uri 'none'; " *
                                 "form-action 'none'; frame-ancestors 'none'",
    "X-Content-Type-Options" => "nosniff",
    "X-Frame-Options" => "DENY",
    "Referrer-Policy" => "no-referrer",
    "Cross-Origin-Opener-Policy" => "same-origin",
    "Cross-Origin-Resource-Policy" => "same-origin",
]

const _STATUS_TEXT = Dict(400 => "bad request", 404 => "not found", 405 => "method not allowed",
                          413 => "request too large", 415 => "unsupported media type",
                          421 => "unknown host", 500 => "internal error")

struct _RequestError <: Exception
    status::Int
    msg::String
end

_response(status, body, type; extra = Pair{String,String}[], cache = "no-store") =
    HTTP.Response(status, vcat(["Content-Type" => type, "Cache-Control" => cache], extra), body)
_json(status, x) = _response(status, JSON3.write(x), "application/json; charset=utf-8")
_error(status, msg) = _json(status, (error = msg,))

function _harden!(res::HTTP.Response)
    if res.status >= 400 && !startswith(HTTP.header(res, "Content-Type"), "application/json")
        res.body = Vector{UInt8}(JSON3.write((error = get(_STATUS_TEXT, res.status, "error"),)))
        HTTP.setheader(res, "Content-Type" => "application/json; charset=utf-8")
        HTTP.setheader(res, "Cache-Control" => "no-store")
    end
    for h in _SECURITY_HEADERS
        HTTP.setheader(res, h)
    end
    HTTP.setheader(res, "Content-Length" => string(length(res.body)))
    return res
end

mutable struct _Job
    id::String
    status::Symbol              # :queued, :running, :done, :failed
    messages::Vector{String}
    error::Union{Nothing,String}
    analysis::Union{Nothing,Analysis}
    recorded::Union{Nothing,NamedTuple}
    request::Dict{String,Any}
    input::InputSet
end

_finished(job::_Job) = job.status in (:done, :failed, :cancelled)

# The worker task owns the store: it closes it when it stops, after the last
# job, so the database is never closed under a running analysis.
mutable struct _App
    store::UsageStore
    engine::Function            # () -> IDEngine
    checkpoint::Union{Nothing,String}
    portal_url::String
    uploads::Dict{String,InputSet}
    upload_order::Vector{String}
    jobs::Dict{String,_Job}
    job_order::Vector{String}
    queue::Channel{_Job}
    lock::ReentrantLock
    max_records::Int
    max_kept::Int
    max_pending::Int
    defaults::LCAParams
    accept_rank::Symbol
    stopping::Threads.Atomic{Bool}
    worker::Union{Nothing,Task}
    listener::Any
end

const _APP = Ref{Union{Nothing,_App}}(nothing)

# Runs the queued jobs one at a time. After `shutdown` the jobs still queued
# are cancelled; the running one completes and is recorded. The store is
# closed when the loop ends.
function _worker_loop(app::_App)
    try
        for job in app.queue
            if app.stopping[]
                lock(app.lock) do
                    job.status = :cancelled
                    push!(job.messages, "Cancelled: the server was stopped")
                end
                continue
            end
            _run_job!(app, job)
        end
    finally
        close(app.store)
    end
end

_query() = HTTP.queryparams(HTTP.URI(Genie.Requests.request().target))
_payload() = JSON3.read(String(copy(Genie.Requests.request().body)))

function _is_json(req)
    t = lowercase(HTTP.header(req, "Content-Type"))
    return startswith(strip(t), "application/json")
end

function _keep!(dict, order, key, value, max_kept)
    dict[key] = value
    push!(order, key)
    while length(order) > max_kept
        delete!(dict, popfirst!(order))
    end
end

# Jobs are evicted only once finished (oldest first): queued and running jobs
# stay reachable. Call with `app.lock` held.
function _keep_job!(app::_App, job::_Job)
    app.jobs[job.id] = job
    push!(app.job_order, job.id)
    excess = length(app.job_order) - app.max_kept
    excess > 0 || return
    for id in copy(app.job_order)
        excess > 0 || break
        if _finished(app.jobs[id])
            delete!(app.jobs, id)
            filter!(!=(id), app.job_order)
            excess -= 1
        end
    end
end

_n_pending(app::_App) = count(j -> j.status in (:queued, :running), values(app.jobs))

function _get_job(app)
    id = get(_query(), "id", "")
    job = lock(() -> get(app.jobs, id, nothing), app.lock)
    job === nothing && throw(_RequestError(404, "unknown analysis"))
    return job
end

function _done_analysis(app)
    job = _get_job(app)
    job.status === :done || throw(_RequestError(409, "analysis not finished"))
    return job, job.analysis
end

_num(x, name, lo, hi) = (x isa Real && !(x isa Bool) && isfinite(x) && lo <= x <= hi) ? x :
    throw(_RequestError(400, "$name must be a number within [$lo, $hi]"))

# LCA parameters of a request: the server defaults, overridden by the fields
# present in `params`.
function _params_from(p, d::LCAParams)
    kw = Dict{Symbol,Any}(:min_identity => d.min_identity, :min_coverage => d.min_coverage,
                          :max_delta => d.max_delta, :max_matches => d.max_matches,
                          :min_align_length => d.min_align_length,
                          :min_agreement => d.min_agreement, :rank_cap => d.rank_cap,
                          :rank_thresholds => copy(d.rank_thresholds), :ranks => d.ranks)
    q = get(p, :params, nothing)
    q === nothing && return LCAParams(; kw...)
    q isa JSON3.Object || throw(_RequestError(400, "params must be an object"))
    haskey(q, :min_identity) && (kw[:min_identity] = _num(q[:min_identity], "min_identity", 0, 100))
    haskey(q, :min_coverage) && (kw[:min_coverage] = _num(q[:min_coverage], "min_coverage", 0, 100))
    haskey(q, :max_delta) && (kw[:max_delta] = _num(q[:max_delta], "max_delta", 0, 100))
    haskey(q, :max_matches) && (kw[:max_matches] = Int(_num(q[:max_matches], "max_matches", 1, 1000)))
    haskey(q, :min_align_length) &&
        (kw[:min_align_length] = Int(_num(q[:min_align_length], "min_align_length", 0, 100_000)))
    if haskey(q, :min_agreement)
        v = _num(q[:min_agreement], "min_agreement", 0, 100)
        v > 50 || throw(_RequestError(400, "min_agreement must be above 50 (a strict majority)"))
        kw[:min_agreement] = v
    end
    if haskey(q, :rank_cap)
        q[:rank_cap] isa Bool || throw(_RequestError(400, "rank_cap must be true or false"))
        kw[:rank_cap] = q[:rank_cap]
    end
    if haskey(q, :rank_thresholds)
        t = q[:rank_thresholds]
        t isa JSON3.Object || throw(_RequestError(400, "rank_thresholds must be an object"))
        limits = copy(kw[:rank_thresholds])
        for (k, v) in pairs(t)
            r = Symbol(k)
            haskey(limits, r) || throw(_RequestError(400, "unknown rank in rank_thresholds: $k"))
            limits[r] = _num(v, "rank threshold for $k", 0, 100)
        end
        kw[:rank_thresholds] = limits
    end
    return LCAParams(; kw...)
end

function _text_param(p, key, default, maxlen)
    v = get(p, key, default)
    v isa AbstractString || throw(_RequestError(400, "$key must be text"))
    s = strip(v)
    length(s) <= maxlen || throw(_RequestError(400, "$key is longer than $maxlen characters"))
    any(iscntrl, s) && throw(_RequestError(400, "$key contains control characters"))
    return String(s)
end

function _options(app)
    p = app.defaults
    return (databases = sort!([(key = String(k), code = v) for (k, v) in BOLDassignment.BOLD_DATABASES];
                              by = x -> x.key),
            modes = sort!([(key = String(k), mi = v.mi, maxh = v.maxh, max_sequences = v.max_sequences)
                           for (k, v) in BOLDassignment.BOLD_MODES]; by = x -> x.mi, rev = true),
            defaults = params_dict(p), ranks = String.(p.ranks), accept_rank = String(app.accept_rank),
            lca_version = LCA_VERSION,
            global_note = GLOBAL_STATS_NOTE)
end

function _handle_upload(app)
    p = _payload()
    p isa JSON3.Object || throw(_RequestError(400, "expected a JSON object"))
    fasta = get(p, :fasta, nothing)
    fasta isa AbstractString || throw(_RequestError(400, "fasta must be text"))
    name = _text_param(p, :filename, "input.fasta", 200)
    input = try
        read_asvs(Vector{UInt8}(codeunits(fasta)); name, max_records = app.max_records)
    catch err
        err isa ArgumentError || rethrow()
        throw(_RequestError(400, err.msg))
    end
    token = string(uuid4())
    lock(() -> _keep!(app.uploads, app.upload_order, token, input, app.max_kept), app.lock)
    return _json(200, merge((token = token,), input_summary(input)))
end

function _handle_analyze(app)
    p = _payload()
    p isa JSON3.Object || throw(_RequestError(400, "expected a JSON object"))
    token = _text_param(p, :token, "", 64)
    input = lock(() -> get(app.uploads, token, nothing), app.lock)
    input === nothing && throw(_RequestError(404, "unknown upload; upload the FASTA again"))
    isempty(input.unique) && throw(_RequestError(400, "no valid sequence to analyse"))
    db = _text_param(p, :db, "animal_public", 40)
    haskey(BOLDassignment.BOLD_DATABASES, Symbol(db)) || throw(_RequestError(400, "unknown database"))
    mode = _text_param(p, :mode, "rapid", 40)
    haskey(BOLDassignment.BOLD_MODES, Symbol(mode)) || throw(_RequestError(400, "unknown mode"))
    params = try
        _params_from(p, app.defaults)
    catch err
        err isa ArgumentError || rethrow()
        throw(_RequestError(400, err.msg))
    end
    accept = Symbol(_text_param(p, :accept_rank, String(app.accept_rank), 20))
    accept in params.ranks || throw(_RequestError(400, "unknown accept_rank"))
    request = Dict{String,Any}("project" => _text_param(p, :project, "default", 200),
                               "user" => _text_param(p, :user, "", 200),
                               "db" => Symbol(db), "mode" => Symbol(mode), "params" => params,
                               "accept_rank" => accept)
    job = _Job(string(uuid4()), :queued, ["Queued"], nothing, nothing, nothing, request, input)
    lock(app.lock) do
        app.stopping[] && throw(_RequestError(503, "the server is stopping"))
        _n_pending(app) >= app.max_pending &&
            throw(_RequestError(503, "too many analyses waiting ($(app.max_pending)); try again later"))
        _keep_job!(app, job)
        put!(app.queue, job)
    end
    return _json(200, (id = job.id, status = job.status))
end

function _run_job!(app, job)
    job.status = :running
    note(msg) = lock(() -> push!(job.messages, msg), app.lock)
    try
        r = job.request
        a = run_analysis(job.input; project = r["project"], user = r["user"], params = r["params"],
                         accept_rank = r["accept_rank"], db = r["db"], mode = r["mode"],
                         engine = app.engine(), checkpoint = app.checkpoint,
                         portal_url = app.portal_url, progress = note)
        rec = record_analysis!(app.store, a)
        rec.recorded || note("Same analysis already recorded ($(rec.analysis_id)): usage counters not incremented")
        foreach(note, a.errors)
        job.analysis = a
        job.recorded = rec
        job.status = :done
        note("Finished")
    catch err
        err isa InterruptException && rethrow()
        @error "Analysis failed" exception = (err, catch_backtrace())
        job.error = first(sprint(showerror, err), 2000)
        job.status = :failed
    end
end

function _job_json(app, job)
    msgs = lock(() -> copy(job.messages), app.lock)
    base = (id = job.id, status = job.status, messages = msgs, error = job.error,
            input = input_summary(job.input), project = job.request["project"],
            user = job.request["user"])
    job.status === :done || return base
    return merge(base, (summary = analysis_summary(job.analysis), recorded = job.recorded,
                        provenance_id = job.analysis.id))
end

function _matches_json(a::Analysis, hash::AbstractString)
    rows = [r for r in match_table(a) if r.sequence_hash == hash]
    isempty(rows) && !any(u -> u.hash == hash, a.input.unique) &&
        throw(_RequestError(404, "unknown sequence"))
    return rows
end

function _register_routes(app::_App)
    R = Genie.Router
    for (path, (_, type)) in _STATIC
        body = _STATIC_BODY[path]
        R.route(path; named = Symbol(:boldtrace_static, path)) do
            _response(200, body, type; cache = endswith(path, ".png") ? "public, max-age=86400" : "no-store")
        end
    end
    wrap(f) = () -> try
        f()
    catch err
        err isa _RequestError && return _error(err.status, err.msg)
        err isa ArgumentError && return _error(400, err.msg)
        @error "Request failed" exception = (err, catch_backtrace())
        _error(500, "internal error, see the server log")
    end
    post(path, f, name) = R.route(path; method = R.POST, named = name) do
        _is_json(Genie.Requests.request()) || return _error(415, "send application/json")
        wrap(f)()
    end
    get_(path, f, name) = R.route(wrap(f), path; named = name)

    get_("/health", () -> _json(200, (status = "ok",)), :boldtrace_health)
    get_("/api/options", () -> _json(200, _options(app)), :boldtrace_options)
    post("/api/upload", () -> _handle_upload(app), :boldtrace_upload)
    post("/api/analyze", () -> _handle_analyze(app), :boldtrace_analyze)
    get_("/api/job", () -> _json(200, _job_json(app, _get_job(app))), :boldtrace_job)
    get_("/api/assignments", () -> _json(200, assignment_table(_done_analysis(app)[2])),
         :boldtrace_assignments)
    get_("/api/matches", () -> begin
             _, a = _done_analysis(app)
             _json(200, _matches_json(a, get(_query(), "hash", "")))
         end, :boldtrace_matches)
    get_("/api/references", () -> _json(200, reference_usage(_done_analysis(app)[2])),
         :boldtrace_references)
    get_("/api/attribution", () -> _json(200, taxonomic_attribution(_done_analysis(app)[2])),
         :boldtrace_attribution)
    get_("/api/export", () -> begin
             _, a = _done_analysis(app)
             name = get(_query(), "file", "")
             name in EXPORT_FILES || throw(_RequestError(404, "unknown export file"))
             hash = get(_query(), "hash", "")
             if !isempty(hash)
                 name == "bold_matches.csv" || throw(_RequestError(400, "hash applies only to bold_matches.csv"))
                 k = findfirst(u -> u.hash == hash, a.input.unique)
                 k === nothing && throw(_RequestError(404, "unknown sequence"))
                 # File name from the first input ID, reduced to safe characters.
                 id = replace(first(a.input.unique[k].ids), r"[^A-Za-z0-9._-]" => "_")
                 name = "bold_matches_" * first(id, 80) * ".csv"
             end
             pids = unique(r.process_id for r in values(a.references))
             body = export_file(a, startswith(name, "bold_matches") ? "bold_matches.csv" : name;
                                citations = citation_links(app.store, pids),
                                sequence_hash = isempty(hash) ? nothing : hash)
             type = endswith(name, ".json") ? "application/json; charset=utf-8" :
                    endswith(name, ".tsv") ? "text/tab-separated-values; charset=utf-8" :
                    "text/csv; charset=utf-8"
             _response(200, body, type;
                       extra = ["Content-Disposition" => "attachment; filename=\"$name\""])
         end, :boldtrace_export)
    get_("/api/stats", () -> begin
             q = _query()
             project = get(q, "project", "")
             user = get(q, "user", "")
             personal = usage_stats(app.store; project = isempty(project) ? nothing : project,
                                    user = isempty(user) ? nothing : user)
             _json(200, (personal = personal, global_ = usage_stats(app.store)))
         end, :boldtrace_stats)
    return nothing
end

function _check_request(req::HTTP.Request, allowed_hosts, max_body)
    if !isempty(allowed_hosts)
        h = lowercase(strip(HTTP.header(req, "Host")))
        name = startswith(h, '[') ? first(split(h[2:end], ']')) : first(split(h, ':'))
        (h in allowed_hosts || name in allowed_hosts) || throw(_RequestError(421, "unknown host"))
    end
    cl = HTTP.header(req, "Content-Length")
    if !isempty(cl)
        n = tryparse(Int, strip(cl))
        (n === nothing || n < 0) && throw(_RequestError(400, "invalid Content-Length"))
        n > max_body && throw(_RequestError(413, "request larger than $max_body bytes"))
    end
    return nothing
end

function _read_body(stream::HTTP.Stream, max_body::Int)
    buf = IOBuffer()
    while !eof(stream)
        write(buf, readavailable(stream))
        position(buf) > max_body && throw(_RequestError(413, "request larger than $max_body bytes"))
    end
    return take!(buf)
end

# Checks host and body size before Genie reads the request (Genie would read
# bodies of any size in memory), then hardens every response.
function _stream_handler(allowed_hosts::Set{String}, max_body::Int)
    return function (stream::HTTP.Stream)
        req = stream.message
        close_after = false
        res = try
            _check_request(req, allowed_hosts, max_body)
            req.body = _read_body(stream, max_body)
            HTTP.closeread(stream)
            Genie.Server.setup_http_listener(req, HTTP.Response(); stream)
        catch err
            if err isa _RequestError
                close_after = err.status in (400, 413)
                _error(err.status, err.msg)
            else
                @error "Request failed" exception = (err, catch_backtrace())
                _error(500, "internal error, see the server log")
            end
        end
        _harden!(res)
        close_after && HTTP.setheader(res, "Connection" => "close")
        req.response = res
        res.request = req
        HTTP.startwrite(stream)
        write(stream, res.body)
        if close_after
            close(stream.stream)
            stream.ntoread = 0
        end
        return
    end
end

"""
    serve(; database = "boldtrace.sqlite", host = "127.0.0.1", port = 8080, async = false,
          engine = () -> IDEngine(), checkpoint = "bold_checkpoint",
          portal_url = "https://portal.boldsystems.org", max_body = 64 * 2^20,
          max_records = 100_000, max_connections = 64, allowed_hosts = String[],
          max_kept = 50, max_pending = 20, defaults = LCAParams(), accept_rank = :phylum)

Start the web application (Genie.jl) on `http://host:port` with the usage
database at `database` (SQLite, created if needed). Sections: Upload,
Analysis, Assignments, Reference usage, Taxonomic attribution, Statistics
and Export.

Analyses run one at a time in a background task and query the official BOLD
v5 identification engine through `engine()` (rate limits of the API are
respected by `BOLDassignment.IDEngine`). Raw BOLD results are stored in
`checkpoint` (a directory, `nothing` to disable) and reused. `defaults` (an
[`LCAParams`](@ref)) and `accept_rank` are the settings proposed by the page;
every analysis can change them. At most
`max_pending` analyses may be queued or running (further requests get HTTP
503). Uploads and finished analyses are kept in memory (the last
`max_kept`); queued and running analyses are never evicted. Analyses and
events are stored in the database.

With `async = false` the call blocks until Ctrl-C; with `async = true` it
returns the server, stopped by [`shutdown`](@ref). The default `host`
accepts only local connections. There is no authentication: on a shared
server use a reverse proxy with HTTPS and set `allowed_hosts`.
"""
function serve(; database::AbstractString = "boldtrace.sqlite", host::AbstractString = "127.0.0.1",
               port::Integer = 8080, async::Bool = false, engine::Function = () -> IDEngine(),
               checkpoint::Union{Nothing,AbstractString} = "bold_checkpoint",
               portal_url::AbstractString = PORTAL_URL, max_body::Integer = 64 * 2^20,
               max_records::Integer = 100_000, max_connections::Integer = 64,
               allowed_hosts = String[], max_kept::Integer = 50, max_pending::Integer = 20,
               defaults::LCAParams = LCAParams(), accept_rank::Symbol = :phylum)
    accept_rank in defaults.ranks || throw(ArgumentError("accept_rank :$accept_rank is not a rank"))
    max_kept >= 1 || throw(ArgumentError("max_kept must be at least 1"))
    max_pending >= 1 || throw(ArgumentError("max_pending must be at least 1"))
    old = _APP[]
    old === nothing || old.stopping[] ||
        throw(ArgumentError("a BOLDTrace server is already running in this session; call shutdown first"))
    ip = string(Sockets.getaddrinfo(host))
    _check_port(ip, Int(port))
    store = open_store(database)
    app = _App(store, engine, checkpoint === nothing ? nothing : abspath(checkpoint),
               String(portal_url), Dict(), String[], Dict(), String[], Channel{_Job}(Inf),
               ReentrantLock(), Int(max_records), Int(max_kept), Int(max_pending), defaults,
               accept_rank,
               Threads.Atomic{Bool}(false), nothing, nothing)
    app.worker = errormonitor(Threads.@spawn _worker_loop(app))
    _APP[] = app
    Genie.config.app_env = "prod"
    Genie.config.server_document_root = mktempdir()
    Genie.config.server_handle_static_files = false
    Genie.config.server_signature = ""
    Genie.config.log_requests = false
    Genie.config.websockets_server = false
    _register_routes(app)
    hosts = Set{String}(lowercase(strip(String(h))) for h in allowed_hosts)
    listener = try
        HTTP.listen!(_stream_handler(hosts, Int(max_body)), parse(Sockets.IPAddr, ip), Int(port);
                     max_connections = Int(max_connections))
    catch
        _stop!(app; timeout = Inf)
        rethrow()
    end
    app.listener = listener
    server = Genie.Server.ServersCollection(webserver = listener)
    @info "BOLDTrace at http://$(host == "0.0.0.0" ? "localhost" : host):$port (database $(abspath(database)))"
    async && return server
    try
        _wait(listener)
    finally
        _stop!(app; timeout = Inf)
    end
    return server
end

# Stops accepting analyses, cancels the queued ones and waits (at most
# `timeout` seconds) for the running analysis; the worker then closes the
# store. Returns whether the worker has finished.
function _stop!(app::_App; timeout::Real = Inf)
    lock(app.lock) do
        app.stopping[] = true
        close(app.queue)
    end
    task = app.worker
    task === nothing && return true
    if !istaskdone(task)
        running = lock(() -> any(j -> j.status === :running, values(app.jobs)), app.lock)
        running && @info "Waiting for the running analysis to finish before closing the database"
        if isfinite(timeout)
            timedwait(() -> istaskdone(task), Float64(timeout); pollint = 0.1)
        else
            wait(task)
        end
    end
    if istaskdone(task)
        _APP[] === app && (_APP[] = nothing)
        return true
    end
    @warn "The running analysis has not finished after $timeout s; the database will be " *
          "closed when it ends"
    return false
end

function _wait(listener)
    interactive = isinteractive()
    interactive || Base.exit_on_sigint(false)
    try
        while isopen(listener)
            sleep(0.5)
        end
    catch err
        err isa InterruptException || rethrow()
    finally
        close(listener)
        interactive || Base.exit_on_sigint(true)
    end
    return nothing
end

function _check_port(ip::AbstractString, port::Int)
    sock = try
        Sockets.listen(parse(Sockets.IPAddr, ip), port)
    catch err
        err isa Base.IOError || rethrow()
        error("cannot listen on $ip:$port ($(err.msg)): is another server using the port?")
    end
    close(sock)
    return nothing
end

"""
    shutdown(server; timeout = Inf) -> Bool

Stop a server started with `serve(...; async = true)`. New requests are
refused at once; analyses still queued are cancelled; the running analysis,
if any, is completed and recorded, then the database is closed. Waits at
most `timeout` seconds for the running analysis and returns whether
everything has stopped; after a timeout the database is closed as soon as
the analysis ends.
"""
function shutdown(server::Genie.Server.ServersCollection; timeout::Real = Inf)
    Genie.Server.down(server)
    app = _APP[]
    app === nothing && return true
    return _stop!(app; timeout)
end
