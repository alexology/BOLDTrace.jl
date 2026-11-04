# Persistent usage events in SQLite. The SQL is kept portable to PostgreSQL:
# TEXT/INTEGER/DOUBLE PRECISION columns, ISO-8601 UTC timestamps as TEXT, no
# AUTOINCREMENT, upserts with ON CONFLICT. Moving to PostgreSQL needs another
# DBInterface driver and `$n` placeholders instead of `?`.

const SCHEMA_VERSION = 1

"""
Note shown with every global statistic.
"""
const GLOBAL_STATS_NOTE = "Global statistics count only the identifications recorded by this " *
                          "application, not all uses of BOLD worldwide."

const _SCHEMA = [
    """CREATE TABLE IF NOT EXISTS schema_info (
        version INTEGER NOT NULL)""",
    """CREATE TABLE IF NOT EXISTS analyses (
        analysis_id TEXT PRIMARY KEY,
        fingerprint TEXT NOT NULL UNIQUE,
        created_at TEXT NOT NULL,
        project TEXT NOT NULL,
        user_name TEXT NOT NULL,
        input_name TEXT,
        input_file_sha256 TEXT NOT NULL,
        input_set_hash TEXT NOT NULL,
        n_records INTEGER NOT NULL,
        n_unique INTEGER NOT NULL,
        database_code TEXT,
        mode TEXT,
        lca_version TEXT NOT NULL,
        params_json TEXT NOT NULL,
        provenance_json TEXT NOT NULL)""",
    """CREATE TABLE IF NOT EXISTS asv_results (
        analysis_id TEXT NOT NULL REFERENCES analyses(analysis_id),
        sequence_hash TEXT NOT NULL,
        asv_ids TEXT NOT NULL,
        status TEXT NOT NULL,
        taxon TEXT,
        taxon_rank TEXT,
        accepted INTEGER NOT NULL,
        best_identity DOUBLE PRECISION,
        n_selected INTEGER NOT NULL,
        PRIMARY KEY (analysis_id, sequence_hash))""",
    """CREATE TABLE IF NOT EXISTS bold_references (
        record_key TEXT PRIMARY KEY,
        process_id TEXT NOT NULL,
        marker TEXT,
        record_id TEXT,
        sequence_hash TEXT,
        taxon TEXT,
        taxon_rank TEXT,
        identified_by TEXT,
        identification_method TEXT,
        bin_uri TEXT,
        recordsets TEXT,
        first_seen TEXT NOT NULL,
        last_seen TEXT NOT NULL)""",
    """CREATE TABLE IF NOT EXISTS identification_events (
        event_key TEXT PRIMARY KEY,
        analysis_id TEXT NOT NULL REFERENCES analyses(analysis_id),
        project TEXT NOT NULL,
        user_name TEXT NOT NULL,
        created_at TEXT NOT NULL,
        sequence_hash TEXT NOT NULL,
        record_key TEXT NOT NULL,
        process_id TEXT NOT NULL,
        identity DOUBLE PRECISION NOT NULL,
        coverage DOUBLE PRECISION,
        align_length INTEGER,
        included_in_lca INTEGER NOT NULL,
        decision TEXT NOT NULL,
        accepted INTEGER NOT NULL,
        final_taxon TEXT,
        final_rank TEXT)""",
    """CREATE TABLE IF NOT EXISTS citation_links (
        process_id TEXT NOT NULL,
        doi TEXT NOT NULL,
        relation TEXT NOT NULL,
        evidence_source TEXT NOT NULL,
        verification_status TEXT NOT NULL,
        PRIMARY KEY (process_id, doi, relation))""",
    "CREATE INDEX IF NOT EXISTS events_project ON identification_events(project)",
    "CREATE INDEX IF NOT EXISTS events_user ON identification_events(user_name)",
    "CREATE INDEX IF NOT EXISTS events_record ON identification_events(record_key)",
]

"""
    UsageStore

Database of recorded analyses and identification events. Open with
[`open_store`](@ref), close with `close(store)`.
"""
mutable struct UsageStore
    db::SQLite.DB
    path::String
    lock::ReentrantLock
end

"""
    open_store(path = ":memory:") -> UsageStore

Open (and create if needed) the SQLite database at `path`.
"""
function open_store(path::AbstractString = ":memory:")
    db = SQLite.DB(path)
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    for sql in _SCHEMA
        DBInterface.execute(db, sql)
    end
    v = _scalar(db, "SELECT MAX(version) FROM schema_info")
    if v === missing || v === nothing
        DBInterface.execute(db, "INSERT INTO schema_info (version) VALUES (?)", (SCHEMA_VERSION,))
    elseif v > SCHEMA_VERSION
        throw(ArgumentError("database schema version $v is newer than this package ($SCHEMA_VERSION)"))
    end
    return UsageStore(db, String(path), ReentrantLock())
end

Base.close(s::UsageStore) = (SQLite.close(s.db); nothing)

function _scalar(db, sql, params = ())
    for row in DBInterface.execute(db, sql, params)
        return row[1]
    end
    return missing
end

_db(x) = x === nothing ? missing : x
_db(x::Bool) = x ? 1 : 0

"""
    record_analysis!(store, a::Analysis) -> NamedTuple

Store the analysis, its per-sequence results, the reference records and one
identification event per (unique sequence, reference record). An analysis
with the same fingerprint (same project, user, sequence set, BOLD database
and mode, LCA version and parameters) is not recorded again: the result has
`recorded = false` and the identifier of the existing analysis. Event keys
are unique, so no event can be counted twice.
"""
function record_analysis!(store::UsageStore, a::Analysis)
    lock(store.lock) do
        existing = _scalar(store.db, "SELECT analysis_id FROM analyses WHERE fingerprint = ?",
                           (a.fingerprint,))
        existing === missing ||
            return (analysis_id = String(existing), recorded = false, n_events = 0)
        now = _iso(a.created_at)
        n_events = 0
        SQLite.transaction(store.db) do
            src = a.source
            DBInterface.execute(store.db,
                """INSERT INTO analyses (analysis_id, fingerprint, created_at, project, user_name,
                   input_name, input_file_sha256, input_set_hash, n_records, n_unique,
                   database_code, mode, lca_version, params_json, provenance_json)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (a.id, a.fingerprint, now, a.project, a.user, a.input.name, a.input.file_sha256,
                 a.input.set_hash, length(a.input.records), length(a.input.unique),
                 string(get(src, "database_code", "")), string(get(src, "mode", "")),
                 LCA_VERSION, JSON3.write(params_dict(a.params)), JSON3.write(provenance(a))))
            for (i, u) in enumerate(a.input.unique)
                res = a.results[i]
                ok = accepted(a, i)
                DBInterface.execute(store.db,
                    """INSERT INTO asv_results (analysis_id, sequence_hash, asv_ids, status, taxon,
                       taxon_rank, accepted, best_identity, n_selected)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (a.id, u.hash, join(u.ids, ";"), String(res.status), _db(res.taxon),
                     res.rank === nothing ? missing : String(res.rank), _db(ok),
                     res.best_identity, res.n_selected))
                for m in a.matches[i]
                    ref = a.references[m.key]
                    key = bytes2hex(SHA.sha256(join((a.fingerprint, u.hash, m.key), '\n')))
                    DBInterface.execute(store.db,
                        """INSERT INTO identification_events (event_key, analysis_id, project,
                           user_name, created_at, sequence_hash, record_key, process_id, identity,
                           coverage, align_length, included_in_lca, decision, accepted,
                           final_taxon, final_rank)
                           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                           ON CONFLICT (event_key) DO NOTHING""",
                        (key, a.id, a.project, a.user, now, u.hash, m.key, ref.process_id,
                         m.identity, m.coverage, m.align_length, _db(m.included),
                         m.included && !supports(a, i, m) ? "dissent:$(m.dissent)" : String(m.reason),
                         _db(ok), _db(res.taxon),
                         res.rank === nothing ? missing : String(res.rank)))
                    n_events += 1
                end
            end
            for ref in values(a.references)
                taxon, rank = reference_taxon(ref, a.params.ranks)
                DBInterface.execute(store.db,
                    """INSERT INTO bold_references (record_key, process_id, marker, record_id,
                       sequence_hash, taxon, taxon_rank, identified_by, identification_method,
                       bin_uri, recordsets, first_seen, last_seen)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                       ON CONFLICT (record_key) DO UPDATE SET
                         record_id = COALESCE(excluded.record_id, bold_references.record_id),
                         sequence_hash = COALESCE(excluded.sequence_hash, bold_references.sequence_hash),
                         taxon = excluded.taxon, taxon_rank = excluded.taxon_rank,
                         identified_by = COALESCE(excluded.identified_by, bold_references.identified_by),
                         identification_method = COALESCE(excluded.identification_method,
                                                          bold_references.identification_method),
                         bin_uri = excluded.bin_uri,
                         recordsets = COALESCE(excluded.recordsets, bold_references.recordsets),
                         last_seen = excluded.last_seen""",
                    (ref.key, ref.process_id, ref.marker, ref.record_id, ref.sequence_hash, taxon,
                     rank, ref.identified_by, ref.identification_method, ref.bin_uri,
                     ref.recordsets, now, now))
            end
        end
        return (analysis_id = a.id, recorded = true, n_events = n_events)
    end
end

"""
    usage_stats(store; project = nothing, user = nothing, analysis_id = nothing) -> NamedTuple

Counts over the recorded identification events, optionally restricted to a
project, a user or one analysis (personal counter); without filters, global
counter of the application:

- `n_analyses`, `n_projects`: recorded analyses and their projects, including
  analyses without any BOLD match;
- `n_events` (unique sequence × reference record pairs retrieved),
  `n_lca_events` (pairs in an LCA set), `n_contribution_events` (pairs in an
  LCA set of an accepted assignment);
- `n_asvs` (distinct sequences, by hash), `n_process_ids`,
  `n_reference_sequences` (distinct reference sequence hashes; records whose
  sequence was not retrieved count by record key), `n_process_ids_in_lca`;
- `note`: [`GLOBAL_STATS_NOTE`](@ref).
"""
function usage_stats(store::UsageStore; project = nothing, user = nothing, analysis_id = nothing)
    where_, params = _filters(project, user, analysis_id)
    q(sql) = coalesce(_scalar(store.db, sql, params), 0)
    # Analyses and projects come from `analyses`, so that analyses without any
    # BOLD match count too; the other figures describe reference usage.
    an = "FROM analyses e $where_"   # same column names as the events
    lock(store.lock) do
        ev = "FROM identification_events e $where_"
        (scope = project === nothing && user === nothing && analysis_id === nothing ? "global" : "personal",
         project = project, user = user, analysis_id = analysis_id,
         n_analyses = q("SELECT COUNT(*) $an"),
         n_projects = q("SELECT COUNT(DISTINCT e.project) $an"),
         n_events = q("SELECT COUNT(*) $ev"),
         n_lca_events = q("SELECT COUNT(*) $ev $(_and(where_)) e.included_in_lca = 1"),
         n_contribution_events = q("SELECT COUNT(*) $ev $(_and(where_)) e.included_in_lca = 1 AND e.accepted = 1 AND e.decision = 'selected'"),
         n_asvs = q("SELECT COUNT(DISTINCT e.sequence_hash) $ev"),
         n_process_ids = q("SELECT COUNT(DISTINCT e.process_id) $ev"),
         n_process_ids_in_lca = q("SELECT COUNT(DISTINCT e.process_id) $ev $(_and(where_)) e.included_in_lca = 1"),
         n_reference_sequences = q("""SELECT COUNT(DISTINCT COALESCE(r.sequence_hash, 'record:' || e.record_key))
                                      FROM identification_events e
                                      JOIN bold_references r ON r.record_key = e.record_key $where_"""),
         note = GLOBAL_STATS_NOTE)
    end
end

_and(where_) = isempty(where_) ? "WHERE" : "AND"

function _filters(project, user, analysis_id)
    conds = String[]
    params = Any[]
    project === nothing || (push!(conds, "e.project = ?"); push!(params, project))
    user === nothing || (push!(conds, "e.user_name = ?"); push!(params, user))
    analysis_id === nothing || (push!(conds, "e.analysis_id = ?"); push!(params, analysis_id))
    return (isempty(conds) ? "" : "WHERE " * join(conds, " AND ")), params
end

"""
    list_analyses(store; project = nothing, user = nothing) -> Vector{NamedTuple}

Recorded analyses, most recent first.
"""
function list_analyses(store::UsageStore; project = nothing, user = nothing)
    where_, params = _filters(project, user, nothing)
    where_ = replace(where_, "e.user_name" => "user_name", "e.project" => "project")
    lock(store.lock) do
        rows = NamedTuple[]
        for r in DBInterface.execute(store.db,
                """SELECT analysis_id, created_at, project, user_name, input_name, n_records,
                   n_unique, database_code FROM analyses $where_ ORDER BY created_at DESC""", params)
            push!(rows, (analysis_id = r.analysis_id, created_at = r.created_at, project = r.project,
                         user = r.user_name, input_name = r.input_name, n_records = r.n_records,
                         n_unique = r.n_unique, database_code = r.database_code))
        end
        rows
    end
end

"""
    add_citation!(store, link::CitationLink)

Store a bibliographic link (replaces the evidence and status of an existing
link with the same process ID, DOI and relation).
"""
function add_citation!(store::UsageStore, c::CitationLink)
    lock(store.lock) do
        DBInterface.execute(store.db,
            """INSERT INTO citation_links (process_id, doi, relation, evidence_source,
               verification_status) VALUES (?, ?, ?, ?, ?)
               ON CONFLICT (process_id, doi, relation) DO UPDATE SET
                 evidence_source = excluded.evidence_source,
                 verification_status = excluded.verification_status""",
            (c.process_id, c.doi, c.relation, c.evidence_source, c.verification_status))
    end
    return nothing
end

"""
    citation_links(store, process_ids = nothing) -> Vector{CitationLink}

Stored bibliographic links, optionally only for `process_ids`.
"""
function citation_links(store::UsageStore, process_ids = nothing)
    wanted = process_ids === nothing ? nothing : Set(String.(process_ids))
    links = CitationLink[]
    lock(store.lock) do
        for r in DBInterface.execute(store.db, "SELECT * FROM citation_links ORDER BY process_id, doi")
            wanted === nothing || r.process_id in wanted || continue
            push!(links, CitationLink(r.process_id, r.doi, r.relation, r.evidence_source,
                                      r.verification_status))
        end
    end
    return links
end
