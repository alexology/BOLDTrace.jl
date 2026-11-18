# BOLDTrace.jl

Identification of amplicon sequence variants (ASVs) through
[BOLD Systems v5](https://id.boldsystems.org), lowest common ancestor (LCA)
assignment, traceability of the BOLD reference records used and
documentation of the people who identified the reference specimens. It
includes a web interface built on [Genie.jl](https://genieframework.com).

This package is part of an experimental effort to use AI to port and develop
bioinformatics and quantitative ecology tools in Julia, with the goal of
creating a unified environment for high-performance scientific computing.
High performance is a goal, not a measured guarantee.

```
FASTA → BOLD identification engine → matches → LCA → assignment → usage counters → export
```

## Installation

This package is not registered in the Julia General registry, and neither
is its dependency
[BOLDassignment.jl](https://github.com/alexology/BOLDassignment.jl). Install
both from GitHub, BOLDassignment.jl first:

```julia
using Pkg
Pkg.add(url = "https://github.com/alexology/BOLDassignment.jl")
Pkg.add(url = "https://github.com/alexology/BOLDTrace.jl")
```

or, from local copies, `Pkg.develop(path = "path/to/BOLDassignment.jl")` and
then `Pkg.develop(path = "path/to/BOLDTrace.jl")`. With Julia ≥ 1.11
`Project.toml` also points to `../BOLDassignment.jl` (`[sources]`), so a
checkout next to BOLDassignment.jl works directly.

## Web application

```
julia --threads=2 examples/run_app.jl
```

opens `http://127.0.0.1:8080` with seven sections:

- **Upload**: multi-FASTA, validation, number of sequences and of unique
  sequences, errors and warnings, file and sequence-set hashes;
- **Analysis**: BOLD library and search mode, every LCA parameter (identity,
  coverage, distance from the best match, number of matches, alignment
  length, minimum agreement, rank cap and its threshold for each rank),
  progress and errors. After an analysis, **Recompute LCA** applies new
  parameters to the BOLD hits and metadata already downloaded, in seconds
  and without any BOLD request; recomputations are recorded in the usage
  statistics only when "Record in usage statistics" is ticked;
- **Assignments**: ASV ID, taxon, rank, best identity, coverage, matches in
  the LCA, agreement and dissenting references; click a row to see every
  reference, included or excluded, with the reason. Downloads: assignments,
  all hits, and the hits of the selected ASV;
- **Reference usage**: process ID, taxon, identifier and the three usage
  counts (download);
- **Taxonomic attribution**: statistics per taxonomic identifier (download);
- **Statistics**: personal (project/user) and global counters;
- **Export**: the six files listed below.

`serve(; database, checkpoint, host, port, allowed_hosts, ...)` runs the
server; `serve(...; async = true)` returns it for `shutdown`. Analyses run
one at a time (the BOLD rate limits apply to the whole client); at most
`max_pending` (20) may wait, and queued or running analyses are never
dropped from memory; the last `max_kept` (50) finished analyses are kept,
always including the one just returned. `shutdown(server; timeout = Inf)`
cancels the queued analyses, lets the running one finish and be recorded,
then closes the database. It acts only on the server it was given and can be
called again safely. `serve(...; defaults = LCAParams(...), accept_rank = :genus)` sets
the parameters proposed by the page; each analysis can still change them.
The default host accepts only local connections and there is no authentication: on a
shared server use a reverse proxy with HTTPS and `allowed_hosts`. Request
bodies above `max_body` (64 MiB) are refused before being read; responses
carry a restrictive Content Security Policy.

## Use from Julia

```julia
using BOLDTrace

input = read_asvs("asvs.fasta")                    # validation, duplicates, hashes
a = run_analysis(input; project = "lake2026", user = "AL",
                 params = LCAParams(min_identity = 90, min_coverage = 80, max_delta = 1),
                 db = :animal_public, mode = :rapid, checkpoint = "bold_checkpoint")

b = reanalyze(a; params = LCAParams(min_agreement = 95))   # no BOLD request

assignment_table(a)          # one row per input record
match_table(a)               # every (ASV, reference) match with the LCA decision
reference_usage(a)           # retrieved / LCA / contribution counts per record
taxonomic_attribution(a)     # per identified_by value

store = open_store("boldtrace.sqlite")
record_analysis!(store, a)   # not recorded twice
usage_stats(store; project = "lake2026")   # personal
usage_stats(store)                         # global
export_analysis(a, "results/")
```

## Input

`read_asvs` reads FASTA with FASTX. Identifiers (first word of the header)
are kept. Sequences are normalized (no gaps or whitespace, upper case,
`U → T`) and hashed with SHA-256; identical sequences are analysed once and
mapped back to every input ID. Records with invalid characters, empty
sequences or an identifier reused for a different sequence are reported and
skipped without stopping the analysis; short sequences (< 100 bases, the
minimum overlap of the engine) and ambiguous ones are warnings. A file that
is not valid FASTA is rejected as a whole.

## Tests

```
julia --project=. --threads=4 test/runtests.jl
```

Simulated identification engine and portal (formats observed on BOLD),
the real BOLD response of BOLDassignment.jl, SQLite and the HTTP server; no
network access.

## License

MIT (see `LICENSE`). BOLD data retrieved with this package are subject to the
terms of use of BOLD Systems. The dnambientale.it logo in `src/assets` is the
one used by ASVAtlas.jl.
