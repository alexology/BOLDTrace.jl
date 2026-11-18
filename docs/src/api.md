# API

```@meta
CurrentModule = BOLDTrace
```

```@docs
BOLDTrace
```

## Analysis

```@docs
Analysis
MatchRow
ReferenceRecord
accepted
analyze_hits
is_incomplete
reanalyze
reference_conflicts
reference_taxon
run_analysis
supports
```

## Input

```@docs
InputRecord
InputSet
UniqueSequence
input_summary
read_asvs
sequence_hash
```

## LCA (BOLDTrace.LCA)

```@docs
LCA
```

```@docs
LCA.LCAParams
LCA.LCAResult
LCA.LCA_VERSION
LCA.Match
LCA.MatchDecision
LCA.lca
LCA.DEFAULT_RANKS
LCA.DEFAULT_RANK_THRESHOLDS
```

## Tables

```@docs
analysis_summary
assignment_table
match_table
reference_usage
taxonomic_attribution
```

## Export and provenance

```@docs
CITATION_RELATIONS
CITATION_STATUSES
CitationLink
EXPORT_FILES
export_analysis
export_file
provenance
used_process_ids
write_table
```

## Usage counters and citations

```@docs
GLOBAL_STATS_NOTE
UsageStore
add_citation!
citation_links
list_analyses
open_store
record_analysis!
usage_stats
```

## Web application

```@docs
serve
shutdown
```
