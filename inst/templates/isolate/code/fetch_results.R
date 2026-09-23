# Copy the pipeline outputs gmaio reads from `pipeline_remote` (e.g. a server) into
# `pipeline_results`, skipping reads, BAMs, assemblies and bins.
# Safe to rerun while the pipeline is still going: only new or changed files are copied.
# Run from the project directory in a terminal: Rscript code/fetch_results.R

gmaio::fetch_pipeline_results("config.yml")
