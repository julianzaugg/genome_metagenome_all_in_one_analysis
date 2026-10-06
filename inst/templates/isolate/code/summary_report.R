# One-page HTML summary of the project: headline numbers and the results the analysis scripts recorded,
# with links to the full tables and figures. Optional: set outputs: report: true in config.yml, run the
# analysis scripts, then this one (last; it can be rerun at any time). Needs rmarkdown and pandoc.

project.l <- gmaio::load_processed()
config.l <- project.l$config
if (!gmaio::report_enabled(config.l)){
  gmaio::skip_analysis("The summary report is off; set outputs: report: true in config.yml and rerun the analysis scripts")
}
gmaio::write_summary_report(project.l)
