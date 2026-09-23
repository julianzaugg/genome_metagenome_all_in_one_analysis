# Bulk processing: load every available pipeline output, link samples to metadata,
# write the standard tables and save the processed project for the analysis scripts.
# Run from the project directory: Rscript code/main.R

config.l <- gmaio::read_config("config.yml")

project.l <- gmaio::start_project(config.l)
project.l <- gmaio::add_sylph(project.l)
project.l <- gmaio::add_singlem(project.l)
project.l <- gmaio::add_mags(project.l)
project.l <- gmaio::add_gene_catalogue_functions(project.l)
project.l <- gmaio::add_nonpareil(project.l)
project.l <- gmaio::add_genomespot(project.l)
project.l <- gmaio::add_mobile_elements(project.l)
project.l <- gmaio::add_palettes(project.l)

gmaio::write_project_tables(project.l)
gmaio::save_processed(project.l)
