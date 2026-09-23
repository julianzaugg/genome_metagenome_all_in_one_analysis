# Bulk processing for isolate genomes: load every available pipeline output, link genomes
# to metadata, write the standard tables and save the processed project.
# Run from the project directory: Rscript code/main.R

config.l <- gmaio::read_config("config.yml")

project.l <- gmaio::start_project(config.l)
project.l <- gmaio::add_isolate_genomes(project.l)
project.l <- gmaio::add_comparisons(project.l)
project.l <- gmaio::add_mobile_elements(project.l)
project.l <- gmaio::add_genomespot(project.l)
project.l <- gmaio::add_palettes(project.l)

gmaio::write_isolate_tables(project.l)
gmaio::save_processed(project.l)
