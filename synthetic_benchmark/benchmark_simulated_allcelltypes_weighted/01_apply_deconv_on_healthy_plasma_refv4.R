source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions
library(rstan)

read_dir_plasma <-  here('data', 'interim', 'case_studies', 'shen2018', '01_get_read_counts_refv2')
write_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '01_apply_deconv_on_healthy_plasma_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)


# ---- important arguments ----

ref_cts <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v4.rds'))
ref_anc <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v4.rds'))

# ---- load sample read counts ----

md_samples <- fread(here('data', 'metadata', 'shen2018', 'sample_metadata_processed_shen2018.csv'))

# ---- apply deconvolution on 3 healthy plasma samples with low, medium, high coverage ----

plasma_samples <- c('Control_24', 'Control_3', 'Control_14') 
# Read coverage respectively: 1. low (21 mil reads), 2. medium (44 mil reads), 3. high (58 mil reads) respectively

for (COVERAGE in 1) {
# for (COVERAGE in 2) {
# for (COVERAGE in 3) {
  
  plasma_sample <- switch (as.character(COVERAGE),
                           '1' = plasma_samples[1],
                           '2' = plasma_samples[2],
                           '3' = plasma_samples[3])
  
  output <- decemedip(sample_bam_file = md_samples$bamdir[md_samples$sample == plasma_sample],
                      paired_end = FALSE,
                      diagnostics = TRUE,
                      ref_cts = ref_cts,
                      ref_anc = ref_anc
  )
  fit <- output$posterior
  
  saveRDS(fit, here(write_dir, paste0('stanfit_healthy_plasma_cov', COVERAGE, '.rds')))
}

