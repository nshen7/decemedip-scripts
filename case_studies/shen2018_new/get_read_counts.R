source('code/SETPATHS.R')
library(argparse)
devtools::load_all('../decemedip/')

write_dir <- here('data', 'interim', 'case_studies', 'shen2018_new', 'get_read_counts')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

#### ONLY NEED HEALTHY CONTROL SAMPLES

# ---- load metadata ----

md_samples <- fread(here('data', 'metadata', 'shen2018', 'sample_metadata_processed_shen2018.csv')) |>
  filter(group == 'Control')

## Cell type specific reference
ref_cts.se <- decemedip::hg19.ref.cts.se

## Anchor region reference
ref_anc.se <- decemedip::hg19.ref.anc.se

# ---- obtain read count matrix ----

## Read counts from cell-free samples
samples_cts.se <- getRoiReadCount(sample_bam_files = md_samples$bamdir,
                                 sample_names = md_samples$sample,
                                 sample_paired = FALSE,
                                 roi = granges(ref_cts.se),
                                 col_data = md_samples)
saveRDS(samples_cts.se, here(write_dir, paste0('summarizedexperiment_healthy_samples_cts_regions.rds')))

samples_anc.se <- getRoiReadCount(sample_bam_files = md_samples$bamdir,
                                 sample_names = md_samples$sample,
                                 sample_paired = FALSE,
                                 roi = granges(ref_anc.se),
                                 col_data = md_samples)
saveRDS(samples_anc.se, here(write_dir, paste0('summarizedexperiment_healthy_samples_anc_regions.rds')))


