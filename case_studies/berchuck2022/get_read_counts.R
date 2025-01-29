source('code/SETPATHS.R')
library(argparse)
devtools::load_all('../decemedip/', compile = FALSE)

write_dir <- here('data', 'interim', 'case_studies', 'berchuck2022', 'get_read_counts')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

# ---- load metadata ----

md_tissue <- fread(here('data', 'metadata', 'Berchuck2022_LuCap_PDX_MeDIP', 
                        'sample_metadata_processed_Berchuck2022_LuCap_PDX_MeDIP.csv'))

## Cell type specific reference (hg19)
ref_cts.se <- decemedip::hg19.ref.cts.se

## Anchor region reference (hg19)
ref_anc.se <- decemedip::hg19.ref.anc.se


# ---- obtain read count matrix ----

## Read counts from RCC tissues
tissue_cts.se <- getRoiReadCount(sample_bam_files = md_tissue$bam_dir,
                                 sample_names = md_tissue$Sample_Name,
                                 sample_paired = TRUE,
                                 roi = granges(ref_cts.se),
                                 col_data = md_tissue)
saveRDS(tissue_cts.se, here(write_dir, paste0('summarizedexperiment_LuCap_PDX_samples_cts_regions.rds')))

tissue_anc.se <- getRoiReadCount(sample_bam_files = md_tissue$bam_dir,
                                 sample_names = md_tissue$Sample_Name,
                                 sample_paired = TRUE,
                                 roi = granges(ref_anc.se),
                                 col_data = md_tissue)
saveRDS(tissue_anc.se, here(write_dir, paste0('summarizedexperiment_LuCap_PDX_samples_anc_regions.rds')))


