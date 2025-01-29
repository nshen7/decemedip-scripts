source('code/SETPATHS.R')
library(argparse)
devtools::load_all('../decemedip/', compile = FALSE)

## Bash Command Line Argument Parsing 
parser <- ArgumentParser()
parser$add_argument("-b", "--binwidth", type = "integer", default = 1,
                    help = "Binwidth of reference regions")

args <- parser$parse_args()
BINWIDTH <- args$binwidth

write_dir <- here('data', 'interim', 'case_studies', 'shen2018', '01_get_read_counts_refv2')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

# ---- load metadata ----

md_samples <- fread(here('data', 'metadata', 'shen2018', 'sample_metadata_processed_shen2018.csv'))

## Cell type specific reference
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v2.rds'))
## Anchor region reference
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v2.rds'))


# ---- obtain read count matrix ----

## Read counts from cell-free samples
samples_cts.se <- getRoiReadCount(sample_bam_files = md_samples$bamdir,
                                 sample_names = md_samples$sample,
                                 sample_paired = FALSE,
                                 roi = granges(ref_cts.se) |> resize(width = BINWIDTH, fix = 'center'),
                                 col_data = md_samples)
saveRDS(samples_cts.se, here(write_dir, paste0('summarizedexperiment_samples_cts_regions_bw', BINWIDTH, '.rds')))

samples_anc.se <- getRoiReadCount(sample_bam_files = md_samples$bamdir,
                                 sample_names = md_samples$sample,
                                 sample_paired = FALSE,
                                 roi = granges(ref_anc.se) |> resize(width = BINWIDTH, fix = 'center'),
                                 col_data = md_samples)
saveRDS(samples_anc.se, here(write_dir, paste0('summarizedexperiment_samples_anc_regions_bw', BINWIDTH, '.rds')))


