source('code/SETPATHS.R')
library(argparse)
devtools::load_all('../decemedip/', compile = FALSE)

## Bash Command Line Argument Parsing 
parser <- ArgumentParser()
parser$add_argument("-b", "--binwidth", type = "integer", default = 1,
                    help = "Binwidth of reference regions")

args <- parser$parse_args()
BINWIDTH <- args$binwidth

write_dir <- here('data', 'interim', 'process_tissue_samples', 'berchuck2022', '04_get_read_counts_refv2')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

# ---- load metadata ----

md_tissue <- fread(here('data', 'metadata', 'Berchuck2022_LuCap_PDX_MeDIP', 'sample_metadata_processed_Berchuck2022_LuCap_PDX_MeDIP.csv'))

## Cell type specific reference (hg19)
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v2.rds'))
## Anchor region reference (hg19)
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v2.rds'))


# ---- obtain read count matrix ----

## Read counts from RCC tissues
tissue_cts.se <- getRoiReadCount(sample_bam_files = md_tissue$bam_dir,
                                 sample_names = md_tissue$Sample_Name,
                                 sample_paired = TRUE,
                                 roi = granges(ref_cts.se) |> resize(width = BINWIDTH, fix = 'center'),
                                 col_data = md_tissue)
saveRDS(tissue_cts.se, here(write_dir, paste0('summarizedexperiment_LuCap_PDX_samples_cts_regions_bw', BINWIDTH, '.rds')))

tissue_anc.se <- getRoiReadCount(sample_bam_files = md_tissue$bam_dir,
                                 sample_names = md_tissue$Sample_Name,
                                 sample_paired = TRUE,
                                 roi = granges(ref_anc.se) |> resize(width = BINWIDTH, fix = 'center'),
                                 col_data = md_tissue)
saveRDS(tissue_anc.se, here(write_dir, paste0('summarizedexperiment_LuCap_PDX_samples_anc_regions_bw', BINWIDTH, '.rds')))


