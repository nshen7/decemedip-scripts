source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/', compile = FALSE)
library(argparse)

md_dir <- here('data', 'metadata', 'medip_vs_wgbs')
read_dir_medip <- here('data', 'raw', 'medip_vs_wgbs')

write_dir <- here('data', 'interim', 'medip_vs_wgbs_ref_regions_v4', '01_get_counts_medip')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

## Bash Command Line Argument Parsing
parser <- ArgumentParser()
parser$add_argument("-b", "--binwidth", type = "integer", default = 1,
                    help = "Binwidth of reference regions")

args <- parser$parse_args()
BINWIDTH <- args$binwidth

# ---- load data ----

md_medip <- read_tsv(here(md_dir, 'metadata_medip.tsv')) |>
  mutate(bam_file = here(read_dir_medip, 'medip', paste0(`File accession`, '.bam')))

## Cell type specific reference
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v4.rds'))
## Anchor region reference
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v4.rds'))

# ---- Use MEDIPS to process MeDIP-seq profile (hg19) of K562 and GM12878 cell lines ----

wrapper <- function(cell_line, marker_type) {
  
  target_regions.gr <- switch (
    marker_type,
    'cts' = granges(ref_cts.se) |> resize(width = BINWIDTH, fix = 'center'),
    'anc' = granges(ref_anc.se) |> resize(width = BINWIDTH, fix = 'center')
  )
  
  idx <- which(md_medip$`Biosample term name` == cell_line)
  medip_count.se <- getRoiReadCount(sample_bam_files = md_medip$bam_file[idx],
                                    sample_names = md_medip$`Biosample term name`[idx],
                                    sample_paired = FALSE,
                                    roi = target_regions.gr)
  
  medip_count.df <- assays(medip_count.se)[[1]] |> as.data.frame()
  names(medip_count.df) <- 'medip_read'
  medip_count.df <- cbind(as.data.frame(target_regions.gr), medip_count.df)
  
  fwrite(medip_count.df, here(write_dir, paste0('medip_', cell_line, '_', marker_type, '_bw', BINWIDTH, '.txt.gz')))
}

wrapper(cell_line = 'K562', marker_type = 'cts')
wrapper(cell_line = 'GM12878', marker_type = 'cts')
wrapper(cell_line = 'K562', marker_type = 'anc')
wrapper(cell_line = 'GM12878', marker_type = 'anc')