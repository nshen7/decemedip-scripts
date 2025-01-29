source('code/SETPATHS.R')
library(parallel)
library(argparse)
mc.cores <- 14

md_dir <- here('data', 'metadata', 'medip_vs_wgbs')
read_dir_wgbs <- here('data', 'interim', 'medip_vs_wgbs', '02_lift_over_wgbs')

write_dir <- here('data', 'interim', 'medip_vs_wgbs_ref_regions_v4', '02_get_counts_wgbs')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

## Bash Command Line Argument Parsing
parser <- ArgumentParser()
parser$add_argument("-b", "--binwidth", type = "integer", default = 1,
                    help = "Binwidth of reference regions")

args <- parser$parse_args()
BINWIDTH <- args$binwidth

# ---- Fixed arguments ----
min_total_read <- 0 # min total read count per CpG in WGBS

# ---- load in reference regions and metadata ----
## Cell type specific reference
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v4.rds'))
## Anchor region reference
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v4.rds'))

md_wgbs <- read_tsv(here(md_dir, 'metadata_wgbs.tsv')) %>% filter(`File type` == 'bed')


# ---- Bin counts for WGBS profile of K562 and GM12878 cell lines ----
wrapper <- function(cell_line, marker_type) {
  
  target_regions.gr <- switch (
    marker_type,
    'cts' = granges(ref_cts.se) |> resize(width = BINWIDTH, fix = 'center'),
    'anc' = granges(ref_anc.se) |> resize(width = BINWIDTH, fix = 'center')
  )
  
  wgbs.df <- fread(here(read_dir_wgbs, paste0('WGBS_', cell_line, '_hg19.csv.gz'))) %>%
    filter(total_read >= min_total_read)
  wgbs.gr <- makeGRangesFromDataFrame(wgbs.df, keep.extra.columns = T)
  print('finished read wgbs')
  
  ## Find indices for WGBS CpGs for each reference site
  hits <- findOverlaps(wgbs.gr, target_regions.gr)
  bin_inds <- split(queryHits(hits), subjectHits(hits))
  print('finished hits')
  
  countByBin <- function(j) {
    ind <- bin_inds[[j]]
    gr <- wgbs.gr[ind]
    return(data.frame(total_read = sum(gr$total_read),
                      meth_read  = sum(gr$meth_read)))
  }
  bin_counts.df <- do.call(rbind, mclapply(1:length(bin_inds), function(j) countByBin(j), mc.cores = mc.cores))
  print('finished counting')
  
  all_bin_counts.df <- data.frame(total_read = rep(0, length(target_regions.gr)),
                                  meth_read  = rep(0, length(target_regions.gr)))
  all_bin_counts.df[as.integer(names(bin_inds)), ] <- bin_counts.df
  wgbs_bin.df <- cbind(as.data.frame(target_regions.gr), all_bin_counts.df) |>
    mutate(methyl_frac = meth_read / total_read)
  print('finished formatting')
  
  fwrite(wgbs_bin.df, here(write_dir, paste0('wgbs_', cell_line, '_', marker_type, '_bw', BINWIDTH, '.txt.gz')))
}

wrapper(cell_line = 'K562', marker_type = 'cts')
wrapper(cell_line = 'GM12878', marker_type = 'cts')
wrapper(cell_line = 'K562', marker_type = 'anc')
wrapper(cell_line = 'GM12878', marker_type = 'anc')
