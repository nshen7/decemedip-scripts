source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions
library(argparse)
library(rstan)

# ---- Bash Command Line Argument Parsing ----

# ## for testing
# COVERAGE   <- 1
# PROPORTION <- 0.3
# SEED       <- 2


parser <- ArgumentParser()
parser$add_argument("-c", "--coverage", type = "integer", help = "read coverage")
parser$add_argument("-p", "--proportion", type = "double", help = "targeted proportion of cancer tissue")
parser$add_argument("-s", "--seed", type = "integer", help = "read coverage")

args <- parser$parse_args()
COVERAGE   <- args$coverage # Read coverage: 1, 2 representing medium, high respectively
PROPORTION <- args$proportion
SEED       <- args$seed

# ---- set file paths ----

read_dir_tissue <- here('data', 'interim', 'case_studies', 'berchuck2022', 'get_read_counts')
read_dir_plasma <-  here('data', 'interim', 'case_studies', 'shen2018_new', 'get_read_counts')
write_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_pdx_mixture', '01_benchmark_on_synthetic_mixtures_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

# ---- load processed data ----
tissue_samples <- c('LuCaP_78CR', 'LuCaP_81')
plasma_samples <- c('Control_7', 'Control_14') # Read coverage: ~48 mil reads, ~58 mil reads
name_seg <- paste0('_cov', COVERAGE, '_prop', PROPORTION, '_seed', SEED)

md_plasma <- fread(here('data', 'metadata', 'shen2018', 'sample_metadata_processed_shen2018.csv')) |>
  filter(group == 'Control')

tissue_cts.se <- readRDS(here(read_dir_tissue, paste0('summarizedexperiment_LuCap_PDX_samples_cts_regions.rds')))
tissue_anc.se <- readRDS(here(read_dir_tissue, paste0('summarizedexperiment_LuCap_PDX_samples_anc_regions.rds')))
plasma_cts.se <- readRDS(here(read_dir_plasma, paste0('summarizedexperiment_healthy_samples_cts_regions.rds')))
plasma_anc.se <- readRDS(here(read_dir_plasma, paste0('summarizedexperiment_healthy_samples_anc_regions.rds')))

plasma_sample <- switch (as.character(COVERAGE),
                         '1' = plasma_samples[1],
                         '2' = plasma_samples[2])
tissue_sample <- switch (as.character(COVERAGE),
                         '1' = tissue_samples[1],
                         '2' = tissue_samples[2])

counts_cts_plasma <- assays(plasma_cts.se)$counts[,plasma_sample] |> unlist()
counts_anc_plasma <- assays(plasma_anc.se)$counts[,plasma_sample] |> unlist()
counts_cts_tissue <- assays(tissue_cts.se)$counts[,tissue_sample] |> unlist()
counts_anc_tissue <- assays(tissue_anc.se)$counts[,tissue_sample] |> unlist()

if (sum(counts_cts_tissue) < 0.5 * sum(counts_cts_plasma)) stop('Not enough CTS reads from tissue!')
if (sum(counts_anc_tissue) < 0.5 * sum(counts_anc_plasma)) stop('Not enough anchor reads from tissue!')


# ---- utils ----

generateSyntheticMixture <- function(
    target_proportion, # targeted proportion of reads from tissue that is spiked in to healthy plasma
    reads_plasma, # vector of read counts of reference regions from a healthy plasma sample
    reads_tissue, # vector of read counts of reference regions from a tissue sample
    seed # random seed
) {
  
  if (length(reads_plasma) != length(reads_tissue)) stop('Lengths of input vectors not aligned.')
  if (class(reads_plasma) != "integer" | class(reads_tissue) != "integer") stop('Input vectors contain non-integer values.')
  
  set.seed(seed)
  
  ## Total number of reads wish to sample from tissue
  target_count <- round(target_proportion * sum(reads_plasma))

  ## Subsample read counts from tissue sample
  sample_prob_tissue <- target_count / sum(reads_tissue)
  sample_prob_plasma <- 1 - target_proportion
  reads_tissue_sampled <- map_int(reads_tissue, ~ rbinom(1, size = .x, prob = sample_prob_tissue))
  reads_plasma_sampled <- map_int(reads_plasma, ~ rbinom(1, size = .x, prob = sample_prob_plasma))

  reads_mix <- reads_plasma_sampled + reads_tissue_sampled
}

# ---- generate mixture and apply decemedip ----

counts_cts <- generateSyntheticMixture(target_proportion = PROPORTION, 
                                       reads_plasma = counts_cts_plasma, 
                                       reads_tissue = counts_cts_tissue, 
                                       seed = SEED)
counts_anc <- generateSyntheticMixture(target_proportion = PROPORTION, 
                                       reads_plasma = counts_anc_plasma, 
                                       reads_tissue = counts_anc_tissue, 
                                       seed = SEED)

output <- decemedip(counts_cts = counts_cts, 
                    counts_anc = counts_anc, 
                    diagnostics = FALSE)
fit <- output$posterior

smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> 
  as.data.frame() |>
  mutate(coverage = COVERAGE, proportion = PROPORTION, seed = SEED)
fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu', name_seg,'.csv.gz')))

smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> 
  as.data.frame() |>
  mutate(coverage = COVERAGE, proportion = PROPORTION, seed = SEED)
fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma', name_seg,'.csv.gz')))

cell_types <- colnames(decemedip::hg19.ref.cts.se)
smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
  as.data.frame() |>
  mutate(cell_type = factor(cell_types, levels = cell_types)) |>
  relocate(cell_type) |>
  mutate(coverage = COVERAGE, proportion = PROPORTION, seed = SEED)
fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi', name_seg,'.csv.gz')))


# plotDiagnostics(fit, plot_type = 'y_fit')

