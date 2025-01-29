source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions
library(rstan)
library(argparse)
# source('code/synthetic_benchmark/util_functions.R')

# ---- Directories ----

read_dir_plasma <-  here('data', 'interim', 'case_studies', 'shen2018', '01_get_read_counts_refv2')

code_dir <- here('code', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted')
read_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '02_extract_parameter_posterior_refv4')

write_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '03_benchmark_on_synthetic_mixtures_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)


# ---- Parsed arguments from bash ----
# ## For testing
# COVERAGE <- 1 # Read coverage: low (1), medium (2), high (3) respectively
# PROPORTION <- 0.35 # Proportion of the targeted cell type
# SEED <- 9 # Index of random sample from regression parameter posterior

parser <- ArgumentParser()
parser$add_argument("-c", "--coverage", type = "integer", help = "read coverage")
parser$add_argument("-p", "--proportion", type = "double", help = "targeted proportion of cancer tissue")
parser$add_argument("-s", "--seed", type = "integer", help = "read coverage")

args <- parser$parse_args()
COVERAGE   <- args$coverage # Read coverage: low (1), medium (2), high (3) respectively
PROPORTION <- args$proportion # Proportion of the targeted cell type
SEED       <- args$seed

# ---- important arguments to decemedip function ----

ref_cts <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                        'summarizedexperiment_cell_type_specific_atlas_v4.rds'))
ref_anc <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                        'summarizedexperiment_all_tissue_um_atlas_v4.rds'))


# ---- Generate read counts of the targeted cell type ----

# for (CELLTYPE in 24:25) {
for (CELLTYPE in 7:25) {
  
  name_seg <- paste0('_cellType', CELLTYPE, '_cov', COVERAGE, '_prop', PROPORTION, '_seed', SEED)
  
  ## Prepare regression parameters
  w_mu    <- fread(here(read_dir, paste0('param_posterior_w_mu_cov', COVERAGE, '.csv')))$mean
  w_sigma <- fread(here(read_dir, paste0('param_posterior_w_sigma_cov', COVERAGE, '.csv')))$mean
  pi    <- fread(here(read_dir, paste0('param_posterior_pi_cov', COVERAGE, '.csv')))$mean
  stopifnot(length(w_mu) == length(w_sigma))
  
  ## Prepare model variables
  X <- rbind(assays(ref_cts)[[1]], assays(ref_anc)[[1]])
  z <- c(rowData(ref_cts)$n_cpgs_100bp, rowData(ref_anc)$n_cpgs_100bp) |> log1p()
  
  ## Spike in targeted cell type proportion
  pi <- pi * (1 - PROPORTION)
  pi[CELLTYPE] <- pi[CELLTYPE] + PROPORTION
  
  ## Prepare input for stan model
  knots_z <- quantile(z, probs = c(0, 1))
  data_list <- list('N' = nrow(X), 'K' = ncol(X), 'X' = X, 'z' = z, 
                    'knots_z' = knots_z, 'n_knot_z' = length(knots_z), 'degree_z'  = 3,
                    'w_mu' = w_mu, 'w_sigma' = w_sigma, 'L' = length(w_sigma),
                    'pi' = pi)
  
  ## Generate the read counts
  model_dir <- here(code_dir, 'read_simulator.rds')
  if (file.exists(model_dir)) {
    model <- readRDS(model_dir)
  } else{
    model <- stan_model(here(code_dir, 'read_simulator.stan'))
    saveRDS(model, model_dir)
  }
  
  samp <- sampling(model,
                   data = data_list,
                   iter = 1, chain = 1,
                   algorithm = "Fixed_param",
                   seed = SEED)
  counts <- as.matrix(samp, pars = "y_sim") |> c() |> as.integer()
  counts_cts <- counts[1:nrow(ref_cts)] 
  counts_anc <- counts[(nrow(ref_cts)+1):length(counts)]
  
  output <- decemedip(counts_cts = counts_cts, counts_anc = counts_anc, 
                      diagnostics = FALSE, 
                      seed = 2023,
                      ref_cts = ref_cts,
                      ref_anc = ref_anc)
  fit <- output$posterior
  # plotDiagnostics(output, 'model_fit')
  
  smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> 
    as.data.frame() |>
    mutate(targeted_cell_type = CELLTYPE, coverage = COVERAGE, proportion = PROPORTION, seed = SEED)
  fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu', name_seg,'.csv.gz')))
  
  smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> 
    as.data.frame() |>
    mutate(targeted_cell_type = CELLTYPE, coverage = COVERAGE, proportion = PROPORTION, seed = SEED)
  fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma', name_seg,'.csv.gz')))
  
  cell_types <- colnames(ref_cts)
  smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
    as.data.frame() |>
    mutate(cell_type = factor(cell_types, levels = cell_types)) |>
    relocate(cell_type) |>
    mutate(targeted_cell_type = CELLTYPE, coverage = COVERAGE, proportion = PROPORTION, seed = SEED)
  fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi', name_seg,'.csv.gz')))
  
}
