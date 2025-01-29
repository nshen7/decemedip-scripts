source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions
library(rstan)

read_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '01_apply_deconv_on_healthy_plasma')
write_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '02_extract_parameter_posterior')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)


for (COVERAGE in 1:3) {
  
  name_seg <- paste0('_cov', COVERAGE)
  
  fit <- readRDS(here(read_dir, paste0('stanfit_healthy_plasma_cov', COVERAGE, '.rds')))
  
  w_mu    <- monitor(extract(fit, pars=c("w_mu"),    permuted = FALSE), print = FALSE) |>
    as.data.frame() |>
    mutate(param_name = paste0('w_mu_', 1:5)) |>
    relocate(param_name)
  fwrite(w_mu, here(write_dir, paste0('param_posterior_w_mu', name_seg, '.csv')))
  
  w_sigma <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), print = FALSE) |>
    as.data.frame() |>
    mutate(param_name = paste0('w_sigma_', 1:5)) |>
    relocate(param_name)
  fwrite(w_sigma, here(write_dir, paste0('param_posterior_w_sigma', name_seg, '.csv')))
  
  cell_types <- colnames(decemedip::hg19.ref.cts.se)
  smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), print = FALSE) |>
    as.data.frame() |>
    mutate(cell_type = factor(cell_types, levels = cell_types)) |>
    relocate(cell_type)
  fwrite(smr_pi.df, here(write_dir, paste0('param_posterior_pi', name_seg,'.csv')))

  print(COVERAGE)
}
