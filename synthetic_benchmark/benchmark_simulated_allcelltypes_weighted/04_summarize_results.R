source('code/SETPATHS.R')

# ---- set file paths ----

read_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '03_benchmark_on_synthetic_mixtures')
write_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '04_summarize_results')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

# ---- summarize results ----

smr_pi.df          <- data.frame()
smr_pi_nonblood.df <- data.frame()
smr_w_mu.df        <- data.frame()
smr_w_sigma.df     <- data.frame()

for (CELLTYPE in 7:25) {
  for (COVERAGE in 1:3) {
    for (PROPORTION in c(0.01, seq(0.05, 0.5, 0.05))) {
      for (SEED in 1:10) {
        
        name_seg <- paste0('_cellType', CELLTYPE, '_cov', COVERAGE, '_prop', PROPORTION, '_seed', SEED)
        
        smr_pi_temp.df     <- fread(here(read_dir, paste0('fitted_pi', name_seg,'.csv.gz'))) |> 
          dplyr::slice(targeted_cell_type[1])
        smr_pi_nonblood_temp.df <- fread(here(read_dir, paste0('fitted_pi', name_seg,'.csv.gz'))) |>
          dplyr::slice(8:n()) |> # remove blood cells
          group_by(targeted_cell_type, coverage, proportion, seed) |>
          summarise(total = sum(mean)) |>
          ungroup()
          
        smr_w_mu_temp.df    <- fread(here(read_dir, paste0('fitted_w_mu', name_seg,'.csv.gz'))) |>
          mutate(param_name = paste0('w_mu_', 1:5)) |>
          relocate(param_name)
        smr_w_sigma_temp.df <- fread(here(read_dir, paste0('fitted_w_sigma', name_seg,'.csv.gz'))) |>
          mutate(param_name = paste0('w_sigma_', 1:5)) |>
          relocate(param_name)
        
        smr_pi.df          <- rbind(smr_pi.df, smr_pi_temp.df)
        smr_pi_nonblood.df <- rbind(smr_pi_nonblood.df, smr_pi_nonblood_temp.df)
        smr_w_mu.df        <- rbind(smr_w_mu.df, smr_w_mu_temp.df)
        smr_w_sigma.df     <- rbind(smr_w_sigma.df, smr_w_sigma_temp.df)
      }
    }
  }
  cat(CELLTYPE, ' ')
}

fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_summary.csv.gz')))
fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_summary.csv.gz')))
fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_summary.csv.gz')))
fwrite(smr_pi_nonblood.df, here(write_dir, paste0('fitted_pi_nonblood_summary.csv.gz')))

