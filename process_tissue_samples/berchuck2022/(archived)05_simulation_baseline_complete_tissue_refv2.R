source('code/SETPATHS.R')
library(rstan)
library(bayesplot)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = 6)

code_dir <- here('code', 'process_tissue_samples')
read_dir <- here('data', 'interim', 'process_tissue_samples', 'berchuck2022', '04_get_read_counts_refv2')
write_dir <- here('data', 'interim', 'process_tissue_samples', 'berchuck2022', '05_simulation_baseline_complete_tissue_refv2')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'process_tissue_samples', 'berchuck2022', '05_simulation_baseline_complete_tissue_refv2')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

# ---- load in counts ----
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v2.rds'))
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v2.rds'))

BINWIDTH <- 1
tissue_cts.se <- readRDS(here(read_dir, paste0('summarizedexperiment_LuCap_PDX_samples_cts_regions_bw', BINWIDTH, '.rds')))
tissue_anc.se <- readRDS(here(read_dir, paste0('summarizedexperiment_LuCap_PDX_samples_anc_regions_bw', BINWIDTH, '.rds')))

## Sanity check
all(rowData(tissue_cts.se)$probe == rowData(ref_cts.se)$probe)
all(rowData(tissue_anc.se)$probe == rowData(ref_anc.se)$probe)

table(rowData(tissue_cts.se)$n_cpgs_100bp)
table(rowData(tissue_anc.se)$n_cpgs_100bp)

colQuantiles(rbind(assays(tissue_cts.se)[[1]], assays(tissue_anc.se)[[1]]), probs = c(seq(0, 1, 0.25), 0.999)) 
colSums(rbind(assays(tissue_cts.se)[[1]], assays(tissue_anc.se)[[1]]))

## Empirical corr matrix of reference cell types
mat <- as.matrix(assays(ref_cts.se)[[1]])
Xi <- cor(mat)

# ---- utils ----
prepareInputWrapper <- function(
    sample_name_tissue, # sample name of reference regions from a cancer tissue sample
    seed # random seed for generating read mixture
) {
  
  ## ! Always CTS before anchor 
  reads_tissue <- c(assays(tissue_cts.se)$counts[,sample_name_tissue], assays(tissue_anc.se)$counts[,sample_name_tissue]) 
  y <- reads_tissue
  X <- rbind(assays(ref_cts.se)$beta, assays(ref_anc.se)$beta)
  z <- c(rowData(ref_cts.se)$n_cpgs_100bp, rowData(ref_anc.se)$n_cpgs_100bp) |> log1p()
  
  # if (!all(names(y) == rownames(X))) stop('names of y does not match rownames of X')
  # if (!all(names(y) == names(z))) stop('names of y does not match names of z')
  
  return(list('N' = nrow(X), 'K' = ncol(X), 'X' = X, 'y' = y, 'z' = z))
}

################
##### TEST #####
################

# ---- Fit model with test data ----

# for (sample_name_tissue in colnames(tissue_cts.se)) {
for (sample_name_tissue in c('LuCaP_96CR', 'LuCaP_208.2', 'LuCaP_173.2')) {
  
  # ---- Prepare model input for testing ----
  
  seed <- 1818
  
  input_list <- prepareInputWrapper(sample_name_tissue = sample_name_tissue,
                                    seed = seed)
  
  prior_list <- list(
    'alpha'     = rep(1, input_list$K),
    's_mu'      = 1,
    's_sigma'   = 1,
    'n_knot_x'  = 2,
    'knots_x'   = c(0, 1),
    'degree_x'  = 3,
    # 'n_knot_z'  = 3,
    # 'knots_z'   = quantile(input_list$z, probs = c(0, 0.5, 1)),
    'n_knot_z'  = 2,
    'knots_z'   = quantile(input_list$z, probs = c(0, 1)),
    'degree_z'  = 3,
    'Xi' = Xi, 's_theta' = 1, 's_tau' = 1
  )
  
  data_list <- c(input_list, prior_list)
  
  # ---- Fit decemedip-model-2 model with test data ----
  
  # model <- stan_model(here(code_dir, 'decemedip-model.stan'))
  # model <- stan_model(here(code_dir, 'decemedip-model-2-no-intercept.stan'))
  # model <- stan_model(here(code_dir, 'decemedip-model-2.stan'))
  model <- stan_model(here(code_dir, 'decemedip-model-3.stan'))
  fit <- sampling(model,
                  data = data_list,
                  iter = 1000, chain = 4,
                  control = list(adapt_delta = 0.95,
                                 max_treedepth = 15))
  # plot(fit, pars=c("w_mu"))
  # plot(fit, pars=c("w_sigma"))
  # plot(fit, pars=c("tau"))
  # plot(fit, pars=c("theta"))
  # plot(fit, pars=c("eta"))
  # plot(fit, pars=c("pi"))

  smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
  # fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_decemedip_model_', sample_name_tissue,'.csv')))
  # fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_decemedip_model_2_noIntercept_', sample_name_tissue,'.csv')))
  # fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_decemedip_model_2_', sample_name_tissue,'.csv')))
  fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_decemedip_model_3_', sample_name_tissue,'.csv')))

  smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
  # fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_decemedip_model_', sample_name_tissue,'.csv')))
  # fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_decemedip_model_2_noIntercept_', sample_name_tissue,'.csv')))
  # fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_decemedip_model_2_', sample_name_tissue,'.csv')))
  fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_decemedip_model_3_', sample_name_tissue,'.csv')))

  smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
    as.data.frame() |>
    mutate(cell_type = factor(colnames(data_list$X), levels = rev(colnames(data_list$X)))) |>
    relocate(cell_type)
  # fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_decemedip_model_', sample_name_tissue,'.csv')))
  # fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_decemedip_model_2_noIntercept_', sample_name_tissue,'.csv')))
  # fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_decemedip_model_2_', sample_name_tissue,'.csv')))
  fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_decemedip_model_3_', sample_name_tissue,'.csv')))
  
  smr_pi.df |>
    ggplot(aes(y = cell_type)) +
    geom_point(aes(x = mean)) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`)) +
    ggtitle(paste0('RCC tissue sample: ', sample_name_tissue))
  # ggsave(here(plot_dir, paste0('interval_fitted_pi_decemedip_model_', sample_name_tissue,'.png')))
  # ggsave(here(plot_dir, paste0('interval_fitted_pi_decemedip_model_2_noIntercept_', sample_name_tissue,'.png')))
  # ggsave(here(plot_dir, paste0('interval_fitted_pi_decemedip_model_2_', sample_name_tissue,'.png')))
  ggsave(here(plot_dir, paste0('interval_fitted_pi_decemedip_model_3_', sample_name_tissue,'.png')))
  
  # traceplot(fit, pars = c("y_sim[1]"), inc_warmup = TRUE)
  #
  # y <- data_list$y
  y_sim <- as.matrix(fit, pars = "y_sim")
  # idx <- sample(nrow(y_sim), 100)
  # ppc_rootogram(y, y_sim[idx, ]) + xlim(-.5, 100)
  # ppc_dens_overlay(y, y_sim[idx, ]) + xlim(0, 20)
  # ppc_hist(y, y_sim[sample(nrow(y_sim), 8), ], binwidth = 1) + xlim(0, 20)
  # ppc_stat(y, y_sim[idx, ], stat = "max")
  # prop_zero <- function(x) mean(x == 0)
  # ppc_stat(y, y_sim[idx, ], stat = "prop_zero")
  # ppc_stat_2d(y, y_sim[idx, ], stat = c("mean", "sd"))
  # ppc_error_hist(y, y_sim[sample(nrow(y_sim), 8), ], binwidth = 1) + xlim(-30, 30)


  plot.df <- data.frame(
    y = data_list$y,
    x = as.matrix(data_list$X) %*% matrix(smr_pi.df$mean, ncol = 1),
    z = data_list$z,
    y_pred = colMeans(y_sim),
    y_pred_2.5 = colQuantiles(y_sim, probs = 0.025),
    y_pred_97.5 = colQuantiles(y_sim, probs = 0.975)
  )

  plot.df |>
    ggplot(aes(x = x)) +
    geom_ribbon(aes(ymin = y_pred_2.5, ymax = y_pred_97.5), fill = 'lightgrey') +
    geom_point(aes(y = y), size = 0.5, color = 'orange2') +
    geom_point(aes(y = y_pred), size = 0.5) +
    facet_wrap(~ exp(z), scales = 'free') +
    theme_classic()
  # ggsave(here(plot_dir, paste0('point_pred_y_vs_x_decemedip_model_', sample_name_tissue,'.png')))
  # ggsave(here(plot_dir, paste0('point_pred_y_vs_x_decemedip_model_2_noIntercept_', sample_name_tissue,'.png')))
  # ggsave(here(plot_dir, paste0('point_pred_y_vs_x_decemedip_model_2_', sample_name_tissue,'.png')))
  ggsave(here(plot_dir, paste0('point_pred_y_vs_x_decemedip_model_3_', sample_name_tissue,'.png')))
  
}

