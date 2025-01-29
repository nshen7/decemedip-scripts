source('code/SETPATHS.R')
library(rstan)
library(bayesplot)
library(argparse)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = 6)
source('code/sim_studies/util_functions.R')


## Bash Command Line Argument Parsing 
parser <- ArgumentParser()
parser$add_argument("-g", "--group", type = "character", 
                    help = "Sample group")

args <- parser$parse_args()
GROUP <- args$group


code_dir <- here('code', 'sim_studies')
read_dir <-  here('data', 'interim', 'case_studies', 'shen2018', '01_get_read_counts_refv2')
write_dir <- here('data', 'interim', 'case_studies', 'shen2018', '02_apply_deconv')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'case_studies', 'shen2018', '02_apply_deconv')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

# ---- load in counts ----
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v2.rds'))
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v2.rds'))

BINWIDTH <- 1
samples_cts.se <- readRDS(here(read_dir, paste0('summarizedexperiment_samples_cts_regions_bw', BINWIDTH, '.rds')))
samples_anc.se <- readRDS(here(read_dir, paste0('summarizedexperiment_samples_anc_regions_bw', BINWIDTH, '.rds')))

## Sanity check
all(rowData(samples_cts.se)$probe == rowData(ref_cts.se)$probe)
all(rowData(samples_anc.se)$probe == rowData(ref_anc.se)$probe)

table(rowData(samples_cts.se)$n_cpgs_100bp)
table(rowData(samples_anc.se)$n_cpgs_100bp)

colSums(rbind(assays(samples_cts.se)[[1]], assays(samples_anc.se)[[1]]))
colQuantiles(rbind(assays(samples_cts.se)[[1]], assays(samples_anc.se)[[1]]), probs = c(seq(0, 1, 0.25), 0.999)) 


# --- plot PCA of the samples ----
cts_mat <- assays(samples_cts.se)[[1]] |> subset(rowVars(assays(samples_cts.se)[[1]]) > 0)
anc_mat <- assays(samples_anc.se)[[1]] |> subset(rowVars(assays(samples_anc.se)[[1]]) > 0)
pca_result <- prcomp(rbind(cts_mat, anc_mat) |> t() |> log1p(), scale. = TRUE)  # Scale the data
summary(pca_result)
data.frame(colData(samples_cts.se), pca_result$x[, 1:2]) |>
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(color = group))

# ---- util function ----

prepareInputWrapper <- function(
    sample_name, # sample name of reference regions from a cancer tissue sample
    seed # random seed for generating read mixture
) {
  
  ## ! Always CTS before anchor 
  reads <- c(assays(samples_cts.se)$counts[,sample_name], 
             assays(samples_anc.se)$counts[,sample_name]) 
  idx_incl <- which(reads < 400)
  y <- reads[idx_incl]
  X <- rbind(assays(ref_cts.se)$beta, assays(ref_anc.se)$beta)[idx_incl,]
  z <- c(rowData(ref_cts.se)$n_cpgs_100bp, rowData(ref_anc.se)$n_cpgs_100bp)[idx_incl] |> log1p()
  
  # if (!all(names(y) == rownames(X))) stop('names of y does not match rownames of X')
  # if (!all(names(y) == names(z))) stop('names of y does not match names of z')
  
  return(list('N' = nrow(X), 'K' = ncol(X), 'X' = X, 'y' = y, 'z' = z))
}

################
##### TEST #####
################


# ---- Fit decemedip-model model with test data ----

samples_in_group <- colnames(samples_cts.se)[grep(GROUP, colnames(samples_cts.se))]

for (sample_name in samples_in_group) {
  
  # ---- Prepare model input for testing ----
  
  seed <- 1818
  
  input_list <- prepareInputWrapper(sample_name = sample_name,
                                    seed = seed)
  
  
  ## Empirical corr matrix of reference cell types
  mat <- as.matrix(assays(ref_cts.se)[[1]])
  Xi <- cor(mat)
  
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
  
  # model <- stan_model(here(code_dir, 'decemedip-model-2.stan'))
  model <- stan_model(here(code_dir, 'decemedip-model-3.stan'))
  fit <- sampling(model,
                  data = data_list,
                  iter = 1000, chain = 4,
                  control = list(adapt_delta=0.95),
                  seed = 10)
  # plot(fit, pars=c("theta"))
  # plot(fit, pars=c("w_mu"))
  # plot(fit, pars=c("w_sigma"))
  # plot(fit, pars=c("pi"))
  
  smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
  # fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_decemedip_model_2_', sample_name,'.csv')))
  fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_decemedip_model_3_', sample_name,'.csv')))
  
  smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
  # fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_decemedip_model_2_', sample_name,'.csv')))
  fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_decemedip_model_3_', sample_name,'.csv')))
  
  smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
    as.data.frame() |>
    mutate(cell_type = factor(colnames(data_list$X), levels = rev(colnames(data_list$X)))) |>
    relocate(cell_type)
  # fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_decemedip_model_2_', sample_name,'.csv')))
  fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_decemedip_model_3_', sample_name,'.csv')))
  
  smr_pi.df |>
    ggplot(aes(y = cell_type)) +
    geom_point(aes(x = mean)) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`)) +
    ggtitle(paste0('Lung tissue sample: ', sample_name))
  # ggsave(here(plot_dir, paste0('interval_fitted_pi_decemedip_model_2_', sample_name,'.png')), width = 9, height = 7)
  ggsave(here(plot_dir, paste0('interval_fitted_pi_decemedip_model_3_', sample_name,'.png')), width = 9, height = 7)
  
  # traceplot(fit, pars = c("y_sim[1]"), inc_warmup = TRUE)
  #
  y <- data_list$y
  y_sim <- as.matrix(fit, pars = "y_sim")
  idx <- sample(nrow(y_sim), 100)
  # ppc_rootogram(y, y_sim[idx, ]) + xlim(-.5, 100)
  # ppc_dens_overlay(y, y_sim[idx, ]) + xlim(0, 20)
  # ppc_hist(y, y_sim[sample(nrow(y_sim), 8), ], binwidth = 1) + xlim(0, 20)
  # ppc_stat(y, y_sim[idx, ], stat = "max")
  # prop_zero <- function(x) mean(x == 0)
  # ppc_stat(y, y_sim[idx, ], stat = "prop_zero")
  # ppc_error_hist(y, y_sim[sample(nrow(y_sim), 8), ], binwidth = 1) + xlim(-30, 30)
  
  p1 <- ppc_dens_overlay(y, y_sim[idx, ]) + xlim(0, 20)
  p2 <- ppc_stat_2d(y, y_sim[idx, ], stat = c("mean", "sd"))
  p3 <- ppc_stat(y, y_sim[idx, ], stat = "max")
  prop_zero <- function(x) mean(x == 0)
  p4 <- ppc_stat(y, y_sim[idx, ], stat = "prop_zero")
  cowplot::plot_grid(p1, p2, p3, p4, label_size = 12)
  # ggsave(here(plot_dir, paste0('bayesplot_diagnostics_decemedip_model_2_', sample_name,'.png')), width = 12, height = 10)
  ggsave(here(plot_dir, paste0('bayesplot_diagnostics_decemedip_model_3_', sample_name,'.png')), width = 12, height = 10)

  smr_y_sim.df <- monitor(extract(fit, pars=c("y_sim"), permuted = FALSE), print = FALSE, digits_summary = 5) |>
    as.data.frame()
  plot.df <- data.frame(
    y = data_list$y,
    x = as.matrix(data_list$X) %*% matrix(smr_pi.df$mean, ncol = 1),
    z = data_list$z,
    y_pred = smr_y_sim.df$mean,
    y_pred_2.5 = smr_y_sim.df$`2.5%`,
    y_pred_97.5 = smr_y_sim.df$`97.5%`
  )

  plot.df |>
    ggplot(aes(x = x)) +
    geom_ribbon(aes(ymin = y_pred_2.5, ymax = y_pred_97.5), fill = 'lightgrey') +
    geom_point(aes(y = y), size = 0.5, color = 'orange2') +
    geom_point(aes(y = y_pred), size = 0.5) +
    facet_wrap(~ exp(z), scales = 'free_y') +
    theme_classic()
  # ggsave(here(plot_dir, paste0('point_pred_y_vs_x_decemedip_model_2_', sample_name,'.png')), width = 9, height = 7)
  ggsave(here(plot_dir, paste0('point_pred_y_vs_x_decemedip_model_3_', sample_name,'.png')), width = 9, height = 7)
  
}

