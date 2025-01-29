source('code/SETPATHS.R')
library(rstan)
library(bayesplot)
library(decemedip)
devtools::load_all('../decemedip/', compile = FALSE)

rstan_options(auto_write = TRUE)
options(mc.cores = 6)

code_dir <- here('code', 'medip_vs_wgbs_ref_regions_v4')
md_dir <- here('data', 'metadata', 'medip_vs_wgbs')
read_dir_medip <- here('data', 'interim', 'medip_vs_wgbs_ref_regions_v4', '01_get_counts_medip')
read_dir_wgbs <- here('data', 'interim', 'medip_vs_wgbs_ref_regions_v4', '02_get_counts_wgbs')

write_dir <- here('data', 'interim', 'medip_vs_wgbs_ref_regions_v4', '03_fit_model_and_plot')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'medip_vs_wgbs_ref_regions_v4', '03_fit_model_and_plot')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

## Arguments
BINWIDTH <- 1
# cell_line <- 'K562'
cell_line <- 'GM12878'

# ---- Read in MeDIP-seq and WGBS counts ----

medip_read <- c(
  fread(here(read_dir_medip, paste0('medip_', cell_line, '_cts_bw', BINWIDTH, '.txt.gz')))$medip_read,
  fread(here(read_dir_medip, paste0('medip_', cell_line, '_anc_bw', BINWIDTH, '.txt.gz')))$medip_read
)

wgbs_mf <- c(
  fread(here(read_dir_wgbs, paste0('wgbs_', cell_line, '_cts_bw', BINWIDTH, '.txt.gz')))$methyl_frac,
  fread(here(read_dir_wgbs, paste0('wgbs_', cell_line, '_anc_bw', BINWIDTH, '.txt.gz')))$methyl_frac
)
  
n_cpgs_100bp <- c(
  fread(here(read_dir_medip, paste0('medip_', cell_line, '_cts_bw', BINWIDTH, '.txt.gz')))$n_cpgs_100bp,
  fread(here(read_dir_medip, paste0('medip_', cell_line, '_anc_bw', BINWIDTH, '.txt.gz')))$n_cpgs_100bp
)

df <- data.frame(medip_read, wgbs_mf, n_cpgs_100bp) |> 
  mutate(weights_raw = c(rep(1,   SummarizedExperiment::nrow(decemedip::hg19.ref.cts.se)), 
                         rep(0.5, SummarizedExperiment::nrow(decemedip::hg19.ref.anc.se)))) |> 
  na.omit() ### REMOVE sites with NA methyl fractions in WGBS

nrow(df) # = 3341 for K562; 3406 for GM12878

# ---- Fit same GLM as decemedip model ----

xtilde      <- df$wgbs_mf
y           <- df$medip_read
z           <- df$n_cpgs_100bp |> log1p()
weights_raw <- df$weights_raw

input_list <- list('N' = length(xtilde), 
                   'y' = y, 
                   'z' = z, 
                   'xtilde' = xtilde,
                   'weights_raw' = weights_raw)

### DEFAULT settings in decemedip
prior_list <- list(
  's_mu'      = 3,
  's_sigma'   = 3,
  'n_knot_z'  = 2,
  'knots_z'   = quantile(input_list$z, probs = c(0, 1)),
  'degree_z'  = 3
)

data_list <- c(input_list, prior_list)
saveRDS(data_list, here(write_dir, paste0('datalist_decemedip_glm_', BINWIDTH, 'bp_', cell_line, '.rds')))
fit <- stan(
  file = here(code_dir, 'decemedip_glm.stan'), 
  data = data_list,
  iter = 2000, chains = 4, seed = 123
)
saveRDS(fit, here(write_dir, paste0('stanfit_decemedip_glm_', BINWIDTH, 'bp_', cell_line, '.rds')))


# ----- Plot predicted vs. true val -----

fit <- readRDS(here(write_dir, paste0('stanfit_decemedip_glm_', BINWIDTH, 'bp_', cell_line, '.rds')))

y_sim <- as.matrix(fit, pars = "y_sim")
plot.df <- data.frame(
  y = data_list$y,
  x = data_list$xtilde,
  z = data_list$z,
  y_pred = colMeans(y_sim),
  y_pred_2.5 = colQuantiles(y_sim, probs = 0.025),
  y_pred_97.5 = colQuantiles(y_sim, probs = 0.975)
) |>
  mutate(n_cpgs_100bp = factor(exp(z)-1))

plot.df |>
  ggplot(aes(x = x)) +
  geom_ribbon(aes(ymin = y_pred_2.5, ymax = y_pred_97.5), fill = 'lightgrey') +
  geom_linerange(aes(ymin = y, ymax = y_pred), size = 0.5, color = 'darkgrey') +
  geom_point(aes(y = y), size = 0.5, color = 'orange2') +
  geom_point(aes(y = y_pred), size = 0.5) +
  facet_wrap(~ n_cpgs_100bp, scales = 'free_y') +
  theme_classic()
ggsave(here(plot_dir, paste0('point_y&yhat_vs_x_', BINWIDTH, 'bp_', cell_line, '.png')),
       height = 6, width = 8)
       

# smr.df <- plot.df |>
#   mutate(bin = cut(y_pred, breaks = seq(0, max(y_pred), 1))) |>
#   filter(!is.na(bin)) |>
#   group_by(bin) |>
#   # group_by(bin, n_cpgs_100bp) |>
#   summarise(n_bins             = n(),
#             bin_start          = min(y_pred),
#             y_mean   = mean(y),
#             y_median = median(y),
#             y_q0.025     = quantile(y, 0.025),
#             y_q0.25     = quantile(y, 0.25),
#             y_q0.75     = quantile(y, 0.75),
#             y_q0.975     = quantile(y, 0.975)) |> 
#   ungroup()
# 
# smr.df |>
#   ggplot(aes(x = bin_start)) +
#   geom_errorbar(aes(ymin = pmax(y_q0.025, 0),
#                     ymax = pmin(y_q0.975, max(y_median))),
#                 color = 'darkgrey',
#                 alpha = 1) +
#   geom_errorbar(aes(ymin = pmax(y_q0.25, 0),
#                     ymax = pmin(y_q0.75, max(y_median))),
#                 color = 'lightpink2', size = 1) +
#   geom_point(aes(y = y_mean, size = n_bins), color = 'dodgerblue4') +
#   geom_abline(intercept = 0, slope = 1, color = 'black', linetype = 'dashed', size = 1) +
#   scale_size_continuous(range = c(0.5, 2)) + 
#   # facet_wrap(~ n_cpgs_100bp, scales = 'fixed') +
#   ylim(0, max(smr.df$bin_start)) +
#   theme_classic() 
