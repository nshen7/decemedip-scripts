source('code/SETPATHS.R')

# ---- set file paths ----
read_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_pdx_mixture', '01_benchmark_on_synthetic_mixtures')
write_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_pdx_mixture', '02_summarize_results')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
plot_dir <- here('plots', 'synthetic_benchmark_new', 'benchmark_pdx_mixture', '02_summarize_results')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

## Sample pairs used in simultion
# tissue_samples <- c('LuCaP_78CR', 'LuCaP_81')
# plasma_samples <- c('Control_7', 'Control_14') # Read coverage: ~48 mil reads, ~58 mil reads

# ---- get results ----
# ## RUN ONCE
# #### joint the tables of results from synthetic benchmark
# smr_pi.df      <- data.frame()
# smr_w_mu.df    <- data.frame()
# smr_w_sigma.df <- data.frame()
# for (COVERAGE in 1:2) {
#   for (PROPORTION in c(0.01, seq(0.05, 0.5, 0.05))) {
#     for (SEED in 1:10) {
# 
#       name_seg <- paste0('_cov', COVERAGE, '_prop', PROPORTION, '_seed', SEED)
#       
#       smr_pi_temp.df      <- fread(here(read_dir, paste0('fitted_pi', name_seg,'.csv.gz')))
#       smr_w_mu_temp.df    <- fread(here(read_dir, paste0('fitted_w_mu', name_seg,'.csv.gz')))
#       smr_w_sigma_temp.df <- fread(here(read_dir, paste0('fitted_w_sigma', name_seg,'.csv.gz')))
# 
#       smr_pi.df      <- rbind(smr_pi.df, smr_pi_temp.df)
#       smr_w_mu.df    <- rbind(smr_w_mu.df, smr_w_mu_temp.df)
#       smr_w_sigma.df <- rbind(smr_w_sigma.df, smr_w_sigma_temp.df)
# 
#     }
#     cat(PROPORTION, ' ')
#   }
# }
# 
# fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_summary.csv.gz')))
# fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_summary.csv.gz')))
# fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_summary.csv.gz')))

# ---- plots on pi ----

smr_pi.df <- fread(here(write_dir, paste0('fitted_pi_summary.csv.gz'))) |>
  mutate(coverage = factor(coverage, labels = c('LuCaP_78CR & Control_7 (library size ~ 48 million reads)',
                                                'LuCaP_81 & Control_14 (library size ~ 58 million reads)')))

smr_pi.df |>
  filter(cell_type == 'Prostate') |>
  ggplot(aes(proportion, mean)) +
  geom_linerange(aes(group = proportion, ymin = `2.5%`, ymax = `97.5%`), 
                 position = position_dodge2(width = 0.035),
                 linewidth = 0.5, color="lightskyblue2", alpha = 0.5) +
  geom_linerange(aes(group = proportion, ymin = `25%`, ymax = `75%`), 
                 position = position_dodge2(width = 0.035),
                 linewidth = 0.5, color="deepskyblue4", alpha = 0.5) +
  geom_linerange(aes(group = proportion, x = proportion, y = proportion,
                     xmin = proportion - 0.02, xmax = proportion + 0.02), 
                 linewidth = 1, color="darkorange2") +
  # geom_linerange(data = pdx_pi.df, aes(group = proportion, x = proportion, y = proportion,
  #                    xmin = proportion - 0.02, xmax = proportion + 0.02), 
  #                linewidth = 1, color="darkorange2") +
  geom_point(aes(group = proportion), 
             position = position_dodge2(width = 0.035),
             size = 0.5, color = 'grey20', shape = 10) +
  facet_wrap(~ coverage) +
  ylab('Estimated proportion of prostate') +
  xlab('Mixed-in proportion of PRAD PDX reads') +
  theme_classic()
ggsave(here(plot_dir, paste0('interval_estimProp_vs_mixedProp.png')), width = 8, height = 4)

