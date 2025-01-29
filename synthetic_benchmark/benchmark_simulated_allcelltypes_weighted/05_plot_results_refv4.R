source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions

# ---- set file paths ----
read_dir <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '04_summarize_results_refv4')
plot_dir <- here('plots', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '05_plot_results_refv4')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---- utils ----
CTNAMES <- colnames(decemedip::hg19.ref.cts.se)
getCellTypeName <- function(i) CTNAMES[i]

# ---- estimation of pi ----
smr_pi.df <- fread(here(read_dir, paste0('fitted_pi_summary.csv.gz'))) |>
  mutate(cell_type = gsub('_', ' ', cell_type)) |>
  mutate(cell_type = factor(cell_type, levels = unique(cell_type))) |>
  filter(targeted_cell_type != 7)

quantile(smr_pi.df$Rhat)
quantile(smr_pi.df$Bulk_ESS)

for (COVERAGE in 1:3) {
  
  smr_pi.df |>
    filter(coverage == COVERAGE) |>
    ggplot(aes(proportion, mean)) +
    geom_boxplot(aes(group = factor(proportion)), width = 0.5) +
    geom_jitter(size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = 'dashed', color = 'darkgrey') +
    facet_wrap(~ cell_type, nrow = 3) +
    ylab('Estimated cell type proportion from deconvolution') +
    xlab('Mixed-in proportion of reads from targeted cell type') +
    theme_classic()
  ggsave(here(plot_dir, paste0('boxplot_estimProp_vs_mixedProp_cov', COVERAGE, '.png')), width = 12, height = 6)
  
  smr_pi.df |>
    filter(coverage == COVERAGE) |>
    ggplot(aes(x = proportion, y = mean)) +
    geom_linerange(aes(group = proportion, ymin = `2.5%`, ymax = `97.5%`), 
                   position = position_dodge2(width = 0.035),
                   linewidth = 0.5, color="lightskyblue2", alpha = 0.5) +
    geom_linerange(aes(group = proportion, ymin = `25%`, ymax = `75%`), 
                   position = position_dodge2(width = 0.035),
                   linewidth = 0.5, color="deepskyblue4", alpha = 0.5) +
    geom_linerange(aes(group = proportion, x = proportion, y = proportion,
                       xmin = proportion - 0.02, xmax = proportion + 0.02), 
                   linewidth = 1, color="darkorange2") +
    geom_point(aes(group = proportion), 
               position = position_dodge2(width = 0.035),
               size = 0.5, color = 'grey20', shape = 10) +
    # geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = 'dashed', color = 'darkgrey') +
    facet_wrap(~ cell_type, nrow = 3, scales = 'free_y') +
    ylab('Estimated cell type proportion from deconvolution') +
    xlab('Mixed-in proportion of reads from targeted cell type') +
    theme_classic()
  ggsave(here(plot_dir, paste0('interval_estimProp_vs_mixedProp_cov', COVERAGE, '.png')), width = 12, height = 6)
  
}



COVERAGECOLORS <- grafify:::graf_palettes$PrGn_div[3:1]
names(COVERAGECOLORS) <- c('Low', 'Medium', 'High')
smr_pi.df |>
  mutate(coverage = factor(coverage, levels = c(1, 2, 3), labels = c('Low', 'Medium', 'High'))) |>
  group_by(cell_type, coverage, seed) |>
  summarise(correlation = cor(mean, proportion)) |>
  ungroup() |>
  ggplot(aes(cell_type, correlation, group = interaction(coverage, cell_type), color = coverage)) +
  geom_line(position = position_dodge(width = 0.6)) +
  geom_point(size = 1, position = position_dodge(width = 0.6)) +
  theme_classic() +
  # cowplot::theme_cowplot() +
  scale_color_manual(values = COVERAGECOLORS, name = 'Coverage') +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Pearson correlation') + xlab('Cell type')
ggsave(here(plot_dir, paste0('point_correlation_vs_cellype.png')), width = 6, height = 4)

# ---- estimated total non-blood proportions ----

smr_pi_nonblood.df <- fread(here(read_dir, paste0('fitted_pi_nonblood_summary.csv.gz'))) |>
  filter(targeted_cell_type != 7) |>
  mutate(targeted_cell_type = map_chr(targeted_cell_type, getCellTypeName))

for (COVERAGE in 1:3) {
    
  smr_pi_nonblood.df |>
    filter(coverage == COVERAGE) |>
    ggplot(aes(proportion, total)) +
    geom_boxplot(aes(group = factor(proportion)), width = 0.5, color = 'grey20') +
    geom_jitter(size = 0.5, color = 'grey20') +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = 'dashed', color = 'darkgrey') +
    facet_wrap(~ targeted_cell_type, nrow = 3) +
    ylab('Estimated proportion of prostate from deconvolution') +
    xlab('Mixed-in proportion of reads from targeted cell type') +
    theme_classic()
  ggsave(here(plot_dir, paste0('boxplot_estimProp_vs_mixedProp_cov', COVERAGE, '.png')), width = 12, height = 6)
  
  smr_pi_nonblood.df |>
    filter(coverage == COVERAGE) |>
    ggplot(aes(x = proportion, y = total)) +
    geom_linerange(aes(group = proportion, x = proportion, y = proportion,
                       xmin = proportion - 0.02, xmax = proportion + 0.02), 
                   linewidth = 1, color="darkorange2") +
    geom_point(aes(group = proportion),
               position = position_dodge2(width = 0.035),
               size = 0.3, color = 'grey20', shape = 10) +
    facet_wrap(~ targeted_cell_type, nrow = 3, scales = 'free_y') +
    ylab('Estimated proportion of prostate from deconvolution') +
    xlab('Mixed-in proportion of reads from targeted cell type') +
    theme_classic()
  ggsave(here(plot_dir, paste0('interval_estimProp_vs_mixedProp_cov', COVERAGE, '.png')), width = 12, height = 6)
  
}

# ---- estimation of regression parameters ----


for (COVERAGE in 1:3) {
  
  smr_w_mu.df <- fread(here(read_dir, paste0('fitted_w_mu_summary.csv.gz'))) |>
    filter(coverage == COVERAGE)
  smr_w_sigma.df <- fread(here(read_dir, paste0('fitted_w_sigma_summary.csv.gz'))) |>
    filter(coverage == COVERAGE)

  read_dir_2 <- here('data', 'interim', 'synthetic_benchmark_new', 'benchmark_simulated_allcelltypes_weighted', '02_extract_parameter_posterior_refv4')
  true_w_mu.df    <- fread(here(read_dir_2, paste0('param_posterior_w_mu_cov', COVERAGE, '.csv'))) |>
    dplyr::select(param_name, mean, `2.5%`, `97.5%`) |>
    dplyr::rename('true_mean' = mean, 'true_2.5%' = `2.5%`, 'true_97.5%' = `97.5%`)
  true_w_sigma.df <- fread(here(read_dir_2, paste0('param_posterior_w_sigma_cov', COVERAGE, '.csv'))) |>
    dplyr::select(param_name, mean, `2.5%`, `97.5%`)  |>
    dplyr::rename('true_mean' = mean, 'true_2.5%' = `2.5%`, 'true_97.5%' = `97.5%`)

  w_mu_avg.df <- smr_w_mu.df |>
    group_by(param_name, proportion) |>
    summarise(avg_across_celltypes = mean(mean))
  w_sigma_avg.df <- smr_w_sigma.df |>
    group_by(param_name, proportion) |>
    summarise(avg_across_celltypes = mean(mean))
    
  smr_w_mu.df |>
    left_join(true_w_mu.df, by = 'param_name') |>
    ggplot(aes(proportion, mean)) +
    geom_ribbon(aes(ymin = `true_2.5%`, ymax = `true_97.5%`), fill = 'grey') + 
    geom_path(aes(group = factor(targeted_cell_type):factor(seed), 
                  color = factor(targeted_cell_type)), size = 0.5) +
    geom_hline(aes(yintercept = true_mean), linetype = 'dashed', size = 1) + 
    geom_path(data = w_mu_avg.df, aes(proportion, avg_across_celltypes), color = 'green', size = 1) +
    facet_wrap(~ param_name, nrow = 2, scale = 'free_y') +
    ylab('Estimated regression parameters for mu from synthetic data') +
    xlab('Mixed-in proportion of reads from targeted cell type') +
    theme_classic()
  ggsave(here(plot_dir, paste0('lines_estimMuParams_vs_mixedProp_cov', COVERAGE, '.png')), width = 12, height = 6)
  
  smr_w_sigma.df |>
    left_join(true_w_sigma.df, by = 'param_name') |>
    ggplot(aes(proportion, mean)) +
    geom_ribbon(aes(ymin = `true_2.5%`, ymax = `true_97.5%`), fill = 'grey') + 
    geom_path(aes(group = factor(targeted_cell_type):factor(seed), 
                  color = factor(targeted_cell_type)), size = 0.5) +
    geom_hline(aes(yintercept = true_mean), linetype = 'dashed', size = 1) + 
    geom_path(data = w_sigma_avg.df, aes(proportion, avg_across_celltypes), color = 'green', size = 1) +
    facet_wrap(~ param_name, nrow = 2, scale = 'free_y') +
    ylab('Estimated regression parameters for sigma from synthetic data') +
    xlab('Mixed-in proportion of reads from targeted cell type') +
    theme_classic()
  ggsave(here(plot_dir, paste0('lines_estimSigmaParams_vs_mixedProp_cov', COVERAGE, '.png')), width = 12, height = 6)
  
}

