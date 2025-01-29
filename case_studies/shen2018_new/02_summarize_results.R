source('code/SETPATHS.R')
library(rstan)
library(bayesplot)
library(argparse)
library(Rsamtools)

code_dir <- here('code', 'sim_studies')
read_dir <-  here('data', 'interim', 'case_studies', 'shen2018_new', '01_apply_deconv')
write_dir <- here('data', 'interim', 'case_studies', 'shen2018_new', '02_summarize_results')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'case_studies', 'shen2018_new', '02_summarize_results')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

md_samples <- fread(here('data', 'metadata', 'shen2018', 'sample_metadata_processed_shen2018.csv'))

######################################
# ---- plot estimated proportions ----
######################################

# ## Run once
# all_samples_pi.df <- data.frame()
# all_samples_pi_nonblood.df <- data.frame()
# 
# for (i in seq_len(nrow(md_samples))) {
# 
#   sample <- md_samples$sample[i]
#   group <- md_samples$group[i]
# 
#   sample_pi.df <- fread(here(read_dir, paste0('fitted_pi_', sample,'.csv'))) |>
#     select(cell_type, mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
#     mutate(sample = sample, group = group) |>
#     relocate(sample, group)
#   
#   sample_pi_nonblood.df <- sample_pi.df |>
#     dplyr::slice(8:n()) |> # remove blood cells
#     group_by(sample, group) |>
#     summarise(total = sum(mean)) |>
#     ungroup()
# 
#   all_samples_pi.df <- all_samples_pi.df |> rbind(sample_pi.df)
#   all_samples_pi_nonblood.df <- all_samples_pi_nonblood.df |> rbind(sample_pi_nonblood.df)
#   
#   cat(sample, ' ')
# }
# fwrite(all_samples_pi.df, here(write_dir, paste0('estimated_pi.csv.gz')))
# fwrite(all_samples_pi_nonblood.df, here(write_dir, paste0('estimated_pi_nonblood.csv.gz')))

all_samples_pi.df <- fread(here(write_dir, paste0('estimated_pi.csv.gz'))) |>
  left_join(md_samples, by = c('sample', 'group')) |>
  mutate(group = relevel(factor(group), ref = 'Control'), 
         cell_type = factor(cell_type, levels = unique(cell_type))) 

all_samples_pi.df |>
  ggplot(aes(group, mean)) +
  geom_jitter(size = 0.5, color = 'steelblue') +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~ cell_type, scales = 'free_y') +
  xlab('Patient group') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Estimated cell type proportions')
ggsave(here(plot_dir, paste0('boxplot_deconvProp_vs_patientGroup.png')), width = 15, height = 12)


all_samples_pi.df |>
  mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type)))) |>
  ggplot(aes(mean, cell_type)) +
  geom_jitter(size = 0.5, color = 'steelblue') +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~ group, scales = 'free_y') +
  ylab('Patient group') +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Estimated cell type proportions')
ggsave(here(plot_dir, paste0('boxplot_cellType_vs_deconvProp.png')), width = 12, height = 15)


## Try remove the low coverage samples
filtered_samples_pi.df <- all_samples_pi.df |>
  filter(million_reads >= 10)
  # filter(cts_total_reads >= 0.5*nrow(ref_cts.se))

filtered_samples_pi.df %>%
  ggplot(aes(group, mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`)) +
  geom_jitter(size = 0.5, color = 'steelblue') +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~ cell_type, scales = 'free_y') +
  xlab('Patient group') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Estimated cell type proportions')
ggsave(here(plot_dir, paste0('boxplot_deconvProp_vs_patientGroup_filteredLowCoverage.png')), width = 15, height = 12)


filtered_samples_pi.df %>%
  mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type)))) |>
  ggplot(aes(mean, cell_type)) +
  geom_jitter(size = 0.5, color = 'steelblue') +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~ group, scales = 'free_y') +
  ylab('Patient group') +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Estimated cell type proportions')
ggsave(here(plot_dir, paste0('boxplot_cellType_vs_deconvProp_filteredLowCoverage.png')), width = 12, height = 15)


### Proportion of non-blood tissues

all_samples_pi_nonblood.df <- fread(here(write_dir, paste0('estimated_pi_nonblood.csv.gz')))
all_samples_pi_nonblood.df |> 
  mutate(group = relevel(factor(group), 'Control')) |>
  ggplot(aes(group, total)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(size = 0.5, color = 'grey20') +
  xlab('Patient group') +
  theme_classic() +
  ylab('Estimated non-blood cell type proportions')
ggsave(here(plot_dir, paste0('violin_deconvNonBloodProp_vs_cellType.png')), width = 6, height = 5)


######################################
# ---- plot estimated w_mu ----
######################################

# ### RUN ONCE
# all_samples_w_mu.df <- data.frame(sample = NULL, group = NULL, million_reads = NULL, mean = NULL,
#                                   `2.5%` = NULL, `25%` = NULL, `50%` = NULL, `75%` = NULL, `97.5%` = NULL)
# 
# for (i in seq_len(nrow(md_samples))) {
# 
#   sample <- md_samples$sample[i]
#   group <- md_samples$group[i]
#   bam_file <- md_samples$bamdir[i]
#   million_reads <- md_samples$million_reads[i]
# 
#   sample_w_mu.df <- fread(here(read_dir, paste0('fitted_w_mu_', sample,'.csv'))) |>
#     select(mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
#     mutate(sample = sample, group = group, million_reads = million_reads, parameter = paste0('w_mu_', c(1:(n()-1), 0))) |>
#     relocate(sample, group, million_reads, parameter)
# 
#   all_samples_w_mu.df <- all_samples_w_mu.df |> rbind(sample_w_mu.df)
# 
#   cat(sample, ' ')
# }
# fwrite(all_samples_w_mu.df, here(write_dir, paste0('estimated_w_mu.csv.gz')))


all_samples_w_mu.df <- fread(here(write_dir, paste0('estimated_w_mu.csv.gz')))
new_labels <- c('w_mu_0' = expression(w[0]^mu),
                'w_mu_1' = expression(w[1]^mu),
                'w_mu_2' = expression(w[2]^mu),
                'w_mu_3' = expression(w[3]^mu),
                'w_mu_4' = expression(w[4]^mu)
)
all_samples_w_mu.df |>
  mutate(group = relevel(factor(group), ref = 'Control'),
         parameter = factor(parameter, labels = new_labels)) |> 
  ggplot(aes(million_reads, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = group), linewidth = 0.5, alpha = 0.5) + 
  geom_point(aes(color = group), size = 1) +
  geom_smooth(method = 'lm', color = 'tomato4', linetype = 'dashed') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  facet_wrap(~ parameter, labeller = label_parsed) +
  ggsci::scale_color_simpsons() +
  theme_classic() +
  xlab('Total coverage (million reads)') +
  ylab('Estimate')
ggsave(here(plot_dir, paste0('point&interval_w_mu_vs_coverage.png')))

######################################
# ---- plot estimated w_sigma ----
######################################

### RUN ONCE
all_samples_w_sigma.df <- data.frame(sample = NULL, group = NULL, million_reads = NULL, mean = NULL,
                                  `2.5%` = NULL, `25%` = NULL, `50%` = NULL, `75%` = NULL, `97.5%` = NULL)

for (i in seq_len(nrow(md_samples))) {

  sample <- md_samples$sample[i]
  group <- md_samples$group[i]
  bam_file <- md_samples$bamdir[i]
  million_reads <- md_samples$million_reads[i]
  
  sample_w_sigma.df <- fread(here(read_dir, paste0('fitted_w_sigma_', sample,'.csv'))) |>
    select(mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group, million_reads = million_reads, parameter = paste0('w_sigma_', c(1:(n()-1), 0))) |>
    relocate(sample, group, million_reads, parameter)

  all_samples_w_sigma.df <- all_samples_w_sigma.df |> rbind(sample_w_sigma.df)

  cat(sample, ' ')
}
fwrite(all_samples_w_sigma.df, here(write_dir, paste0('estimated_w_sigma.csv.gz')))


all_samples_w_sigma.df <- fread(here(write_dir, paste0('estimated_w_sigma.csv.gz')))
new_labels <- c('w_sigma_0' = expression(w[0]^sigma),
                'w_sigma_1' = expression(w[1]^sigma),
                'w_sigma_2' = expression(w[2]^sigma),
                'w_sigma_3' = expression(w[3]^sigma),
                'w_sigma_4' = expression(w[4]^sigma)
)
all_samples_w_sigma.df |>
  mutate(group = relevel(factor(group), ref = 'Control'),
         parameter = factor(parameter, labels = new_labels)) |> 
  ggplot(aes(million_reads, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = group), linewidth = 0.5, alpha = 0.5) + 
  geom_point(aes(color = group), size = 1) +
  geom_smooth(method = 'lm', color = 'tomato4', linetype = 'dashed') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  facet_wrap(~ parameter, labeller = label_parsed) +
  ggsci::scale_color_simpsons() +
  theme_classic() +
  xlab('Total coverage (million reads)') +
  ylab('Estimate')
ggsave(here(plot_dir, paste0('point&interval_w_sigma_vs_coverage.png')), width = 8, height = 7)

