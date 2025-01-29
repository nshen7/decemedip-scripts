source('code/SETPATHS.R')
.libPaths('/home/nshen7/R/cfmedip_extra_libs')
library(argparse)

read_dir <- here('data', 'interim', 'case_studies', 'baca2023', '01_summarize_results_refv4')
plot_dir <- here('plots', 'case_studies', 'baca2023', '02_make_plots_refv4')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

md_samples <- fread(here(read_dir, 'processed_metadata.csv'))
GROUPS <- unique(md_samples$group)
COLORS <- colorRampPalette(grafify:::graf_palettes$okabe_ito[1:8])(length(GROUPS))
names(COLORS) <- GROUPS

## Cancer types that has corresponded cell types in reference atlas:
## Lung, Colorectal, Renal, Breast, Esophageal, Prostate, Bladder, Hepatocellular
CANCERTYPES <- c('Lung cancer', 'Colorectal cancer',     'Breast cancer', 'Prostate cancer', 'Hepatocellular cancer')
RELEVANTCTS <- c('Lung_cells',  'Colon_epithelial_cells','Breast',        'Prostate',        'Hepatocytes')


######################################
# ---- estimated proportions ----
######################################


all_samples_pi.df <- fread(here(read_dir, paste0('estimated_pi.csv.gz'))) |>
  left_join(md_samples, by = c('sample', 'group')) |>
  mutate(group = relevel(factor(group), ref = 'Healthy'), 
         cell_type = factor(cell_type, levels = unique(cell_type))) |>
  filter(! cell_type %in% unique(cell_type)[1:7])  # remove blood cells


all_samples_pi.df |>
  ggplot(aes(group, mean)) +
  geom_jitter(size = 0.5, color = 'steelblue') +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~ cell_type) +
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
  xlab('Estimated cell type proportions')
ggsave(here(plot_dir, paste0('boxplot_cellType_vs_deconvProp.png')), width = 12, height = 15)


for (i in seq_along(CANCERTYPES)) {
  # this_cancer_type <- CANCERTYPES[i]
  this_cell_type   <- RELEVANTCTS[i]
  
  cell_type_name <- gsub('_', ' ', this_cell_type)
  if (cell_type_name %in% c('Breast', 'Prostate')) cell_type_name <- paste(cell_type_name, 'cells')
  
  plot.df <- all_samples_pi.df |>
    filter(group %in% c('Healthy', CANCERTYPES)) |>
    filter(cell_type == this_cell_type) |>
    # arrange(desc(mean)) |>
    mutate(group = recode(group, 'Hepatocellular cancer' = 'Hepato- cellular cancer')) |> 
    mutate(sample = factor(sample, levels = sample))
  
  COLORS2 <- COLORS
  names(COLORS2)[which(names(COLORS2) == "Hepatocellular cancer")] <- 'Hepato- cellular cancer'
  
  plot.df |>
    ggplot(aes(x = sample, y = mean, color = group)) +
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`),
                   position = position_dodge2(width = 0.3),
                   linewidth = 1, alpha = 0.3) +
    geom_linerange(aes(ymin = `25%`, ymax = `75%`),
                   position = position_dodge2(width = 0.3),
                   linewidth = 1, alpha = 1) +
    geom_point(position = position_dodge2(width = 0.3),
               size = 1, color = 'grey20', shape = 19) +
    scale_color_manual(values = COLORS2) +
    theme_classic() +
    facet_grid(~ group, scales = "free_x", space = "free", switch = "x", labeller = label_wrap_gen(width=10)) +  # Split x-axis by groups
    scale_x_discrete(labels = NULL) +  # Remove sample names from the x-axis
    ggtitle(cell_type_name) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      strip.text.x = element_text(size = 8, face = "bold"),  # Style for group labels
      strip.text.y = element_text(size = 8),  # Style for group labels
      # strip.placement = "outside",  # Move facet title below the x-axis
      strip.background = element_blank(),  # Remove box around facet title
      panel.spacing = unit(0.5, "lines"),  # Reduce space between facets
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.line.x = element_blank(),   # Remove x-axis line
      axis.line.y = element_blank(),   # Remove x-axis line
      legend.position = "none"         # Remove legend
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1)) +
    xlab(NULL) +
    ylab('Estimated proportions\nwith credible intervals')
  ggsave(here(plot_dir, paste0('interval_deconvProp_vs_patientGroup_cellTypeIs', this_cell_type, '.png')),
         width = 10, height = 5)
}


# 
# ## Try remove the low ctDNA samples
# filtered_samples_pi.df <- all_samples_pi.df |>
#   filter(ctDNA >= 0.1 | group == 'Healthy') |>
#   group_by(group) 
# 
# filtered_samples_pi.df %>%
#   ggplot(aes(group, mean)) +
#   geom_jitter(size = 0.5, color = 'steelblue') +
#   geom_boxplot(alpha = 0.5) +
#   facet_wrap(~ cell_type) +
#   xlab('Patient group') +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ylab('Estimated cell type proportions')
# ggsave(here(plot_dir, paste0('boxplot_deconvProp_vs_patientGroup_filteredCtDNA.png')), width = 15, height = 12)
# 
# all_samples_pi.df |>
#   mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type)))) |>
#   ggplot(aes(mean, cell_type)) +
#   geom_jitter(size = 0.5, color = 'steelblue') +
#   geom_boxplot(alpha = 0.5) +
#   facet_wrap(~ group, scales = 'free_y') +
#   ylab('Patient group') +
#   theme_classic() +
#   xlab('Estimated cell type proportions')
# ggsave(here(plot_dir, paste0('boxplot_cellType_vs_deconvProp_filteredCtDNA.png')), width = 12, height = 15)
# 
# for (i in seq_along(CANCERTYPES)) {
#   this_cancer_type <- CANCERTYPES[i]
#   this_cell_type   <- RELEVANTCTS[i]
#   
#   cell_type_name <- gsub('_', ' ', this_cell_type)
#   if (cell_type_name %in% c('Breast', 'Prostate')) cell_type_name <- paste(cell_type_name, 'cells')
#   
#   filtered_samples_pi.df |>
#     filter(group %in% c('Healthy', this_cancer_type)) |>
#     filter(cell_type == this_cell_type) |>
#     arrange(desc(mean)) |>
#     mutate(sample = factor(sample, levels = sample)) |>
#     ggplot(aes(x = sample, y = mean, color = group)) +
#     geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`),
#                    position = position_dodge2(width = 0.3),
#                    linewidth = 1, alpha = 0.3) +
#     geom_linerange(aes(ymin = `25%`, ymax = `75%`),
#                    position = position_dodge2(width = 0.3),
#                    linewidth = 1, alpha = 1) +
#     geom_point(position = position_dodge2(width = 0.3),
#                size = 1, color = 'grey20', shape = 19) +
#     scale_color_manual(values = COLORS) +
#     theme_classic() +
#     facet_grid(~ group, scales = "free_x", space = "free", switch = "x", labeller = label_wrap_gen(width=20)) +  # Split x-axis by groups
#     scale_x_discrete(labels = NULL) +  # Remove sample names from the x-axis
#     ggtitle(cell_type_name) +
#     theme(
#       plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
#       strip.text.x = element_text(size = 8, face = "bold"),  # Style for group labels
#       strip.text.y = element_text(size = 8),  # Style for group labels
#       # strip.placement = "outside",  # Move facet title below the x-axis
#       strip.background = element_blank(),  # Remove box around facet title
#       panel.spacing = unit(0.5, "lines"),  # Reduce space between facets
#       axis.text.x = element_blank(),  # Remove x-axis text
#       axis.ticks.x = element_blank(),  # Remove x-axis ticks
#       axis.line.x = element_blank(),   # Remove x-axis line
#       axis.line.y = element_blank(),   # Remove x-axis line
#       legend.position = "none"         # Remove legend
#     ) + 
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
#                        limits = c(0, 1),
#                        breaks = c(0, 0.5, 1)) +
#     xlab(NULL) +
#     ylab('Estimated proportions\nwith credible intervals')
#   ggsave(here(plot_dir, paste0('interval_deconvProp_vs_patientGroup_cellTypeIs', this_cell_type, '_filteredCtDNA.png')),
#          width = 5, height = 2.5)
# }

### -------- Proportion of non-blood tissues ------
## Using all samples
all_samples_pi_nonblood.df <- fread(here(read_dir, paste0('estimated_pi_nonblood.csv.gz'))) |>
  left_join(md_samples |> select(sample, ctDNA), by = 'sample')  |>
  mutate(group = relevel(factor(group), 'Healthy'))

all_samples_pi_nonblood.df |> 
  ggplot(aes(group, total, color = group)) +
  geom_violin(aes(fill = group), alpha = 0.5, width = 2) +
  geom_jitter() +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  xlab('Patient group') +
  ylab('Estimated non-blood cell type proportions') +
  theme_classic() +
  theme(legend.position = "none")
ggsave(here(plot_dir, paste0('violin_deconvNonBloodProp_vs_cellType.png')), width = 7.5, height = 5)

# ## Try remove the low ctDNA samples
# filtered_samples_pi_nonblood.df <- all_samples_pi_nonblood.df |>
#   filter(ctDNA >= 0.1 | group == 'Healthy') 
# 
# filtered_samples_pi_nonblood.df |> 
#   ggplot(aes(group, total, color = group)) +
#   geom_violin(aes(fill = group), alpha = 0.5, width = 2) +
#   geom_jitter() +
#   scale_fill_manual(values = COLORS) +
#   scale_color_manual(values = COLORS) +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
#   xlab('Patient group') +
#   ylab('Estimated non-blood cell type proportions') +
#   theme_classic() +
#   theme(legend.position = "none")
# ggsave(here(plot_dir, paste0('violin_deconvNonBloodProp_vs_cellType_filteredCtDNA.png')), width = 7.5, height = 5)

# ---- ctDNA vs. respective tissue proportions ----

all_samples_pi.df |> 
  filter(group %in% CANCERTYPES) |>
  filter(paste(group, cell_type) %in% paste(CANCERTYPES, RELEVANTCTS)) |>
  filter(ctDNA > 0) |>
  ggplot(aes(ctDNA, mean, color = group)) + 
  geom_abline(intercept = 0, slope = 0.01, linetype = 'dashed', color = 'grey10') +
  geom_smooth(aes(group = group), method = 'lm', size = 1, se = FALSE) +
  geom_point() +
  scale_color_manual(values = COLORS) +
  xlab('ichorCNA ctDNA (%)') +
  ylab('Estimated non-blood cell type proportions') +
  labs(color = "Patient group") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  theme_classic() 
ggsave(here(plot_dir, paste0('point_deconvRespectiveProp_vs_ctDNA.png')), width = 7.5, height = 5)

  
# ---- ctDNA vs. non-blood proportions ----

all_samples_pi_nonblood.df |> 
  filter(group %in% CANCERTYPES) |>
  filter(ctDNA > 0) |>
  ggplot(aes(ctDNA, total, color = group)) + 
  geom_abline(intercept = 0, slope = 0.01, linetype = 'dashed', color = 'grey10') +
  geom_smooth(aes(group = group), method = 'lm', size = 1, se = FALSE) +
  geom_point() +
  scale_color_manual(values = COLORS) +
  xlab('ichorCNA ctDNA (%)') +
  ylab('Estimated non-blood cell type proportions') +
  labs(color = "Patient group") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  theme_classic() 
ggsave(here(plot_dir, paste0('point_deconvNonBloodProp_vs_ctDNA.png')), width = 7.5, height = 5)




######################################
# ---- estimated w_mu ----
######################################

all_samples_w_mu.df <- fread(here(write_dir, paste0('estimated_w_mu.csv.gz')))
new_labels <- c('w_mu_0' = expression(w[0]^mu),
                'w_mu_1' = expression(w[1]^mu),
                'w_mu_2' = expression(w[2]^mu),
                'w_mu_3' = expression(w[3]^mu),
                'w_mu_4' = expression(w[4]^mu)
)
all_samples_w_mu.df |>
  mutate(group = relevel(factor(group), ref = 'Healthy'),
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
# ---- estimated w_sigma ----
######################################

all_samples_w_sigma.df <- fread(here(write_dir, paste0('estimated_w_sigma.csv.gz')))
new_labels <- c('w_sigma_0' = expression(w[0]^sigma),
                'w_sigma_1' = expression(w[1]^sigma),
                'w_sigma_2' = expression(w[2]^sigma),
                'w_sigma_3' = expression(w[3]^sigma),
                'w_sigma_4' = expression(w[4]^sigma)
)
all_samples_w_sigma.df |>
  mutate(group = relevel(factor(group), ref = 'Healthy'),
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

