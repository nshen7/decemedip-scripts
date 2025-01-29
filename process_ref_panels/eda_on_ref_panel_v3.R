source('code/SETPATHS.R')

plot_dir <- here('plots', 'metadata', 'references', 'MethAtlas_refv3')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v3.rds'))
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v3.rds'))

dev.off()
png(here(plot_dir, 'heatmap_beta_on_cts_sites.png'), height = 600, width = 1200, res = 100)
pheatmap::pheatmap(t(assays(ref_cts.se)$beta), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)
dev.off()

png(here(plot_dir, 'heatmap_beta_on_anc_sites.png'), height = 600, width = 1200, res = 100)
pheatmap::pheatmap(t(assays(ref_anc.se)$beta), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)
dev.off()

assays(ref_cts.se)$beta |>
  pivot_longer(cols = everything(), names_to = 'Cell type', values_to = 'Beta') |>
  mutate(`Cell type` = factor(`Cell type`, levels = rev(colnames(ref_cts.se)))) |>
  group_by(`Cell type`) |>
  mutate(`Average beta` = mean(Beta)) |>
  ggplot(aes(Beta, `Cell type`)) +
  geom_jitter(size = 0.3, color = 'steelblue') +
  geom_boxplot(alpha = 0.5) +
  geom_point(aes(`Average beta`, color = 'Mean')) + 
  theme_classic()
ggsave(here(plot_dir, 'boxplot_beta_on_cts_sites.png'), width = 7, height = 6)


###################################
####### Correlation heatmap #######
###################################
## ---- util ----

plotCorr <- function(corr.mat) {
  p <- corr.mat |>
    as.matrix() |>
    as.data.frame() |>
    mutate(cell_type_1 = colnames(ref_cts.se)) |>
    pivot_longer(cols = -c('cell_type_1'), values_to = 'correlation', names_to = 'cell_type_2') |>
    mutate(cell_type_1 = factor(cell_type_1, levels = rev(colnames(ref_cts.se))),
           cell_type_2 = factor(cell_type_2, levels = rev(colnames(ref_cts.se)))
    ) |>
    ggplot(aes(cell_type_1, cell_type_2, fill = correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(-0.25, 1), 
                         space = "Lab", 
                         name = "Correlation\nbetween\ncell types") +
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    xlab(NULL) + ylab(NULL) + 
    coord_fixed()
  return(p)
}


## ---- CTS only ----

corr.mat <- cor(assays(ref_cts.se)[[1]], method = 'pearson') 

plotCorr(corr.mat)
ggsave(here(plot_dir, 'heatmap_sample_corr_on_cts_sites.png'), height = 10, width = 12)

## ---- CTS + anchor ----

corr.mat <- cor(rbind(assays(ref_cts.se)[[1]], assays(ref_anc.se)[[1]]), method = 'pearson') 

plotCorr(corr.mat)
ggsave(here(plot_dir, 'heatmap_sample_corr_on_cts&anc_sites.png'), height = 10, width = 12)


## Compare distribution of correlation with ref panel v2
corr_v3 <- reshape2::melt(corr.mat)

ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v2.rds'))
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v2.rds'))
corr.mat_v2 <- cor(rbind(assays(ref_cts.se)[[1]], assays(ref_anc.se)[[1]]), method = 'pearson') 
corr_v2 <- reshape2::melt(corr.mat_v2)

ggplot() +
  geom_point(aes(corr_v2$value, corr_v3$value)) +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic()


###########################################
###### Overlap with cancer biomarker ######
###########################################

# dmrs.gr <- readxl::read_xlsx(here('data', 'metadata', 'Ibrahim2022_TCGA_tumor_biomarkers', 'SuppTable5_all_DMB_info.xlsx')) |>
dmrs.gr <- readxl::read_xlsx(here('data', 'metadata', 'Ibrahim2022_TCGA_tumor_biomarkers', 'SuppTable4_all_DMR_info.xlsx')) |>
  select(-c(`...1`)) |>
  rename(V1 = 'dmrName') |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
seqlevelsStyle(dmrs.gr) <- "UCSC"

## check which CTS regions overlaps with hyper-methylated cancer biomarkers 
dmrs_hyper.gr <- dmrs.gr |> subset(value > 0)
hits <- findOverlaps(granges(ref_cts.se), dmrs_hyper.gr)
ref_cts_olap.gr <- granges(ref_cts.se)[queryHits(hits)]
dmrs_hyper_olap.gr <- dmrs_hyper.gr[subjectHits(hits)]

combined.df <- cbind(data.frame(ref_cts_olap.gr), mcols(dmrs_hyper_olap.gr))
table(combined.df$label)

## Effect size of the tumor DMRs that overlap with CTS sites
combined.df |>
  mutate(label = gsub('(.*) hypermethylated', '\\1', label)) |>
  mutate(is_prostate = ifelse(label == 'Prostate', TRUE, FALSE)) |>
  ggplot(aes(label, value, color = is_prostate)) +
  geom_boxplot() +
  geom_jitter(size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(~tumor_type) +
  ylab('Effect size') +
  xlab('Cell type')

## Number of CTS sites that overlap with tumor DMRs
olap.df <- combined.df |>
  # right_join(as.data.frame(mcols(ref_cts.se)), by = c('probe', 'label', 'pos', 'n_cpgs_100bp')) |>
  mutate(label = gsub('(.*) hypermethylated', '\\1', label)) |>
  group_by(probe, label) |>
  summarise(n_olap_dmrs = sum(!is.na(dmrName)))
table(olap.df$label)
counts.df <- olap.df |> group_by(label) |> summarise(n = n())
olap.df |>
  ggplot(aes(label, n_olap_dmrs)) +
  geom_violin(width = 4, scale = 'count') +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_text(data = counts.df, aes(x = label, y = 15, label = paste0('n=', n))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab('# cancer types that has DMRs overlapping with CTS site') +
  xlab('Cell type') 

