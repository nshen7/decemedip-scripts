source('code/SETPATHS.R')

plot_dir <- here('plots', 'metadata', 'references', 'MethAtlas_refv2')

ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v2.rds'))
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v2.rds'))

dev.off()
png(here(plot_dir, 'heatmap_beta_on_cts_sites.png'), height = 600, width = 1200, res = 100)
pheatmap::pheatmap(t(assays(ref_cts.se)$beta), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)
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


# ###########################################
# ###### Overlap with cancer biomarker ######
# ###########################################
# 
# # dmrs.gr <- readxl::read_xlsx(here('data', 'metadata', 'Ibrahim2022_TCGA_tumor_biomarkers', 'SuppTable5_all_DMB_info.xlsx')) |>
# dmrs.gr <- readxl::read_xlsx(here('data', 'metadata', 'Ibrahim2022_TCGA_tumor_biomarkers', 'SuppTable4_all_DMR_info.xlsx')) |>
#   select(-c(`...1`)) |>
#   rename(V1 = 'dmrName') |>
#   makeGRangesFromDataFrame(keep.extra.columns = TRUE)
# seqlevelsStyle(dmrs.gr) <- "UCSC"
# 
# ## check which CTS regions overlaps with RCC-specific biomarkers
# dmrs_rcc_hyper.gr <- dmrs.gr |> subset(tumor_type == 'KIRC' & value > 0)
# 
# idx_olap_cts <- which(countOverlaps(granges(ref_cts.se), dmrs_rcc_hyper.gr) > 0)
# 
# ## Type of CTS regions overlap with KIRC tumour specific hyper-DMRs
# table(granges(ref_cts.se)$label[idx_olap_cts]) |> sort(decreasing = TRUE)
# 
# ## Heatmap of the reference regions that overlap with KIRC tumour specific hyper-DMRs
# beta.mat <- assays(ref_cts.se[idx_olap_cts, ])$beta
# png(here(plot_dir, 'heatmap_beta_on_cts_sites_olap_KIRC_DMRs.png'), height = 500, width = 700, res = 100)
# pheatmap::pheatmap(t(beta.mat), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE)
# dev.off()
# 
# beta.mat |>
#   pivot_longer(cols = everything(), names_to = 'Cell type', values_to = 'Beta') |>
#   mutate(`Cell type` = factor(`Cell type`, levels = rev(colnames(beta.mat)))) |>
#   group_by(`Cell type`) |>
#   mutate(`Average beta` = mean(Beta)) |>
#   ggplot(aes(Beta, `Cell type`)) +
#   geom_jitter(size = 0.3, color = 'steelblue') +
#   geom_boxplot(alpha = 0.5) +
#   geom_point(aes(`Average beta`, color = 'Mean')) + 
#   geom_vline(xintercept = 0.25, color = 'red', linetype = 'dashed') + 
#   theme_classic()
# ggsave(here(plot_dir, 'boxplot_beta_on_cts_sites_olap_KIRC_DMRs.png'), width = 7, height = 6)
# 
# dmrs_olap.gr <- GRanges()
# for (this_label in unique(granges(ref_cts.se)$label)) {
#   idx_olap_dmr <- which(countOverlaps(dmrs_rcc_hyper.gr, 
#                                       granges(ref_cts.se) |> subset(label == this_label)
#                                       ) > 0)
#   temp <- dmrs_rcc_hyper.gr[idx_olap_dmr]
#   temp$overlapLabel <- this_label
#   dmrs_olap.gr <- c(dmrs_olap.gr, temp)
# }
# 
# table(dmrs_olap.gr$overlapLabel)
# dmrs_olap.gr |>
#   as.data.frame() |>
#   ggplot(aes(value, overlapLabel)) +
#   geom_point(color = 'steelblue')
