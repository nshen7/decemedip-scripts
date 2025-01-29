## What's changed in V3 from V4:
##  1. Removed probes that are in tumor-related DMRs

source('code/SETPATHS.R')
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(rtracklayer)

read_dir <- write_dir <- here('data', 'metadata', 'references', 'MethAtlas')

## 450K array probe annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## CpG coordinates in hg19
if (! 'CPG_COORDS' %in% ls())
  CPG_COORDS <- fread(here('data', 'interim', 'cpg_coords_in_hg19', 'cpg_coords_in_hg19.csv.gz')) |> makeGRangesFromDataFrame()

# ---- Load pan-cancer DMR atlas ----

tumor_dmrs.gr <- readxl::read_xlsx(here('data', 'metadata', 'Ibrahim2022_TCGA_tumor_biomarkers', 'SuppTable4_all_DMR_info.xlsx')) |>
  select(-c(`...1`)) |>
  rename(V1 = 'dmrName') |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
seqlevelsStyle(tumor_dmrs.gr) <- "UCSC"

# ---- Load atlas ----

ref_atlas.df <- fread(here(read_dir, 'reference_atlas.csv')) |> ## only for obtaining colnames
  rename('CpGs' = 'probe')
full_atlas_0.df <- fread(here(read_dir, 'full_atlas.csv.gz'), header = F, col.names = colnames(ref_atlas.df)) |>
  na.omit() |> 
  as.data.frame()

## Remove probes that locates in cancer-associated DMRs
probe_coords.gr <- anno450k[match(full_atlas_0.df$probe, rownames(anno450k)), c("chr", "pos")] |> 
  as.data.frame() |>
  mutate(start = pos, end = pos) |>
  makeGRangesFromDataFrame()
idx_olap_with_tumor_dmrs <- (countOverlaps(probe_coords.gr, tumor_dmrs.gr) > 0)
full_atlas.df <- full_atlas_0.df[!idx_olap_with_tumor_dmrs, ]

## Extract the beta matrix
full_beta.df <- full_atlas.df[, -1]
rownames(full_beta.df) <- full_atlas.df$probe

## Obtain diff of target and MAX background beta for each cell type
full_diff_max.df <- full_atlas.df |> select('probe')
for (i in seq_len(ncol(full_beta.df)))  {
  label <- colnames(full_beta.df)[i]
  beta_tg <- full_beta.df[[ i]]
  beta_bg <- full_beta.df[, -i] |> as.matrix() |> rowMaxs()
  full_diff_max.df[[label]] <- beta_tg - beta_bg
  print(i)
}

## Obtain diff of target and MEAN background beta for each cell type
full_diff_mean.df <- full_atlas.df |> select('probe')
for (i in seq_len(ncol(full_beta.df)))  {
  label <- colnames(full_beta.df)[i]
  beta_tg <- full_beta.df[[ i]]
  beta_bg <- full_beta.df[, -i] |> as.matrix() |> rowMeans()
  full_diff_mean.df[[label]] <- beta_tg - beta_bg
  print(i)
}

## Check number of probes hypermethylated in targeted cell type
colSums(full_diff_max.df[, -1] > 0) |> quantile()
 #   0%   25%   50%   75%  100% 
 # 2285  6173 11112 16163 90992 

colSums(full_diff_max.df[, -1] > 0 & full_diff_mean.df[, -1] > 0.18) |> sort(decreasing = TRUE)
   # Erythrocyte_progenitors     Colon_epithelial_cells                Hepatocytes           Cortical_neurons 
   #                   23251                      15667                       8870                       7539 
   #         CD8T-cells_EPIC                    Thyroid      Pancreatic_duct_cells    Pancreatic_acinar_cells 
   #                    7024                       3657                       3485                       3232 
   #            B-cells_EPIC      Pancreatic_beta_cells             Monocytes_EPIC                     Breast 
   #                    2466                       2182                       1969                       1561 
   #                Prostate            CD4T-cells_EPIC                Left_atrium              NK-cells_EPIC 
   #                    1498                       1320                        947                        893 
   #              Adipocytes           Neutrophils_EPIC                     Kidney              Uterus_cervix 
   #                     767                        717                        706                        569 
   #                 Bladder                 Lung_cells       Head_and_neck_larynx Vascular_endothelial_cells 
   #                     533                        490                        411                        340 
   #                Upper_GI 
   #                     105


# ---- assemble ref atlas ----
K_PER_CT <- 100

## Select cell type specific marker probes that are top-ranked in diff between 
## target cell type and average background
getCTSHyperProbes <- function(label, top_k) {
  idx <- which(full_diff_max.df[[label]] > 0)
  df <- data.frame(probe = full_atlas.df[['probe']][idx],
                   label = paste0(label, ' hypermethylated'),
                   diff_mean = full_diff_mean.df[[label]][idx],
                   diff_max  = full_diff_max.df[[label]][idx]) |>
    group_by(label) |>
    slice_max(order_by = diff_mean, n = top_k)
  return(df)
}

cts_probes.temp <- map_dfr(colnames(full_diff_max.df)[-1],
                           getCTSHyperProbes, top_k = K_PER_CT) |>
  select(probe, label)
ref_cts_beta.temp <- full_beta.df[match(cts_probes.temp$probe, full_atlas.df$probe), ]

## The complete reference panel
ref_cts.df <- cbind(
  cts_probes.temp, 
  as.data.frame(anno450k[match(cts_probes.temp$probe, rownames(anno450k)), c("chr", "pos")])
) |>
  mutate(start = pos, end = pos) |>
  as.data.frame()
rownames(ref_cts.df) <- ref_cts.df$probe

ref_cts_beta.df <- full_beta.df[match(cts_probes.temp$probe, rownames(full_beta.df)), ]

## Sanity check
sum(duplicated(ref_cts.df)) # = 0; no duplicated probes
all(rownames(ref_cts.df) == ref_cts.df$probe) # = TRUE; all probes correctly mapped
all(rownames(ref_cts_beta.df) == ref_cts.df$probe) # = TRUE; all probes correctly mapped

## Write out
ref_cts.gr <- makeGRangesFromDataFrame(ref_cts.df, keep.extra.columns = T)
ref_cts.gr$n_cpgs_100bp <- countOverlaps(ref_cts.gr |> resize(width = 100, fix = 'center'), CPG_COORDS)
ref_cts.se <- SummarizedExperiment(rowData = ref_cts.gr,
                                   assays = list(beta = ref_cts_beta.df))
saveRDS(ref_cts.se, here(write_dir, 'summarizedexperiment_cell_type_specific_atlas_v4.rds'))

table(ref_cts.gr$n_cpgs_100bp)
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
# 333 333 257 237 239 221 184 167 134 126  80  72  49  32  24   5   4   3 

# ---- Anchor sites (i.e., all-tissue M/U sites) from MethAtlas ----

K_PER_DIRECTION <- 500

# Find blocks where WGBS samples are all methylated (beta > 1-margin) or all unmethylated (beta < margin)
margin <- 0.1
idx_all_m <- which(rowMeans(full_beta.df > 1-margin) >= 1)
idx_all_u <- which(rowMeans(full_beta.df < margin) >= 1)

buildSE <- function(idx, label) {
  
  subset.df <- full_atlas.df[idx, ]
  beta.df   <- subset.df[, -1]
  
  ## Compute average beta values per marker
  avg_beta  <- rowMeans(beta.df)
  if (label == 'All-tissue U') avg_beta_rank <- rank(avg_beta) else avg_beta_rank <- rank(-avg_beta)
  
  probe_list <- unlist(subset.df[, 1])
  coords     <- anno450k[match(probe_list, rownames(anno450k)), c("chr", "pos")]
  md.df      <- data.frame(probe = probe_list, coords) %>%
    mutate(label = label, margin = margin, 
           start = pos, end = pos,
           avg_beta = avg_beta,
           avg_beta_rank = avg_beta_rank) 
  rownames(md.df) <- NULL
  
  se <- SummarizedExperiment(rowData = makeGRangesFromDataFrame(md.df, keep.extra.columns = TRUE),
                             assays = list(beta = beta.df))
  return(se)
}

u_atlas.se <- buildSE(idx = idx_all_u, label = 'All-tissue U')
rowData(u_atlas.se)$n_cpgs_100bp <- countOverlaps(granges(u_atlas.se) |> resize(width = 100, fix = 'center'), CPG_COORDS)
m_atlas.se <- buildSE(idx = idx_all_m, label = 'All-tissue M')
rowData(m_atlas.se)$n_cpgs_100bp <- countOverlaps(granges(m_atlas.se) |> resize(width = 100, fix = 'center'), CPG_COORDS)

table(rowData(u_atlas.se)$n_cpgs_100bp)
table(rowData(m_atlas.se)$n_cpgs_100bp)

## Sample the anchor probes to follow similar distribution of CpG density as the CTS probes
target_freqs_0 <- data.frame(table(ref_cts.gr$n_cpgs_100bp) / length(ref_cts.gr)) |> 
  rename('Var1' = 'n_cpgs_100bp', 'Freq' = 'prob') |> 
  mutate(n_cpgs_100bp = as.integer(n_cpgs_100bp))

sampleRowsAccordingToCGDensity <- function(se, k, seed = 2022) {
  
  set.seed(seed)
  
  # Freqs in u_atlas.se/m_atlas.se region set
  original_freqs <- data.frame(table(rowData(se)$n_cpgs_100bp)) |> 
    rename('Var1' = 'n_cpgs_100bp', 'Freq' = 'count') |> 
    mutate(n_cpgs_100bp = as.integer(n_cpgs_100bp))
  
  target_freqs <- target_freqs_0 |> 
    left_join(original_freqs, by = 'n_cpgs_100bp') |>
    mutate(adjusted_prob = prob / count)
  sampling_prob <- rowData(se) |> as.data.frame() |> left_join(target_freqs, by = 'n_cpgs_100bp') |> pull(adjusted_prob)
  sampling_prob[is.na(sampling_prob)] <- 0
  sampled.se <- se[sample(x = 1:nrow(se), size = k, prob = sampling_prob, replace = FALSE), ]
  
  return(sampled.se)
}

u_atlas_sampled.se <- sampleRowsAccordingToCGDensity(se = u_atlas.se, k = K_PER_DIRECTION)
m_atlas_sampled.se <- sampleRowsAccordingToCGDensity(se = m_atlas.se, k = K_PER_DIRECTION)

table(rowData(u_atlas_sampled.se)$n_cpgs_100bp)
table(rowData(m_atlas_sampled.se)$n_cpgs_100bp)


ref_anc.se <- rbind(u_atlas_sampled.se, m_atlas_sampled.se)
saveRDS(ref_anc.se, here(write_dir, 'summarizedexperiment_all_tissue_um_atlas_v4.rds'))

table(granges(ref_anc.se)$label)
# All-tissue M All-tissue U 
#          500          500



# ---- lift over the reference probes from hg19 to hg38 ----

## Cell type specific reference
ref_cts.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_cell_type_specific_atlas_v4.rds'))
## Anchor region reference
ref_anc.se <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                           'summarizedexperiment_all_tissue_um_atlas_v4.rds'))

# Download and import the chain file
chain_file <- here('data', 'raw', 'hg19ToHg38.over.chain')
chain <- import.chain(chain_file)

# Lift over the reference probes
ref_cts_lifted.list <- liftOver(granges(ref_cts.se), chain)
idx_unmapped_cts <- which(elementNROWS(ref_cts_lifted.list) == 0)

ref_anc_lifted.list <- liftOver(granges(ref_anc.se), chain)
idx_unmapped_anc <- which(elementNROWS(ref_anc_lifted.list) == 0)

# Remove unmapped probes
ref_cts_2.se <- ref_cts.se; ref_anc_2.se <- ref_anc.se
if (length(idx_unmapped_cts) > 0) ref_cts_2.se <- ref_cts_2.se[-idx_unmapped_cts, ]
if (length(idx_unmapped_anc) > 0) ref_anc_2.se <- ref_anc_2.se[-idx_unmapped_ans, ]

# Output lifted reference panel
ref_cts_lifted.se <- SummarizedExperiment(
  rowData = unlist(ref_cts_lifted.list),
  assays = list(beta = assays(ref_cts_2.se)$beta)
)
genome(ref_cts_lifted.se) <- "hg38"
saveRDS(ref_cts_lifted.se, here(write_dir, 'summarizedexperiment_cell_type_specific_atlas_hg38_v4.rds'))

ref_anc_lifted.se <- SummarizedExperiment(
  rowData = unlist(ref_anc_lifted.list),
  assays = list(beta = assays(ref_anc_2.se)$beta)
)
genome(ref_anc_lifted.se) <- "hg38"
saveRDS(ref_anc_lifted.se, here(write_dir, 'summarizedexperiment_all_tissue_um_atlas_hg38_v4.rds'))


