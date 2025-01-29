source('code/SETPATHS.R')
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

read_dir <- write_dir <- here('data', 'metadata', 'references', 'MethAtlas')

## 450K array probe annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## CpG coordinates in hg19
if (! 'CPG_COORDS' %in% ls())
  CPG_COORDS <- fread(here('data', 'interim', 'cpg_coords_in_hg19', 'cpg_coords_in_hg19.csv.gz')) |> makeGRangesFromDataFrame()

# ---- Load atlas ----
ref_atlas.df <- fread(here(read_dir, 'reference_atlas.csv')) ## only to obtain colnames
full_atlas.df <- fread(here(read_dir, 'full_atlas.csv.gz'), header = F, col.names = colnames(ref_atlas.df)) |>
  na.omit() |> 
  as.data.frame() 
full_beta.df <- full_atlas.df[, -1]
rownames(full_beta.df) <- full_atlas.df$CpGs

# ## Select cell type specific hyper-methylated probes (according to method described in MethAtlas paper)
# full_beta_normed.df <- map_dfc(full_beta.df, ~ .x / sum(.x))
# getCTSHyperProbes <- function(col, top_k) {
#   idx <- order(full_beta_normed.df[[col]], decreasing = T)[seq_len(top_k)]
#   return(data.frame(probe = full_atlas.df[['CpGs']][idx], rank = seq_len(top_k)))
# }
# cts_hyper_md.df <- map_dfr(colnames(full_beta.df), getCTSHyperProbes, top_k = 200)


## Obtain diff of target and MAX background beta for each cell type
full_diff_max.df <- full_atlas.df |> select('CpGs')
for (i in seq_len(ncol(full_beta.df)))  {
  label <- colnames(full_beta.df)[i]
  beta_tg <- full_beta.df[[ i]]
  beta_bg <- full_beta.df[, -i] |> as.matrix() |> rowMaxs()
  full_diff_max.df[[label]] <- beta_tg - beta_bg
  print(i)
}

## Obtain diff of target and MEAN background beta for each cell type
full_diff_mean.df <- full_atlas.df |> select('CpGs')
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
 # 2586  6675 12716 16670 94073

colSums(full_diff_mean.df[, -1] > 0.2) |> quantile()
 #   0%   25%   50%   75%  100% 
 # 1273  4834 12346 18253 28950 

colSums(full_diff_max.df[, -1] > 0 & full_diff_mean.df[, -1] > 0.2) |> quantile()
  #  0%   25%   50%   75%  100% 
  # 119   626  1563  3362 20014 

## Select cell type specific marker probes that are top-ranked in diff_max 
getCTSHyperProbes <- function(label, top_k) {
  idx <- order(full_diff_max.df[[label]], decreasing = T)[seq_len(top_k)]
  return(data.frame(CpGs = full_atlas.df[['CpGs']][idx],
                    diff_mean = full_diff_mean.df[[label]][idx],
                    diff_max  = full_diff_max.df[[label]][idx],
                    diff_max_rank  = seq_len(top_k), 
                    label = label))
}
cts_probes.df <- map_dfr(colnames(full_diff_max.df)[-1], getCTSHyperProbes, top_k = 2000) |>
  group_by(label) |>
  mutate(diff_mean_rank = rank(-diff_mean))
cts_atlas.df <- cbind(cts_probes.df, 
                      as.data.frame(anno450k[match(cts_probes.df$CpGs, rownames(anno450k)), c("chr", "pos")])) |>
  mutate(start = pos, end = pos)
cts_atlas_beta.df <- full_beta.df[match(cts_atlas.df$CpGs, full_atlas.df$CpGs), ]

## Sanity check
sum(duplicated(cts_atlas.df)) # = 0; no duplicated probes
all(rownames(cts_atlas.df) == cts_atlas.df$CpGs) # = TRUE; all probes correctly mapped
all(rownames(cts_atlas_beta.df) == cts_atlas.df$CpGs) # = TRUE; all probes correctly mapped

cts_atlas.gr <- makeGRangesFromDataFrame(cts_atlas.df, keep.extra.columns = T)
cts_atlas.gr$n_cpgs_100bp <- countOverlaps(cts_atlas.gr |> resize(width = 100, fix = 'center'), CPG_COORDS)
cts_atlas.se <- SummarizedExperiment(rowData = cts_atlas.gr,
                                     assays = list(beta = cts_atlas_beta.df))
saveRDS(cts_atlas.se, here(write_dir, 'summarizedexperiment_cell_type_specific_atlas.rds'))

# ---- All-tissue M/U regions from MethAtlas ----

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
  md.df      <- data.frame(CpGs = probe_list, coords) %>%
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
m_atlas.se <- buildSE(idx = idx_all_m, label = 'All-tissue M')

all_tissue_atlas.se <- rbind(u_atlas.se, m_atlas.se)
rowData(all_tissue_atlas.se)$n_cpgs_100bp <- countOverlaps(granges(all_tissue_atlas.se) |> resize(width = 100, fix = 'center'), CPG_COORDS)
saveRDS(all_tissue_atlas.se, here(write_dir, 'summarizedexperiment_all_tissue_um_atlas.rds'))

table(granges(all_tissue_atlas.se)$label)
# All-tissue M All-tissue U 
#         4944        78343 
