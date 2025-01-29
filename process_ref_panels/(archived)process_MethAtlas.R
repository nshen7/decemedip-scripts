source('code/SETPATHS.R')
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

read_dir <- write_dir <- here('data', 'metadata', 'references', 'MethAtlas')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

## 450K array probe annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ---- Cell-type specific regions from MethAtlas ----

cts_atlas.df <- readxl::read_xlsx(here('data', 'metadata', 'references', 'MethAtlas', 'supp_data_1.xlsx'), 
                                  sheet = 'Table 1') %>%
  dplyr::rename(label = name, start = from, end = to) %>%
  group_by(label) %>% 
  mutate(name = paste(label, 1:n())) %>%
  ungroup()

md_cols <- c('acc', 'chr', 'pos', 'start', 'end', 'group', 'label', 'name')
cts_atlas_beta.df <- cts_atlas.df %>% select(-md_cols)
cts_atlas_md.df <- cts_atlas.df %>% 
  select(md_cols) %>% 
  mutate(direction = ifelse(group > 0, '+', '-'))

## Obtain target and mean background beta
getTargetAndBackgroundBeta <- function(i, md = cts_atlas_md.df, betas = cts_atlas_beta.df) {
  
  ct_label <- md$label[i]
  cn <- colnames(betas)
  target <- betas[i, which(cn == ct_label)] %>% unlist
  background <- betas[i, -which(cn == ct_label)] %>% unlist
  
  return(data.frame(beta_target = target,
                    beta_backgound_mean = mean(background)))
}

beta_bg_and_tg.df <- map_dfr(1:nrow(cts_atlas_md.df), getTargetAndBackgroundBeta)
cts_atlas_md.df <- cbind(cts_atlas_md.df, beta_bg_and_tg.df)
rownames(cts_atlas_md.df) <- NULL

cts_atlas.se <- SummarizedExperiment(rowData = makeGRangesFromDataFrame(cts_atlas_md.df, keep.extra.columns = T),
                                     assays = list(beta = cts_atlas_beta.df))
saveRDS(cts_atlas.se, here(write_dir, 'summarizedexperiment_cell_type_specific_atlas.rds'))

# ---- All-tissue M/U regions from MethAtlas ----

full_atlas.df <- fread(here(read_dir, 'full_atlas.csv.gz'), header = F)

# Find blocks where WGBS samples are all methylated (beta > 1-margin) or all unmethylated (beta < margin)
margin <- 0.05
idx_all_m <- which(rowMeans(full_atlas_beta.df > 1-margin) >= 1)
idx_all_u <- which(rowMeans(full_atlas_beta.df < margin) >= 1)

buildSE <- function(idx, label, bin_size = 100) {
  
  subset.df <- full_atlas.df[idx, ]
  
  probe_list <- unlist(subset.df[, 1])
  coords     <- anno450k[rownames(anno450k) %in% probe_list, c("chr", "pos")]
  md.df      <- data.frame(acc = probe_list, coords) %>%
    mutate(label = label, margin = margin, start = pos - bin_size/2, end = pos + bin_size/2) 
  rownames(md.df) <- NULL
    
  beta.df <- subset.df[, -1]
  colnames(beta.df) <- colnames(cts_atlas_beta.df)
  
  se <- SummarizedExperiment(rowData = makeGRangesFromDataFrame(md.df, keep.extra.columns = T),
                             assays = list(beta = beta.df))
  return(se)
}

u_atlas.se <- buildSE(idx = idx_all_u, label = 'All-tissue U')
m_atlas.se <- buildSE(idx = idx_all_m, label = 'All-tissue M')

all_tissue_atlas.se <- rbind(u_atlas.se, m_atlas.se)
saveRDS(all_tissue_atlas.se, here(write_dir, 'summarizedexperiment_all_tissue_um_atlas.rds'))

