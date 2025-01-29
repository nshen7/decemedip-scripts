source('code/SETPATHS.R')
library(SummarizedExperiment)

read_dir <- write_dir <- here('data', 'metadata', 'references', 'MethAtlas')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

# ---- Load cell-type specific regions from MethAtlas ----

cts_all.df <- readxl::read_xlsx(here('data', 'metadata', 'references', 'MethAtlas', 'supp_data_1.xlsx'), sheet = 'Table 1') %>%
  dplyr::rename(label = name, start = from, end = to) %>%
  group_by(label) %>% 
  mutate(name = paste(label, 1:n())) %>%
  ungroup()
  
md_cols <- c('acc', 'chr', 'pos', 'start', 'end', 'group', 'label', 'name')
cts_all_md.df <- cts_all.df %>% select(md_cols) %>%
  mutate(direction = ifelse(group > 0, '+', '-')) 

cts_all_beta.df <- cts_all.df %>% select(-md_cols)

cts_all_beta_with_name.df <- cbind(name = cts_all_md.df$name, cts_all_beta.df)
fwrite(cts_all_beta_with_name.df, here(write_dir, 'CTS_regions_beta.csv'))

cts.se <- SummarizedExperiment(
  rowData = makeGRangesFromDataFrame(cts_all_md.df, keep.extra.columns = T),
  assays = list(beta = cts_all_beta.df)
)


# ---- Compute difference between target cell type and background cell types ----

getTargetAndBackgroundBeta <- function(i, region_md = cts_all_md.df, region_betas = cts_all_beta.df) {
  
  ct_label <- region_md$label[i]
  cn <- colnames(region_betas)
  target <- region_betas[i, which(cn == ct_label)] %>% unlist
  background <- region_betas[i, -which(cn == ct_label)] %>% unlist
  
  return(data.frame(beta_target = target,
                    beta_backgound_mean = mean(background)))
}

beta_bg_and_tg.df <- do.call(rbind, lapply(1:nrow(cts_all_md.df), getTargetAndBackgroundBeta))

cts_all_md_processed.df <- cbind(cts_all_md.df, beta_bg_and_tg.df) %>%
  mutate(beta_diff = beta_target - beta_backgound_mean)
fwrite(cts_all_md_processed.df, here(write_dir, 'processed_CTS_regions_metadata.csv'))


cts_all_md_processed.df %>%
  ggplot(aes(label, beta_diff)) + 
  geom_boxplot() +
  facet_wrap(~ direction, ncol = 1, scales = 'free_y') + 
  theme_classic()

