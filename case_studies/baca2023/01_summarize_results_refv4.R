source('code/SETPATHS.R')
library(rstan)
library(bayesplot)
library(argparse)
library(Rsamtools)

read_dir <-  here('data', 'interim', 'case_studies', 'baca2023', 'from_ze')
write_dir <- here('data', 'interim', 'case_studies', 'baca2023', '01_summarize_results_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

load(here(read_dir, 'Decemedip_NM_Results.RDATA')) ## loads object `NM_pheno`
cell_types <- colnames(decemedip::hg19.ref.cts.se)
md_samples <- NM_pheno |>
  select(-any_of(cell_types)) |>
  mutate(sample = paste0(study_name, '_unique.sorted.dedup'),
         group = cancer_type,
         million_reads = fragments / 1e6) |> 
  ## Remove patients with library size < 5 million reads
  filter(million_reads >= 5) |>
  select(-cancer_type) |>
  # Only records major cancer types:
  mutate(group = ifelse(group == 'NEPC', 'Prostate', group),
         group = ifelse(group %in% c('NSCLC', 'SCLC'), 'Lung', group),
         group = ifelse(group == 'Small cell bladder', 'Bladder', group)
  ) |>
  # Remove cancer types that does not have obvious corresponding cell type present in reference panel
  filter(! group %in% c('Glioma', 'Melanoma', 'Merkel cell', 'Ovarian', 'Thymic')) |>
  # Change cancer type names 
  mutate(group = ifelse(group == 'Healthy', group, paste(group, 'cancer'))) |>
  # Remove cancer types that contains less than 5 patients
  group_by(group) |>
  mutate(n = n()) |>
  filter(n >= 5) |>
  ungroup()
fwrite(md_samples, here(write_dir, 'processed_metadata.csv'))

######################################
# ---- estimated proportions ----
######################################

all_samples_pi.df <- data.frame()
all_samples_pi_nonblood.df <- data.frame()

for (i in seq_len(nrow(md_samples))) {

  sample <- md_samples$sample[i]
  group <- md_samples$group[i]

  sample_pi.df <- fread(here(read_dir, 'NM_output_01022025', paste0('fitted_pi_', sample,'.csv.gz'))) |>
    select(cell_type, mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group) |>
    relocate(sample, group)

  sample_pi_nonblood.df <- sample_pi.df |>
    dplyr::slice(8:n()) |> # remove blood cells
    group_by(sample, group) |>
    summarise(total = sum(mean)) |>
    ungroup()

  all_samples_pi.df <- all_samples_pi.df |> rbind(sample_pi.df)
  all_samples_pi_nonblood.df <- all_samples_pi_nonblood.df |> rbind(sample_pi_nonblood.df)

  cat(sample, ' ')
}
fwrite(all_samples_pi.df, here(write_dir, paste0('estimated_pi.csv.gz')))
fwrite(all_samples_pi_nonblood.df, here(write_dir, paste0('estimated_pi_nonblood.csv.gz')))

######################################
# ---- estimated w_mu ----
######################################

### RUN ONCE
all_samples_w_mu.df <- data.frame(sample = NULL, group = NULL, million_reads = NULL, mean = NULL,
                                  `2.5%` = NULL, `25%` = NULL, `50%` = NULL, `75%` = NULL, `97.5%` = NULL)

for (i in seq_len(nrow(md_samples))) {

  sample <- md_samples$sample[i]
  group <- md_samples$group[i]
  million_reads <- md_samples$million_reads[i]

  sample_w_mu.df <- fread(here(read_dir, 'NM_output_01022025', paste0('fitted_w_mu_', sample,'.csv.gz'))) |>
    select(mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group, million_reads = million_reads, parameter = paste0('w_mu_', c(1:(n()-1), 0))) |>
    relocate(sample, group, million_reads, parameter)

  all_samples_w_mu.df <- all_samples_w_mu.df |> rbind(sample_w_mu.df)

  cat(sample, ' ')
}
fwrite(all_samples_w_mu.df, here(write_dir, paste0('estimated_w_mu.csv.gz')))



######################################
# ---- estimated w_sigma ----
######################################

### RUN ONCE
all_samples_w_sigma.df <- data.frame(sample = NULL, group = NULL, million_reads = NULL, mean = NULL,
                                  `2.5%` = NULL, `25%` = NULL, `50%` = NULL, `75%` = NULL, `97.5%` = NULL)

for (i in seq_len(nrow(md_samples))) {

  sample <- md_samples$sample[i]
  group <- md_samples$group[i]
  million_reads <- md_samples$million_reads[i]
  
  sample_w_sigma.df <- fread(here(read_dir, 'NM_output_01022025', paste0('fitted_w_sigma_', sample,'.csv.gz'))) |>
    select(mean, `2.5%`, `25%`, `50%`, `75%`, `97.5%`) |>
    mutate(sample = sample, group = group, million_reads = million_reads, parameter = paste0('w_sigma_', c(1:(n()-1), 0))) |>
    relocate(sample, group, million_reads, parameter)

  all_samples_w_sigma.df <- all_samples_w_sigma.df |> rbind(sample_w_sigma.df)

  cat(sample, ' ')
}
fwrite(all_samples_w_sigma.df, here(write_dir, paste0('estimated_w_sigma.csv.gz')))


