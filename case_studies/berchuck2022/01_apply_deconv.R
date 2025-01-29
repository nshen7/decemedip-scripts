source('code/SETPATHS.R')
library(decemedip)
devtools::load_all('../decemedip/') # to load new changes in functions
library(argparse)
library(rstan)

## Bash Command Line Argument Parsing 
parser <- ArgumentParser()
parser$add_argument("-i", "--index", type = "integer", 
                    help = "Index of the sample")

args <- parser$parse_args()
INDEX <- args$index


code_dir <- here('code', 'sim_studies')
write_dir <- here('data', 'interim', 'case_studies', 'berchuck2022', '01_apply_deconv')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'case_studies', 'berchuck2022', '01_apply_deconv')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

# ---- Load data ----
md_samples <- fread(here('data', 'metadata', 'Berchuck2022_LuCap_PDX_MeDIP', 'sample_metadata_processed_Berchuck2022_LuCap_PDX_MeDIP.csv'))

##### THESE ARGS HAS BEEN SET TO DEFAULT IN PACKAGE
# # ---- important arguments to decemedip function ----
# 
# ref_cts <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
#                         'summarizedexperiment_cell_type_specific_atlas_v3.rds'))
# ref_anc <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
#                         'summarizedexperiment_all_tissue_um_atlas_v3.rds'))
# 
# stan_input_params <- list(
#   'alpha'     = rep(1, SummarizedExperiment::ncol(ref_cts)),
#   's_mu'      = 3,
#   's_sigma'   = 3,
#   'n_knot_z'  = 0,
#   'degree_z'  = 3,
#   'Xi'        = cor(as.matrix(SummarizedExperiment::assays(ref_cts)[[1]])),
#   's_theta'   = 3,
#   's_tau'     = 3
# )
# 
# weight_cts <- 1
# weight_anc <- 0.5

# ---- Fit decemedip-model model ----

sample_name <- md_samples$Sample_Name[INDEX]

## !! MCMC chains get stuck for some samples, in that case, change a seed
output <- decemedip(sample_bam_file = md_samples$bam_dir[md_samples$Sample_Name == sample_name],
                    paired_end = TRUE,
                    diagnostics = TRUE)
fit <- output$posterior
# plotDiagnostics(output, 'model_fit')

smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_', sample_name,'.csv')))

smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_', sample_name,'.csv')))

smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
  as.data.frame() |>
  mutate(cell_type = factor(colnames(decemedip::hg19.ref.cts.se), levels = rev(colnames(decemedip::hg19.ref.cts.se)))) |>
  relocate(cell_type)
fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_', sample_name,'.csv')))

plotDiagnostics(decemedip_output = output, plot_type = 'model_fit')
ggsave(here(plot_dir, paste0('bayesplot_diagnostics_', sample_name,'.png')), width = 12, height = 10)

plotDiagnostics(decemedip_output = output, plot_type = 'y_fit')
ggsave(here(plot_dir, paste0('point_pred_y_vs_x_', sample_name,'.png')), width = 12, height = 10)
