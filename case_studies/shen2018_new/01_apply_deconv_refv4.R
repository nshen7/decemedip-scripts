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
write_dir <- here('data', 'interim', 'case_studies', 'shen2018_new', '01_apply_deconv_refv4')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
plot_dir <- here('plots', 'case_studies', 'shen2018_new', '01_apply_deconv_refv4')
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

md_samples <- fread(here('data', 'metadata', 'shen2018', 'sample_metadata_processed_shen2018.csv'))

# ---- important arguments to decemedip function ----

ref_cts <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                        'summarizedexperiment_cell_type_specific_atlas_v4.rds'))
ref_anc <- readRDS(here('data', 'metadata', 'references', 'MethAtlas', 
                        'summarizedexperiment_all_tissue_um_atlas_v4.rds'))

# ---- Fit decemedip-model model ----

sample_name <- md_samples$sample[INDEX]

output <- decemedip(sample_bam_file = md_samples$bamdir[md_samples$sample == sample_name],
                    paired_end = FALSE,
                    diagnostics = TRUE, 
                    seed = 2025,
                    ref_cts = ref_cts,
                    ref_anc = ref_anc)
fit <- output$posterior
# plotDiagnostics(output, 'model_fit')

smr_w_mu.df <- monitor(extract(fit, pars=c("w_mu"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
fwrite(smr_w_mu.df, here(write_dir, paste0('fitted_w_mu_', sample_name,'.csv')))

smr_w_sigma.df <- monitor(extract(fit, pars=c("w_sigma"), permuted = FALSE), digits_summary = 5) |> as.data.frame()
fwrite(smr_w_sigma.df, here(write_dir, paste0('fitted_w_sigma_', sample_name,'.csv')))

smr_pi.df <- monitor(extract(fit, pars=c("pi"), permuted = FALSE), digits_summary = 5) |>
  as.data.frame() |>
  mutate(cell_type = factor(colnames(ref_cts), levels = rev(colnames(ref_cts)))) |>
  relocate(cell_type)
fwrite(smr_pi.df, here(write_dir, paste0('fitted_pi_', sample_name,'.csv')))


plotDiagnostics(decemedip_output = output, plot_type = 'model_fit')
ggsave(here(plot_dir, paste0('bayesplot_diagnostics_', sample_name,'.png')), width = 12, height = 10)

plotDiagnostics(decemedip_output = output, plot_type = 'y_fit')
ggsave(here(plot_dir, paste0('point_pred_y_vs_x_', sample_name,'.png')), width = 12, height = 10)
