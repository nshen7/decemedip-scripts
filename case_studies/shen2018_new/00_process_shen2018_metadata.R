source('code/SETPATHS.R')

write_dir <- here('data', 'metadata', 'shen2018')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

md_old <- fread('/scratch/st-kdkortha-1/kdkortha/shen2018_raw/cfMeDIP_Transfer_Korthauer_md5sum.csv', header = F)
md_new <- md_old |> 
  dplyr::select(V1) |>
  tidyr::separate(V1, into = c('group', 'sample'), sep = '/') |>
  dplyr::filter(!grepl('\\.bai$', sample)) |>
  dplyr::mutate(bamdir = paste0('/scratch/st-kdkortha-1/kdkortha/shen2018_raw/cfMeDIP_Transfer_Korthauer/', group, '/', sample)) |>
  dplyr::mutate(sample = gsub('(.*).bam', '\\1', x = sample),
                group = gsub('cfMeDIP-(.*)', '\\1', x = group)) |>
  dplyr::mutate(million_reads = map_int(bamdir, ~countBam(BamFile(.x))$records)/1e6)

fwrite(md_new, here(write_dir, 'sample_metadata_processed_shen2018.csv'))
