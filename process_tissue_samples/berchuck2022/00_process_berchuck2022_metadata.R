source('code/SETPATHS.R')

# Publication: Berchucket al. Detecting Neuroendocrine Prostate Cancer Through 
#              Tissue-Informed Cell-Free DNA Methylation Analysis. Clin Cancer 
#              Res 1 March 2022; 28 (5): 928–938. 
#              https://doi.org/10.1158/1078-0432.CCR-21-3762

read_dir <- write_dir <- here('data', 'metadata', 'Berchuck2022_LuCap_PDX_MeDIP')
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

sample_md <- readxl::read_xlsx(here(read_dir, 'NIHMS1766178-supplement-1.xlsx'), col_types = c('text', 'text'))
sample_md$`LuCaP PDX`[which(sample_md$`LuCaP PDX` == "145.19999999999999")] = "145.2" # fix bug produced by xlsx
sample_md$`LuCaP PDX`[which(sample_md$`LuCaP PDX` == "86CR")] = "96CR" # type from the publication supp table (prompted by Sylvan)

bam_md <- fread(here(read_dir, '01LuCaP_sample_sheet.csv')) 

join_md <- bam_md |>
  mutate(Sample_Name_2 = gsub(pattern = '(LuCaP_.*)_.*', replacement = '\\1', x = Sample_Name)) |>
  mutate(`LuCaP PDX` = gsub(pattern = 'LuCa(P|p)_(.*)$', replacement = '\\2', x = Sample_Name_2)) |>
  mutate(file_name = gsub('(.*)_1.fq.gz', '\\1', fq1)) |> 
  mutate(bam_dir = paste0('/scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/data/raw/Berchuck2022_LuCap_PDX_MeDIP/bam/', file_name, '_sorted.bam')) |>
  left_join(sample_md, by = join_by(`LuCaP PDX`)) |>
  select(-Sample_Name_2)

# These two samples were not included in the publication
join_md$`Tumor Type`[join_md$`LuCaP PDX` == '208.2'] = 'NEPC' 
join_md$`Tumor Type`[join_md$`LuCaP PDX` == '173.2'] = 'double-negative'

for (i in seq_len(nrow(join_md))) {
  bam_file <- join_md$bam_dir[i]
  if (file.exists(bam_file))
    join_md$million_reads[i] <- countBam(BamFile(bam_file))$records/1e6
  print(i)
}

fwrite(join_md, here(write_dir, 'sample_metadata_processed_Berchuck2022_LuCap_PDX_MeDIP.csv'))

## Write sample name list for fastq processing

write.table(join_md |> select(file_name), 
            here('data', 'raw', 'Berchuck2022_LuCap_PDX_MeDIP', 'sample_name_list.txt'),
            row.names = F, col.names = F, quote = F)
write.table(join_md |> filter(`LuCaP PDX` %in% c('96CR', '208.2', '173.2')) |> select(file_name), 
            here('data', 'raw', 'Berchuck2022_LuCap_PDX_MeDIP', 'sample_name_list_2.txt'),
            row.names = F, col.names = F, quote = F)
