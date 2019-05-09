library(ogbox)
library(dplyr)
library(magrittr)
library(readr)
library(tibble)
library(readxl)

softDown('GSE79238',file = 'data-raw/GSE79238.soft',overwrite = TRUE)


neuroSeqMeta = softParser('data-raw/GSE79238.soft')


neuroSeqMeta %<>% select(`!Sample_characteristics_ch1 = age (postnatal day)`,
                         `!Sample_characteristics_ch1 = ercc_mix (ercc mix used)`,
                         `!Sample_characteristics_ch1 = ercc5 (10^-5 dilution of ercc in ul)`,
                         `!Sample_characteristics_ch1 = region`,
                         `!Sample_characteristics_ch1 = Sex`,
                         `!Sample_characteristics_ch1 = strain`,
                         `!Sample_characteristics_ch1 = tissue`,
                         `!Sample_characteristics_ch1 = weight (g)`,
                         `!Sample_characteristics_ch1 = num_cells (# of sorted cells in the sample)`,
                         `!Sample_geo_accession`,
                         `!Sample_title`,
                         `!Sample_description`,
                         `!Sample_source_name_ch1`)

names(neuroSeqMeta) = c('Age','Ercc Mix', 'Ercc dilution','Region','Sex','Strain','Tissue','Weight','Cell Count','GSM','name','sample_label','source')


regionMap = list(Cerebellum = '^Cerebellum',
                 Hippocampus = '^Hippocampus',
                 Hypothalamus = '^Hypothalamus',
                 Cortex = '^Isocortex',
                 Medulla = '^Medulla',
                 Midbrain = '^Midbrain',
                 Olfactory = '^Olfactory',
                 `Spinal Cord` = '^Spinal Cord',
                 Striatum = '^Striatum',
                 Thalamus = '^Thalamus',
                 Pons = '^Pons')

neuroSeqHierarchy = list(All =
                             list(Cerebrum = list(Cortex = '',
                                                  Hippocampus = '',
                                                  Olfactory = '',
                                                  Striatum = ''),
                                  `Brain stem` = list(
                                      Interbrain = list(Thalamus = '',
                                                        Hypothalamus = ''),
                                      Midbrain = '',
                                      Hindbrain = list(Medulla = '',
                                                       Pons = '')
                                  ),
                                  Cerebellum = '',
                                  `Spinal Cord` = ''
                             ))

neuroSeqMeta$ogregion = neuroSeqMeta$source %>% sapply(function(x){
    out = regionMap %>% sapply(function(y){
        grepl(pattern = y, x,perl = TRUE,ignore.case = TRUE)
    }) %>% which %>% names
    if(length(out) == 0 | length(out)>1){
        out = ''
    }
    return(out)
})



neuroSeqMeta$cellType = neuroSeqMeta$name %>% gsub(' rep[0-9]','',.)


# unasigned cells
neuroSeqMeta$source[neuroSeqMeta$ogregion == ""]  %>% table %>% sort 

neuroSeqSupplement = read_xlsx('data-raw/elife-38619-supp2-v2.xlsx')

neuroSeqMeta$sample_id = neuroSeqMeta$sample_label %>% stringr::str_extract('[0-9]*?$') %>% as.integer()
# missing_id = neuroSeqMeta$sample_label[!neuroSeqMeta$sample_id %in% neuroSeqSupplement$sample_id]
neuroSeqMeta %<>% filter(sample_id %in% neuroSeqSupplement$sample_id)
neuroSeqMeta %<>% arrange(sample_id)

neuroSeqSupplement %<>% filter(sample_id %in% neuroSeqMeta$sample_id) %>% arrange(sample_id)

neuroSeqSupplement$sample_label %<>% stringr::str_replace_all('\\.up_','_')
neuroSeqSupplement$sample_label %<>% stringr::str_replace_all('\\.low_','_')

neuroSeqSupplement$sample_label %<>% stringr::str_replace_all('\\.\\.','/')

neuroSeqMeta$sample_label[!neuroSeqMeta$sample_label %in% neuroSeqSupplement$sample_label]
neuroSeqSupplement$sample_label[!neuroSeqSupplement$sample_label %in% neuroSeqMeta$sample_label]


neuroSeqMeta$transmitter = neuroSeqSupplement$transmitter

# read counts in ------
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79238/suppl/GSE79238_htseq-counts.txt.gz',destfile = 'data-raw/GSE79238_counts.gz')


neuroSeqCounts = read_tsv('data-raw/GSE79238_counts.gz')

neuroSeqCounts %<>% column_to_rownames("symbol")

neuroSeqMeta %<>% arrange(transmitter, cellType)

neuroSeqCounts = neuroSeqCounts[,neuroSeqMeta$sample_label]






neuroSeqMeta %>% 
    filter(ogregion =='Cortex') %$%
    sample_label %>%
    {neuroSeqCounts['Cox6a2',.]} %>% unlist %>% {.>2000} %>% which


neuroSeqMeta %>% 
    filter(ogregion =='Cortex') %$%
    sample_label %>%
    {neuroSeqCounts['Sst',.]} %>% unlist %>% plot




neuroSeqMeta %>% 
    filter(region =='Cortex') %$% cellType %>% table



