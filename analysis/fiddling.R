library(ogbox)
library(dplyr)
library(magrittr)
library(readr)
library(tibble)

# softDown('GSE79238',file = 'data-raw/GSE79238.soft')

neuroSeqMeta = softParser('data-raw/GSE79238.soft')


neuroSeqMeta %<>% select(`!Sample_characteristics_ch1 = age (postnatal day)`,
                         `!Sample_characteristics_ch1 = ercc_mix (ercc mix used)`,
                         `!Sample_characteristics_ch1 = ercc5 (10^-5 dilution of ercc in ul)`,
                         `!Sample_characteristics_ch1 = region`,
                         `!Sample_characteristics_ch1 = Sex`,
                         `!Sample_characteristics_ch1 = strain`,
                         `!Sample_characteristics_ch1 = tissue`,
                         `!Sample_characteristics_ch1 = weight (g)`,
                         `!Sample_geo_accession`,
                         `!Sample_description`,
                         `!Sample_source_name_ch1`)

names(neuroSeqMeta) = c('Age','Ercc Mix', 'Ercc dilution','Region','Sex','Strain','Tissue','Weight','GSE','sample_name','source')

download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79238/suppl/GSE79238_htseq-counts.txt.gz',destfile = 'data-raw/GSE79238_counts.gz')


neuroSeqCounts = read_tsv('data-raw/GSE79238_counts.gz')

neuroSeqCounts %<>% column_to_rownames("symbol")

neuroSeqCounts = neuroSeqCounts[,neuroSeqMeta$sample_name]


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

neuroSeqMeta$broadSource = neuroSeqMeta$source %>% sapply(function(x){
    out = regionMap %>% sapply(function(y){
        grepl(pattern = y, x,perl = TRUE,ignore.case = TRUE)
    }) %>% which %>% names
    if(length(out) == 0 | length(out)>1){
        out = ''
    }
    return(out)
})

neuroSeqMeta$source[neuroSeqMeta$broadSource == ""]  %>% table %>% sort 

neuroSeqMeta$sample_name %>% stringr::str_split('\\.|_') %>% sapply(function(x){
    x[-((length(x)-1):length(x))] %>% paste(collapse = '_')
}) %>% table

neuroSeqMeta$source %>% sapply(assignRegion,regionMap)
