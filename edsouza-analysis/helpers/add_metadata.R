
# Celltype annotations
annotations <- readRDS('/home/ubuntu/data/cxcr4-pdac/combine_samples/annotations.RDS')
annotations$broad_annot <- annotations$Primary
epi_cells <- rownames(annotations[annotations$Primary=='Epithelial',])
malig_cells <- rownames(annotations[annotations$Primary=='Malignant',])
annotations[epi_cells, 'broad_annot'] <- paste0('Epithelial::', gsub('\\:.*$', '', annotations[epi_cells, 'Secondary']))
annotations$granular_annot <- paste0(annotations$Primary, '::', annotations$Secondary)
annotations[malig_cells, 'granular_annot'] <- paste0('_', annotations[malig_cells, 'granular_annot'])  # Literally just make the color pop more


# Patient-level metadata
metadata <- data.frame(
        orig.ident = paste0('GM', seq(16, 40)),
        pt.id = c('001-01', '001-03', '001-04', '001-05', '001-06', '001-08', '001-12', '001-18', '001-19', '002-01', '002-02',
                  '001-01', '001-03', '001-04', '001-05', '001-06', '001-08', '001-12', '001-18', '001-19', '002-01', '002-02',
                  '001-03', '001-04', '001-06'),
        pct.change = c(-31, -14, -50, 38, -66, -100, -53, -75, -30, -18, -21,
                       -31, -14, -50, 38, -66, -100, -53, -75, -30, -18, -21,
                       -14, -50, -66),
        treatment.stage = c(rep('PRE', 11), rep('ON_TRT', 11), rep('ON_PROG', 3)),
        age = c(62, 58, 58, 58, 46, 68, 58, 69, 63, 68, 47,
                62, 58, 58, 58, 46, 68, 58, 69, 63, 68, 47,
                58, 58, 46),
        race = c('White', 'White', 'White', 'Black/AA', 'Black/AA', 'White', 'White', 'Black/AA', 'White', 'White', 'White',
                 'White', 'White', 'White', 'Black/AA', 'Black/AA', 'White', 'White', 'Black/AA', 'White', 'White', 'White',
                 'White', 'White', 'Black/AA'),
        sex = c('M', 'F', 'M', 'M', 'M', 'F', rep('M', 5),
                'M', 'F', 'M', 'M', 'M', 'F', rep('M', 5),
                'F', 'M', 'M'),
        ethnicity = c('Not Hispanic', 'Not Hispanic', 'Dominican', 'Not Hispanic', 'Dominican', rep('Not Hispanic', 6),
                        'Not Hispanic', 'Not Hispanic', 'Dominican', 'Not Hispanic', 'Dominican', rep('Not Hispanic', 6),
                        'Not Hispanic', 'Dominican', 'Dominican'),
        response = c('PR', 'SD', 'PR', 'PD', rep('PR', 5), 'SD', 'SD',
                        'PR', 'SD', 'PR', 'PD', rep('PR', 5), 'SD', 'SD',
                        'SD', 'PR', 'PR')
        )

metadata <- metadata %>% mutate(biopsy.site = case_when(
    orig.ident %in% c("GM22", "GM33") ~ "pancreas",
    orig.ident == "GM36" ~ "lung",
    orig.ident == "GM40" ~ "lymph node",
    TRUE ~ "liver"  # default
  ))


metadata$effective.treatment.stage <- ifelse(metadata$treatment.stage=='PRE', 'Pre-treatment',
                                        ifelse(metadata$treatment.stage=='ON_TRT', 'Treatment\nw/o progression',
                                        ifelse(metadata$treatment.stage=='ON_PROG', 'Treatment\nw/ progression', NA)))
metadata[(metadata$response=='PD') & (metadata$treatment.stage=='ON_TRT'), 'effective.treatment.stage'] <- 'Treatment\nw/ progression'
metadata$effective.treatment.stage <- factor(metadata$effective.treatment.stage, levels=c('Pre-treatment', 'Treatment\nw/o progression', 'Treatment\nw/ progression'))    


# Add celltype annotations and patient metadata
add_metadata <- function(obj){
        if (!('pt.id' %in% colnames(obj@meta.data))) {                
                md_backup <- obj@meta.data
                md_new <- md_backup %>% rownames_to_column('TEMP_BARCODE') %>%
                        merge(metadata, by.x='orig.ident', by.y='orig.ident') %>%
                        column_to_rownames('TEMP_BARCODE') 
                obj@meta.data <- md_new
                }
        return(obj %>% AddMetaData(annotations))
}


make_granular_labels <- function() {
  annotations <- readRDS(objects[['annotations']])
  annotations$broad_annot <- annotations$Primary
  epi_cells <- rownames(annotations[annotations$Primary=='Epithelial',])
  malig_cells <- rownames(annotations[annotations$Primary=='Malignant',])
  annotations[epi_cells, 'broad_annot'] <- paste0('Epithelial::', gsub('\\:.*$', '', annotations[epi_cells, 'Secondary']))
  annotations$granular_annot <- paste0(annotations$Primary, '::', annotations$Secondary)
  annotations[malig_cells, 'granular_annot'] <- paste0('_', annotations[malig_cells, 'granular_annot'])  # Literally just make the color pop more
  return(annotations)
}