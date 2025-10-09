library(RColorBrewer)
library(circlize)

source('/home/ubuntu/projects/edsouza-summer2023/cxcr4-pdac/helpers/add_metadata.R')

colors_only_malig <- c('Group_1'='#E7298A', 'Unselected'='#666666')

# RColorBrewer::brewer.pal(8, 'Dark2')
colors_primary <- c('Myeloid'='#1B9E77',
                    'T cell'='#D95F02',
                    'B cell'='#7570B3',
                    'Malignant'='#E7298A',
                    'Fibroblast'='#66A61E',
                    'Hepatocyte'='#E6AB02',
                    'Endothelial'='#A6761D',
                    'Epithelial'='#666666')

# From colorblind-friendly RColorBrewer palettes
colors_treatment_stage <- c('ON_PROG'='#1B9E77',
                            'ON_TRT'='#7570B3',
                            'PRE'='#E6AB02')
colors_effective_treatment_stage <- c('Pre-treatment' = '#E6AB02', 
                                      'Treatment\nw/o progression' = '#bdb9f1', 
                                      'Treatment\nw/ progression' = '#117356')

colors_treatment_response <- c('PR'='#4daf4a', 
                              'PD'='#e41a1c', 
                              'SD'='#377eb8')

no_ticks <- theme(axis.text.x=element_blank(), 
                  axis.ticks.x=element_blank(), 
                  axis.text.y=element_blank(), 
                  axis.ticks.y=element_blank()) 



b_pal <- scales::hue_pal(h = c(180, 240))(4)
t_pal <- scales::hue_pal(h = c(300, 360))(7)
myel_pal <- scales::hue_pal(h = c(60, 120))(9)
colors_immune_comp <- c('Cycling B cell'    =b_pal[1],
                        'Memory B cell'     =b_pal[2],
                        'Naive B'           =b_pal[3],
                        'Plasma cell'       =b_pal[4],
                        'cDC'              =myel_pal[1],
                        'Cycling myeloid'  =myel_pal[2],
                        'HSC'              =myel_pal[3],
                        'Langerhans'       =myel_pal[4],
                        'M1'               =myel_pal[5],
                        'M2'               =myel_pal[6],
                        'Mast cell'        =myel_pal[7],
                        'Monocyte'         =myel_pal[8],
                        'pDC'              =myel_pal[9],
                        'CD4 Effector'      =t_pal[1],
                        'CD4 Memory'        =t_pal[2],
                        'CD8 Tem'           =t_pal[3],
                        'CD8 Tex'           =t_pal[4],
                        'MAIT cells'        =t_pal[5],
                        'NK'                =t_pal[6],
                        'Tregs'             =t_pal[7])

# RColorBrewer::brewer.pal(8, 'Dark2')
colors_fib_comp <- c('Myofibr.'='#D95F02',
                    'Hepatocyte trans. fibr.'='#E6AB02',
                    'Mesothelial fibr.'='#A6761D',
                    'LN fibr. retic.'='#1B9E77',
                    'Pericyte fibr.'='#E7298A',
                    'CXCL12+ fibr.'='#7570B3',
                    'CXCL14+ fibr.'='#66A61E',
                    'Alveolar fibr.'='#666666')

biopsy_site_cols <- c('pancreas'='#ff7f00', 'lung'='#377eb8', 'lymph node'='#4daf4a', 'liver'='#e41a1c')

##########################################################
# NMF Colors
# Helper funcs for factor metadata

# From ColorBrewer
pt_id_cols <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', 
                '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', 
                '#cab2d6', '#6a3d9a', '#ffff99')
pt_id_cols <- pt_id_cols[1:length(unique(metadata$pt.id))]
names(pt_id_cols) <- unique(metadata$pt.id)

# From https://mokole.com/palette.html
orig_ident_cols <- c('#2f4f4f', '#8b4513', '#191970', '#006400', '#808000', 
                      '#4682b4', '#9acd32', '#20b2aa', '#8b008b', '#ff0000', 
                      '#ffa500', '#ffff00', '#00ff00', '#00fa9a', '#8a2be2', 
                      '#4169e1', '#dc143c', '#0000ff', '#d8bfd8', '#ff7f50', 
                      '#ff00ff', '#db7093', '#ff1493', '#f5deb3', '#ee82ee')
orig_ident_cols <- orig_ident_cols[1:length(unique(metadata$orig.ident))]
names(orig_ident_cols) <- unique(metadata$orig.ident)



get_col_meta <- function(factorname, desired_var, metadata) {
  samplename <- strsplit(factorname, '_')[[1]][1]
  metadata[metadata$orig.ident==samplename, desired_var]
}

make_col_colors <- function(metadata) {

  age_pal <- RColorBrewer::brewer.pal(n = 8, name = "PuRd")
  col_annot_colors <- list(
      sex=c('M'='#377eb8', 'F'='#e41a1c'),
      pt.id=pt_id_cols,
      orig.ident=orig_ident_cols,
      biopsy.site=biopsy_site_cols,
      age=colorRamp2(c(min(metadata$age), max(metadata$age)), c(age_pal[2], age_pal[8])),
      pct.change=colorRamp2(c(-100, 0, 100), c('#4dac26', '#f7f7f7', '#d01c8b')),
      response=colors_treatment_response,  
      treatment.stage=colors_treatment_stage
    )

    return(col_annot_colors)
}

make_column_ha <- function(mat, metadata) {
  factors <- colnames(mat)

  HeatmapAnnotation(
    pt.id = sapply(factors, function(f) get_col_meta(f, 'pt.id', metadata)),
    orig.ident = sapply(factors, function(f) get_col_meta(f, 'orig.ident', metadata)),
    sex = sapply(factors, function(f) get_col_meta(f, 'sex', metadata)),
    age = sapply(factors, function(f) get_col_meta(f, 'age', metadata)),
    empty1 = anno_empty(border = FALSE, height=unit(1, "mm", metadata)),
    biopsy.site = sapply(factors, function(f) get_col_meta(f, 'biopsy.site', metadata)),
    treatment.stage = sapply(factors, function(f) get_col_meta(f, 'treatment.stage', metadata)),
    pct.change = sapply(factors, function(f) get_col_meta(f, 'pct.change', metadata)),
    response = sapply(factors, function(f) get_col_meta(f, 'response', metadata)),
    empty2 = anno_empty(border = FALSE, height=unit(1, "mm", metadata)),
    col=make_col_colors(metadata)
  )
}

dotplot_theme <- theme(legend.position='right',
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.title = element_blank())
dotplot_colors <- colorRampPalette((RColorBrewer::brewer.pal(9, "Greens")))(100)


diversity_plot_theme <- theme(legend.position='none', 
                    axis.text=element_text(size=14, color='black'),
                    axis.title=element_text(size=16, color='black'))
lighten <- function(x) adjustcolor(x, offset = c(0.15, 0.15, 0.15, 0))



swarm_plot_colors <- scale_color_gradient2(low=scales::muted('blue'), mid='white', high=scales::muted('red')) 
swarm_plot_theme <-  theme(legend.position='none', 
                        axis.text=element_text(color='black'),
                        plot.title=element_text(hjust=1)) 