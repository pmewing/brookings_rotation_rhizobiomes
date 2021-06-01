libs = c('here',
         'ggplot2',
         'ggpmisc',
         'data.table', 
         'magrittr',
         'grid',
         'gridExtra',
         'lavaan')

for (i in libs) {
  if (!(require(i, character.only=TRUE))) {
    install.packages(i, Ncpu=4)
    library(i, character.only=TRUE)
  }
}


load_pcwOrd = function() {
  here('Scripts', 
       'pcwOrd_0.2.1.Rdata') %>% 
    load(envir=.GlobalEnv)
}
load_pcwOrd()


load(here('Data', 
          'Microbes as Mechanism Results.Rdata'))  # loads ord_ls, among other objects


ICs = lapply(sem_ls, function(x) {
  x[['df']] = NULL
  IC = sapply(x, AIC)
  min_IC = IC[which.min(IC)]
  del_IC = min_IC - IC['no_mediation']
  chi_sq = names(min_IC) %>% 
    x[[.]] %>% 
    anova %>% 
    .[2, c('Chisq', 'Pr(>Chisq)')]
  out = c(min_IC, 
          del_IC = del_IC, 
          chi_sq) %>% 
    as.data.frame
  names(out)[4] = 'Chisq_p_val'
  return(out)
}) #%>% 
  # do.call(cbind, .) %>% 
  # set_colnames(names(sem_ls))

plot_these = c('corn_seedling_FOLLOW_SOY',
               'corn_flowering_ROTATION_LENGTH',
               'soy_seedling_FOLLOW_CORN',
               'soy_flowering_ROTATION_LENGTH')

# change order and colors. matches factor levels using names. 
palette = c(CS = '#A8AA39', # MEDIUM YELLOW,
            COWwS = '#6F1C1A', # DARK RED
            CPWwS = '#E56764', # LIGHT RED
            CSSwP = '#4C639A', # LIGHT BLUE
            CSSwSf = '#17264A') # DARK BLUE

# for the title of each sub-plot. Names are specified in names(ord_ls)
name_parser = c(seedling = 'c) Seedling',
                     flowering = 'd) Flowering',
                     ITS = 'Fungal Community',
                     `16S` = 'Bacterial Community',
                     FOLLOW_SOY = 'a Soybean Legacy',
                     FOLLOW_CORN = 'a Maize Legacy',
                     ROTATION_LENGTH = 'a Four Species Rotation')

# for the title of each figure (group of plots). Names are crop specified in names(ord_ls).
crop_name = c(corn = 'Corn',
              soy = 'Soybean')

# for saving. Names are specified in names(ord_ls)
figure_filename = c(corn = 'FIG 4 - Corn Biomass from Microbes.jpeg',
                    soy = 'FIG 7 - Soybean Biomass from Microbes.jpeg')

save = FALSE
save = TRUE


plot_axis_effect = function(barcode, response, partial, constraint, contrast, data, partial_response=TRUE, title, palette=NULL,
                            autoflip_x=TRUE) {
  
  label.x=0.9
  
  yax_name = expression(paste('Biomass [g plant' ^-1, ']', sep=''))
  xax_name = name_parser[barcode] %>% 
    paste('Response to', name_parser[contrast])
  xax = paste0('AXIS_', barcode)
  if (!(xax %in% names(data))) xax %<>% paste(contrast, sep='_')
  
  pltdf = data
  if (partial_response) {
    pltdf[, response] = paste(response, 
                              paste(partial, collapse='+'), 
                              sep='~') %>% 
      formula %>% 
      lm(data=pltdf) %>% 
      residuals
    yax_name = expression(paste('Marginal Biomass [g plant' ^-1, ']', sep=''))
  }

  flip_x = contrast=='ROTATION_LENGTH'  # flip x axis for rotation length (or, make sure 2-year is negative)  
  if (flip_x & autoflip_x) {
    mean_CS = tapply(pltdf[, xax], pltdf$ROTATION, mean) %>% 
      .['CS']
    if (mean_CS > 0) {
      pltdf[, xax] %<>% multiply_by(-1)
      label.x=0.1
      
    }
  }
  
  
  
  plt_out = 
    ggplot(pltdf, aes_string(x=xax,
                             y=response))
  
  if (is.null(palette)) {
    plt_out = plt_out + scale_color_brewer(palette='Dark2')
  } else {
    plt_out = plt_out + scale_color_manual(values=palette)
  }
  
  plt_out = plt_out +
    geom_smooth(aes_string(x=xax,
                           y=response),
                color='darkgray',
                fill='lightgray',
                linetype='dashed',
                formula='y~x',
                method='lm',
                inherit.aes=FALSE) +
    geom_point(aes_string(shape=partial[1], # if partialing out multiple factors, keep the first
                          color=constraint),
               size=2,
               alpha=0.9) +
    stat_poly_eq(aes(label=stat(eq.label)),
                 formula=y~x,
                 npcx = label.x,
                 size=3,
                 parse=TRUE) +
    stat_poly_eq(aes(label=stat(rr.label)),
                 formula=y~x,
                 npcx = label.x,
                 npcy = 0.9,
                 size=3,
                 parse=TRUE) +
    xlab(xax_name) +
    ylab(yax_name) +
    ggtitle(title) +
    theme_minimal() +
    theme(panel.grid=element_line(size=0.2),
          legend.title=element_blank(),
          title=element_text(size=8),
          axis.title=element_text(size=8),
          axis.text=element_text(size=8),
          legend.text=element_text(size=8))
  return(plt_out)
}


grid_arrange_shared_legend = function(plots, title=NULL, position='bottom') {
  # plots <- list(...)
  nplots = length(plots)
  ncol = floor(sqrt(nplots))
  nrow = ceiling(sqrt(nplots))
  if (position %in% c('left', 'right')) {
    ncol = ncol + 1
    nrow = nrow
    lheight_adjust = 1
    lwidth_adjust = 0.95
  } else {
    ncol = ncol
    nrow = nrow + 1
    lheight_adjust = 0.95
    lwidth_adjust = 1
  }
  
  g <- ggplotGrob(plots[[1]] + theme(legend.position=position, 
                                     legend.box='vertical', 
                                     legend.spacing=unit(0, 'in'),
                                     legend.margin=margin()))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth = sum(legend$width)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = ncol,
    nrow = nrow,
    heights = unit.c(unit(lheight_adjust, "npc") - lheight, rep(lheight, nrow-1)),
    # widths = unit.c(unit(lwidth_adjust, 'npc') - lwidth, rep(lwidth, ncol-1)),
    top=textGrob(title, x=unit(0.05, 'npc'), just='left'))
}

# plot corn
corn_plots = grepl('corn', plot_these) %>% 
  plot_these[.] %>% 
  lapply(function(x) {
  df = sem_ls[[x]] %>% 
    extract2('df')
  barcode = ICs[[x]] %>% 
    names %>% 
    extract(1) %>% 
    strsplit('_') %>% 
    extract2(1) %>% 
    extract(3)
  contrast_name = strsplit(x, '_') %>% 
    extract2(1) %>% 
    extract(c(3,4)) %>% 
    paste(collapse='_')
  title = strsplit(x, '_') %>% 
    extract2(1) %>% 
    extract(2) %>% 
    name_parser[.]
  
  df[, 'ROTATION'] %<>%
    as.character %>% 
    gsub('CSCS', 'CS', .) %>% 
    gsub('W', 'Ww', .) %>% 
    factor(levels=names(palette))
  
  if (is.na(barcode)) {
    out = NULL
  } else {
    out = plot_axis_effect(barcode, 
                      'AV_DW',
                      'YEAR',
                      'ROTATION',
                      contrast_name,
                      df,
                      title=title,
                      palette=palette)
  }
  
  return(out)
})

if (save) {
  here('Results', figure_filename['corn']) %>% 
     jpeg(height=7.5, width=4, units='in', res=300, quality=100)
  grid_arrange_shared_legend(corn_plots, position='bottom')
  dev.off()
}


# plot soybean
soy_plots = corn_plots = grepl('soy', plot_these) %>% 
  plot_these[.] %>% 
  lapply(function(x) {
    df = sem_ls[[x]] %>% 
      extract2('df')
    barcode = ICs[[x]] %>% 
      names %>% 
      extract(1) %>% 
      strsplit('_') %>% 
      extract2(1) %>% 
      extract(3)
    contrast_name = strsplit(x, '_') %>% 
      extract2(1) %>% 
      extract(c(3,4)) %>% 
      paste(collapse='_')
    title = strsplit(x, '_') %>% 
      extract2(1) %>% 
      extract(2) %>% 
      name_parser[.]
    df[, 'ROTATION'] %<>%
      as.character %>% 
      gsub('CSCS', 'CS', .) %>% 
      gsub('W', 'Ww', .) %>% 
      factor(levels=names(palette))
    
    if (is.na(barcode)) {
      out = NULL
    } else {
      out = plot_axis_effect(barcode, 
                             'AV_DW',
                             'YEAR',
                             'ROTATION',
                             contrast_name,
                             df,
                             title=title,
                             palette=palette)
    }
    
    return(out)
  })
for (i in seq_len(length(soy_plots))) {
  if (is.null(soy_plots[[i]])) soy_plots[[i]] = NULL
}

if (save) {
  here('Results', figure_filename['soy']) %>% 
    jpeg(height=4, width=4, units='in', res=300, quality=100)
  grid_arrange_shared_legend(soy_plots, position='bottom')
  dev.off()
}
