####
# constrained ordinations of each community * crop * time point
# Components:
#   1. sites and labeled centeroids
#   2. convex hull?
#   3. axes - labeled with % (constrained) variance
#   4. annotations: Rsq, p-value
####

libs = c('here', 
              'ggplot2',
              'data.table',
              'magrittr',
              'grid',
              'gridExtra')

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

functions = c('Convenience Functions.R')#,
# 'Ordination Functions.R')
source_functions = function() {
  if (!require(here)) install.packages('here'); library(here)
  for (i in functions) source(here('Scripts', i))
}
source_functions()

load(here('Data', 
          'Microbes as Mechanism Results.Rdata'))  # loads ord_ls, among other objects

plot_these = 'None' %>% 
  grepl(., names(ord_ls)) %>%
  ord_ls[.]

# change order and colors. matches factor levels using names. 
palette = c(CS = '#A8AA39', # MEDIUM YELLOW,
            COWwS = '#6F1C1A', # DARK RED
            CPWwS = '#E56764', # LIGHT RED
            CSSwP = '#4C639A', # LIGHT BLUE
            CSSwSf = '#17264A') # DARK BLUE

# for the title of each sub-plot. Names are specified in names(ord_ls)
name_parser = c(seedling = 'Seedling',
                flowering = 'Flowering',
                ITS = 'Fungi',
                `16S` = 'Bacteria')
letter_parser = list(
  seedling = c(`16S` = 'a)',
               ITS = 'b)'),
  flowering = c(`16S` = 'c)',
                ITS = 'd)')
)

# for the title of each figure (group of plots). Names are crop specified in names(ord_ls).
crop_name = c(corn = 'Corn',
              soy = 'Soybean')

title = c()
for (i in 1:2){
  for (j in 3:4) {
    title = c(title,
              paste(name_parser[i], 
                    name_parser[j],
                    sep=', '))
  }
}
title %<>% paste(letters[1:4], ., sep=') ')

# for saving. Names are specified in names(ord_ls)
figure_filename = c(corn = 'FIG 2 - Corn Ordinations.jpeg',
                    soy = 'FIG 5 - Soybean Ordinations.jpeg')

#### Change below for different axes, etc ####


# title for each sub-plot, given the name of the object in ord_ls
generate_title = function(name, name_parser, letter) {
  name %<>% as.character
  nn = strsplit(name, '_')[[1]]
  out = name_parser[nn] %>% 
    na.omit
  letter_id = names(out)
  lets = letter_parser[[letter_id[1]]][letter_id[2]]
  
  out %<>% paste(collapse=', ') %>% 
    paste(lets, .)
  return(out)
}

# generate R-squared annotation
.annotate_rsq = function(x) {
  rsq = x$results$model_F['Rsq'] %>% 
    round(2) %>% 
    sprintf("%.2f", .)
  
  p_val = x$results$model_F['p_val'] %>% 
    as.numeric %>%
    signif(2)
  
  out = substitute(italic(R)^2 == a,#'\n'italic(p) == b,
                    list(a=rsq))
  return(out)
}

# generate p-value annotation
.annotate_pval = function(x) {
  p_val = x$results$model_F['p_val'] %>% 
    as.numeric %>%
    signif(2) %>% 
    round(4)
  
  oper = '='
  if (p_val < 0.0001) {
    p_val = 0.0001
    oper = '<'
  } else if (p_val < 0.01) {
    p_val %<>% sprintf('%.4f', .)
  } else if (p_val <= 0.1) {
    p_val %<>% sprintf('%.3f', .)
  }
  
  out = substitute(italic(p)~a~b,
                   list(a=oper,
                        b=p_val))
  return(out)
}

# calculate a convex hull for each group based on axes in coordinates
.chull_coords = function(x, coordinates=c(1,2)) {
  xy_values = x[, coordinates]
  vertices = chull(xy_values)
  out = x[vertices, ]
  return(out)
}

# generate axis labels with percent variance explained
.axis_label = function(x) {
  x %<>% ord_variance %>% 
    extract2('constrained')
  nn = colnames(x) %>% 
    strsplit('Axis') %>% 
    sapply(extract, 2)
  vv = x['pct_group', ] %>% 
    round(1)
  out = paste('Axis ', nn, ' [', vv, '%]', sep='')
  return(out)
}

# grid_arrange_shared_legend = function(plots, title=NULL) {
#   # plots <- list(...)
#   g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   grid.arrange(
#     do.call(arrangeGrob, lapply(plots, function(x)
#       x + theme(legend.position="none"))),
#     legend,
#     ncol = 1,
#     heights = unit.c(unit(0.95, "npc") - lheight, lheight), top=textGrob(title, x=unit(0.05, 'npc'), just='left'))
# }

# plot
fig2_plot = function(obj, 
                     axes=c(1,2), 
                     color='ROTATION', 
                     shape='YEAR',
                     title=NULL, 
                     palette=NA) {
  
  # Make dataframes:
  #   1. point data
  #   2. centeroid data
  #   3. axis variances
  obs = obj$observed
  spl = ord_scores(obs, 'row', 'principle', axes=axes, add_X=TRUE, add_Z=TRUE)
  
  if (!is.null(names(palette))) {
    spl[, color] %<>% 
      as.character %>% 
      gsub('CSCS', 'CS', .) %>% 
      gsub('W', 'Ww', .) %>% 
      factor(levels=names(palette))
  }

  ord_sum = ord_variance(obs)
  
  # axes
  X = colnames(spl)[axes[1]]
  Y = colnames(spl)[axes[2]]
  
  # calculate convex hull
  hulls = split(spl, spl[, color]) %>% 
    lapply(.chull_coords, axes) %>% 
    do.call(rbind, .)
  
  # make axis labs
  ll = .axis_label(obs)[axes]
  xlab = ll[1]
  ylab = ll[2]
  
  # x and y limits
  pltlims = c(
    min = min(spl$pcAxis1, spl$pcAxis2),
    max = max(spl$pcAxis1, spl$pcAxis2)
  )
  
  # annotation
  pltrange = diff(pltlims)
  y0 = pltlims[2] + pltrange*0.02
  x0 = pltlims[2] - pltrange*0.13
  y1 = y0 - pltrange*0.07
  
  rsq = .annotate_rsq(obj)
  pval = .annotate_pval(obj)

  # plot
  plt = ggplot() +
    coord_equal(xlim=pltlims,
                ylim=pltlims) +
    scale_fill_manual(values=palette) +
    scale_color_manual(values=palette) +
    geom_hline(aes(yintercept=0),
               size=0.4,
               color='gray') +
    geom_vline(aes(xintercept=0),
               size=0.4,
               color='gray') +
    geom_polygon(data=hulls,
                 aes_string(x=X,
                            y=Y,
                            fill=color,
                            color=NULL),
                 alpha=0.2,
                 show.legend=FALSE) +
    geom_point(data=spl,
               aes_string(x=X,
                          y=Y,
                          color=color,
                          shape=shape),
               size=2,
               alpha=0.8) +
    labs(x=xlab,
         y=ylab,
         title=title) +
    annotate('text', 
             x=x0, 
             y=y0, 
             label=rsq, 
             size=3, 
             color='#333333') +
    annotate('text', 
             x=x0, 
             y=y1, 
             label=pval, 
             size=3, 
             color='#333333') +
    theme_minimal() +
    theme(panel.grid=element_line(size=0.2),
          legend.title=element_blank(),
          legend.text=element_text(size=8),
          plot.title=element_text(size=8, hjust=0),
          axis.title=element_text(size=8),
          axis.text=element_text(size=8))
  return(plt)
}


# Generate Plots
tt = names(plot_these) %>% 
  lapply(function(x) {
    obj = plot_these[[x]]
    title = generate_title(x, name_parser)
    
    out = fig2_plot(obj, title=title, palette=palette)
    return(out)
  }) %>% 
  set_names(names(plot_these))


## Save Plots
for (i in names(figure_filename)) {
  to_plot = grepl(i, names(tt)) %>% 
    which %>% 
    tt[.]
  top_title = crop_name[i]
  
  # save to pdf
  here('Results', figure_filename[i]) %>% 
    jpeg(height=7.5, width=7.5, units='in', res=300, quality=100)
  grid_arrange_shared_legend(to_plot)#, top_title)
  dev.off()
  # grid.arrange(grobs=to_plot, top=textGrob(top_title))
}

