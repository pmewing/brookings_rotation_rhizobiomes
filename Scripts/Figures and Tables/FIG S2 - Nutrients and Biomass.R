functions = c('Convenience Functions.R')#,
# 'Ordination Functions.R')
source_functions = function() {
  require(here)
  for (i in functions) source(here('Scripts', i))
}
source_functions()

# Libraries
libs = c('magrittr', 
         'here',
         'ggplot2',
         'gridExtra',
         'reshape2',
         'car')
load_libs(libs)

# Data
plot = 'Master Plot Data.csv' %>% 
  here('Data', .) %>% 
  read.csv
rownames(plot) = plot$SAMPLE_ID
plot[, 1:14] %<>% lapply(as.factor)

# Names and such
prettycrop = c(corn = 'Corn',
               soy = 'Soybean')
prettytime = c(seedling = 'Seedling',
               flowering = 'Flowering')
prettyvar = c(AV_DW = 'Biomass',
              P = 'Phosphorus',
              K = 'Potassium',
              N = 'Nitrogen')

palette = c('#E56764',  # dark red
            'white',  # light gray
            '#4C639A')  # dark blue

# Output
plot_out = 'FIG S2 - Biomass with Nutrients.jpeg'

#### Begin ####
crops = c('corn', 'soy')
times = c('seedling', 'flowering')

#outliers: that over-leverage a linear model of biomass ~ nutrients
is_outlier = c('corn_2016_118_t1',  # this is the second-highest biomass
               'corn_2016_336_t1')   # this has the highest biomass and highest potassium
is_outlier = rownames(plot) %in% is_outlier
plot %<>% .[!is_outlier, ]

pltdf = lapply(crops, function(cc) {
  tmp = lapply(times, function(tt) {
    df = subset(plot, CROP == cc & SAMPLING == tt, select = c(names(prettyvar), 'YEAR'))
    
    # remove year effects by centering/scaling within year
    df %<>% 
      split(df$YEAR) %>% 
      lapply(extract, names(prettyvar)) %>% 
      lapply(scale) %>% 
      do.call(rbind, .) %>% 
      cor(use='complete.obs')
    
    is_lower = lower.tri(df, diag=FALSE)
    df[is_lower] = 0
    df %<>% melt
    
    df$CROP = cc
    df$SAMPLING = tt
    return(df)
  })
  tmp = do.call(rbind, tmp)
}) %>% 
  do.call(rbind, .)

    
pltdf[, 'Var1'] %<>% prettyvar[.] %>% 
  factor(levels=prettyvar)
pltdf[, 'Var2'] %<>% prettyvar[.] %>% 
  factor(levels=rev(prettyvar))
pltdf[, 'CROP'] %<>% prettycrop[.] %>% 
  factor(levels=prettycrop)
pltdf[, 'SAMPLING'] %<>% prettytime[.] %>% 
  factor(levels=prettytime)
pltdf$LABEL = sapply(pltdf$value, function(x) ifelse(x==0, '', sprintf('%.2f', round(x, 2))))

plt = ggplot(data=pltdf) +
  facet_grid(CROP ~ SAMPLING, switch='y') +
  coord_equal() +
  scale_fill_gradient2(low=palette[1],
                       mid=palette[2],
                       high=palette[3],
                       limits=c(-1, 1)) +
  geom_tile(aes(x=Var1,
                  y=Var2,
                  fill=value),
            show.legend=FALSE) +
  geom_text(aes(x=Var1,
                 y=Var2,
                 label=LABEL),
             size=2.5) +
  labs(x=NULL,
       y=NULL) +
  theme_minimal() +
  theme(panel.grid=element_line(size=0.2),
        # legend.title=element_blank(),
        title=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text=element_text(size=8),
        axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        strip.placement='outside'
  )

here('Results', plot_out) %>% 
  jpeg(height=3.5, width=4, units='in', res=300, quality=100)
plt
dev.off()
