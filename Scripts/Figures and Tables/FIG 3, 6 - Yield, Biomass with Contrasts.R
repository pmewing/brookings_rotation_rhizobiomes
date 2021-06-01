#### Plot Figure 2: Biomass and Yield with Contrasts

functions = c('Convenience Functions.R')#,
# 'Ordination Functions.R')
source_functions = function() {
  if (!require(here)) install.packages('here'); library(here)
  for (i in functions) source(here('Scripts', i))
}
source_functions()

# Libraries
libs = c('magrittr', 
         'here',
         'ggplot2',
         'grid',
         'gridExtra')
load_libs(libs)

# Data
plot = 'Master Plot Data.csv' %>% 
  here('Data', .) %>% 
  read.csv
rownames(plot) = plot$SAMPLE_ID
plot[, 1:14] %<>% lapply(as.factor)

contrast_results = 'Contrast Results CROP.csv' %>% 
  here('Results', .) %>% 
  read.csv

# Names and such
pretty_crop = c(corn = 'Corn',
               soy = 'Soybean')
pretty_measure = c(seedling = 'Biomass, g/plant',
               flowering = 'Biomass, g/plant',
               yield = 'kg/ha')
pretty_cont = c(ROTATION_LENGTH = 'Rotation Diversity:\nFour vs Two Species',
               FOLLOW_CORN = 'Previous Crop Legacy:\nCorn vs Wheat',
               FOLLOW_SOY = 'Previous Crop Legacy:\nSoybean vs Other',
               WITH_PEA = 'Rotation With Pea:\nYes vs No',
               FOLLOW_SUNFLOWER = 'Previous Crop Legacy:\nSunflower vs Other')
pretty_title = c('Seedling', 
                  'Flowering',
                  'Yield,')  # yield
pretty_bartitle = paste(c('a)', 'c)', 'e)'), pretty_title, pretty_measure) %>% 
  set_names(names(pretty_measure))
pretty_efftitle = paste(c('b)', 'd)', 'f)'), pretty_title, pretty_measure) %>% 
  set_names(names(pretty_measure))

palette = c(CS = '#A8AA39', # MEDIUM YELLOW,
            COWwS = '#6F1C1A', # DARK RED
            CPWwS = '#E56764', # LIGHT RED
            CSSwP = '#4C639A', # LIGHT BLUE
            CSSwSf = '#17264A') # DARK BLUE

# Output
plot_out = c(corn = 'FIG 3 - Corn Yield, Biomass with Contrasts.jpeg',
             soy = 'FIG 6 - Soybean Yield, Biomass with Contrasts.jpeg')

#### Begin yield/bimass dataframe ####
facts = c('CROP', 'ROTATION', 'YEAR', 'SAMPLING')
# pltdf = subset(plot, CROP==crop, select=c('AV_DW', 'ROTATION', 'YEAR', 'SAMPLING', 'CROPYIELD'))
pltdf = plot[, c(facts, 'AV_DW')]
yld = plot[, c(facts, 'CROPYIELD')]
colnames(yld) %<>% gsub('CROPYIELD', 'AV_DW', .)
yld$SAMPLING = 'yield'
pltdf$CROPYIELD = NULL
pltdf = rbind(pltdf, yld) %>% 
  as.data.frame
pltdf %<>% na.omit

# Calculate median, min, max, and standard error. ggplot will calculate mean
ff = paste(facts, collapse='+') %>% 
  paste('AV_DW', ., sep='~') %>% 
  as.formula
center = pltdf %>%
  aggregate(ff, data=., median)
low = pltdf %>%
  aggregate(ff, data=., min) %>%
  set_colnames(c(facts, 'LOW'))
high = pltdf %>%
  aggregate(ff, data=., max) %>%
  set_colnames(c(facts, 'HIGH'))
err = pltdf %>% 
  aggregate(ff, 
            data=., 
            function(x) {
              sd(x)/sqrt(length(x)-1)
            }) %>%
  set_colnames(c(facts, 'STERR'))

# merge into a single dataframe (for barplots)
pltdf = merge(center, low, by=facts) %>% 
  merge(high, by=facts) %>% 
  merge(err, by=facts)

pltdf[, 'ROTATION'] %<>% 
  as.character %>% 
  gsub('CSCS', 'CS', .) %>% 
  gsub('W', 'Ww', .) %>% 
  factor(levels=names(palette))

# Format names
pltdf$MEASURE = 
  pltdf$SAMPLING %>% 
  as.character %>% 
  pretty_measure[.] %>% 
  factor(levels=unique(pretty_measure))
pltdf$TITLE = 
  pltdf$SAMPLING %>% 
  as.character %>% 
  pretty_bartitle[.] %>% 
  factor(levels=pretty_bartitle)

#### Format names for contrasts plotting dataframe ####
eff_df = contrast_results
eff_df$MEASURE =
  eff_df$SAMPLING %>% 
  as.character %>% 
  pretty_measure[.]
eff_df$CONTRAST %<>% 
  as.character %>% 
  pretty_cont[.] %>% 
  factor(levels=rev(pretty_cont))
eff_df$SIGNIF %<>% gsub('\\.', '^', .)
eff_df$TITLE =
  eff_df$SAMPLING %>% 
  as.character %>% 
  pretty_efftitle[.] %>% 
  factor(levels=pretty_efftitle)

# Plot
dark_gray = '#444444'
light_gray = '#AAAAAA'
signif_offset = -0.4
plt = lapply(names(pretty_crop), function(x) {
  bars = subset(pltdf, CROP==x) %>% 
    ggplot(aes(x=YEAR, 
               y=AV_DW, 
               fill=ROTATION)) + 
    scale_fill_manual(values=palette) +
    scale_color_manual(values=rep(dark_gray, 5)) +
    facet_wrap(~ TITLE,
               scales='free_y',
               ncol=1) +
    geom_col(position='dodge',
             width=0.8) +
    geom_errorbar(aes(x=YEAR,
                      y=AV_DW,
                      color=ROTATION,
                      ymin=AV_DW - STERR,
                      ymax=AV_DW + STERR),
                  width=0,
                  position=position_dodge(width=0.8),
                  show.legend=FALSE,
                  inherit.aes=FALSE) +
    labs(x=NULL,
         y=NULL,
         fill='Rotation') +
    # ggtitle(prettycrop[x]) +
    theme_minimal() +
    theme(legend.position='left',
          legend.title=element_blank(),
          legend.text=element_text(size=8),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_line(size=0.2),
          axis.title=element_text(size=8),
          axis.text=element_text(size=8),
          strip.text=element_text(hjust=0,
                                  size=8),
          panel.spacing.y=unit(2, 'lines'))
  
  effs = subset(eff_df, CROP==x) %>% 
      ggplot(aes(x=coefficients,
                 y=CONTRAST)) +
      facet_wrap(~TITLE,
                 scales='free_x',
                 ncol=1) +
      geom_vline(aes(xintercept=0),
                 color=light_gray,
                 size=0.4) +
      geom_point(size=4,
                 shape=18,
                 color=dark_gray) +
      geom_errorbarh(aes(xmin=coefficients-sigma,
                         xmax=coefficients+sigma),
                     height=0,
                     color=dark_gray) +
      geom_text(aes(label=SIGNIF),
                position=position_nudge(y=signif_offset),
                size=3) +
      xlab(NULL) +
      ylab(NULL) +
      # ggtitle(' ') +
      theme_minimal() +
      theme(panel.grid=element_line(size=0.2),
            legend.title=element_blank(),
            title=element_text(size=8),
            axis.title=element_text(size=8),
            axis.text=element_text(size=8),
            strip.text=element_text(hjust=0,
                                    size=8),
            panel.spacing.y=unit(2, 'lines')
      )
  
  out = list(bars, effs)
}) %>% 
  set_names(names(pretty_crop))

#### Export ####
for (i in names(plt)) {
  title = pretty_crop[i]
  here('Results', plot_out[i]) %>% 
  jpeg(width=7.5,
      height=5,
      units='in',
      res=300,
      quality=100)
  grid.arrange(grobs=plt[[i]], nrow=1)
  dev.off()
}
