libs = c('here', 
         'magrittr', 
         'reshape2',
         'lavaan')

for (i in libs) {
  if (!(require(i, character.only=TRUE))) {
    install.packages(i, Ncpu=4)
    library(i, character.only=TRUE)
  }
}


load(here('Data', 
          'Microbes as Mechanism Results.Rdata'))  # loads ord_ls, among other objects

load(here('Data',
          'Contrasts.Rdata'))

load(here('Data',
          'Crop ANOVAs.Rdata'))

load(here('Data',
          'ASV_Table_Filtering.Rdata'))

crop_contrast = read.csv(here('Results',
                              'Contrast Results CROP.csv'))

taxa_results %<>% subset(!(CROP=='soy' & SAMPLING=='flowering'))

# Table S3 - Microbial Community Data
here('Results', 'SUPPL TAB 3 - ASV Table Stats.txt') %>% 
  sink
cat('16S\n')
filtering_results$bac
cat('\nITS\n')
filtering_results$fun
sink()

# Table S5: Contrast Values of Treatments
contrasts %<>% 
  names %>% 
  lapply(function(x) {
    tt = contrasts[[x]] %>% 
      round(3) %>% 
      as.data.frame
    tt$CONTRAST = rownames(tt)
    tt$CROP = x
    out_cols = c('CROP', 'CONTRAST', 'CSCS', 'COWS', 'CPWS', 'CSSwP', 'CSSwSf')
    
    return(tt[, out_cols])
    
  }) %>% 
  do.call(rbind, .) %>% 
  set_rownames(NULL)

contrasts$CONTRAST %<>% gsub('_', ' ', .) %>% 
  tolower
names(contrasts) %<>% 
  gsub('CROP', 'Crop', .) %>% 
  gsub('CONTRAST', 'Contrast', .)

contrasts[, 3:7] %<>%
  apply(1, function(x) {
    mn = abs(x[x!=0]) %>% 
      min
    x = x / mn
    return(x)
}) %>% 
  t %>% 
  round

'SUPPL TAB 5 - Contrasts.csv' %>% 
  here('Results', .) %>% 
  write.csv(contrasts, ., row.names=FALSE)



# Table S6 - PERMANOVA results
# all results are in ord_ls
make_permanova = function(ord_ls) {
  permanova= sapply(ord_ls, function(x) x$results$model_F) %>%   # all permanove results
    round(3)
  Df = sapply(ord_ls, function(x) paste(x$Df[c('model_total', 'free')], collapse=','))  # degrees of freedom
  permanova = rbind(permanova, Df = Df) %>%
    t %>% 
    as.data.frame %>% 
    type.convert
  
  unpack_identifiers = function(id) {
    x = strsplit(id, '_')[[1]]
    
    nsplt = length(x)
    
    out = c(Barcode = x[3],
            Crop = x[1],
            Sampling = x[2],
            Constraint = paste(x[4:nsplt], collapse=' ')
    )
    return(out)
  }
  
  model_ID = rownames(permanova) %>% 
    sapply(unpack_identifiers) %>% 
    t
  
  permanova %<>% cbind(model_ID, .) %>% 
    set_rownames(NULL)
  
  permanova$Constraint %<>% gsub('None', 'ROTATION', .) %>% 
    tolower %>% 
    factor(levels=tolower(c('ROTATION', 'ROTATION LENGTH', 'FOLLOW SOY', 'FOLLOW CORN')))
  
  sorter = with(permanova, paste(Barcode, Crop, Sampling, as.numeric(Constraint))) %>% 
    order
  permanova %<>% .[sorter, ]
  
  drop_cols = sapply(colnames(permanova), function(x) x %in% c('var_fitted', 'var_resid', 'Rsq_adj'))
  permanova %<>% .[!drop_cols]
  
  return(permanova)
}
permanova_results = make_permanova(ord_ls)

"SUPPL TAB 6 - PERMANOVA Results.csv" %>% 
  here('Results', .) %>% 
  write.csv(permanova_results, ., row.names=FALSE)

# Table S7: Crop ANOVA Statistics
names(crop_anova) %<>% 
  strsplit('_') %>% 
  lapply(function(x) {
    paste(x[1], x[3]) %>% 
      gsub('seedling', 'seedling biomass', .) %>% 
      gsub('flowering', 'flowering biomass', .) %>% 
      tolower
  })

reorder = sort(names(crop_anova))
crop_anova %<>% .[reorder]
 
'SUPPL TAB 7 - Crop ANOVA results.txt' %>% 
  here('Results', .) %>% 
  sink
for (i in names(crop_anova)) {
  cat(i)
  cat('\n\n')
  print(as.data.frame(crop_anova[[i]]))
  cat('\n\n')
}
sink()

# Table S8: Crop Contrast Results
crop_contrast$CONTRAST %<>% 
  gsub('_', ' ', .) %>% 
  tolower
crop_contrast %<>% 
  lapply(function(x) {
    if (is.numeric(x)) x %<>% round(3)
    return(x)
  }) %>% 
  as.data.frame

names(crop_contrast) %<>%
  tolower %>% 
  sapply(function(x) {
    substr(x, 1, 1) %<>% toupper
    return(x)
  }) %>% 
  gsub('Tstat', 't_Stat', .) %>% 
  gsub('Pvalues', 'p_Value', .) %>% 
  gsub('Coefficients', 'Difference', .) %>% 
  gsub('Sigma', 'Std_Error', .)

crop_contrast$Signif %<>% gsub('\\.', '^', .)

reorder = with(crop_contrast, paste(Crop, Sampling, Contrast)) %>% 
  order
crop_contrast %<>% .[reorder, ]

'SUPPL TAB 8 - Crop Contrast Results.csv' %>% 
  here('Results', .) %>% 
  write.csv(crop_contrast, ., row.names=FALSE)

# Table S10: SEM Results
measures = c(`Log-likelihood` = 'logl',
             Parameters = 'npar',
             Observations = 'ntotal',
             AIC = 'aic',
             BIC = 'bic',
             CFI = 'cfi',
             `Chi-squared` = 'chisq',
             `Deg freedom` = 'df',
             `p value` = 'pvalue')
sem_tests = lapply(sem_ls, function(x) {
  out = sapply(x, function(y) {
    if (class(y) == 'lavaan') {
      vec = fitmeasures(y, fit.measures=measures) %>% 
        set_names(names(measures))
      return(vec)
    }
  }) %>% 
    do.call(rbind, .) %>% 
    as.data.frame
  out = data.frame(model=rownames(out),
                   out) %>% 
    set_rownames(NULL)
  return(out)
})
sem_tests %<>% 
  names %>% 
  lapply(function(x) {
    vals = sem_tests[[x]]
    out = data.frame(model_group = x,
                     vals)
    return(out)
  }) %>% 
  do.call(rbind, .) %>% 
  set_rownames(NULL) %>% 
  lapply(unlist) %>% 
  as.data.frame(stringsAsFactors=FALSE)

sem_tests[, c('Log.likelihood', 'AIC', 'BIC')] %<>% lapply(function(x) as.numeric(x) %>% round(2))
sem_tests[, c('Chi.squared', 'p.value', 'CFI')] %<>% lapply(function(x) as.numeric(x) %>% round(3))

sem_tests$crop = sem_tests$model_group %>% 
  as.character %>% 
  strsplit('_') %>% 
  sapply(extract, 1)
sem_tests$sampling = sem_tests$model_group %>% 
  as.character %>% 
  strsplit('_') %>% 
  sapply(extract, 2)
sem_tests$contrast = sem_tests$model_group %>% 
  as.character %>% 
  strsplit('_') %>% 
  sapply(extract, -c(1, 2)) %>% 
  sapply(paste, collapse=' ') %>% 
  tolower

drop_rows = sem_tests$contrast=='its combined'
sem_tests %<>% .[!drop_rows, ]
sem_tests$contrast %<>%
  sapply(function(x) ifelse(nchar(x)==0, '--', x))

#
mediation_parser = list(
  no_mediation = c('None', 'None'),
  partial_mediation_ITS = c('Partial', 'None'),
  partial_mediation_16S = c('None', 'Partial'),
  partial_mediation_full = c('Partial', 'Partial'),
  full_mediation_ITS = c('Full', 'None'),
  full_mediation_16S = c('None', 'Full'),
  full_mediation_full = c('Full', 'Full'),
  no_microbe = c('None', 'None'),
  full_microbe = c('Full, Seedling', 'Full, Flowering'),
  biomass_only = c('None', 'None')
) %>% 
  do.call(rbind, .) %>% 
  set_colnames(c('mediation by ITS', 'mediation by 16S'))

mediation_parser %<>% .[as.character(sem_tests$model), ]


# alter value for 16S in soy - not tested. 
rn = which(sem_tests$crop == 'soy')
mediation_parser[rn, 'mediation by 16S'] = '--'

# alter value for 16s in corn seedling follow soy - not tested
rn = with(sem_tests, crop=='corn' & sampling=='seedling' & contrast=='follow soy') %>% 
  which
mediation_parser[rn, 'mediation by 16S'] = '--'

sem_tests %<>% 
  data.frame(mediation_parser,
             stringsAsFactors=FALSE)

sem_tests = apply(sem_tests[, c('crop', 'sampling', 'contrast')], 1, paste, collapse='') %>% 
  split(sem_tests, .) %>% 
  lapply(function(x) {
    val = rep('', nrow(x))
    mn = which.min(x$AIC)
    val[mn] = '*'
    x$selected = val
    return(x)
  }) %>% 
  do.call(rbind, .) %>% 
  set_rownames(NULL)

keep_cols = c('crop', 'sampling', 'contrast', 'mediation.by.ITS', 'mediation.by.16S', 
              'selected', 
              'Log.likelihood', 'Parameters', 'Observations', 'AIC', 'BIC', 'CFI', 'Chi.squared', 'Deg.freedom', 'p.value')

sem_tests %<>% .[, keep_cols]
names(sem_tests) %<>%
  gsub('\\.', ' ', .) %>% 
  gsub('Log likelihood', 'log-likelihood', .) %>% 
  # gsub('npar', 'parameters', .) %>% 
  # gsub('ntotal', 'observations', .) %>% 
  gsub('Chi squared', 'chi-squared', .) %>% 
  gsub('Deg freedom', 'deg freedom', .) %>% 
  sapply(function(x) {
    substr(x, 1, 1) %<>% toupper
    return(x)
  }) %>% 
  gsub('P value', 'p value', .)

'SUPPL TAB 10 - SEM Statistics.csv' %>% 
  here('Results', .) %>% 
  write.csv(sem_tests, ., row.names=FALSE, na='')

# Table S11: Candidate Taxa
"SUPPL TAB 11 - Candidate Taxa.csv" %>% 
  here('Results', .) %>% 
  write.csv(taxa_results, ., row.names=FALSE)
