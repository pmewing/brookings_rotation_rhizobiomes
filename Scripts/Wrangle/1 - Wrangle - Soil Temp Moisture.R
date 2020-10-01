# Processing 2016 field data (soil temp/moisture) - corn/soybean, T1/T2


# Levels (so far):
# Crop: Corn, Soybean
# Time: T1, T2
# Year: 2016

#### Setup ####
{
  libs = c('lattice',
           'magrittr', 'here', 'reshape2')  # workflow
  
  load_libs = function() {
    for(i in libs) library(i, character.only=TRUE)
  }
  load_libs()
  
  # Directories
  dat_dir = file.path("Data", "Raw", "Temp Moist")
  gen_df_filename = function(crop, year, time) {
    require(here)
    here(dat_dir, 
         paste(crop, year, time, "Field.csv"))
  }
  
  results_dir = 'Data'
   
 
  # Load a dataframe for column names
  master_dir = file.path('Data', 'Raw')
  master_df = '2016_plot_data.csv'
  
  # Convenience functions
  rename_columns = function(data, x) {
    # data = data.frame or matrix
    # x = data.frame or matrix. 
    # First column = old name.
    # Second column = new name.
    
    x = matrix(x, ncol=2)
    nn = colnames(data)
    
    for (i in 1:nrow(x)) {
      nn = gsub(x[i, 1], x[i, 2], nn)
    }
    # return(nn)
    
    colnames(data) = nn
    return(data)
  }
  
  mean_na = function(x) mean(x, na.rm=TRUE)
  var_na = function(x) var(x, na.rm=TRUE)
  
  # Load master df for names
  df = here(master_dir, master_df) %>%
    read.csv(stringsAsFactors=FALSE)
  names(df) %<>% toupper


  # convert filename to df standard
  sampling_key = c(T1 = 'seedling',
                   T2 = 'flowering')

  # identifying variables
  keep_factors = c('CROP', 'PLOT', 'REP', 'YEAR', 'SAMPLING')

  df_list = list()  # for conveniently saving individual dataframes
}

#### Corn T1 2016 ####
{
  crop = 'Corn'
  time = 'T1'
  year = '2016'
  
  # Load
  dd = gen_df_filename(crop, year, time) %>% 
    read.csv(stringsAsFactors=FALSE)
  
  colnames(dd) %<>%
    toupper
  
  # ID Name Changes
  names(dd)
  names(df)
  new_names = rbind(
    c('PLOT_NUMBER', 'PLOT'),
    c('TEMP', 'TEMPERATURE'),
    c('SOIL_MOISTURE_RWC75', 'VMOISTURE'),
    c('SHOOT_P_MOISTURE', 'SHOOT_MOISTURE')
  )
  
  dd %<>% rename_columns(new_names)
  
  # Make data changes (ex. renaming reps)
  dd[, "REP"] %<>% paste0('R', .)
  dd[, 'YEAR'] = year
  dd[, 'SAMPLING'] = sampling_key[time]
  dd[, "CROP"] %<>% tolower
  dd[, "SHOOT_MOISTURE"] = (dd[, "SHOOT_FW"] - dd[, "SHOOT_DW"]) / dd[, "SHOOT_FW"]
  
  extra_factors = c('PLANT_NUMBER')
  keep_variables = c('TEMPERATURE', 'VMOISTURE', 'SHOOT_HEIGHT', 'SHOOT_MOISTURE')
  dd[, c(keep_factors, extra_factors)] %<>% lapply(as.factor)
  
  # for reshaping the dataframe
  cast_formula = paste(keep_factors, collapse="+") %>% 
    paste0("~variable") %>% 
    formula
  
  means = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=mean_na)
  
  variances = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=var_na)
  names(variances) %<>% sapply(function(x) {
      if (is.numeric(variances[, x])) {
        x = paste0(x, "_VAR")
      }
      return(x)
    })
  
  # Add to list of dataframes
  nn = paste(crop, time, year)
  df_list[[nn]] = merge(means, variances, by=keep_factors)
}

#### Corn T2 2016 ####
{
  crop = 'Corn'
  time = 'T2'
  year = '2016'
  
  # Load
  dd = gen_df_filename(crop, year, time) %>% 
    read.csv(stringsAsFactors=FALSE, skip=1)
  
  colnames(dd) %<>%
    toupper
  
  # ID Name Changes
  head(dd)
  names(df)
  new_names = rbind(
    # c('PLOT_NUMBER', 'PLOT'),
    c('TEMP.C', 'TEMPERATURE'),
    c('STD.VWC', 'VMOISTURE'),
    c('SHOOT.FW', 'SHOOT_FW')
  )
  
  dd %<>% rename_columns(new_names)
  names(dd) %<>% gsub('\\.', '', .)
  
  # Make data changes
  dd[, 'REP'] = dd[, "PLOT"] %>% 
    substr(1, 1) %>% 
    paste0('R', .)
  dd[, 'YEAR'] = year
  dd[, 'SAMPLING'] = sampling_key[time]
  dd[, 'CROP'] = tolower(crop)
  
  
  extra_factors = c('PLANT')
  keep_variables = c('TEMPERATURE', 'VMOISTURE', 'SHOOT_FW')
  dd[, c(keep_factors, extra_factors)] %<>% lapply(as.factor)
  
  # for reshaping the dataframe
  cast_formula = paste(keep_factors, collapse="+") %>% 
    paste0("~variable") %>% 
    formula
  
  means = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=mean_na)
  
  variances = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=var_na)
  names(variances) %<>% sapply(function(x) {
    if (is.numeric(variances[, x])) {
      x = paste0(x, "_VAR")
    }
    return(x)
  })
  
  tmp = df[, c(keep_factors, 'AV_DW')]
  tmp[, keep_factors] %<>% lapply(as.factor)
  means %<>% merge(tmp, by=keep_factors)
  means[, "SHOOT_MOISTURE"] = (means[, "SHOOT_FW"] - means[, "AV_DW"]) / means[, "SHOOT_FW"]
  means[, c('AV_DW', 'SHOOT_FW')] = NULL
  
  variances[, "SHOOT_FW_VAR"] = NULL
  
  # Add to list of dataframes
  nn = paste(crop, time, year)
  df_list[[nn]] = merge(means, variances, by=keep_factors)
}

#### Soybean T1 2016 ####
{
  crop = 'Soybean'
  time = 'T1'
  year = '2016'
  
  # Load
  dd = gen_df_filename(crop, year, time) %>% 
    read.csv(stringsAsFactors=FALSE)
  colnames(dd) %<>%
    toupper
  
  # ID Name Changes
  names(dd)
  names(df)
  new_names = rbind(
    c('PLOT.NUMBER', 'PLOT'),
    c('TEMP.', 'TEMPERATURE'),
    c('STD_VWC', 'VMOISTURE'),
    c('SHOOT_P_MOISTURE', 'SHOOT_MOISTURE'),
    c('PLANT.NUMBER', 'PLANT_NUMBER')
  )
  
  dd %<>% rename_columns(new_names)
  
  # Make data changes (ex. renaming reps)
  dd[, "REP"] %<>% paste0('R', .)
  dd[, 'YEAR'] = year
  dd[, 'SAMPLING'] = sampling_key[time]
  dd[, "CROP"] %<>% tolower
  dd[, "SHOOT_MOISTURE"] = (dd[, "SHOOT_FW"] - dd[, "SHOOT_DW"]) / dd[, "SHOOT_FW"]
  dd[, "NODULATION"] %<>% 
    gsub(" ", "", .) %>% 
    equals('Y') %>% 
    as.numeric
  
  extra_factors = c('PLANT_NUMBER')
  keep_variables = c('TEMPERATURE', 'VMOISTURE', 'SHOOT_HEIGHT', 'SHOOT_MOISTURE', 'NODULATION')
  dd[, c(keep_factors, extra_factors)] %<>% lapply(as.factor)
  
  # for reshaping the dataframe
  cast_formula = paste(keep_factors, collapse="+") %>% 
    paste0("~variable") %>% 
    formula
  
  means = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=mean_na)
  
  variances = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=var_na)
  names(variances) %<>% sapply(function(x) {
    if (is.numeric(variances[, x])) {
      x = paste0(x, "_VAR")
    }
    return(x)
  })
  variances[, "NODULATION_VAR"] = NULL
  
  # Add to list of dataframes
  nn = paste(crop, time, year)
  df_list[[nn]] = merge(means, variances, by=keep_factors)
}

#### Soybean T2 2016 ####
{
  crop = 'Soybean'
  time = 'T2'
  year = '2016'

  
  dd = gen_df_filename(crop, year, time) %>% 
    read.csv(stringsAsFactors=FALSE, skip=1)
  
  colnames(dd) %<>%
    toupper
  
  # ID Name Changes
  names(dd)
  names(df)
  new_names = rbind(
    c('PLOT.NUMBER', 'PLOT'),
    c('TEMP.', 'TEMPERATURE'),
    c('STD_VWC', 'VMOISTURE'),
    c('SHOOT_P_MOISTURE', 'SHOOT_MOISTURE'),
    c('PLANT.NUMBER', 'PLANT_NUMBER')
  )
  
  dd %<>% rename_columns(new_names)
  
  # Make data changes (ex. renaming reps)
  dd[, "REP"] %<>% paste0('R', .)
  dd[, 'YEAR'] = year
  dd[, 'SAMPLING'] = sampling_key[time]
  dd[, "CROP"] %<>% tolower
  dd[, "SHOOT_MOISTURE"] = (dd[, "SHOOT_FW"] - dd[, "SHOOT_DW"]) / dd[, "SHOOT_FW"]
  dd[, "NODULATION"] %<>% 
    gsub(" ", "", .) %>% 
    equals('Y') %>% 
    as.numeric
  
  extra_factors = c('PLANT_NUMBER')
  keep_variables = c('TEMPERATURE', 'VMOISTURE', 'SHOOT_MOISTURE', 'NODULATION')
  dd[, c(keep_factors, extra_factors)] %<>% lapply(as.factor)
  
  # for reshaping the dataframe
  cast_formula = paste(keep_factors, collapse="+") %>% 
    paste0("~variable") %>% 
    formula
  
  means = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=mean_na)
  
  variances = dd[, c(keep_factors, extra_factors, keep_variables)] %>% 
    melt %>% 
    dcast(cast_formula, fun.aggregate=var_na)
  names(variances) %<>% sapply(function(x) {
    if (is.numeric(variances[, x])) {
      x = paste0(x, "_VAR")
    }
    return(x)
  })
  variances[, "NODULATION_VAR"] = NULL
  
  # Add to list of dataframes
  nn = paste(crop, time, year)
  df_list[[nn]] = merge(means, variances, by=keep_factors)
}

#### Combine ####

# Ensure each df has the same data
ls_data = sapply(df_list, function(x) {
  tt = sapply(x, is.numeric)
  tt = names(x)[tt]
  return(tt)
}) %>% 
  unlist %>% 
  unique %T>%
  print

df_list %<>% lapply(function(x) {
  for (i in ls_data) {
    tt = i %in% names(x)
    if (!tt) x[, i] = NA
  }
  return(x)
})

# Reorder each df
df_list %<>% lapply(function(x) {
  out = x[, c(keep_factors, ls_data)]
  return(out)
})

# Combine
out = do.call(rbind, df_list) %>% 
  data.frame %>% 
  set_rownames(NULL)
out[, keep_factors] %<>% lapply(as.character)
out[, "CROP"] %<>% gsub('soybean', 'soy', .)

#### Export ####

ww = FALSE
if (ww) {
  fileout = '2016 Soil Temperature Moisture.csv'
  write.csv(out,
            here(results_dir, fileout),
            row.names=FALSE)
}


# out = merge(df, out, by=keep_factors, all=TRUE)
# 
# #### Export ####
# fileout = "Master Plot Data.csv"
# 
# ww = FALSE
# if (ww) write.csv(out, 
#                   here(results_dir, fileout),
#                   row.names=FALSE)
