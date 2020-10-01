# Munge 2016, 2017 data

##############
#### Load ####
##############

# Libraries
libs = c('here', 'lattice', 'magrittr')

for (i in libs) {
  library(i, character.only=TRUE)
}

# Helper Functions
functions = here("Scripts", 'Wrangle', "Wrangle - Helper Functions.R")
source(functions)

# Directories
dir = 'Data'
dir_raw = 'Raw'

# Plot data - soils, nutrients, biomass, descriptors, etc.
in_plot17 = '2017_plot_data.txt'                  # Raw upload from Soledad
in_plot16 = '2016_plot_data.csv'                  # Raw upload, converted to CSV manually

# Temp/moisture data, including plant height, nodulation
in_temp16 = '2016 Soil Temperature Moisture.csv'  # from Wrangle - Soil Temp Moisture.R

# Fungal taxa tables, for the fungal sample key. Old data, so we're only after the sample keys.
dir_fun = file.path('Fungi', 'OLD Fungi')             # Retained for sample keys. The OTU table is from an older pipeline (pre-removal of ITSExpress step)
in_fun16 = 'all_data_2016_Fungi_Jan2020.txt'          # Same as 2016 plot data, but containing fungal sample keys (and an old, CSS-normalized version of the OTU table)
in_fun17 = 'all_data_2017_Fungi_Jan2020.txt'          # Same as 2017 plot data, but containing fungal sample keys (and a old, CSS-normalized version of the OTU table)


# Bacterial taxa tables, for the bacterial sample key
dir_bac = 'Bacteria'
in_bac16 = 'all_data_2016_2019analysis.txt'           # Same as 2016 plot data, but containing 16S sample keys (and a CSS-normalized version of the OTU table)
in_bac17 = 'all_data_2017_2019analysis.txt'           # Same as 2017 plot data, but containing 16S sample keys (and a CSS-normalized version of the OTU table)

# all counts
in_counts = "masterdatasheet_counts.txt"  # contains raw counts of the (old) ITS and current (only) 16S, but with different ASV names to convert later

# export?
ww = FALSE

###################
#### Plot Data ####
###################

plot16 = here(dir, dir_raw, in_plot16) %>% 
  read.csv(stringsAsFactors=FALSE) %>% 
  format_names

# 2017 doesn't play nicely with read.table, so will try manual assembly
plot17 = here(dir, dir_raw, in_plot17) %>% 
  read_tough_files %>% 
  format_names

# Find mismatches between names
mismatches = mismatch_names(plot16, plot17) %T>% print

# Changes to make:
# 1. Remove "X" from a number of soil measures in master df
# 2. Change stages to match physiology ex. STAGE_1 becomes STAGE_V2)?
# 3. Convert spaces and periods to underscores

names(plot16) %<>% 
  sapply(function(x) {
    if (substr(x, 1, 1) == 'X') {
      x = substr(x, 2, 1E6)
    }
    return(x)
  }) %>%
  gsub('STAGE_1', 'STAGE_V2', .) 

names(plot17) %<>% 
  gsub('STAGE_3', 'STAGE_V7', .) %>% 
  gsub('STAGE_4', 'STAGE_R3', .) %>% 
  gsub('CURRENTCROP', 'CROP', .)

mismatches = mismatch_names(plot16, plot17) %T>% print


#############################
#### Check factor levels ####
#############################

factors = names(plot16)[1:12] %T>% print

# Set factors all as character
plot16[factors] %<>% lapply(as.character)
plot17[factors] %<>% lapply(as.character)


factor_levels = 
  lapply(factors, function(x) {
    in_plot16 = unique(plot16[, x])
    in_plot17 = unique(plot17[, x])
    out = mismatch_names(in_plot16, in_plot17)
    return(out)
  }) %>% 
  set_names(factors) %T>%
  print

# Changes to make:
# 1. sampleid in plot16 is bacteria id
# 3. update sample_id to reflect crop, year, plot, time

plot16[, 'BACTERIA_ID'] = plot16[, 'SAMPLEID']
plot17[, 'BACTERIA_ID'] = NA

plots = rbind(plot16,
              plot17) %>% 
  as.data.frame %>% 
  make_sampleID


# Ensure all rotation flags are correct - one COWS is labeled as 2-year, for example.
plots[, c('CROP', 'ROTATION', 'PREVIOUS_CROP', 'YEARS', 'TYPE', 'ROT')] %>% 
  unique

plots[plots$ROTATION=='COWS', 'YEARS'] = 'Four'

plots[, c('CROP', 'ROTATION', 'PREVIOUS_CROP', 'YEARS', 'TYPE', 'ROT')] %>% 
  unique

#################################
#### Format temperature data ####
#################################

temp16 = here(dir, dir_raw, in_temp16) %>% 
  read.csv(stringsAsFactors=FALSE) %>% 
  make_sampleID

temp16[, 1:5] %<>% lapply(as.character)

temp = temp16  # for easy incorporation later

# Reduce to sample id and data
temp %<>% data_and_key

#############################################
#### Load and format bacteria sample IDs ####
#############################################

bac16_raw = here(dir, dir_raw, dir_bac, in_bac16) %>%
  read.table(header=TRUE, stringsAsFactors=FALSE) %>%
  format_names

bac16 = bac16_raw[, factors]

bac16[, 'BACTERIA_ID'] = bac16[, 'SAMPLEID']
bac16 %<>%
  make_sampleID(cleanup=FALSE) %>%  # will need 'SAMPLEID' later
  extract(c('SAMPLEID', 'SAMPLE_ID', 'BACTERIA_ID'))

bac17_raw = here(dir, dir_raw, dir_bac, in_bac17) %>%
  read_tough_files %>%
  format_names

bac17 = bac17_raw[, factors]

tag = sapply(bac17[, 'SAMPLEID'], function(x) {
  ss = strsplit(x, "_")[[1]]
  paste('B',
        substr(ss[1], 1, 1),
        substr(ss[2], 2, 10000),
        '.',
        sep="") %>%
    toupper
})

bac17[, 'BACTERIA_ID'] = paste0(tag, bac17[, 'SAMPLEID'])
bac17 %<>%
  make_sampleID(cleanup=FALSE) %>%
  extract(c('SAMPLEID', 'SAMPLE_ID', 'BACTERIA_ID'))

bac = rbind(bac16, bac17) %>%
  as.data.frame


##########################################
#### Load and format fungi sample IDs ####
##########################################

fun16_raw = here(dir, dir_raw, dir_fun, in_fun16) %>%
  read_tough_files %>%
  format_names

fun16 = fun16_raw[, factors]

fun16[, 'FUNGI_ID'] = fun16[, 'SAMPLEID'] %>%
  gsub("_", ".", .)

fun16 %<>%
  make_sampleID(cleanup=FALSE) %>%
  extract(c('SAMPLEID', 'SAMPLE_ID', 'FUNGI_ID'))

fun17_raw = here(dir, dir_raw, dir_fun, in_fun17) %>%
  read_tough_files %>%
  format_names

fun17 = fun17_raw[, factors]

tag = sapply(fun17[, 'SAMPLEID'], function(x) {
  ss = strsplit(x, '_')[[1]]
  paste('F',
        substr(ss[1], 1, 1),
        substr(ss[2], 2, 1000),
        ".",
        sep="") %>%
    toupper
})

fun17[, 'FUNGI_ID'] = paste0(tag, fun17[, 'SAMPLEID'])
fun17 %<>%
  make_sampleID(cleanup=FALSE) %>%
  extract(c('SAMPLEID', 'SAMPLE_ID', 'FUNGI_ID'))

fungi = rbind(fun16, fun17) %>%
  as.data.frame


###############
#### Merge ####
###############

# plots + temp. Drops one row from temp (missing in plot data)
df = merge(plots, temp, by='SAMPLE_ID', all.x=TRUE)
dim(df)

# data + bacteria
df = merge(bac, df, by='SAMPLE_ID', all.y=TRUE)
dim(df)
# ensure the bacteria_id info from 2016 plot data matches the OTU table.
tt = cbind(df[, 'BACTERIA_ID.x'], df[, 'BACTERIA_ID.y']) %>% 
  na.omit
sum(tt[, 1] != tt[, 2])
# clean up
df[, 'BACTERIA_ID.y'] = NULL
names(df) %<>% gsub('ID.x', 'ID', .)

# data + fungi
df = merge(fungi, df, by='SAMPLE_ID', all.y=TRUE)

# clean up
tt = grepl('SAMPLEID', colnames(df))
df[, tt] = NULL


#################################
#### Export Master Dataframe ####
#################################

file_out = 'Master Plot Data.csv'


if (ww) {
  write.csv(df,
            here(dir, file_out),
            row.names=FALSE)
}


############################################################
#### Munge microbials into multi-year tables by kingdom ####
############################################################

counts = here(dir, dir_raw, in_counts) %>% 
  read_tough_files %>% 
  format_names %>% 
  make_sampleID

counts = merge(df[, c('SAMPLE_ID', 'BACTERIA_ID', 'FUNGI_ID')], 
               counts, 
               by='SAMPLE_ID')

# subset bacterial ASVs
bac = counts
rownames(bac) = bac[, 'BACTERIA_ID']

tt = grepl('BASV', colnames(bac))
bac = bac[, tt] %>% 
  as.matrix

colnames(bac) %<>% gsub('BASV', 'ASV_', .)

# won't repeat for fungi because this is the old data


####################################
#### Export Sequencing Matrices ####
####################################

out_bac = 'Master 16S Counts.csv'


if (ww) {
  here(dir, out_bac) %>% 
    write.csv(bac, .)
}

