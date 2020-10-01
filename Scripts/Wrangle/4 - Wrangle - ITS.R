# Soledad re-processed the ITS data to squeeze a few more reads out (by dropping an ITSExpress step)
# This script will:
#   1. Load the new count and taxonomy tables
#   2. Harmonize row names with FUNGI_ID
#   3. Process a new funguild table
#   4. Export



library(here)
library(magrittr)

here('Scripts', 'Convenience Functions.R') %>% 
  source

in_dir = file.path('Raw', 'Fungi', 'NEW fungi')
in_count = 'fungi_2016_2017_counts.txt'
in_tax = 'taxfungi_2016_2017_counts.txt'
in_envi = 'Master Plot Data.csv'

raw_count = here('Data', in_dir, in_count) %>% 
  read.table(row.names=1, stringsAsFactors=FALSE)

raw_tax = here('Data', in_dir, in_tax) %>% 
  read.table(row.names=1, stringsAsFactors=FALSE)

envi = here('Data', in_envi) %>% 
  read.csv(stringsAsFactors=FALSE)


# The new data has two formats:
#   - crop_sampling_plot_primer_barcode (2017 format)
#   - PrimerCropSampling-plot_barcode_lane (as AA#-###_S###_L001) (2016 format)
#
# The old data's versions of this are:
#   - PrimerCropSampling.crop_sampling_plot (2017 format)
#   - PrimerCropSampling.plot (2016 format)
#
# Additionally, the new data is more complete, so will need to merge this in to the data.frame.
#
# Approach:
#   1. Split new data into the two formats, based on whether it contains 'FC' or 'FS', or 'corn' or 'soy'.
#   2. For the 2017 format:
#     - split on '_', 
#     - keep the first 3 terms
#     - concatenate with '_'
#     - Construct a prefix based on crop and sampling
#     - concatenate with '.'
#     - use as rownames for 2017 format
#     - Make a data.frame including year, sampling, crop, and plot information.
#   3. For the 2016 format:
#     - Split on '_'
#     - Keep only the first term
#     - Replace '-' with '.'
#     - use as rownames for 2016 format
#     - Make a data.frame including year, sampling, crop, and plot info
#   4. Rbind the 2016 and 2017 rowname data.frames
#   5. Merge into envi as a new column
#   6. Check results.

# split
is_2017 = grepl('corn', rownames(raw_count)) | grepl('soy', rownames(raw_count))
ct17 = raw_count[is_2017, ]
ct16 = raw_count[!is_2017, ]

# 2017 data
rn = rownames(ct17) %>% 
  strsplit('_') %>%
  lapply(extract, 1:3) %>% 
  sapply(paste, collapse='_')

prefix = data.frame(
  primer = 'F',
  crop = grepl('corn', rn) %>% 
    sapply(ifelse, 'C', 'S'),
  time = grepl('t1', rn) %>% 
    sapply(ifelse, '1', '2')
) %>% 
  apply(1, paste, collapse='')

# sample info
rn17 = strsplit(rn, '_') %>% 
  do.call(rbind, .) %>% 
  as.data.frame %>% 
  set_colnames(c('CROP', 'SAMPLING', 'PLOT'))
sampling_replace = c(t1 = 'seedling',
                     t2 = 'flowering')
rn17$SAMPLING %<>% sampling_replace[.]
rn17$YEAR = 2017  
rn17$NEW_ID = paste(prefix, rn, sep='.')
rn17 %<>% sapply(as.character) %>% 
  as.data.frame(stringsAsFactors=FALSE)

# Check
sum(!(rn17$NEW_ID %in% envi$FUNGI_ID))  # looks good

# Set as rownames
rownames(ct17) = rn17$NEW_ID

# 2016 data
rn = rownames(ct16) %>% 
  strsplit('_') %>% 
  sapply(extract, 1) %>% 
  gsub('-', '.', .)

# Sample info
CROP = grepl('FS', rn) %>% 
  sapply(ifelse, 'soy', 'corn')
SAMPLING = c('seedling', 'flowering') %>% 
  .[substr(rn, 3, 3) %>% as.numeric]
PLOT = strsplit(rn, '\\.') %>% 
  sapply(extract, 2) %>% 
  as.character
rn16 = data.frame(CROP,
                  SAMPLING,
                  PLOT,
                  YEAR = 2016,
                  NEW_ID = rn) %>% 
  sapply(as.character) %>% 
  as.data.frame(stringsAsFactors=FALSE)

# check
sum(!(rn16$NEW_ID %in% envi$FUNGI_ID))

rownames(ct16) = rn16$NEW_ID

# combine new count matrices
out_count = rbind(ct16, ct17)
colnames(out_count) %<>% gsub('ASV', 'ASV_', .)

# update FUNGI_ID
new_rn = rbind(rn16, rn17) %>% 
  as.data.frame(stringsAsFactors=FALSE)
out_envi = merge(envi, new_rn, by=c('CROP', 'SAMPLING', 'YEAR', 'PLOT'), all.x=TRUE)
is_mismatched = out_envi$FUNGI_ID != out_envi$NEW_ID
out_envi[is_mismatched, ]

out_envi$FUNGI_ID = out_envi$NEW_ID
out_envi$NEW_ID = NULL

#### reformat taxonomy
out_tax = sapply(raw_tax, 
             function(x) {
               strsplit(x, '__') %>% 
                 sapply(extract, 2)
             }) %>% 
  set_rownames(rownames(raw_tax))
rownames(out_tax) %<>% gsub('ASV', 'ASV_', .)

#### Update funguild
funguild_script = file.path('FUNGuild', 'Guilds_v1.1.py')
python_call = 'python3'  # replace with the correct path as needed
out_fg_tax = 'FUNGuild Taxonomy.txt'

fg_tax = out_tax[, c('Genus', 'Species')] %>% 
  apply(1, paste, collapse=';') 
fg_tax = data.frame(OTU_ID = names(fg_tax),
                    taxonomy = paste0(';', fg_tax))
here('Data', in_dir, out_fg_tax) %>% 
  write.table(fg_tax, 
              ., 
              row.names=FALSE,
              sep='\t',
              quote=FALSE)

in_quotes = function(x) paste("'", x, "'", sep="")

pypath = here(funguild_script)
path_to_fg_tax = here('Data', in_dir, out_fg_tax)

funguild_call = paste(python_call,
                      in_quotes(pypath),
                      '-otu',
                      in_quotes(path_to_fg_tax),
                      '-db fungi -m -u')

system(funguild_call)


# Wrangle messy funguild output to a .csv
raw_fg = 
  tools::file_path_sans_ext(out_fg_tax) %>% 
  paste0('.guilds') %>% 
  paste(tools::file_ext(out_fg_tax), sep=".") %>% 
  here('Data', in_dir, .) %>% 
  read.table(sep='\t', stringsAsFactors=FALSE, fill=TRUE, header=TRUE) %>% 
  extract(, 1:9)
colnames(raw_fg) %<>% gsub('\\.', ' ', .)


# Export
# Funguild
tools::file_path_sans_ext(out_fg_tax) %>% 
  paste0('.guilds') %>% 
  paste0('.csv') %>% 
  here('Data', .) %>% 
  write.csv(raw_fg, ., row.names=FALSE)

# counts 
here('Data', 'Master ITS Counts.csv') %>% 
  write.csv(out_count, ., row.names=TRUE)

# taxonomy
here('Data', 'Master ITS Taxonomy.csv') %>% 
  write.csv(out_tax, ., row.names=TRUE)

# environmental
here('Data', 'Master Plot Data.csv') %>% 
  write.csv(out_envi, ., row.names=FALSE)
