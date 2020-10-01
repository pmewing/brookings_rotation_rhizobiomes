# Process bacteria ASVs

library(here)
library(magrittr)

in_dir = 'Data'
in_bac_dir = file.path('Raw', 'Bacteria')

# load raw 16S tables
in_tag = c('ASV_taxonomy_', '_16S.txt')
years = c('2016', '2017')
in_bac = sapply(years, function(x) paste(in_tag[1], x, in_tag[2], sep=""))

# Write?
ww = FALSE

## Bacteria Tables
bac = list()
for (i in years) {
  tt = here(in_dir, in_bac_dir, in_bac[i]) %>% 
    read.table(row.names=1, stringsAsFactors=FALSE, sep=" ")
  tt['ASV'] = rownames(tt)
  bac[[i]] = tt
}

# ensure matches between ASVs
matches = rownames(bac[[1]]) %in% rownames(bac[[2]])
matches = rownames(bac[[1]])[matches]

match_tab = lapply(bac, function(x) x[matches, ]) %>% 
  lapply(function(x) x[is.na(x)] = '')
sum(match_tab[[1]] != match_tab[[2]]) # Zero. good. 

# Combine
bac = do.call(rbind, bac)
dd = duplicated(bac$ASV)
bac = bac[!dd, ]

# Rename ASV to include underscore (for compatibility with read tables)
bac['ASV'] = with(bac, paste(
  substr(ASV, 1, 3), 
  '_',
  substr(ASV, 4, 9999),
  sep='')
)
rownames(bac) = bac[,'ASV']
bac['ASV'] = NULL

# Export

if (ww) {
  here(in_dir,
       'Master 16S Taxonomy.csv') %>% 
    write.csv(bac, .)
}
