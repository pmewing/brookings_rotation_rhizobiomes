# Functions for wrangling the data

read_tough_files = function(directory, split_char='\t') {
  # For files that read.table has trouble with, but are otherwise fine.
  # Returns a data.frame. Does not convert characters to factors. 
  
  require(magrittr)
  
  raw = readLines(directory) %>% 
    strsplit(split=split_char)
  
  nc = max(sapply(raw, length))
  nr = length(raw)
  out = matrix(ncol=nc, nrow=nr-1) %>% 
    as.data.frame %>% 
    set_names(raw[[1]])
  
  for (i in 2:nr) out[i-1, ] = raw[[i]]
  
  out %<>% type.convert(as.is=TRUE)
  
  return(out)
}

format_names = function(x) {
  # x is a dataframe
  # replace periods and spaces with underscores
  # all caps names
  
  require(magrittr)
  
  colnames(x) %<>% 
    gsub('\\.', '_', .) %>% 
    gsub(' ', '_', .) %>% 
    toupper
  return(x)
}

mismatch_names = function(a, b) {
  # a and b are dataframes or vectors
  # find levels of a not in b, and b not in a
  # returns a well-named list
  
  require(magrittr)
  
  a_lvls = a
  if (is.data.frame(a)) a_lvls %<>% names
  b_lvls = b
  if (is.data.frame(b)) b_lvls %<>% names
  
  n = sapply(a_lvls, function(x) !(x %in% b_lvls)) %>% 
    a_lvls[.]
  
  p = sapply(b_lvls, function(x) !(x %in% a_lvls)) %>% 
    b_lvls[.]
  
  out = list(n, p)
  names(out) = c(deparse(substitute(a)),
                 deparse(substitute(b)))
  return(out)
}

make_sampleID = function(x, cleanup=TRUE) {
  # make a column called SAMPLE_ID from CROP, YEAR, PLOT, and a SAMPLING time flag
  # x is a data.frame
  # If cleanup = TRUE, then deletes a column called 'SAMPLEID'
  
  tt = ifelse(x[, 'SAMPLING']=='seedling',
              't1', 
              't2')
  x[, 'SAMPLE_ID'] = with(x, 
                          paste(CROP, YEAR, PLOT, tt, 
                                sep="_"))
  
  t1 = length(unique(x[, 'SAMPLE_ID']))
  t2 = nrow(x)
  
  if (t1 != t2) warning('SAMPLE_ID is not unique')
  
  if (cleanup) x[, "SAMPLEID"] = NULL
  
  return(x)
}

data_and_key = function(x, key='SAMPLE_ID', also_keep=NULL) {
  # retain numeric columns plus a sample key from data.frame x
  # will also keep a list of column names in also_keep
  x = as.data.frame(x)
  
  tt = sapply(x, is.numeric)
  tt = colnames(x)[tt]
  tt = c('SAMPLE_ID', also_keep, tt)
  
  return(x[, tt])
}

corner = function(x, position='topleft', nr=6, nc=6) {
  # Shows nr rows and nc columns of data.frame/matrix x. Position defines which corner.
  
  stopifnot(position %in% c('topleft', 'topright', 'bottomleft', 'bottomright'))
  
  rows = c(1:nr)
  if (!grepl('top', position)) {
    rows = sort(nrow(x)-rows+1)
  }
  
  cols = c(1:nc)
  if (!grepl('left', position)) {
    cols = sort(ncol(x)-cols+1)
  }
  
  return(x[rows, cols])
}

make_taxa_table = function(x, row_names='SAMPLEID', asv_tag='ASV') {
  out = x
  
  if (!is.null(row_names)) rownames(out) = x[, row_names]
  
  cols = grepl('ASV', colnames(x))
  
  out = out[, cols]
  out = as.matrix(out)
  return(out)
}


merge_taxa_tables = function(..., row_as_sample=TRUE, na_val=0) {
  # combine OTU or ASV tables into one, such that mismatching columns are propogated.
  
  tables = list(...)
  
  if (row_as_sample) tables = lapply(tables, t) # columns are sample units
  
  tables = lapply(tables, function(x) { # transpose and make a column with ASV identifiers
    x = as.data.frame(x)
    x[, 'COLNAME'] = rownames(x)
    return(x)
  })
  
  # merge
  out = tables[[1]]
  for (i in 2:length(tables)) {
    out = merge(out, tables[[i]], by='COLNAME', all=TRUE)
  }
  
  # clean up and reconvert to standard table.
  rownames(out) = out[, 'COLNAME']
  out[, 'COLNAME'] = NULL
  out = as.matrix(out)
  
  out[is.na(out)] = na_val
  
  if (row_as_sample) out = t(out)
  
  return(out)
}