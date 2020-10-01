# Contents:
#   1. corner(x, n=6, loc='topleft'): Show n*n cells of a matrix or dataframe from stated corner.
#   2. subset_matrix(data, metadata, column, value): select rows of data corresponding to value in metadata.

#### Corner ####
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

subset_matrix = function(data, metadata, column, value, key) {
  # subset data matrix based on column value in metadata
  # ensures metadata rownames match those of data.
  
  keep_rows = which(metadata[, column] == value)
  keep_rows = metadata[keep_rows, key]
  
  avail_rows = rownames(data)
  tt = keep_rows %in% avail_rows
  keep_rows = keep_rows[tt]
  
  out = data[keep_rows, ]
  out = as.matrix(out)
  
  return(out)
}

clean_matrix = function(X, min_occurrances=2, min_taxa=1, min_reads=100) {
  
  low_reads = rowSums(X) < min_reads
  X = X[!low_reads, ]
  
  dim_diff = 1
  while (dim_diff != 0) {
    in_dim = prod(dim(X))
    
    pa = X > 0
    
    X = X[rowSums(pa) >= min_taxa, 
            colSums(pa) >= min_occurrances]
    out_dim = prod(dim(X))
    dim_diff = out_dim - in_dim
  }
  return(X)
}

make_factor = function(X, grouping_vars, key) {
  # return a matrix of grouping_vars 
  # with names given by key.
  
  X = na.omit(X[, c(grouping_vars, key)])
  
  out = as.matrix(X[, grouping_vars])
  rownames(out) = X[, key]
  
  return(out)
}

id_taxa = function(taxa, table, phylo='Genus', depth=2, return_phylum=TRUE) {
  
  if (class(table) == 'phyloseq') {
    table %<>%
      tax_table %>% 
      data.frame
  }
  table = as.matrix(table)
  
  out = lapply(taxa, function(x) {
  tax_id = x
  
  taxa = table[x, ]
  
  phy_depths = colnames(table)
  
  tt = TRUE
  ix_start = which(phy_depths == phylo)
  ix = ix_start + 1
  while (tt) {
    ix = ix - 1
    tt = is.na(taxa[ix])
  }

  if (ix == ix_start) {
    indices = seq(ix-depth+1, ix)
    out = paste(taxa[indices], collapse=" ")
  } else {
    indices = seq(ix-depth+2, ix)
    out = paste(taxa[indices], collapse=" ")
    out = paste('Unknown', out)
  }
  out = matrix(out)
  colnames(out) = 'Taxa'

  if (return_phylum) {
    phy = matrix(taxa['Phylum'])
    colnames(phy) = 'Phylum'
    out = cbind(phy, out)
  } 
  
    rownames(out) = tax_id
    return(out)
  })
  out = do.call(rbind, out)
  
  attr(out, 'Max Depth') = phylo
  
  return(out)
}

make_matrices = function(kingdom='bacteria', 
                         by='CROP',                 # column to select from 
                         choices=c('corn', 'soy'),  # options to match in `by`
                         close_mat=TRUE, 
                         meta=NULL, # starting here, arguments to generically select matrices. 
                         counts_mat=NULL, 
                         key=NULL) {
  
  if (is.null(counts_mat)) {
    mat = bac
    key = 'BACTERIA_ID'
    if (kingdom == 'fungi') {
      key = 'FUNGI_ID'
      mat = fun
    }
  } else {
    mat = counts_mat
    key = key
  }
  
  match_rows = sapply(meta[, key], function(x) x %in% rownames(mat))
  meta = subset(meta, match_rows & meta[, by] %in% choices) %>% 
    lapply(as.factor) %>% 
    as.data.frame
  rownames(meta) = meta[, key]
  mat = mat[rownames(meta), ]
  
  if (close_mat)  mat = mat/rowSums(mat)
  
  return(list(mat=mat,
              env=meta))
}

clr = function(X, base=exp(1)) {
  # Calculate the centered log-ratio of each row in X. The default is a natural log. 
  
  X = as.matrix(X)
  if (any(X==0)) {stop("Your matrix has zeros")}
  if (any(round(rowSums(X), 5) != 1)) warning("Not all rowSums are one. Are these proportions?")
  
  geomean = apply(X, 1, function(y) {
    sum(log(y), base)/length(y)
  })

  out = sweep(log(X, base), 1, geomean, '-')
  
  return(out)
}


get_legend<-function(myggplot){
  # http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
  
  require(ggplot2)
  
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

load_libs = function(libs) {
  for(i in libs) {
    if (!require(i, character.only=TRUE)) {
      options(Ncpus=parallel::detectCores())
      install.packages(i)
      if (!require(i, character.only=TRUE)) {
        if (!require(BiocManager)) install.packages('BiocManager')
        BiocManager::install(i)
        if (!require(i, character.only=TRUE)) {
          options(Ncpus=1)
          stop(paste('Cannot find package:', i))
        }
      }
    }
  }
}

read_matrix = function(path, row.names=1, stringsAsFactors=FALSE, sep=',') {
  dd = read.table(path, 
                  sep=sep, 
                  row.names=row.names, 
                  header=TRUE,
                  stringsAsFactors=stringsAsFactors)
  dd = as.matrix(dd)
  return(dd)
}

close_matrix = function(x, by=c('row', 'column')) {
  switch(match.arg(by),
         'row' = sweep(x, 1, rowSums(x), '/'),
         'column' = sweep(x, 2, colSums(x), '/')
  )
}

grid_arrange_shared_legend = function(plots, title=NULL) {
  require(grid)
  # plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(0.95, "npc") - lheight, lheight), top=textGrob(title, x=unit(0.05, 'npc'), just='left'))
}

# grid_arrange_shared_legend = function(plots, title=NULL, position='bottom') {
#   # plots <- list(...)
#   nplots = length(plots)
#   ncol = floor(sqrt(nplots))
#   nrow = ceiling(sqrt(nplots))
#   if (position %in% c('left', 'right')) {
#     ncol = ncol + 1
#     nrow = nrow
#     lheight_adjust = 1
#     lwidth_adjust = 0.95
#   } else {
#     ncol = ncol
#     nrow = nrow + 1
#     lheight_adjust = 0.95
#     lwidth_adjust = 1
#   }
#   
#   g <- ggplotGrob(plots[[1]] + theme(legend.position=position, 
#                                      legend.box='vertical', 
#                                      legend.spacing=unit(0, 'in'),
#                                      legend.margin=margin()))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   lwidth = sum(legend$width)
#   grid.arrange(
#     do.call(arrangeGrob, lapply(plots, function(x)
#       x + theme(legend.position="none"))),
#     legend,
#     ncol = ncol,
#     nrow = nrow,
#     heights = unit.c(unit(lheight_adjust, "npc") - lheight, rep(lheight, nrow-1)),
#     # widths = unit.c(unit(lwidth_adjust, 'npc') - lwidth, rep(lwidth, ncol-1)),
#     top=textGrob(title, x=unit(0.05, 'npc'), just='left'))
# }