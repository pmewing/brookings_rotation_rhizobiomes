# read.table to a matrix, with row names and preserving strings
read_matrix = function(path, row.names=1, stringsAsFactors=FALSE, sep=',') {
  dd = read.table(path, 
                  sep=sep, 
                  row.names=row.names, 
                  header=TRUE,
                  stringsAsFactors=stringsAsFactors)
  dd = as.matrix(dd)
  return(dd)
}

#' loads in_count and in_tax from here(path). Cleans up the taxonomy matrix. 
#' Employs zero replacement. 
#' Builds a phyloseq object using rownames from id in plot_data. 
#' Adds an attribute called row_weights, based on log10 counts.
load_to_phyloseq = function(id, in_count, in_tax, plot_data, min_tax=3, min_occurrances=5, min_reads=100, path='Data') {
  require(here)
  
  counts = here('Data', in_count) %>% 
    read_matrix %>%    # read.table(row.names=1, header=TRUE, sep=',', stringsAsFactors=FALSE) %>% as.matrix
    clean_matrix(min_tax=min_tax, min_occurrances=min_occurrances, min_reads=min_reads) %>% 
    cmultRepl(output='p-counts')
  weights = rowSums(counts) %>% 
    log10
  weights = weights/sum(weights)
  counts = sweep(counts, 1, rowSums(counts), '/')
  
  taxonomy = here('Data', in_tax) %>% 
    read_matrix %>% 
    .[colnames(counts), ]
  
  # Format environmental data
  samdat = 
    !is.na(
      as.character(plot_data[, id])
    ) 
  samdat = plot_data[samdat, ]
  rownames(samdat) = samdat[, id]
  
  # make phyloseq object
  out = phyloseq(otu_table(counts, taxa_are_rows=FALSE),
                 sample_data(samdat),
                 tax_table(taxonomy))
  
  weights = weights[rownames(otu_table(out))]
  attr(out, 'row_weights') = weights
  return(out)
}

#' Make a vector of contrast scores from the column of ps_obj (phyloseq object). 
#' Contrast is an 1\*p or p\*1 matrix. 
#' Dimensions must be named with a contrast name and levels(sample_data(ps_obj)[, column])
make_contrasts = function(ps_obj, column, contrast) {
  require(phyloseq)
  require(magrittr)
  
  stopifnot(class(ps_obj)=='phyloseq')
  
  stopifnot(class(contrast) == 'matrix')
  if (any(is.null(dimnames(contrast)))) stop('contrast matrix must have dimnames')
  stopifnot(any(dim(contrast)) == 1)
  if (dim(contrast)[2] == 1) contrast %<>% t
  
  sd = sample_data(ps_obj)
  contr = data.frame(sd) %>% 
    extract2(column) %>% 
    as.character
  contr = contrast[, contr] %>% 
    data.frame %>% 
    set_rownames(rownames(sd)) %>% 
    set_colnames(rownames(contrast))
  return(contr)
}

#' Subsets phylo object by crop and sampling. 
#' Builds Y = otu_table(ps_obj), X = sample_data(ps_obj)[constraint], and Z = sample_data(ps_obj)[partial] matrices
#' optionally uses a contrast for X as defined in contrast
#' Performs partial constrained weighted log-ratio analysis with pcwOrd::pcwOrd:
#' - row weights are log10 total reads
#' - column weights are mean relative abundances
#' - otu_table is center-log-ratio-transformed by easyCODA::CLR
#' Produces diagnostic plots and returns the pcwOrd object.
ordinate = function(ps_obj, crop=NULL, sampling=NULL, constraint=NULL, partial=NULL, plot_name=NULL, contrast=NULL, draw_plot=TRUE) {

  ps = ps_obj
  if (!is.null(crop)) {
    ps = switch(crop, 
                'corn' = subset_samples(ps, CROP=='corn'),
                'soy' = subset_samples(ps, CROP=='soy'))
  }
  if (!is.null(sampling)) {
    ps = switch(sampling,
                'flowering' = subset_samples(ps, SAMPLING=='flowering'),
                'seedling' = subset_samples(ps, SAMPLING=='seedling'))
    # ps = subset_samples(ps, SAMPLING == get('sampling'))
  }
  
  Y = otu_table(ps) %>% 
    CLR
  if (!is.null(constraint)) {
    X = sample_data(ps)[, constraint] %>% 
      data.frame
  } else {
    X = NULL
  }
  if (!is.null(partial)) {
    Z = sample_data(ps)[, partial] %>% 
      data.frame
  } else {
    Z = NULL
  }
  rwt = attr(ps_obj, 'row_weights')[rownames(otu_table(ps))]
  rwt = rwt/sum(rwt)
  add_grouping = NULL  # for plotting
  
  # modify components for plotting
  if (!is.null(contrast)) {
    add_grouping = cbind(X, Z)
    
    X = make_contrasts(ps, constraint, contrast) %>% 
      data.frame
  }
  
  ord = pcwOrd(Y, X, Z, weight_rows=rwt)
  
  if (draw_plot) {
    par(mfrow=c(1,2))
    plot_ord(ord, main=plot_name, row_group=add_grouping)
    ord_scree(ord)
    par(mfrow=c(1,1))
  }

  return(ord)
}

#' Regress taxa scores as predictors of response. Automatically select the best model via `step()`` and AIC,
#' or using L1 penalty (lasso, via `glmnet()`).
#' collects data from the taxa-appropriate phyloseq object.
#' 
#' @arg taxa: character vector of predictors matching rownames(otu_table(ps_obj))
#' @arg response: character column name of response variable. Weird order due to assumed default for this purpose.
#' @arg ps_obj: phyloseq object in which to find partial and response data.
#' @arg predictor: additional predictor vectors to add to the model
#' @arg partial: character column name of partialing data in sample_data(ps_obj). Probably don't want to do this. 
#' @arg lasso: logical. Use lasso penalized regression instead of step()?
#' @arg weight_rows: Logical. Weight rows? If so, are found in attr(ps_obj, 'row_weights').
#' @arg verbose: logical. Print iterations of step()?
regress_taxa = function(taxa, response, ps_obj,
                        predictor=NULL, partial=NULL, 
                        lasso=FALSE, weight_rows=TRUE, verbose=FALSE, alpha=1) {
  require(phyloseq)
  require(magrittr)
  
  # row_weights
  if (weight_rows) {
    rwt = attr(ps_obj, 'row_weights')
    rwt = rwt/sum(rwt)
  } else {
    rwt = rep(1, nrow(otu_table(ps_obj)))
  }
  
  # make dataframes
  # response
  y = sample_data(ps_obj)[, response] %>%
    scale
  # predictor, scaled by row weights
  x = otu_table(ps_obj)[, taxa] %>% 
    log %>% 
    scale
  has_na = apply(x, 2, function(x) any(is.na(x)))
  x %<>% .[, !has_na] %>% 
    data.frame


  if (!is.null(partial)) {
    y = sample_data(ps_obj)[, partial] %>%
      data.frame %>%
      lm(y~., data=.) %>%
      residuals
    
    #EXPERIMENTAL
    # x %<>% sapply(function(q) {
    #   dd = sample_data(ps_obj)[, partial] %>% 
    #     data.frame
    #   q = lm(q~., data=dd) %>% 
    #     residuals
    #   return(q)
    # })
  }
  
  if (!is.null(predictor)) {
    predictor = model.matrix(~., predictor)[, -1, drop=FALSE]
    x %<>% data.frame(predictor)
  }
  
  if (lasso) {
    require(glmnet)
    alpha = alpha
    mm = as.matrix(x)
    y = as.matrix(y)
    rwt = as.matrix(rwt)
    cv_mod = cv.glmnet(mm, y, weights=rwt, standardize=FALSE, nfolds=nrow(x), alpha=alpha)
    mod = glmnet(mm, y, weights=rwt, standardize=FALSE, alpha=alpha)
    mod = list(model = mod, model_cv = cv_mod)
    mod = cv_mod
  } else {
    # mm = data.frame(y,
    #                 x)
    mod = lm(y ~ ., data=x, weights=rwt) %>%
      step(direction='both', trace=1*verbose) %>%
      extract2('call') %>%
      eval
  }
  
  return(mod)
}

#' Select the features in feature_F of ord_permutations$results that have an adjusted
#' p-value < max_p_adjust. Name these taxa based on tax_table(plylo_obj).
#' Include additional info like guilds and (contribution, sensu Greenacre 2013) scores. 
id_top_tax = function(ord_permutations, metric, phylo_obj, phylo_depth, max_p_adj, max_p_val=0.05, axes=1, guilds=NULL, auto_invert=FALSE) {
  require(magrittr)
  require(phyloseq)
  
  # Calculate scores
  scores = ord_scores(ord_permutations$observed, 'column', 'contribution', axes=axes)
  scorenames = colnames(scores)
  
  # invert scores?
  if (auto_invert) {
    for (i in scorenames) {
      correlation = cor(ord_scores(ord_permutations$observed, 'row', 'contribution', axes=axes)[, i], 
                        ord_permutations$observed$X)
      if (correlation < 0) {
        paste('Inverting', i) %>% 
          message
        scores[, i] = -1 * scores[, i]
      }
    }
  }
  
  # Select results
  top_tax = ord_permutations$results[[metric]] 
  rownames(top_tax) = top_tax$feature
  top_tax = top_tax[rownames(scores), ]
  
  # ID top taxa
  top_tax %<>% 
    data.frame(scores) %>%   # combine top taxa and scores
    subset(round(p_adj, 1) <= max_p_adj) %>%   # begin subsetting
    subset(p_val <= max_p_val)
  if (nrow(top_tax)==0) {
    message('no taxa retained')
    return()
  }
  
  # subset quantitative information

  if (metric == 'feature_F') {
    top_tax %<>% .[, c('feature', scorenames, 'Rsq', 'F_obs', 'p_val', 'p_adj')]
  } else if (metric == 'feature_scores') {
    top_tax %<>% .[, c(scorenames, 'p_val', 'p_adj')]# %>% 
      # set_rownames(top_tax$feature)
  } else {
    # (metric!='feature_scores') {
    stop('Metric must be feature_F or feature_scores')
  }
  
  # Begin qualitative info: taxa names and guilds
  tt = tax_table(phylo_obj) %>% 
    data.frame
  tax_names = id_taxa(rownames(top_tax), #top_tax$feature,    # this is sllloooowwww
                      tt,
                      phylo=phylo_depth)
  
  # Guilds
  if (!is.null(guilds)) {
    guilds %<>% .[rownames(top_tax), 'Guild']
  } else {
    guilds = NA
  }
  
  out = data.frame(tax_names, guilds, top_tax)
  # sort by p_value
  ss = order(out$p_val, decreasing=FALSE)
  out = out[ss, ]
  
  return(out)
}

#' subset a phyloseq object to retain otus. Optionally retains row weights
#' if present, subsets CROP, and subsets SAMPLING
subset_otu = function(otu, phylo_obj, crop=NULL, sampling=NULL) {
  require(phyloseq)
  require(magrittr)
  
  sub_ps = phylo_obj
  
  has_rw = any(names(attributes(sub_ps)) == 'row_weights')
  if (has_rw) wt = attr(sub_ps, 'row_weights')
  
  otu_table(sub_ps) %<>% 
    .[, otu]
  
  if (!is.null(crop)) {
    sub_ps %<>% subset_samples(CROP==crop)
  }
  
  if (!is.null(sampling)) {
    sub_ps %<>% subset_samples(SAMPLING==sampling)
  }
  
  if (has_rw) {
    wt = sample_data(sub_ps) %>% 
      rownames %>% 
      wt[.]
    wt = wt/sum(wt)
    attr(sub_ps, 'row_weights') = wt
  }
  
  return(sub_ps)
}

#' calculate ord scores from an ordination for axes scores. Then use these 
#' these axes scores to predict response, along with predictor. 
#' Optionally partial out partial from response.
#' Finds data for response and partial in phylo_obj 
predict_from_ordination = function(response, ord, phylo_obj, predictor_data=NULL, partial=NULL, axes=1) {
  score = ord_scores(ord, 'row', 'contribution', axes=axes) %>% scale
  
  if (is.null(predictor_data)) {
    X = score
  } else {
    X = data.frame(predictor_data, 
                   score)
  }
  
  if (!is.null(partial)) {
    X = data.frame(sample_data(phylo_obj)[, partial],
                   X)
  }

  df = sample_data(phylo_obj)[, response] %>% 
    # unlist %>% 
    data.frame(X)
  names(df)[1] = 'y'
  
  init = paste('y~', partial) %>% 
    formula %>% 
    lm(data=df)
  scope = names(df)[-1] %>% 
    paste(collapse='+') %>% 
    paste('y', ., sep='~') %>% 
    formula
  
  # mod = step(init, scope=scope, direction='forward')
  
  mod = lm(y ~ ., data=df)
  return(mod)
}

#' select the lasso model with the error-minimizing lambda
#' and return the coefficients of all non-zeros (excluding the intercept)
summarize_lasso = function(mod, lambda=c('lambda.min', 'lambda.1se')) {
  lambda = match.arg(lambda, several.ok=TRUE)
  cnames = paste0('coef_', lambda)
  
  s = sapply(lambda, function(x) mod[[x]])
  
  mod_s = coef(mod, s=s) %>% 
    as.matrix %>% 
    round(6)
  # mod_s = coef(mod$model, s=mod$model_cv$`lambda`)[, 1]
  is_retained = rowSums(mod_s != 0) > 0
  # is_taxa = names(mod_s) != '(Intercept)'
  mod_s %<>%
    .[is_retained, , drop=FALSE] %>%
    # as.matrix %>% 
    set_colnames(cnames)
  return(mod_s)
}

#' describe the taxa of taxa_coefs - name them, give scores and effects, guild information if available.
describe_predictive_taxa = function(taxa_coefs, phylo_obj, ordination, depth, guild=NULL, round=3) {
  pred_tax = rownames(taxa_coefs) %>% 
    id_taxa(phylo_obj, phylo=depth)
  
  scores = ord_scores(ordination, 'column', 'contribution', axes=1) %>% 
    extract(rownames(pred_tax), ) %>% 
    round(round)
  
  taxa_coefs %<>% round(round)
  
  if (!is.null(guild)) {
    guild = guild[rownames(taxa_coefs), 'Guild']
    out = data.frame(pred_tax, 
                     Guild = guild,
                     Score = scores,
                     taxa_coefs)
  } else {
    out = data.frame(pred_tax, 
                     Score = scores,
                     taxa_coefs)
  }
  
  return(out)
}

#' replace 'LV' with lv_name
#' replace 'CONTRAST' with contrast
#' append tax_lv parameterization string
customize_sem = function(model, tax_lv, lv_name, contrast) {
  gsub('LV', lv_name, model) %>% 
    gsub('CONTRAST', contrast, .) %>% 
    paste(tax_lv, sep=' \n ')
}

#' make/update a data.frame containing relevant information from structure.
#' structure in this case is a pcwOrd_permutations object (from permute_ord()).
#' Gives crop, sampling, barcode, model Rsq, model F, d.f., and pseudo-p. 
update_full_ord_results = function(structure=NULL, crop=NULL, sampling=NULL, barcode=NULL, results=NULL) {
  
  if (is.null(structure)) {  # initiate
    out = matrix(nrow=0, ncol=7) %>% 
      set_colnames(c('CROP', 'SAMPLING', 'BARCODE', 'RSQ', 'F_OBS', 'Df', 'P_VAL')) %>% 
      as.data.frame
  } else {  # update
    out = data.frame(CROP = crop,
                     SAMPLING = sampling,
                     BARCODE = barcode,
                     RSQ = results$results$model_F['Rsq'],
                     F_OBS = results$results$model_F['F_obs'],
                     Df = paste(results$Df['model_total'], results$Df['free'], sep=','),
                     P_VAL = results$results$model_F['p_val']
    ) %>% 
      set_rownames(NULL) %T>%
      print
    
    update_row = cbind(
      structure[, 'CROP']==crop,
      structure[, 'SAMPLING']==sampling,
      structure[, 'BARCODE']==barcode) %>% 
      rowSums %>% 
      equals(3)
    if (any(update_row)) {  # update in-place
      structure[update_row, ] = out
      out = structure
    } else {  # append a new row
      out =
        rbind(structure, 
              out) %>% 
        as.data.frame
    }
  }
  return(out)
}

update_contrast_ord_results = function(structure=NULL, crop=NULL, sampling=NULL, barcode=NULL, contrast=NULL, results=NULL) {
  
  if (is.null(structure)) {  # initiate
    out = matrix(nrow=0, ncol=8) %>% 
      set_colnames(c('CROP', 'SAMPLING', 'BARCODE', 'CONTRAST', 'RSQ', 'F_OBS', 'Df', 'P_VAL')) %>% 
      as.data.frame
  } else {  # update
    out = data.frame(CROP = crop,
                     SAMPLING = sampling,
                     BARCODE = barcode,
                     CONTRAST = contrast,
                     RSQ = results$results$model_F['Rsq'],
                     F_OBS = results$results$model_F['F_obs'],
                     Df = paste(results$Df['model_total'], results$Df['free'], sep=','),
                     P_VAL = results$results$model_F['p_val']
    ) %>% 
      set_rownames(NULL) %T>%
      print
    
    update_row = cbind(
      structure[, 'CROP']==crop,
      structure[, 'SAMPLING']==sampling,
      structure[, 'BARCODE']==barcode,
      structure[, 'CONTRAST']==contrast
    ) %>% 
      rowSums %>% 
      equals(4)
    if (any(update_row)) {  # update in-place
      structure[update_row, ] = out
      out = structure
    } else {  # append a new row
      out =
        rbind(structure, 
              out) %>% 
        as.data.frame
    }
  }
  return(out)
}

update_taxa_results = function(structure=NULL, crop=NULL, sampling=NULL, barcode=NULL, contrast=NULL, results=NULL) {
  if (is.null(structure)) {  # initiate
    out = c('CROP', 'SAMPLING', 'BARCODE', 'CONTRAST', 'ASV', 'PHYLUM', 'TAXA', 'GUILD', 'AXIS', 'RSQ', 'F_OBS', 
            'P_VAL', 'P_ADJ')
    out = matrix(nrow=0, ncol=length(out)) %>% 
      set_colnames(out) %>% 
      as.data.frame
    
  } else {  # update
    out = data.frame(CROP = crop,
                     SAMPLING = sampling,
                     BARCODE = barcode,
                     CONTRAST = contrast,
                     ASV = rownames(results),
                     PHYLUM = results[, 'Phylum'],
                     TAXA = results[, 'Taxa'],
                     GUILD = results[, 'guilds'],
                     AXIS = results[, 'pcAxis1'],
                     RSQ = results[, 'Rsq'],
                     F_OBS = results[, 'F_obs'],
                     P_VAL = results[, 'p_val'],
                     P_ADJ = results[, 'p_adj']
    ) %>% 
      set_rownames(NULL) %T>%
      print
    
    update_row = cbind(
      structure[, 'CROP']==crop,
      structure[, 'SAMPLING']==sampling,
      structure[, 'BARCODE']==barcode,
      structure[, 'CONTRAST']==contrast) %>% 
      rowSums %>% 
      equals(4)
    if (any(update_row)) {  # update in-place
      structure %<>% .[!update_row, ]
    }
    
    out =
      rbind(structure, 
            out) %>% 
      as.data.frame
  }
  return(out)
}

#' combine axis (principle) scores plus response, contrast, and partial data into a single data.frame for modeling.
#' contrast_row is the row of contrasts[[crop]]. Must be a matrix.
generate_prediction_df = function(ordinations, crop, sampling, contrast_row, partial, constraint, response, plot_data,
                                  auto_invert = TRUE) {
  key_key = plot_data[, c('SAMPLE_ID', 'FUNGI_ID', 'BACTERIA_ID')]
  
  find_ord = c(crop, sampling, rownames(contrast_row)) %>% 
    sapply(grepl, names(ordinations)) %>% 
    apply(1, all)

  format_scores = function(x) {            # basically add relevant info for merging and plotting
    fun_key = ifelse(crop=='corn', 'FC', 'FS')
    is_ITS = grepl(fun_key, rownames(x)) %>% 
      any
    
    # Descriptive name
    new_name = ifelse(is_ITS,
                      'ITS', 
                      '16S') %>% 
      paste0('AXIS_', .)
    colnames(x) = new_name
    
    # merge with other ID keys
    merge_on = ifelse(is_ITS,
                      'FUNGI_ID', 
                      'BACTERIA_ID')
    x[, merge_on] = rownames(x)
    x %<>% merge(key_key[, c(merge_on, 'SAMPLE_ID')], 
                 by=merge_on)
    return(x)
  }
    
  df = ordinations[find_ord] %>% 
    lapply(ord_scores, 'row', 'standard', axes=1) %>%   # calculate scores
    lapply(format_scores) %>% 
    Reduce(function(a, b) merge(a, b, by='SAMPLE_ID', all=TRUE), x=.) %>% # combine into single dataframe using full join
    merge(plot_data[, c(partial, constraint, response, 'SAMPLE_ID')], by='SAMPLE_ID', all.x=TRUE) # add plot data using left join
  
  contrast = rownames(contrast_row)
  df[, contrast] = df[, constraint] %>%  # add contrast levels
    as.character %>% 
    contrast_row[, .]
  
  if (auto_invert) {
    # Are axis scores positively correlated with the contrast?
    
    is_axis = grepl("AXIS_", names(df)) %>% 
      names(df)[.]
    for (i in is_axis) {
      correlation = cor(df[, i], 
                        df[, contrast],
                        use='complete.obs')
      if (correlation < 0) {
        paste('Inverting', i) %>% 
          message
        df[, i] = -1 * df[, i]
      }
    }
  }
  
  return(df)
}

#' plot response ~ barcode_axis based on data in mod_df created above.
#' add shapes for partial and colors for constraint
#' response is optionally the residuals of response ~ partial.
#' returns a ggplot object for easy post-hoc formatting
plot_effect = function(barcode, response, partial, constraint, contrast, data, partial_response=TRUE, title) {
  
  yax_name = response
  xax = paste0('AXIS_', barcode)
  if (!(xax %in% names(data))) xax %<>% paste(contrast, sep='_')
  
  pltdf = data
  if (partial_response) {
    pltdf[, response] = paste(response, 
                              paste(partial, collapse='+'), 
                              sep='~') %>% 
      formula %>% 
      lm(data=pltdf) %>% 
      residuals
    yax_name = paste('Marginal', response)
  }
  
  plt_out = 
    ggplot(pltdf, aes_string(x=xax,
                             y=response,
                             shape=partial[1], # if partialing out multiple factors, keep the first
                             color=constraint)) +
    scale_color_brewer(palette='Dark2') +
    geom_smooth(aes_string(x=xax,
                           y=response),
                color='darkgray',
                fill='lightgray',
                linetype='dashed',
                formula='y~x',
                method='lm',
                inherit.aes=FALSE) +
    geom_point(size=2,
               alpha=0.7) +
    xlab(paste(contrast, barcode, 'Axis Score')) +
    ylab(yax_name) +
    ggtitle(title) +
    theme_minimal()
  
  return(plt_out)
}
