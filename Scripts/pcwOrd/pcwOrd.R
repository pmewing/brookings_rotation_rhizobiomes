#### Overview ####
# Partial (weighted) RDA functions
# 
# Motivation: easyCODA performs weighted RDA, but cannot partial out additional effects. 
# Pre-partialing them doesn't play well with either easyCODA or vegan. vegan doesn't 
# perform weighted RDA. Also, easyCODA is s-l-o-w. vegan's code is much faster.
#
# vegan pRDA algorithm: all calls can be inspected with getAnywhere('function_name')
# Call: rda(X, Y, Z) where:
#   X = the community matrix (columns are features)
#   Y = the predictors, ex envi or treatment
#   Z = conditioning/partialing variables, ex. year or other nusiance effects
#
# Steps:
# 1. Convert Y and Z to model matrices (without intercepts)
# 2. Call ordConstrained(X, Y, Z).
#     NOTE: starting with ordConstrained, names for matrices Y and X are switched in the source code (but not in this description)
#   1. Calls initPCA(X), which:
#     - centers X
#     - Divides by root(nrow-1)
#   2. Calls ordPartial(X, Z), which:
#     - Centers Z
#     - (optional) weights rows in Z
#     - Performs linear regression via qr decomposition of Z
#     - Calculates the fitted values of X to qr(Z)
#       - Partialed out variance is the sum-squared of fitted values
#     - Calculates the residuals of X from qr(Z)
#     - Returns X as the residuals, pX
#     - Returns results as pCCA
#   3. Calls ordConstrain(pX, Y, Z), which:
#     - Combines Z and Y into a single matrix Y
#         - Not sure why, as Z is explicitly excluded downstream
#     - Centers columns by row-weighted means, then reweights by root(rowweights)
#     - Performs linear regression via qr decomposition of Y
#     - Calculates fitted values of pX to qr(Y)
#     - Performs an SVD of these fitted values
#     - Calculates the residuals of of pX from qr(Y)
#     - Returns pX as the residuals, rX
#     - Returns results as CCA
#   4. Calls ordResid(rX), which:
#     - Performs an SVD of rX
#     - Returns the results as CA
#
# Scores:
#   - wa (sites) are principle. They can be calculated by multiplying pX %*% v %*% diag(1/d) (d and v from svd(fitted))
#   - bp/cn (constraint/centeroid) variable scores are standard, multiply X %*% u (from svd(fitted))
#
# So overall:
#   - X is weighted by rows and columns at the beginning
#   - Y and Z are weighted by rows and columns as called
#   - The fitted values give centeroid scores (CCA)
#   - The residuals are regressed or decomposed in the next step.
#
# A number of steps are repeated:
#   1. Center and weight a matrix
#   2. Perform linear regression
#   3. Perform SVDs
# These can all be packaged separately for efficiency and clarity.

#### pcwOrd and underlying functions ####
.regress = function(Y, X) {
  # via QR decomposition
  # Y is response, X is predictor
  
  Q = qr(X)
  out = list(
    fitted = qr.fitted(Q, Y),
    resids = qr.resid(Q, Y)
  )
  return(out)
}

.make_dummy = function(X=NULL, nr=NULL) {
  # turn X into a dummy variable
  if (is.null(X)) X = matrix(rep(1, nr))
  
  X = as.data.frame(X)
  X = model.matrix(~., data=X)[, -1, drop=FALSE]
  
  return(X)
}

.prepare_matrix = function(X, rw=NA, cw=NA, as_CA=FALSE) {
  # center and weight a matrix
  
  nr = nrow(X)
  nc = ncol(X)
  
  if (is.na(rw[1])) rw = rep(1, nr)
  if (is.na(cw[1])) cw = rep(1, nc)
  
  if (as_CA) {
    X = X/sum(X)
  } else {
    colmeans = apply(X, 2, weighted.mean, w=rw)
    X = scale(X, center=colmeans, scale=FALSE) # weighted column means?
  }
  
  X = sweep(X, 1, sqrt(rw), '*') # scale rows
  X = sweep(X, 2, sqrt(cw), '*') # scale columns
  
  return(X)
}


.make_weights = function(P, choice, weight=FALSE, as_vegan=FALSE) {
  # Make vectors of matrix weights for either rows or columns: One vector for the response data
  # (wt_init), and another for the predictors. 
  # 
  # as_vegan refers to vegan-style row- and column-weights applied to the response matrix
  # for unweighted ordination, which are 1/(nrow-1) for rows and 1 for columns
  # Otherwise, both rows and columns have easyCODA-style weights of 1/nrow or 1/ncol, respectively.
  # 
  # if weight==TRUE, then weights by marginals (1/rowSums or 1/colSums)
  # 
  # if weight is specified, then those are used
  
  choice = substr(tolower(choice), 1, 1)
  
  n = switch(choice,
             'r'=nrow(P),
             'c'=ncol(P))
  
  if (!weight[1]) {
    wt = rep(1/n, n)
    
    if (as_vegan) {          # vegan style
      wt = switch(choice,
                  'r' = rep(1/(n-1), n), 
                  'c' = rep(1, n))
      
    }
  } else if (is.logical(weight)) {
    wt = switch(choice,
                'r'=1/rowSums(P),
                'c'=1/colSums(P))
  } else {
    wt = weight
  }
  return(wt)
}

.remove_zeros = function(S) {
  # remove columns/eigenvectors from SVD that are (basically) zero
  
  not_zero = S$d > sqrt(.Machine$double.eps)  # from vegan
  
  S$d = S$d[not_zero]
  S$u = S$u[, not_zero, drop=FALSE]
  S$v = S$v[, not_zero, drop=FALSE]
  
  return(S)
}

pcwOrd = function(Y, 
                  X=NULL, 
                  Z=NULL, 
                  weight_rows=FALSE, 
                  weight_columns=FALSE, 
                  as_CA=FALSE, 
                  as_vegan=FALSE) {
  # Conduct a (partial) (constrained) (weighted) ordination. Currently only runs principle component analysis
  # and derivatives. May include correspondence analaysis in the near future.
  #
  # Y: Community matrix. Columns are features, rows are samples. Alternatively, an object from easyCODA::CLR.
  # X: Constraining data.frame. Treatments or environmental variables.
  # Z: Partialing data.frame. Treatments or environmental variables.
  # weight_rows: logical or numeric. If TRUE, will weight by rowSums(Y).
  # weight_columns: logical or nueric. If TRUE, will weight by colSums(Y).
  # as_CA: NOT IMPLEMENTED. run as correspondence analysis? If TRUE, weights of rows and columns are overrridden as TRUE
  #   and Y is double-standardized. If FALSE, Y is centered and weighted as indicated.
  # as_vegan: weight the data matrix Y as in vegan?
  
  method = 'RDA'
  if (as_CA) {
    stop('as_CA is not yet implemented. You might try easyCODA::CCA or vegan::cca.')
    weight_rows = TRUE
    weight_columns = TRUE
    method = 'CCA'
  }
  
  # for easyCODA compatibility
  if (hasName(Y, 'LR.wt') & hasName(Y, 'LR')) {
    weight_columns = Y$LR.wt
    Y = Y$LR
  }
  
  nr = nrow(Y)
  nc = ncol(Y)
  
  # Begin preparing Y
  Y_orig = Y
  Y = as.matrix(Y)
  
  # Define weights
  rw = .make_weights(Y, 'row', weight_rows, as_vegan)
  cw = .make_weights(Y, 'col', weight_columns, as_vegan)
  
  # scale and weight Y
  Y = .prepare_matrix(Y, rw=rw, cw=cw, as_CA)
  
  # Partial Z from Y
  if (is.null(Z)) {
    pY = Y
    partial = NULL
    X_dummy = X
  } else {
    Z_dummy = .make_dummy(Z, nr)
    Z_act = .prepare_matrix(Z_dummy, rw)
    out = .regress(Y, Z_act)
    
    pY = out$resids
    partial = out$fitted
    if (!is.null(X)) {
      X_dummy = cbind(Z, X)
    }
  }
  
  # Constrain pY by X
  if (is.null(X)) {
    rY = pY
    fitted = NULL
    constrained = NULL
    X_weights = NULL
  } else {
    X_dummy = .make_dummy(X_dummy, nr)
    # X_dummy = .make_dummy(X, nr)
    X_act = .prepare_matrix(X_dummy, rw)
    out = .regress(pY, X_act)
    
    rY = out$resids
    fitted = out$fitted
    
    # ordinate
    constrained = .remove_zeros(svd(fitted))
    
    # Calculate group weights
    zcol = ifelse(is.null(Z), 0, ncol(Z_dummy))
    X_weights = 1/sqrt(colSums(X_act^2))[-c(1:zcol)]
  }
  
  # Ordinate residuals
  unconstrained = .remove_zeros(svd(rY))
  
  # To return:
  out = list(
    method = method,
    Y = Y_orig,
    Y_scaled = Y,
    Y_partial = pY,
    Y_residual = rY,
    X = X,
    X_weights = X_weights,
    Z = Z,
    row_weights = rw,
    row_names = rownames(Y_orig),
    col_weights = cw,
    col_names = colnames(Y_orig),
    group_names = colnames(X),
    partial = partial,
    fitted = fitted,
    constrained = constrained,
    unconstrained = unconstrained
  )
  attr(out, 'class') = 'pcwOrd'
  attr(out, 'method') = ifelse(as_CA, 'CA', 'PCA')
  
  return(out)
}
#### Utility Functions ####
ord_scores = function(R, 
                      choice=c('row', 'column', 'centeroid', 'biplot'), 
                      scaling=c('standard', 'principle', 'contribution'), 
                      constrained=TRUE, axes=c(1,2), add_X=FALSE, add_Z=FALSE, 
                      add_grouping=FALSE, 
                      add_label=FALSE,
                      strip_label=NULL,...) {
  # Calculate scores for rows, columns
  # rownames of add_grouping should match the rownames of choice. Can be a vector.
  # ... is generally ignored
  
  stopifnot(class(R) == 'pcwOrd')
  
  choice = match.arg(choice)
  scaling = match.arg(scaling)
  
  is_unconstrained = is.null(R$constrained)
  
  type = ifelse(constrained & !is_unconstrained,
                'constrained',
                'unconstrained')
  
  if (is_unconstrained & (choice %in% c('centeroid', 'biplot'))) {
    msg = paste('cannot plot', choice, 'as this ordination lacks constraints')
    stop(msg)
  }
  
  # need:
  # weights: by choice
  # singular values: by type
  # left- or right-sides: by choice
  # (maybe) Y - for inferring constrained row scores
  # If constrained, pull all unconstrained axes, as well.
  # If these are centeroids, need to calculate the centeroid on the unconstrained axis.
  
  
  
  if (choice == 'centeroid') {
    X_grp = as.data.frame(R$X)
    if (any(sapply(X_grp, is.numeric))) {
      stop('not sure how to calculate centeroids when mixed with continuous variables')
    }
    
    X_grp = apply(X_grp, 1, paste, collapse='_')
  }
  
  if (choice == 'biplot') {
    X_bp = .make_dummy(R$X)
    bp_names = colnames(X_bp)
    # bp_names = sapply(bp_names, function(x) {
    #   for (i in colnames(R$X)) {
    #     x = gsub(i, '', x)
    #   }
    #   return(x)
    # })
    X_bp = .prepare_matrix(X_bp, 
                           rw=R$row_weights)
    X_bp = as.matrix(X_bp)
  }
  
  wt = switch(choice,
              'row' = R$row_weights,
              'column' = R$col_weights,
              'centeroid' = tapply(R$row_weights, X_grp, mean),  # aggregate row_weights by group
              'biplot' = rep(1, ncol(X_bp))
  )
  
  
  
  V = R[[type]]$v
  U = R[[type]]$u
  d = R[[type]]$d
  if (length(R$X_weights)==0) R$X_weights = 1
  
  scores = switch(choice,
                  'row' = U,
                  'column' = V,
                  'centeroid' = apply(as.matrix(U), 2, tapply, X_grp, mean), # aggregate scores by group
                  'biplot' = sweep(as.matrix(crossprod(X_bp, U)), 1, R$X_weights, '*')
  )
  
  # Calculate unconstrained if necessary
  if (type == 'constrained') {
    uV = R$unconstrained$v
    uU = R$unconstrained$u
    ud = R$unconstrained$d
    
    ss = switch(choice, 
                'row' = uU,
                'column' = uV,
                'centeroid' = matrix(rep(0, length(ud)*nrow(scores)),
                                     nrow=nrow(scores)), 
                'biplot' = matrix(rep(0, length(ud)*nrow(scores)), 
                                  nrow=nrow(scores)))
    d = c(d, ud)
    U = cbind(U, uU)
    V = cbind(V, uV)
    scores = cbind(scores, ss)
  }
      
  if (choice=='row' & type=='constrained') {
    # infer point scores based on the column side of svd
    scores = R$Y_partial %*% V
    scores = sweep(scores, 2, d, '/')
  }
  
  out = switch(scaling,
               'contribution' = scores,
               'standard' = sweep(scores, 1, sqrt(wt), '/'),
               'principle' = 
                 sweep(
                   sweep(scores, 1, sqrt(wt), '/'), # scale standard score by SVs
                   2, d, '*')
  )
  
  if (is.null(axes) | is.character(axes)) {
    axes = 1:ncol(out)
  }
  
  # Column/axis names
  n_axes = ncol(out)
  axis_names = rep('Axis', n_axes)
  
  if (type == 'unconstrained') {
    if (!is_unconstrained) {
      axis_names = paste0(rep('u', n_axes), 
                          axis_names,
                          1:n_axes)
    } else {
      axis_names = paste0(axis_names, 1:n_axes)
    }
  }
  if (type == 'constrained') {
    axis_names = 
      paste0(
        c(rep('c', length(R$constrained$d)),
          rep('u', length(R$unconstrained$d))), 
        axis_names,
        c(1:length(R$constrained$d),
          1:length(R$unconstrained$d)))
  }
  if (!is.null(R$partial)) {
    axis_names = paste0('p',
                        axis_names)
  }
  colnames(out) = axis_names
  
  # Row names
  rownames(out) = switch(choice,
                         'row' = R$row_names,
                         'column' = R$col_names,
                         'centeroid' = rownames(out),
                         'biplot' = bp_names
  )
  
  # Select Axes
  out = out[, axes, drop=FALSE]
  out = as.data.frame(out)
  
  # Additional data
  if (add_X & !is.null(R$X)) {
    if (choice == 'column') stop("Cannot add X to column scores")
    if (choice == 'row') {
      grp = as.matrix(R$X)
      rownames(grp) = rownames(out)
      out = data.frame(out, grp)
    } else {
      out = data.frame(
        out,
        GROUP = rownames(out)
      )
    }
  }
  
  if (add_Z & !is.null(R$Z)) {
    if (choice != 'row') stop("Can only add Z to row scores")
    grp = as.matrix(R$Z)
    rownames(grp) = rownames(out)
    out = data.frame(out, grp)
  }
  
  add_grouping = as.data.frame(add_grouping)
  if (add_grouping[1,1] != FALSE | prod(dim(grouping)) > 1) {
    grp = add_grouping[rownames(out), , drop=FALSE]
    out = data.frame(out, grp)
  }
  
  if (add_label==TRUE) add_label = 15
  if (is.numeric(add_label)) {
    constr = type=='constrained'
    top_labels = top_scores(R, n=add_label, choice=choice,
                             scaling=scaling,
                             axes=axes, constrained=constr)
    
    ll = sapply(rownames(out), function(x) {
      ifelse(x %in% rownames(top_labels), x, '.')
    })
    if (!is.null(strip_label)) {
      ll = gsub(strip_label, '', ll)
    }
    ll_levels = ll[!duplicated(ll)] # preserve order
    out$LABEL = factor(ll, levels=ll_levels)
  }
  
  attr(out, 'class') = c('data.frame', 'pcwOrd_scores')
  attr(out, 'scaling') = scaling
  attr(out, 'type') = type
  attr(out, 'nAxes') = length(axes)
  return(out)
}

.axis_variances = function(R) {
  stopifnot(class(R) == 'pcwOrd')
  
  total = c(Total = sum(R[['Y_scaled']]^2))
  
  uc = R[['unconstrained']][['d']]^2
  names(uc) = paste0('uAxis', 1:length(uc))
  
  cc = R[['constrained']]
  if (!is.null(cc)) {
    cc = cc[['d']]^2
    names(cc) = paste0('cAxis', 1:length(cc))
  }
  
  pp = R[['partial']]
  if (!is.null(pp)) {
    pp = c(Partial = sum(pp^2))
  }
  
  out = list(total = total, 
             partialed = pp,
             constrained = cc,
             unconstrained = uc
  )
  
  out = lapply(out, function(vv) {
    rbind(
      value = vv,
      pct_group = 100*vv/sum(vv),
      pct_total = 100*vv/total
    )
  })
  
  return(out)
}

ord_variance = function(R) {
  # Return variances of ordination axes
  
  axis_var = .axis_variances(R)
  
  aggregate_var = sapply(axis_var, rowSums)[-2, ]
  
  out = list(type = switch(attr(R, 'method'),
                           'CA' = 'inertia',
                           'PCA' = 'variance'),
             summary = aggregate_var)
  out = unlist(
    list(out,
         axis_var),
    recursive=FALSE)
  
  return(out)
}

ord_scree = function(R, n=10, as_percent=TRUE, main=NULL, as_grob=FALSE) {
  # Plot variances of ordination axes.
  # n: Number of axes to plot (excluding partial and total, if present and requested)
  # as_percent: logical or 'group'. 
  #   - If FALSE, gives raw variance values and also the total variance.
  #   - If TRUE, gives percent of total variance.
  #   - If 'group', gives percent of variance within the group, and excludes partial and unconstrained axes.
  # as_grob: as a grid object (via ggplot2). Otherwise uses base graphics.
  
  vars = ord_variance(R)
  ylab = vars$type
  ylab = paste0(toupper(substr(ylab, 1, 1)),
                substr(ylab, 2, 9999))
  
  if (as_percent=='group') {
    ylab = paste('Percent of Within-Group', ylab)
  } else if (as_percent) {
    ylab = paste('Percent of Total', ylab)
  }
  
  p = ifelse(is.null(R$partial), 0, 1)
  
  if (as_percent=='group') {
    vars = vars$constrained
    n = min(n, ncol(vars))
    vars = vars['pct_group', 1:n]
    
  } else if (as_percent) {
    vars = with(vars, 
                cbind(partialed, 
                      constrained, 
                      unconstrained)
    )
    
    n = min(n + p, ncol(vars))
    rest = min(n+1, ncol(vars))
    vars = cbind(vars[, 1:n],
                 Remaining = rowSums(vars[, rest:ncol(vars)])
    )
    
    vars = vars['pct_total', ]
    
  } else {
    vars = with(vars, 
                cbind(total,
                      partialed, 
                      constrained, 
                      unconstrained)
    )
    
    n = min(n + p + 1, ncol(vars))
    vars = vars['value', 1:n]
    
  }
  
  colors = c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E')
  color_select = as.factor(substr(names(vars), 1, 2))
  
  if (as_grob) {
    require(ggplot2)
    
    pltdf = data.frame(
      Y = vars,
      X = factor(names(vars), levels=names(vars)),
      grp = color_select
    )
    
    ggplot(pltdf, 
           aes(x=X,
               y=Y,
               fill=grp)) +
      scale_fill_brewer(palette='Dark2',
                        guide=FALSE) +
      geom_col() +
      ylab(ylab) +
      xlab(NULL) +
      ggtitle(main) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=90))
    
  } else {
    barplot(vars, ylab=ylab, horiz=FALSE, col=colors[color_select], border=NA, main=main)
  }
}

permute_ord = function(R,
                       permute_on=c('partial',
                                    'initial',
                                    'fitted'),
                       times=999,
                       return_permutations=FALSE) {
  # Permutational ANOVA of the pcwOrd results.
  # Arguments
  #   R: pcwOrd object
  #   permute_on: level at which to permute.
  #   times: N permutations
  #   return_permutations: whether to return the permuted F values.
  #
  # Notes:
  #   - initial = initial data. Equal to vegan::permutest(model='direct').
  #   - partial = residuals of partialing (default). Equal to vegan::permutest(model='reduced')
  #   - fitted = residuals of constraints. Equal to vegan::permutest(model='full')
  #
  # Permutes the given matrix, runs through the partialing and conditioning steps, and calculates F.
  # Compares observed F to the permuted F. Returns p_val = p(Fobs < Fpermute).
  
  require(permute)
  
  stopifnot('pcwOrd' %in% class(R))
  
  permute_on = match.arg(permute_on)
  
  if (is.null(R$constrained)) stop("Can't test a model without a constraining variable!")
  
  iY = R$Y_scaled    # initial (scaled/weighted) data
  pY = R$Y_partial   # data after partialing Z
  fY = R$fitted      # data fitted to x
  rY = R$Y_residual  # residuals after fitting X
  X = R$X
  Z = R$Z
  rw_orig = R$row_weights
  nr = length(rw_orig)
  
  zcol = ifelse(is.null(Z),
                0,
                ncol(.make_dummy(Z)))
  xcol = ncol(.make_dummy(X))
  
  num_df = xcol
  denom_df = nr - num_df - zcol - 1
  residuals = sum(rY^2)
  fitted = sum(fY^2)
  partial_resids = sum(pY^2)
  total_var = sum(iY^2)
  
  F_stat = (fitted/num_df) / (residuals/denom_df)
  
  # Permute
  
  P = switch(match.arg(permute_on),
             'partial' = pY,
             'initial' = iY,
             'fitted'  = rY
  )
  
  permutations = shuffleSet(nr, times)
  
  permutations = apply(permutations, 1, function(i) {
    
    P = P[i, ]
    rw = rw_orig[i]
    
    # Partial Z from P
    if (is.null(Z)) {
      pP = P
      X_dummy = X
    } else {
      Z_dummy = .make_dummy(Z)
      Z_act = .prepare_matrix(Z_dummy, rw)
      out = .regress(P, Z_act)
      
      pP = out$resids
      X_dummy = cbind(Z, X)
    }
    
    # Constrain pY by X
    X_dummy = .make_dummy(X_dummy)
    X_act = .prepare_matrix(X_dummy, rw)
    out = .regress(pP, X_act)
    
    return(out)
  })
  
  # Taxa Scores
  
  
  fit_perm = sapply(permutations, function(x) {
    sum(x$fitted^2)
  })
  resid_perm = sapply(permutations, function(x) {
    sum(x$resids^2)
  })
  F_perm = (fit_perm / num_df) / (resid_perm / denom_df)
  
  alpha = sum(F_stat >= F_perm)/(times+1)
  
  if (!return_permutations) {
    F_perm = times
  }
  
  out = c(
    total_variance = total_var,
    variance_after_partialing = partial_resids,
    fitted = fitted,
    residuals = residuals,
    num_df = num_df,
    denom_df = denom_df,
    F_stat = F_stat,
    p_val = 1-alpha,
    n = times
  )
  attr(out, 'permutation_model') = permute_on
  
  if (return_permutations) {
    out = list(
      results = out,
      permutations = permutations,
      observed = R
    )
    class(out) = 'pcwOrd_permutations'
  }
  return(out)
}

permute_features = function(P,
                            ...,
                            axes=1,
                            adjustment='BH',
                            tail=c('two', 'upper', 'lower'),
                            calc_F=FALSE,
                            calc_axes=TRUE,
                            return_permutations=FALSE) {
  # if calc_F = TRUE, will return pseudo-F-tests for each feature across specified axes.
  # if calc_axes = TRUE, will return permutation scores of axis loadings.
  # if axes='all', calculates for all axes. 
  # tail refers to whether one- or two-tailed tests should be used for axis scores
  # adjustment is a direct argument for `p.adjust()`
  # ... are arguments for `permute_ord()`. 
  
  #TODO: Change the defaults on this to something more sensible (currently for backward-compatibility)
  
  if (class(P) == 'pcwOrd') {
    message('Permuting the pcwOrd object')
    P = permute_ord(P, ..., return_permutations=TRUE)
  }
  stopifnot(class(P) == 'pcwOrd_permutations')
  
  tail = match.arg(tail)
  
  if (is.character(axes)) {
    if (tolower(axes[1]) != 'all') {
      stop("Don't understand what axes you want. 'all' or `int`")
    }
    axes = seq_len(length(P$observed$constrained$d))
  }
  
  # Extract information from P
  feature_names = P$observed$col_names
  n_perm = P$results['n']
  observed = ord_scores(P$observed,  # same as P$observed$constrained$v, but with row/column names
                        choice='column',
                        axes=axes)
  scorenames = colnames(observed)
  
  
  # Build results stepwise for backward compatibility...
  out = list() 
  
  # F test
  if (calc_F) {
    message("Calculating a pseudo-F test across all (selected) axes")
    if (length(axes) < length(P$observed$constrained$d)) {
      warning('Pseudo-F tests on subsets of axes may not be valid or meaningful. Consider setting `axes="all"`')
    }
    
    # Extract other information
    # col_names = P$observed$col_names
    nr = length(P$observed$row_weights)
    xcol = ncol(.make_dummy(P$observed$X))
    zcol = ifelse(is.null(P$observed$Z),
                  0,
                  ncol(.make_dummy(P$observed$Z)))
    
    num_df = xcol
    denom_df = nr - num_df - zcol - 1
    
    # Extract permuted variances; calculate F scores
    fit_perm = sapply(P$permutations, function(x) {
      colSums(x$fitted^2)
    })
    resid_perm = sapply(P$permutations, function(x) {
      colSums(x$resids^2)
    })
    F_perm = (fit_perm / num_df) / (resid_perm / denom_df)

    
    # Observed fitted and residual variation, and F score
    obs_fitted = colSums(P$observed$fitted^2)
    obs_resid = colSums(P$observed$Y_residual^2)
    F_obs = (obs_fitted/num_df) / (obs_resid/denom_df)
    
    # Percentiles and adjusted P-value
    p_F = sweep(F_perm, 1, F_obs, '<=')  # This is actually alpha
    p_F = rowSums(p_F) / (n_perm + 1)
    p_F = 1 - p_F
    p_adj_F = p.adjust(p_F, method=adjustment)
    
    F_out = cbind(
      fitted = obs_fitted,
      residuals = obs_resid,
      df_num = num_df,
      df_denom = denom_df,
      F_val = F_obs,
      p_val = p_F,
      adj_p_val = p_adj_F
    )
    # Sort by confidence
    sorter = order(F_out[, 'p_val'], decreasing=FALSE)
    F_out = F_out[sorter, ]
    
    out[['F_test']] = F_out
  }
  
  # Axis tests
  if (calc_axes) {
    message('Calculating axis loadings')
    
    # Compare permuted scores with observed scores
    fit_perm = sapply(P$permutations, function(x) {
      scores = svd(x$fitted, nu=0, nv = max(axes))$v
      rownames(scores) = feature_names
      colnames(scores) = scorenames
      return(scores)
    }, simplify='array')  # features * axes * permutations array
    
    # Calculate p_values for each axis
    p_values = sapply(axes, function(x) {
      obs = observed[, x]
      p_val = fit_perm[, x, ] # creates a features * permutations matrix
      p_val = switch(tail, 
                     'upper' = sweep(p_val, 1, obs, '<='),
                     'lower' = sweep(p_val, 1, obs, '>'),
                     'two' = sweep(abs(p_val), 1, abs(obs), '<=')
      )
      p_val = rowSums(p_val) / (n_perm + 1) # this is alpha
      p_val = 1 - p_val
      return(p_val)
    })
    
    # adjust p_values
    p_adjust = apply(p_values, 2, p.adjust, method=adjustment)

    # Calculate observed percent contribution
    pct = observed^2
    
    # transpose to axis-first and sort by confidence
    axes_out = lapply(axes, function(x) {
      summary = data.frame(
        estimate = observed[, x],
        pct_contribution = pct[, x],
        p_val = p_values[, x],
        adj_p_val = p_adjust[, x]
      )
      sort_order = order(summary['p_val'], decreasing=FALSE)
      summary = summary[sort_order, ]
      
      return(summary)
    })
    names(axes_out) = scorenames
    attr(axes_out, 'description') = 'pcwOrd Feature Loading Permutations'
    attr(axes_out, 'test_tails') = tail
    
    out[['axis_tests']] = axes_out
  }

  if (return_permutations) {
    out[['results']] = P$results
    out[['observed']] = P$observed
    out[['permutations']] = P$permutations
    class(out) = class(P)
  }  
  attr(out, 'n_permutations') = n_perm
  attr(out, 'p_adjustment_method') = adjustment
  
  if (length(out) == 1) out = out[[1]]

  return(out)
}

scale_scores = function(score1, score2=NULL, scaling=NULL) {
  # rescale score1 scores. Two methods:
  # 1. Auto - if score2 is given, then score1 will be multiplied by a constant to match score2
  # 2. Defined - if scaling is given as a scaler. Overriden by score2 if given.
  
  stopifnot(any(class(score1) == 'pcwOrd_scores'))
  
  n_axes = attr(score1, 'nAxes')
  
  if (!is.null(score2)) {
    stopifnot(any(class(score2) == 'pcwOrd_scores'))
    
    n2 = attr(score2, 'nAxes')
    n_axes = min(n_axes, n2)
    
    if (attr(score1, 'scaling') == attr(score2, 'scaling')) {
      msg = paste('score1 and score2 both are', attr(score1, 'scaling'), 'coordinates.')
      warning(msg)
    }
    # calculating the scaling factor
    
    ranges = list(score1, score2)
    ranges = sapply(ranges, function(x) {
      x = x[, c(1:n_axes), drop=FALSE]
      x = rbind(rep(0, n_axes),
                x)
      apply(x, 2, function(y) diff(range(y)))
    }) 
    if (n_axes == 1) ranges = matrix(ranges, nrow=1)
    scaling = mean(ranges[, 2] / ranges[, 1])
    message(paste("Rescaling", attr(score1, 'scaling'), "scores by", signif(scaling, 4)))
  }
  

  score1[, c(1:n_axes)] = score1[, c(1:n_axes), drop=FALSE] * scaling
  return(score1)
}

top_scores = function(R, n=15, choice=c('column', 'row', 'centeroid', 'biplot'),
                      scaling=c('contribution', 'standard', 'principle'),
                      axes='all', constrained=TRUE) {
  
  
  scores = ord_scores(R, choice=match.arg(choice), scaling=match.arg(scaling), constrained=constrained, axes=axes)
  if (axes[1] == 'all') axes = 1:ncol(scores)
  
  # if (weighted) {
  #   weights = ord_variance(R)
  #   if (constrained & !is.null(R$constrained)) {
  #     weights = weights[['constrained']]['value', axes]
  #   } else {
  #     weights = weights[['unconstrained']]['value', axes]
  #   }
  # } else {
  #   weights = rep(1, ncol(scores))[axes]
  # }
  # weights = weights[axes]
  
  n = min(n, nrow(scores))
  
  dist = apply(scores, 1, function(x) sqrt(sum(x^2)))
  dist = sort(dist, decreasing=TRUE)
  dist = as.matrix(dist)[1:n, , drop=FALSE]
  colnames(dist) = paste0('mean_', match.arg(scaling))
  
  # attr(dist, 'mean_calculation') = ifelse(weighted, 'weighted', 'unweighted')
  # attr(dist, 'weights') = weights
  attr(dist, 'axes') = axes
  
  return(dist)
  
}

#### Plot (big) ####
plot_ord = function(R, axes=c(1,2), 
                    row_group=NULL, 
                    max_labels=15, 
                    main=NULL,
                    continuous_scale=NULL,
                    discrete_scale=NULL,
                    shape_scale=NULL,
                    col_text='tomato4',
                    col_arrows='black',
                    shape_legend_position=NULL,
                    color_legend_position=NULL,
                    strip_col_label=NULL) {
  
  #### Generate scores ####
  add_X = TRUE
  add_Z = TRUE
  
  # Identify grouping variables
  nX = ncol(data.frame(R$X))
  nZ = ncol(data.frame(R$Z))
  nG = ncol(data.frame(row_group))
  # if (is.null(row_group)) row_group = FALSE # for downstream compatibility
  
  if (sum(nX, nZ, nG) > 2) {
    add_X = FALSE
    add_Z = FALSE
    if (nX > 0) X = as.matrix(R$X)[, 1, drop=FALSE]
    if (nZ > 0) Z = as.matrix(R$Z)[, 1, drop=FALSE]
    if (nG > 0 & !is.null(row_group)) G = as.matrix(row_group)[R$row_names, , drop=FALSE]
    if (nG >= 2) {
      row_group = G  # manually specified takes precedence
    } else if (nG == 1) {
      first_group = data.frame(X, Z)[, 1, drop=FALSE]
      row_group = data.frame(first_group, G)
    } else {
      row_group = data.frame(X, Z)
    }
    rownames(row_group) = R$row_names
  }
  
  if (is.null(row_group)) row_group=FALSE
  row = ord_scores(R, 'row', 'principle', axes=axes, add_X=add_X, add_Z=add_Z, add_grouping=row_group)
  col = ord_scores(R, 'column', 'standard', axes=axes, add_label=max_labels, strip_label=strip_col_label)
  col = scale_scores(col, row)
  
  bip = NULL
  grp = NULL
  # is_triplot = FALSE
  if (is.factor(R$X[[1]])) {
    grp = ord_scores(R, 'centeroid', 'principle', axes=axes, add_label=TRUE)
    # is_triplot = TRUE
  }
  if (!is.null(R$X)){
    if (is.numeric(as.matrix(R$X[1]))) {
      bip = ord_scores(R, 'biplot', 'standard', axes=axes, add_X=TRUE, add_label=TRUE)
      attr(bip, 'nAxes') = min(length(axes), sum(sapply(R$X, is.numeric)))
      bip = scale_scores(bip, row)
      # is_triplot = TRUE
    }
  }
  
  #### Format colors and shapes ####
  # set scales
  if (is.null(continuous_scale)) {
    continuous_scale = colorRampPalette(colors=c('steelblue4', 
                                                 'darkgoldenrod2', 
                                                 'tomato4'))
  }
  
  if (is.null(discrete_scale)) {  # this is c('Dark2', 'Set3') from RColorBrewer
    discrete_scale = c('#1B9E77','#D95F02','#7570B3','#E7298A',
                       '#66A61E','#E6AB02','#A6761D','#666666',
                       '#8DD3C7','#FFFFB3','#BEBADA','#FB8072',
                       '#80B1D3','#FDB462','#B3DE69','#FCCDE5',
                       '#D9D9D9','#BC80BD','#CCEBC5','#FFED6F')
  }
  
  if (is.null(shape_scale)) {
    shape_scale = c(19, 17, 15, 18, 0, 2, 1, 5, 8:14)
  }
  
  # Set defaults to modify if data requires
  continuous_legend = FALSE
  discrete_legend = FALSE
  colors = discrete_scale[1]
  
  shape_legend = FALSE
  color_legend = FALSE
  shapes = 19
  
  # Select columns (of rowscores) to use for formatting
  # Main difference is continuous vs categorical (factors)
  row_groups = row[, -c(1:2), drop=FALSE] # only two scales currently possible
  is_factor = sapply(row_groups, is.factor)
  is_continuous = sapply(row_groups, is.numeric)
  
  # Continuous colors and categorical shapes
  if (any(is_continuous)) { # continuous takes precedence for colors
    # dinking around to make a continuous colorscale
    row_colors = row_groups[is_continuous][,1]
    breakpoints = findInterval(row_colors,   # via stack exchange
                               sort(row_colors))
    colors = continuous_scale(length(breakpoints))[breakpoints]
    continuous_legend = TRUE
    
    # Shapes reserved for the factors
    if (any(is_factor)) {
      row_shape = row_groups[is_factor][[1]]
      shapes = shape_scale[row_shape]
      shape_legend = TRUE
      col_group = FALSE
    }
    
    # Categorical colors and shapes
  } else if (any(is_factor)) {
    # Categorical colors
    row_colors = row_groups[is_factor][[1]]
    colors = discrete_scale[row_colors]
    
    col_group = TRUE
    
    if (is.null(grp)) {
      color_legend = TRUE
    } else if (any(sapply(levels(grp$LABEL), function(x) !(x %in% levels(row_colors))))) {
      color_legend = TRUE
      col_group = FALSE
    }
    
    if (sum(is_factor) > 1) {
      row_shape = row_groups[is_factor][[2]]
      shapes = shape_scale[row_shape]
      shape_legend=TRUE
    }
  }
  
  # Axis lables 
  variances = ord_variance(R)$unconstrained['pct_total', ]
  # rownames(variances)[2] = 'var'
  if (!is.null(R$constrained)) {
    vv = as.matrix(ord_variance(R)$constrained)
    if (ncol(vv) == 1) {
      vv = vv['pct_total', ]
    } else {
      vv = vv['pct_group', ]
    }
    # rownames(vv)[2] = 'var'
    variances = c(vv, variances)
  }
  variances = variances[axes]
  
  axis_labels = paste(colnames(row), ' [', signif(variances, 3), '%]', sep="")
  
  # Plot....
  all_scores = rbind(row[, 1:2],
                     col[, 1:2],
                     bip[, 1:2],
                     grp[, 1:2])
  plot(all_scores, 
       type='n', asp=1, 
       xlab=axis_labels[1], ylab=axis_labels[2], main=main)
  grid(col='gray80', lty=1)
  axis(1)
  axis(2)
  # axes
  abline(h=0, col='gray50')
  abline(v=0, col='gray50')
  
  text(col[, 1:2], 
       labels=col$LABEL, 
       col=col_text, 
       cex=0.6)
  
  points(row[, 1:2], 
         pch=shapes, 
         col=adjustcolor(colors, alpha.f=0.7), 
         cex=1)
  
  if (!is.null(grp)) {
    if (col_group) {
      col_group = discrete_scale[grp$LABEL]
    } else {
      col_group = '#333333'
    }
    text(grp[, 1:2], 
         label=grp$LABEL, 
         col=adjustcolor(col_group, offset=c(-0.2, -0.2, -0.2, 0)), # darken slightly for easier reading
         font=2)
  }
  
  if (!is.null(bip)) {
    arrows(x0=0, y0=0,
           x1=0.9*bip[, 1],
           y1=0.8*bip[, 2],
           length=0.1,
           col=col_arrows)
    text(bip[, 1:2],
         label=bip$LABEL,
         col=adjustcolor(col_arrows, offset=c(-0.2, -0.2, -0.2, 0)), # darken slightly for easier reading
         cex=0.7)
  }
  
  if (any(shape_legend, continuous_legend, discrete_legend, color_legend)) {
    bounds = c(ymin = min(all_scores[, 2]), # for positioning first legend
               ymax = max(all_scores[, 2]),
               xmin = min(all_scores[, 1]),
               xmax = max(all_scores[, 1]))
    # for positioning second legend relative to first
    max_range = max(c(y = bounds['ymax'] - bounds['ymin'],  # max because aspect=1
                      x = bounds['xmax'] - bounds['xmin']))
    
    base_position = c(bounds['xmin'], bounds['ymax'])
    sl_columns = 0 # shape legend columns. For auto-positioning colors
    sl_rows = 0
    
    if (shape_legend) {
      slp = shape_legend_position
      if (is.null(slp)) slp = base_position
      
      base_position = slp  # update for color legend if needed
      posX = slp[1]
      posY = slp[min(2, length(slp))]
      
      lvls = levels(row_shape)
      nlvls = length(lvls)
      sl_columns = ceiling(nlvls/6)
      sl_rows = ceiling(nlvls/sl_columns)
      
      legend(posX, posY,
             legend=lvls,
             col='black',
             pch=shape_scale[1:nlvls],
             cex=0.7,
             bty='n',
             ncol=sl_columns)
    }
    
    if (color_legend | continuous_legend) {
      clp = color_legend_position
      if (is.null(clp)) clp = c(0, -0.1)
      
      x_offset = max_range * (clp[1] + (clp[1]!=0)*0.1*sl_columns)
      y_offset = max_range * (clp[2] - (clp[2]!=0)*0.025*sl_rows)
      
      posX = base_position[1] + x_offset
      posY = base_position[2] + y_offset
      
      if (is.factor(row_colors)) {
        lvls = levels(row_colors)
        colors = discrete_scale
      } else {
        lvls = round(
          c(min(row_colors, na.rm=TRUE),
            mean(row_colors, na.rm=TRUE),
            max(row_colors, na.rm=TRUE)),
          3)
        colors = continuous_scale(3)
      }
      nlvls = length(lvls)
      nc = ceiling(nlvls/6)
      legend(posX, posY,
             legend=lvls,
             col=colors[1:nlvls],
             pch=19,
             cex=0.7,
             bty='n',
             ncol=nc)
    }
  }
  
}