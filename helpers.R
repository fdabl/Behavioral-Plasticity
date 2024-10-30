#' This function takes a data frame and a BayesFactor fit object (with multiple fits)
#' It then picks the model with the highest marginal likelihood and makes predictions
#' based on its parameters, using the full posterior distribution to quantify uncertainty
#' This does not take the sampling variance into account, which is exactly what we want here,
#' since we are interested in looking at the model-implied values
#' The BayesFactor package doesn't have a built-in function to do that, so we do it manually
predict_bf_model <- function(df, fits) {
  ix <- which.max(fits@bayesFactor[, 1])
  best_fit <- fits[ix]
  post <- posterior(best_fit, iterations = 500, progress = FALSE)
  cnames <- colnames(post)
  
  preds <- c()
  for (i in seq(nrow(df))) {
    income <- df[i, 'combined_income_c']
    country <- as.character(df[i, 'country'])
    
    if (ix == 1) {
      
      # Best model includes only income
      pred <- post[, 'mu'] + post[, 'combined_income_c-combined_income_c'] * income
      
    } else if (ix == 2) {
      
      # Best model includes only country
      pred <- post[, 'mu'] + post[, grepl(country, cnames)]
      
    } else if (ix == 3) {
      
      # Best model has effect of combined income and countries but no interaction
      pred <- post[, 'mu'] + post[, grepl(country, cnames)] + post[, 'combined_income_c-combined_income_c'] * income
      
    } else if (ix == 4) {
      
      # Best model has effect of combined income and countries and interaction
      b_int <- post[, 'combined_income_c-combined_income_c'] + post[, paste0('combined_income_c:country-', country)]
      pred <- post[, 'mu'] + post[, grepl(paste0('^country-', country), cnames)] + b_int * income
    }
    
    preds <- cbind(preds, as.numeric(pred))
  }
  
  preds
}


#' This gets the BayesFactor fit object, using
#' a Cauchy prior with scale sqrt(1/2) for the country effects
#' and a Cauchy prior with scale sqrt(1/2) / 4 for the effect of income
get_bf_fit <- function(
    varname, df,
    rscaleFixed = sqrt(1/2),
    rscaleCont = sqrt(1/2) / 4
) {
  
  form <- as.formula(
    paste0(varname, ' ~ combined_income_c + country + combined_income_c:country')
  )
  
  varname <- sym(varname)
  df_narm <- df %>%
    filter(!is.na(!!varname)) %>% 
    mutate(combined_income_c = combined_income - mean(combined_income))
  
  generalTestBF(
    form, data = df_narm, whichModels = 'withmain',
    rscaleFixed = rscaleFixed, rscaleCont = rscaleCont, progress = FALSE
  )
}


#' We get predictions from the model with the highest marginal likelihood / posterior model probability
#' (assuming uniform model priors) and then adding the mean and the 95% credible interval of the 
#' model-implied values to the data frame
get_df_pred_bf <- function(varname, df) {
  varname <- sym(varname)
  df_narm <- df %>%
    filter(!is.na(!!varname)) %>% 
    mutate(combined_income_c = combined_income - mean(combined_income))
  
  fits <- get_bf_fit(varname, df)
  preds <- predict_bf_model(df_narm, fits)
  
  ypred <- apply(preds, 2, mean)
  ypred_lo <- apply(preds, 2, quantile, 0.025)
  ypred_hi <- apply(preds, 2, quantile, 0.975)
  
  df_pred <- df_narm %>%
    select(!!varname, combined_income, country) %>%
    cbind(data.frame(ypred = ypred, ypred_lo = ypred_lo, ypred_hi = ypred_hi)) %>% 
    group_by(country, combined_income) %>% 
    mutate(
      mean_value = mean(!!varname),
      name = as.character(ensym(varname))
    ) %>% 
    ungroup() %>% 
    select(-1)
  
  df_pred
}


#' Function to compute posterior model probabilities
get_post <- function(fit_all) {
  log_bayes_factors <- c(0, fit_all@bayesFactor[, 1]) # add intercept-only model
  
  # Convert the log Bayes factors to numeric
  log_bf <- as.numeric(log_bayes_factors)
  
  # Find the maximum log Bayes factor for numerical stability
  max_log_bf <- max(log_bf)
  
  # Compute the exponentials and adjust by subtracting the maximum value
  exp_adjusted <- exp(log_bf - max_log_bf)
  
  # Compute the posterior model probabilities
  posterior_probs <- exp_adjusted / sum(exp_adjusted)
  
  posterior_probs
}


normalize_probs <- function(models_post, sel) {
  models_post[sel] / sum(models_post[sel])
}


# This is only for the figure 2 model, assuming a uniform prior
get_bf_inclusion_fig2 <- function(fits) {
  models_post <- Rmpfr::mpfr(get_post(fits), prec = 100)
  
  # Combined_income
  prior_income_included <- 2/4 # since exclude include interaction
  prior_odds_income_included <- prior_income_included / (1 - prior_income_included)
  
  models_post_interaction_rm <- normalize_probs(models_post, seq(4))
  
  posterior_income_included <- sum(models_post_interaction_rm[c(2, 4)])
  posterior_odds_income_included <- posterior_income_included / (1 - posterior_income_included)
  
  bf_inclusion_income <- posterior_odds_income_included / prior_odds_income_included
  
  # Combined_income x Country
  prior_interaction <- 1/2 # is compared only to the model with both main effects
  prior_odds_interaction <- prior_interaction / (1 - prior_interaction)
  
  models_post_interaction_only <- normalize_probs(models_post, c(4, 5))
  posterior_interaction <- models_post_interaction_only[2]
  posterior_odds_interaction <- posterior_interaction / (1 - posterior_interaction)
  
  bf_inclusion_interaction <- posterior_odds_interaction / prior_odds_interaction
  
  as.numeric(c(bf_inclusion_income, bf_inclusion_interaction))
}


#' Compute the inclusion Bayes factors for Figure 2 and store them in a data.frame
get_bf_inclusion_fig2_df <- function(
    actions_order,
    rscaleFixed = sqrt(1/2),
    rscaleCont = sqrt(1/2) / 4
) {
  map_names_rev <- setNames(as.character(names(map_names)), unlist(map_names))
  
  set.seed(1) # as there is always some numerical inaccuracy in the Bayes factor calculation
  fit_all <- lapply(actions_order, function(action) {
    get_bf_fit(map_names_rev[action], df, rscaleFixed = rscaleFixed, rscaleCont = rscaleCont)
  })
  
  df_inclusion_bfs <- round(data.frame(do.call('rbind', lapply(fit_all, get_bf_inclusion_fig2))), 3)
  colnames(df_inclusion_bfs) <- c('income', 'income_x_country')
  df_inclusion_bfs$full_name <- actions_order
  df_inclusion_bfs$rscaleFixed <- rscaleFixed
  df_inclusion_bfs$rscaleCont <- rscaleCont
  df_inclusion_bfs
}


#' This function takes the best fitting model for the policy analysis
#' and computes model-implied values from it, and then merges it with the actual
#' data for plotting later on
get_df_pred_policy <- function(
    fit_best, bp_name, policy_name, df, type = 'individual',
    join = TRUE, use_empirical = TRUE
) {
  
  countries <- c('Denmark', 'United States', 'Nigeria', 'India')
  newdat <- expand.grid(
    country = countries,
    combined_income = seq(15) - mean(seq(15)),
    bp_value = seq(5) - mean(seq(5))
  )
  post <- posterior(fit_best, iterations = 100, progress = FALSE)
  coefficient_names <- colnames(post)
  
  bp_name <- ensym(bp_name)
  policy_name <- ensym(policy_name)
  
  if (use_empirical) {
    
    # We use the empirical data distribution to predict
    newdat <- df %>% select(country, combined_income, !!bp_name)
    colnames(newdat) <- c('country', 'combined_income', 'bp_value')
    
  } else {
    x <- seq(1, 5, 0.25)
    newdat <- expand.grid(
      country = countries,
      combined_income = seq(15),
      bp_value = x
    )
  }
  
  mean_bp_value <- mean(newdat$bp_value)
  newdat$bp_value <- newdat$bp_value - mean_bp_value
  newdat$combined_income <- newdat$combined_income - mean(newdat$combined_income)
  
  df_preds <- c()
  for (i in seq(nrow(newdat))) {
    bp_value <- newdat[i, 'bp_value']
    combined_income <- newdat[i, 'combined_income']
    
    country <- as.character(newdat[i, 'country'])
    country_comb <- paste0('country-', country)
    
    ##############
    # Meat model:
    # curtailment_2 + combined_income + curtailment_2:combined_income + country + curtailment_2:country + combined_income:country
    # Finance / flights model:
    # investment_4 + combined_income + investment_4:combined_income + country + investment_4:country
    ##############
    
    ## Let's write this more generally, so we can also apply it in the loosely matched case
    
    # Main effects
    beta_bp <- 0
    beta_in <- 0
    
    name_bp <- paste0(bp_name, '-', bp_name)
    name_in <- 'combined_income-combined_income'
    
    if (name_bp %in% coefficient_names) {
      beta_bp <- post[, name_bp]
    }
    
    if (name_in %in% coefficient_names) {
      beta_in <- post[, name_in]
    }
    
    # Intercept (effect of country is always there, so merge it in there since it's dummy coded)
    intercept <- post[, 'mu'] + post[, country_comb]
    
    # Pairwise interactions
    beta_bp_in <- 0
    beta_bp_country <- 0
    beta_in_country <- 0
    
    # Three-way interaction
    beta_bp_in_country <- 0
    
    # Let's check for the pairwise interactions
    
    # Interaction of BP x Income
    name_bp_in <- paste0(bp_name, ':combined_income-combined_income')
    if (name_bp_in %in% coefficient_names) {
      beta_bp_in <- post[, name_bp_in]
    }
    
    # Interaction of BP x Country
    name_bp_country <- paste0(bp_name, ':', country_comb)
    
    if (name_bp_country %in% coefficient_names) {
      beta_bp_country <- post[, name_bp_country]
    }
    
    # Interaction of Income x Country
    name_in_country <- paste0('combined_income:', country_comb)
    if (name_in_country %in% coefficient_names) {
      beta_in_country <- post[, name_in_country]
    }
    
    name_bp_in_country <- paste0(bp_name, ':combined_income:', country_comb)
    if (name_bp_in_country %in% coefficient_names) {
      beta_bp_in_country <- post[, name_bp_in_country]
    }
    
    # Let's make predictions
    ypred <- (
      intercept +                                                      # Intercept for the country
        (beta_bp + beta_bp_country) * bp_value +                       # Main effect of Behavioral plasticity + potential interaction with country
        (beta_in + beta_in_country) * combined_income +                # Main effect of income + potential interaction with country
        (beta_bp_in + beta_bp_in_country) * combined_income * bp_value # Potential interaction between bp and income, potentially moderated by country
    )
    
    d <- data.frame(
      times = seq(nrow(post)), ypred = as.numeric(ypred),
      combined_income = combined_income, country = country, bp_value = bp_value
    )
    
    df_preds <- rbind(df_preds, d)
  }
  
  df_preds$name <- as.character(bp_name)
  
  df_pred_avg <- df_preds %>%
    mutate(
      combined_income = as.numeric(factor(combined_income, labels = seq(15))),
      income_bracket = ifelse(
        combined_income <= 5, 'low',
        ifelse(combined_income < 10, 'medium', 'high')
      )
    ) %>% 
    
    # We average predictions across income brackets
    group_by(country, bp_value, income_bracket, times) %>% 
    summarize(ypred_bracket_post = mean(ypred)) %>% 
    
    # And now we average across the posterior distribution
    group_by(country, bp_value, income_bracket) %>% 
    summarize(
      ypred = mean(ypred_bracket_post),
      ypred_lo = quantile(ypred_bracket_post, 0.025),
      ypred_hi = quantile(ypred_bracket_post, 0.975)
    ) %>% 
    mutate(
      name = as.character(ensym(bp_name))
      # country = as.character(country)
    ) %>% 
    ungroup()
  
  # # We look at individual policies, not all averaged together
  if (type == 'individual') {
    df_pred_avg <- df_pred_avg %>% mutate(bp_value = map_back(bp_value))
  }
  # 
  if (type == 'averaged') {
    df_pred_avg <- df_pred_avg %>% mutate(bp_value = bp_value + mean_bp_value)
    return(df_pred_avg)
  }
  
  if (!join) {
    return(df_pred_avg)
  }
  
  bp_name <- ensym(bp_name)
  policy_name <- ensym(policy_name)
  
  df_emp <- df %>% 
    group_by(country, !!bp_name, combined_income) %>%
    summarize(mean_value = mean(!!policy_name)) %>% 
    rename(bp_value = !!bp_name) %>%
    ungroup() %>% 
    mutate(
      combined_income = as.numeric(factor(combined_income, labels = seq(15))),
      income_bracket = ifelse(
        combined_income <= 5, 'low',
        ifelse(combined_income < 10, 'medium', 'high')
      ),
      bp_value = map_back(bp_value)
    )
  
  df_emp %>% 
    left_join(
      df_pred_avg, by = c('country', 'bp_value', 'income_bracket')
    ) %>% 
    ungroup() %>% 
    mutate(
      income_bracket = factor(
        income_bracket,
        levels = c('low', 'medium', 'high'),
        labels = c('Low (1 - 5)', 'Medium (6 - 9)', 'High (10 - 15)')
      )
    )
}


get_size_order <- function(formulas) {
  size <- sapply(formulas, function(form) length(strsplit(form, '\\+')[[1]]))
  order(size)
}


sort_formulas <- function(formulas) {
  formulas[get_size_order(formulas)]
}


#' Returns the index of models that include an interaction with varname
find_interaction <- function(varname, model_formulas, type = 'exclude') {
  if (model_formulas[1] != '1') {
    model_formulas <- c('1', model_formulas)
  }
  
  sel_int_included <- grepl(paste0(varname, ':'), model_formulas) | grepl(paste0(':', varname), model_formulas)
  
  if (type == 'exclude') {
    res <- which(!sel_int_included)
  } else {
    res <- which(sel_int_included)
  }
  
  res
}


#' Returns the index of the models that include varname1 and varname2
find_main <- function(varname1, model_formulas, varname2 = NULL, type = 'include') {
  
  model_formulas <- lapply(
    strsplit(model_formulas, '\\+'), str_trim
  )
  
  sel_incl <- sapply(model_formulas, function(split_formula) {
    included <- varname1 %in% split_formula
    
    if (!is.null(varname2)) {
      included <- included & (varname2 %in% split_formula)
    }
    
    included
  })
  
  if (type == 'include') {
    res <- which(sel_incl)
  } else {
    res <- which(!sel_incl)
  }
  
  res
}


#' Computes the inclusion Bayes factors for Figure 3 (policy analysis)
#' Takes a BayesFactor object that includes all 19 models
get_bf_inclusion_fig3 <- function(
    fits, varname, prior = 'betabinomial', shape1 = 1, shape2 = 1
) {
  
  model_size_prior <- dbetabinom.ab(seq(0, 7), size = 7, shape1 = shape1, shape2 = shape2)
  x <- model_size_prior[1]
  
  models_prior <- c(
    x,
    rep(x / 3, 3),
    rep(x / 3, 3),
    rep(x / 4, 4),
    rep(x / 3, 3),
    rep(x / 3, 3),
    x,
    x
  )
  
  if (prior == 'uniform') {
    models_prior <- rep(1 / 19, 19)
  }
  
  model_formulas <- c('1', rownames(fits@bayesFactor))
  order_size <- get_size_order(model_formulas)
  model_formulas <- sort_formulas(model_formulas)
  models_post <- Rmpfr::mpfr(get_post(fits), prec = 500)[order_size]
  
  ####
  # Main effect: Behavioral plasticity
  ####
  
  ## Find all models that exclude any interaction with bp, renormalize posterior accordingly
  bp_sel_int_excl <- find_interaction(varname, model_formulas, type = 'exclude')
  models_post_bp_int_rm <- normalize_probs(models_post, bp_sel_int_excl)
  models_prior_bp_int_rm <- normalize_probs(models_prior, bp_sel_int_excl)
  
  ## Find all models from those models that include the main effect bp
  bp_sel_incl <- find_main(varname, model_formulas[bp_sel_int_excl], type = 'include')
  
  ## Find all models from those models that exclude the main effect bp
  bp_sel_excl <- find_main(varname, model_formulas[bp_sel_int_excl], type = 'exclude')
  
  ## Calculate prior odds
  prior_bp_included <- sum(models_prior_bp_int_rm[bp_sel_incl])
  prior_odds_bp_included <- prior_bp_included / (1 - prior_bp_included)
  
  ## Sum the posterior of the models that include the main effect of bp 
  posterior_bp_included <- sum(models_post_bp_int_rm[bp_sel_incl])
  posterior_odds_bp_included <- posterior_bp_included / (1 - posterior_bp_included)
  
  ## Calculate the inclusion Bayes factor
  bf_inclusion_bp <- posterior_odds_bp_included / prior_bp_included
  
  ####
  # Interaction effect: Behavioral plasticity x income
  ####
  
  ## Get all models that include both bp and income as a main effect, renormalize posterior accordingly
  bp_income_total <- find_main(varname, model_formulas, varname2 = 'combined_income', type = 'include')
  models_post_both <- normalize_probs(models_post, bp_income_total)
  models_prior_both <- normalize_probs(models_prior, bp_income_total)
  
  ## From that set of models, get those that include the interaction of the two
  bp_income_int_incl <- find_main(paste0(varname, ':combined_income'), model_formulas[bp_income_total], type = 'include')
  bp_income_int_incl <- bp_income_int_incl[-length(bp_income_int_incl)] # need to remove the model with the three way interaction
  
  ## Calculate prior odds
  prior_twoway <- sum(models_prior_both[bp_income_int_incl])
  prior_odds_toway <- prior_twoway / (1 - prior_twoway)
  
  ## Sum the posterior of the models that include the interaction of bp and income
  posterior_twoway <- sum(models_post_both[bp_income_int_incl])
  posterior_odds_twoway <- posterior_twoway / (1 - posterior_twoway)
  
  bf_inclusion_twoway_income <- posterior_odds_twoway / prior_odds_toway
  
  ####
  # Interaction effect: Behavioral plasticity x country
  ####
  
  ## Get all models that include both bp and income as a main effect, renormalize posterior accordingly
  bp_country_total <- find_main(varname, model_formulas, varname2 = 'country', type = 'include')
  models_post_both <- normalize_probs(models_post, bp_country_total)
  models_prior_both <- normalize_probs(models_prior, bp_country_total)
  
  ## From that set of models, get those that include the interaction of the two
  bp_country_int_incl <- find_main(paste0(varname, ':country'), model_formulas[bp_country_total], type = 'include')
  bp_country_int_incl <- bp_country_int_incl[-length(bp_country_int_incl)] # need to remove the model with the three way interaction
  
  ## Calculate prior odds
  prior_twoway <- sum(models_prior_both[bp_country_int_incl])
  prior_odds_toway <- prior_twoway / (1 - prior_twoway)
  
  ## Sum the posterior of the models that include the interaction of bp and income
  posterior_twoway <- sum(models_post_both[bp_country_int_incl])
  posterior_odds_twoway <- posterior_twoway / (1 - posterior_twoway)
  
  bf_inclusion_twoway_country <- posterior_odds_twoway / prior_odds_toway
  
  ####
  # Interaction effect: Behavioral plasticity x income x country
  ####
  
  # We only compare two models here, so that's just the Bayes factor
  bf_inclusion_threeway <- as.numeric((fits[18] / fits[17])@bayesFactor[, 1])
  
  # Behavioral plasticity x income x country
  as.numeric(c(bf_inclusion_bp, bf_inclusion_twoway_country, bf_inclusion_twoway_income, exp(bf_inclusion_threeway)))
}


#' Get BayesFactor fit for loosely domain-matched policy analysis
get_bf_fit_loose <- function(
    df, policy_name, bp_name,
    rscaleFixed = sqrt(1/2),
    rscaleCont = sqrt(1/2) / 4
) {
  
  bp_name <- sym(bp_name)
  policy_name <- sym(policy_name)
  
  df_narm <- df %>% 
    filter(!is.na(!!bp_name), !is.na(!!policy_name)) %>% 
    mutate(
      !!bp_name := !!bp_name - mean(!!bp_name),
      combined_income = combined_income - mean(combined_income)
    )
  
  formula_str <- paste0(deparse(policy_name), ' ~ ', deparse(bp_name), ' * combined_income * country')
  formula <- as.formula(formula_str)
  
  fit <- generalTestBF(
    formula = formula, data = df_narm, whichModels = 'withmain',
    rscaleFixed = rscaleFixed, rscaleCont = rscaleCont, progress = FALSE
  )
  
  fit
}


#' Returns a data.frame of differences in mean behavioral plasticity
#' for each item and country between the top 10% income earners and the rest 
get_top10_difference <- function(varname, df) {
  
  varname <- sym(varname)
  df_narm <- df %>%
    filter(!is.na(!!varname)) %>% 
    mutate(combined_income_c = combined_income - mean(combined_income))
  
  form <- as.formula(
    paste0(varname, ' ~ combined_income_c + country + combined_income_c:country')
  )
  
  newdat <- expand.grid(
    country = countries,
    combined_income_c = seq(15) - mean(seq(15))
  )
  
  fits <- generalTestBF(
    form, data = df_narm, whichModels = 'withmain',
    rscaleFixed = sqrt(1/2), rscaleCont = sqrt(1/2) / 4, progress = FALSE
  )
  
  full_fit <- fits[4]
  post <- posterior(full_fit, iterations = 500, progress = FALSE)
  cnames <- colnames(post)
  
  preds <- c()
  for (i in seq(nrow(newdat))) {
    income <- newdat[i, 'combined_income_c']
    country <- as.character(newdat[i, 'country'])
    
    # Best model has effect of combined income and countries and interaction
    b_int <- post[, 'combined_income_c-combined_income_c'] + post[, paste0('combined_income_c:country-', country)]
    pred <- post[, 'mu'] + post[, grepl(paste0('^country-', country), cnames)] + b_int * income
    
    preds <- cbind(preds, as.numeric(pred))
  }
  
  newdat$combined_income <- newdat$combined_income_c + mean(seq(15))
  
  df_diffs <- do.call('rbind', lapply(countries, function(country) {
    ix <- which(newdat$country == country)
    ix_lo <- ix[seq(1, 9)]
    ix_hi <- ix[seq(10, 15)]
    mean_diffs <- apply(preds[, ix_hi], 1, mean) - apply(preds[, ix_lo], 1, mean)
    mean_diff <- mean(mean_diffs)
    mean_diff_lo <- as.numeric(quantile(mean_diffs, 0.025))
    mean_diff_hi <- as.numeric(quantile(mean_diffs, 0.975))
    
    data.frame(
      country = country, mean_diff = mean_diff,
      mean_diff_lo = mean_diff_lo, mean_diff_hi = mean_diff_hi
    )
  }))
  
  df_diffs$name <- as.character(varname)
  df_diffs
}


wrap_label <- function(label, break_at = c(4, 8)) {
  
  if (label == 'Increasing energy efficiency requirements for buildings') {
    return('Increasing energy efficiency\nrequirements for buildings')
  }
  
  words <- strsplit(label, " ")[[1]]
  if (length(words) > break_at[1]) {
    if (length(words) > break_at[2]) {
      label <- paste(
        paste(words[seq(break_at[1])], collapse = " "), "\n",
        paste(words[seq(break_at[1] + 1, break_at[2])], collapse = " "), "\n",
        paste(words[seq(break_at[2] + 1, length(words))], collapse = " ")
      )
    } else {
      label <- paste(
        paste(words[seq(break_at[1])], collapse = " "), "\n",
        paste(words[seq(break_at[1] + 1, length(words))], collapse = " ")
      )
    }
  }
  
  padded_label <- gsub('\\s*\\n\\s*', '\n', label)
  padded_label
}


wrap_label_vectorized <- function(labels, break_at = c(4, 8)) {
  sapply(labels, wrap_label, break_at)
}


get_cors_all <- function(df_cors, labels, country) {
  dat <- df_cors[df_cors$country == country, ] %>% select(-country)
  
  p <- ncol(dat)
  cormat <- matrix(NA, nrow = p, ncol = p)
  pvals <- cormat
  
  for (i in seq(p)) {
    for (j in seq(p)) {
      res <- cor.test(
        dat[, i], dat[, j],
        use = 'pairwise.complete.obs', method = 'kendall'
      )
      
      cormat[i, j] <- round(res$estimate, 5)
      pvals[i, j] <- round(res$p.value, 5)
    }
  }
  colnames(cormat) <- rownames(cormat) <- labels
  res <- list('country' = country, 'cormat' = cormat, 'pvals' = pvals)
  
  res
}


get_corrplot_all <- function(res) {
  cormat <- res$cormat 
  diag(cormat) <- NA
  
  corrplot(
    cormat, method = 'color', number.cex = 0.60,
    addCoef.col = 'black',
    tl.cex = 0.50, addrect = 20, tl.col = 'black',
    na.label = ' ',
    is.corr = FALSE,
    number.digits = 2, col.lim = c(-0.70, 0.70)
  )
}


#' Run proportion tests
test_prop <- function(df, country, name, type = 'combined') {
  
  df <- df_prep %>% 
    filter(country == !!country, name == !!name)
  
  if (str_starts(name, 'curtailment')) {
    df$value <- df$value + 1
  }
  
  if (type == 'combined') {
    tab <- table(df$value, df$income_bracket)
  } else {
    tab <- table(df$value, df$combined_income)
  }
  
  res <- contingencyTableBF(tab, sampleType = 'indepMulti', fixedMargin = 'cols')
  
  log_bf <- res@bayesFactor$bf
  log_bf
}
