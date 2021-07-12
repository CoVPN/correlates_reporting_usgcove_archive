# this script contains screener definitions

# set max_var to 20 -- should be much larger than total number of variables anyway
max_var <- 20

# SCREENS	In all 4 screens, if total variables passing the screen > maxVar, 
# then the variables are ranked by significance determined from a logistic regression of 
# each variable with endpoint as outcome, and the most significant variables are selected in the model.	

# screen_all No screen. This screen is run only with non-data-adaptive learners, namely, SL.mean, SL.glm, SL.bayesglm, SL.glm.interaction
screen_all <- eval(parse(text = paste0(
  "function(..., Y, X, id, family = binomial(), max_var = ", max_var, "){\n",
    "out <- SuperLearner:::All(X = X, ...)\n",
    "pvals <- univariate_logistic_pval(Y = Y, X = X[,out,drop=FALSE], family = family)\n",
      "ii <- rank(pvals)\n",
      "if(sum(out) > max_var){\n",
      "  out[ii > max_var] <- FALSE\n",
      "}\n",
    "return(out)\n",
  "}"
)))

# screen_glmnet alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100
screen_glmnet <- eval(parse(text = paste0(
  "function(Y, X, id, family = binomial(), max_var = ", max_var, ", ...){\n",
    "set.seed(123)\n",
    "out <- SuperLearner::screen.glmnet(Y = Y, X = X, family = family, ...)\n",
      "pvals <- univariate_logistic_pval(Y = Y, X = X[,out,drop=FALSE], family = family)\n",
      "ii <- rank(pvals)\n",
      "if(sum(out) > max_var){\n",
      "  out[ii > max_var] <- FALSE\n",
      "}\n",
    "return(out)\n",
  "}"
)))

# screen_univariate_logistic_pval minPvalue = 0.1, minscreen = 2; If number of variables with p-value less than minPvalue is less than minscreen, 
# then issue warning and select the two variables with lowest minimum pvalue.
univariate_logistic_pval <- function(Y, X, family = binomial(), obsWeights = rep(1, length(Y)), ...){
  my_family <- if(all(Y %in% c(0,1))){
    binomial()
  }else{
    gaussian()
  }  
  listp <- apply(X, 2, function(x, Y) {
        fit <- glm(Y ~ x, data = data.frame(Y = Y, x), family = my_family, weights = obsWeights)
        summary_fit <- summary(fit)
        pval <- summary_fit$coefficients[2, 4]
        return(pval)
    }, Y = Y)
    return(listp)
}

screen_univariate_logistic_pval_tmp <- function(Y, X,  family = binomial(), obsWeights = rep(1, length(Y)), id, minPvalue = 0.1, minscreen = 2,...){
	 my_family <- if(all(Y %in% c(0,1))){
    binomial()
  }else{
    gaussian()
  }
  listp <- univariate_logistic_pval(Y = Y, X = X, family = my_family, obsWeights = obsWeights)
	whichVariable <- (listp <= minPvalue)
    if (sum(whichVariable) < minscreen) {
        warning("number of variables with p value less than minPvalue is less than minscreen")
        whichVariable[rank(listp) <= minscreen] <- TRUE
    }
    return(whichVariable)   
}

screen_univariate_logistic_pval <- eval(parse(text = paste0(
  "function(Y, X, family = binomial(), max_var = ", max_var, ", ...){\n",
    "out <- screen_univariate_logistic_pval_tmp(X = X, Y = Y, family = binomial(), ...)\n",
    "pvals <- univariate_logistic_pval(Y = Y, X = X[,out,drop=FALSE], family = family)\n",
      "ii <- rank(pvals)\n",
      "if(sum(out) > max_var){\n",
      "  out[ii > max_var] <- FALSE\n",
      "}\n",
    "return(out)\n",
  "}"
)))

# screen_highcor_random if any two variables have Spearman correlation > 0.9, then choose any one of them at random
screen_highcor_random_tmp <- function(Y, X, family, obsWeights, cor_thresh = 0.9, ...){
	cor_mat <- cor(X, method = "spearman")
	correlated_pairs <- list()
	for(i in 1:(nrow(cor_mat) - 1)){
		for(j in (i+1):ncol(cor_mat)){
			this_cor <- cor_mat[i, j]
			if(this_cor > 0.9){
				correlated_pairs <- c(correlated_pairs, list(c(row.names(cor_mat)[i], colnames(cor_mat)[j])))
			}
		}
	}
	include <- rep(TRUE, ncol(X))
	if(length(correlated_pairs) > 0){
		discard_variables <- unlist(lapply(correlated_pairs, function(z){
			keep <- sample(z, 1)
			discard <- z[!(z %in% keep)]
		}), use.names = FALSE)
		include[colnames(X) %in% discard_variables] <- FALSE
	}
	return(include)
}

screen_highcor_random <- eval(parse(text = paste0(
  "function(Y, X, id, family, max_var = ", max_var, ", ...){\n",
    "out <- screen_highcor_random_tmp(Y = Y, X = X, ...)\n",
    "pvals <- univariate_logistic_pval(Y = Y, X = X[,out,drop=FALSE], family = family)\n",
      "ii <- rank(pvals)\n",
      "if(sum(out) > max_var){\n",
      "  out[ii > max_var] <- FALSE\n",
      "}\n",
    "return(out)\n",
  "}"
)))