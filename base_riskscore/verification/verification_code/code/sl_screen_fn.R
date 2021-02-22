# this script contains screener definitions

# read in max_var
max_var <- readRDS(here::here("output", "max_var.rds"))

# SCREENS	In all 4 screens, if total variables passing the screen > maxVar, 
# then the variables are ranked by significance determined from a logistic regression of 
# each variable with endpoint as outcome, and the most significant variables are selected in the model.	

# screen_all No screen. This screen is run only with non-data-adaptive learners, namely, SL.mean, SL.glm, SL.bayesglm, SL.glm.interaction
screen_all <- eval(parse(text = paste0(
  "function(..., X, family, max_var = ", max_var, "){\n",
    "out <- SuperLearner:::All(X = X, ...)\n",
    "if(sum(out) > max_var){\n",
      "pvals <- univariate_logistic_pval(X = X, family = family, ...)\n",
      "ii <- order(pvals, decreasing=FALSE)[1:min(max_var,length(pvals))]\n",
      "idx <- ii[!is.na(pvals[ii])]\n",
      "out <- rep(FALSE, ncol(X))\n",
      "out[idx] <- TRUE\n",
    "}\n",
    "return(out)\n",
  "}"
)))

# screen_glmnet alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100
screen_glmnet <- eval(parse(text = paste0(
  "function(..., X, family, max_var = ", max_var, "){\n",
    "out <- SuperLearner::screen.glmnet(X = X, ...)\n",
    "if(sum(out) > max_var){\n",
      "pvals <- univariate_logistic_pval(X = X, family = family, ...)\n",
      "ii <- order(pvals, decreasing=FALSE)[1:min(max_var,length(pvals))]\n",
      "idx <- ii[!is.na(pvals[ii])]\n",
      "out <- rep(FALSE, ncol(X))\n",
      "out[idx] <- TRUE\n",
    "}\n",
    "return(out)\n",
  "}"
)))

# screen_univariate_logistic_pval minPvalue = 0.1, minscreen = 2; If number of variables with p-value less than minPvalue is less than minscreen, 
# then issue warning and select the two variables with lowest minimum pvalue.
univariate_logistic_pval <- function(Y, X, family, obsWeights = rep(1, length(Y)), ...){
	listp <- apply(X, 2, function(x, Y) {
        fit <- glm(Y ~ x, family = family, weights = obsWeights)
        summary_fit <- summary(fit)
        pval <- summary_fit$coefficients[2, 4]
        return(pval)
    }, Y = Y)
    return(listp)
}

screen_univariate_logistic_pval_tmp <- function(Y, X, family, obsWeights, minPvalue = 0.1, minscreen = 2){
	listp <- univariate_logistic_pval(Y = Y, X = X, family = family, obsWeights = obsWeights)
	whichVariable <- (listp <= minPvalue)
    if (sum(whichVariable) < minscreen) {
        warning("number of variables with p value less than minPvalue is less than minscreen")
        whichVariable[rank(listp) <= minscreen] <- TRUE
    }
    return(whichVariable)   
}

screen_univariate_logistic_pval <- eval(parse(text = paste0(
  "function(..., X, family, max_var = ", max_var, "){\n",
    "out <- screen_univariate_logistic_pval_tmp(X = X, ...)\n",
    "if(sum(out) > max_var){\n",
      "pvals <- univariate_logistic_pval(X = X, family = family, ...)\n",
      "ii <- order(pvals, decreasing=FALSE)[1:min(max_var,length(pvals))]\n",
      "idx <- ii[!is.na(pvals[ii])]\n",
      "out <- rep(FALSE, ncol(X))\n",
      "out[idx] <- TRUE\n",
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