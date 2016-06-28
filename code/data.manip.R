# This file contains functions for processing and manipulating the
# phenotype data.

# This function defines some additional phenotypes.
create.new.phenotypes <- function(pheno, cohort) {
  if (cohort == "tcf7l2") {
    cat("Creating (binary) fear-conditioning freezing variables.\n")
    pheno <- transform(pheno,
                       d1tonebin    = factor(as.integer(d1tone    > 1)),
                       d2contextbin = factor(as.integer(d2context > 1)))
   
    cat("Creating (binary) agouti phenotype.\n")
    pheno <- transform(pheno,
                       agouti = factor(as.integer(coatcolor == "agouti")))
    
    cat("Creating binary phenotypes based on group and box number.\n")
    pheno <-
      with(pheno,
        cbind(pheno,
              binary.from.categorical(group,paste0("group",levels(group))),
              binary.from.categorical(PPIbox,paste0("PPIbox",levels(PPIbox)))))
    
  } else if (cohort == "cacna1c") {
    cat("Creating (binary) agouti phenotype.\n")
    pheno <- transform(pheno,
                       agouti = factor(as.integer(agouti == "a")))

    cat("Creating binary phenotypes based on group and box number.\n")
    pheno <-
      with(pheno,
        cbind(pheno,
              binary.from.categorical(group,paste0("group",levels(group))),
              binary.from.categorical(ppibox,paste0("ppibox",levels(ppibox)))))
  } else
    error("Choice of cohort is not valid.")
  
  return(pheno)
}

# Removes outlying values for a given phenotype based on the given
# outlier function. It removes outliers by assigning "NA" based on a
# given outlier function (applied to residuals) and as specified by
# the sample ids. Returns an updated "pheno" data frame with an NA
# assigned to each outlying observation
remove.outliers <- function (pheno, phenotype, covariates,
                             outlier.function, outliers) {
  
  # Ensure that GENE is included as a covariate if this isn't WT combined data
  # (only the WT combined data will have "study" in the phenotype data frame)
  if (is.null(transformed.pheno$study)) covariates <- union(covariates, GENE)
  
  if (verbose) cp0("\n*Starting remove.outliers on ", phenotype, ".\n")
  
  # Check if there is an outlier function, if there isn't return the original data
  # If there is, fix it to account for missing data
  if (is.null(outlier.function)) {
    if (verbose) cp0("Did not apply an outlier removal function on ", phenotype, " (none given).\n")
  } else { # There is an outlier function
    # Create a function that tells if x is an outlier (takes into account NAs)
    is.outlier <- function(x) {
      y <- outlier.function(x)
      y[is.na(y)] <- FALSE
      return(y)
    }
    
    # Now fit a model and save the residuals as r
    if (length(covariates) > 0) {
      f <- formula(paste(phenotype, "~", paste(covariates, collapse = " + ")))
      r <- resid(lm(f, transformed.pheno, na.action = na.exclude))
    } else {
      r <- transformed.pheno[[phenotype]]
    }
    
    # Check the residuals for outliers and report the number that will be removed
    if (verbose) {
      n <- sum(is.outlier(r))
      if (n == 0) {
        cp0("No outliers for ", phenotype)
      } else {
        cp0("Removed ", n, " outliers for ", phenotype)
      }
      if (length(covariates) > 0) {
        cp0(" conditioned on [", paste(covariates, collapse = " + "), "]")
      }
      cp0(" (outlier removal function).\n")    
    }
    
    # Assign a "NA" to outliers
    transformed.pheno[is.outlier(r), phenotype] <- NA
  }
  
  # Check if there are other outliers specified to be removed
  if (length(outliers) > 0) {
    to.remove <- match(outliers, transformed.pheno$id)
    transformed.pheno[to.remove, phenotype] <- NA
    
    # Report amount removed
    if (verbose) {
      cp0("Removed ", length(to.remove), " more outliers for ", phenotype)
      if (length(covariates) > 0) {
        cp0(" conditioned on [" ,paste(covariates, collapse = " + "), "]")
      }
      cp0(" (ID specified).\n")
    }
    
  } else { # No ID specified outliers to remove
    if (verbose) cp0("Did not remove any ID specified outliers (none given).\n")
  }

  return(transformed.pheno)
  
}
