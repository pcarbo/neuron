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
