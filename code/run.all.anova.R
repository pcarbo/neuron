# Generates the ANOVA results for all phenotypes analyzed in the
# Tcf7l2 and Cacna1c cohorts.
source("defaults.tcf7l2.R")
source("defaults.cacna1c.R")

# Repeat for each cohort.
for (cohort in c("tcf7l2","cacna1c")) {

  # The label we use for the cohorts is also the same as the name of
  # the data table column specifying the gene het/wild-type.
  gene <- toupper(cohort)
  
  # Get the analysis settings for the selected cohort.
  model.info <- list(tcf7l2  = model.info.tcf7l2,
                     cacna1c = model.info.cacna1c)
  model.info <- model.info[[cohort]]

  # Repeat for each of the selected phenotypes.
  for (phenotype in names(model.info)) {
    cat("COHORT: ",toupper(cohort),", PHENOTYPE: ",toupper(phenotype),"\n",
        sep="")
    source("run.anova.analysis.R")
    cat("\n")
  }
}


