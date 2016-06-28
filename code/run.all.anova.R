# TO DO: Explain here what this script does.

# Script parameters.
phenotype <- "pctdurlight"
cohort    <- "tcf7l2"
gene      <- cohort

source("defaults.tcf7l2.R")
source("defaults.cacna1c.R")
model.info <- list(tcf7l2  = model.info.tcf7l2,
                   cacna1c = model.info.cacna1c)

