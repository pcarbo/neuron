# TO DO: Explain here what this script does.

# Script parameters.
phenotype <- "totalactivity"
cohort    <- "tcf7l2"

# TO DO: Explain what this means.
gene <- cohort

# TO DO: Explain here what these lines of code are doing.
source("defaults.tcf7l2.R")
source("defaults.cacna1c.R")
model.info <- list(tcf7l2  = model.info.tcf7l2,
                   cacna1c = model.info.cacna1c)
model.info <- model.info[[cohort]]

