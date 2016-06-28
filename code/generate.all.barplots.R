# TO DO: Explain here what this script does.

###########################################################################################
  # Specify a phenotype, what kind of data to use ("raw" or "residual"), and whether or not to normalize the data (TRUE or FALSE)
  phenotype           <- "totalactivity"
  cohort              <- "cacna1c"
  type.of.data        <- "outrm.raw" #either "outrm.raw" or "outrm.txf.raw" or "outrm.txf.residual"
  type.of.interaction <- "strain"
  normalized          <- FALSE
  pheno.unit          <- "units"
  transform <- FALSE #the only combination of choices that will not work: transform - FALSE, residual - TRUE
  residual <- FALSE #if false, the data used for graphing and exporting will not be regressed and will be "raw"
  export.data <- TRUE
  rungemma <- FALSE
  ###########################################################################################

source("defaults.tcf7l2.R")
source("defaults.cacna1c.R")
model.info <- model.info.tcf7l2

source("generate.barplot.R")
