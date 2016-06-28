# This file contains a few functions to read and process data stored
# in text files.

# Read Tcf7l2 or Cacna1c cohort phenotype data from a CSV file.
read.pheno <- function (cohort)
  if (cohort == "tcf7l2") {

    # Read the data table from the CSV file.
    cat("Reading phenotype data from ../data/pheno.tcf7l2.csv.\n")
    out <-
      read.csv("../data/pheno.tcf7l2.csv",header = TRUE,check.names = FALSE, 
               as.is = c("TCF7L2","FCtimeofday"))

    # Manually convert some of the table columns to factors.
    cat("Converting some phenotype table columns to factors.\n")
    return(transform(out,
                     TCF7L2       = factor(TCF7L2,c("WT","HET")),
                     FCtimeofday  = factor(FCtimeofday,c("AM","PM")),
                     PPIbox       = factor(PPIbox)))
  } else if (cohort == "cacna1c") {
    
    # Read the data table from the CSV file.
    cat("Reading phenotype data from ../data/pheno.cacna1c.csv.\n")
    out <- read.csv("../data/pheno.cacna1c.csv",header = TRUE,
                    check.names = FALSE,as.is = "CACNA1C")
    
    # Manually convert some of the table columns to factors.
    cat("Converting some phenotype table columns to factors.\n")
    return(transform(out,
                     CACNA1C = factor(CACNA1C,c("WT","HET")),
                     ppibox  = factor(ppibox)))
  } else
    error("Choice of cohort is not valid.")
