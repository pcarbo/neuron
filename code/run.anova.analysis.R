# TO DO: Explain here what this script does.

# Three-Way ANOVAs for a Single Phenotype

###########################################################################################
#Parameters
phenotype <- "pctdurlight"
GENE <- "CACNA1C"
###########################################################################################

#Source
# source("Source/data.manip.R")
# source("Source/transformation.functions.R")
# source("Source/misc.R")

#Default
# analysis         <- model.info[[phenotype]]
# transformation   <- analysis$transformation
# covariates       <- unique(analysis$covariates)
# outliers         <- analysis$outliers

#Interactions
# GENE.STRAIN.INTERACTION <- paste0(GENE, ":strain")
# GENE.SEX.INTERACTION    <- paste0(GENE, ":sex") 
# STRAIN.SEX.INTERACTION  <- paste0("strain:sex")
# THREE.WAY.INTERACTION   <- paste0(GENE, ":strain:sex")

# Load the phenotype data.
raw.pheno <- read.pheno(cohort)

# Prepare data by creating new phenotypes
create.new.phenotypes <- function(raw.pheno.data, GENE, verbose = TRUE) {
  if (GENE == "TCF7L2") {
    if (verbose) cp0("Creating weight gain phenotype.\n")
    weightgain <- resid(lm(bw3 ~ bw1, raw.pheno.data, na.action = na.exclude))
    
    if (verbose) cp0("Creating binary FC freezing variables.\n")
    # If some freezing (greater than 1) observation converted to 1
    # If minimal freezing (less than 1) observation converted to 0
    d1pretrainbin <- factor(as.integer(raw.pheno.data$d1pretrain > 1))
    d1tonebin     <- factor(as.integer(raw.pheno.data$d1tone     > 1)) 
    d2contextbin  <- factor(as.integer(raw.pheno.data$d2context  > 1))
    
    if (verbose) cp0("Creating binary agouti phenotype.\n")
    agouti <- factor(as.integer(raw.pheno.data$coatcolor == "agouti"))
    
    if (verbose) cp0("Creating binary phenotypes based on group number and box number (PPI/OFT).\n")
    groupsbin <- binary.from.categorical(raw.pheno.data$group,
                                         col.names = paste0("group", levels(raw.pheno.data$group)))
    ppiboxbin <- binary.from.categorical(raw.pheno.data$PPIbox,
                                         col.names = paste0("PPIbox", 1:5))
    oftboxbin <- binary.from.categorical(raw.pheno.data$oftbox,
                                         col.names = paste0("oftbox", 1:12))
    # Add new phenotypes to the data frame
    new.pheno.data <- cbind(raw.pheno.data, weightgain, d1pretrainbin, d1tonebin, d2contextbin,
                            agouti, groupsbin, ppiboxbin, oftboxbin)
    return(new.pheno.data)
    
  } else if (GENE == "CACNA1C") {
    if (verbose) cp0("Creating weight gain phenotype.\n")
    weightgain <- resid(lm(d100bw ~ d50bw, raw.pheno.data, na.action = na.exclude))
    
    #if (verbose) cp0("Creating D3 total activity phenotype.\n")
    d3totalactivity <- raw.pheno.data[["d2totalactivity"]] + raw.pheno.data[["d3.d2totalactivity"]]
    
    if (verbose) cp0("Creating binary phenotypes based on group, OFT box, and PPI box.\n")
    groupsbin <- binary.from.categorical(raw.pheno.data$group,
                                         col.names = paste0("group", levels(raw.pheno.data$group)))
    ppiboxbin <- binary.from.categorical(raw.pheno.data$ppibox,
                                         col.names = paste0("ppibox", 1:5))
    oftboxbin <- binary.from.categorical(raw.pheno.data$oftbox,
                                         col.names = paste0("oftbox", 1:12))
    
    if (verbose) cp0("Creating binary agouti phenotype.\n")
    agouti <- factor(as.integer(raw.pheno.data$agouti == "a"))
    raw.pheno.data <- raw.pheno.data[, -which(colnames(raw.pheno.data) == "agouti")]
    
    # Add new phenotypes to the data frame
    new.pheno.data <- cbind(raw.pheno.data, groupsbin, ppiboxbin, oftboxbin,
                            weightgain, d3totalactivity, agouti)
    return(new.pheno.data)
  } else {
    cp0("Something went wrong.\n")
    return(invisible(0))
  }
}

prepared.pheno <- create.new.phenotypes(raw.pheno, GENE = GENE)

#Change sex into a binary variable and other columns into numeric
prepared.pheno$sex <- factor(prepared.pheno$sex)
levels(prepared.pheno$sex) <- c(0, 1)
prepared.pheno$sex <- as.numeric(prepared.pheno$sex)

if (GENE == "CACNA1C") {
  a <- c(25, 75:91)
  prepared.pheno[, c(25, 76:92)] <- sapply(prepared.pheno[, a], as.numeric)
} else if (GENE == "TCF7L2") {
  a <- 76:92
  prepared.pheno[, 76:92] <- sapply(prepared.pheno[, a], as.numeric)
}

#Transform
prepared.pheno <- transform.pheno(prepared.pheno, phenotype, transformation)

out <- analysis$outliers
prepared.pheno <- remove.outliers(prepared.pheno, phenotype, covariates, outlier.function = NULL, out, verbose = TRUE)

#Perform 3-Way ANOVA
run.anova <- function(trans.out.rm.pheno, phenotype, covariates) {
  
  # Ensure that GENE is included as a covariate
  covariates <- union(covariates, GENE)
  
  # From the full data set, extract out columns corresponding to ID, strain, phenotype, and covariate(s)
  model.data <- trans.out.rm.pheno[, c("id", "strain", "sex", phenotype, covariates)]
  
  # Fit the model: phenotype ~ covariates + strain + GENE.INTERACTION
  GENE.INTERACTION <- paste(GENE.STRAIN.INTERACTION, "+", GENE.SEX.INTERACTION, "+", STRAIN.SEX.INTERACTION, "+", THREE.WAY.INTERACTION)
  
  string.formula <- 
    paste(phenotype, "~", paste(covariates, collapse = " + "), "+ strain +", "sex +", GENE.INTERACTION)
  f <- as.formula(string.formula)
  g <- lm(f, model.data)
  
  return(anova(g))
}

anova.results <- run.anova(prepared.pheno, phenotype, covariates)
print(anova.results)
anova.results$"Pr(>F)"

