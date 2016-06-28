# Three-Way ANOVAs for a Single Phenotype

source("misc.R")
source("read.data.R")
source("data.manip.R")
# source("Source/transformation.functions.R")

# Retrieve the analysis settings.
analysis       <- model.info[[phenotype]]
transformation <- analysis$transformation
covariates     <- unique(analysis$covariates)
outliers       <- analysis$outliers

#Interactions
# GENE.STRAIN.INTERACTION <- paste0(GENE, ":strain")
# GENE.SEX.INTERACTION    <- paste0(GENE, ":sex") 
# STRAIN.SEX.INTERACTION  <- paste0("strain:sex")
# THREE.WAY.INTERACTION   <- paste0(GENE, ":strain:sex")

# Load the phenotype data.
raw.pheno <- read.pheno(cohort)

# Generate additional phenotypes.
prepared.pheno <- create.new.phenotypes(raw.pheno,cohort)

# Change sex into a binary quantity, and other columns into numeric values.
prepared.pheno$sex         <- factor(prepared.pheno$sex)
levels(prepared.pheno$sex) <- c(0, 1)
prepared.pheno$sex         <- as.numeric(prepared.pheno$sex)

# if (GENE == "CACNA1C") {
#   a <- c(25, 75:91)
#   prepared.pheno[, c(25, 76:92)] <- sapply(prepared.pheno[, a], as.numeric)
# } else if (GENE == "TCF7L2") {
#   a <- 76:92
#   prepared.pheno[, 76:92] <- sapply(prepared.pheno[, a], as.numeric)
# }

# Transform the selected phenotype, if requested.
if (apply.transform & !is.null(transformation)) {
  cat("Applying transformation to ",phenotype,".\n",sep="")
  prepared.pheno[[phenotype]] <- transformation(prepared.pheno[[phenotype]])
}

stop()

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

