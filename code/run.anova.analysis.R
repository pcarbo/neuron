# TO DO: Explain here what this script does.
# Three-Way ANOVAs for a Single Phenotype
source("misc.R")
source("transformation.functions.R")
source("read.data.R")
source("data.manip.R")

# Retrieve the analysis settings.
analysis       <- model.info[[phenotype]]
transformation <- analysis$transformation
covariates     <- unique(analysis$covariates)
outliers       <- analysis$outliers

# Ensure that 'gene' is included as a covariate.
covariates <- union(covariates,gene)

# Interaction terms for linear regression.
gene.strain.interaction <- paste0(gene,":strain")
gene.sex.interaction    <- paste0(gene,":sex") 
strain.sex.interaction  <- paste0("strain:sex")
three.way.interaction   <- paste0(gene,":strain:sex")

# Load the phenotype data.
raw.pheno <- read.pheno(cohort)

# Generate additional phenotypes.
prepared.pheno <- create.new.phenotypes(raw.pheno,cohort)

# Change sex into a binary quantity, and other columns into numeric values.
prepared.pheno$sex         <- factor(prepared.pheno$sex)
levels(prepared.pheno$sex) <- 0:1
prepared.pheno$sex         <- as.numeric(prepared.pheno$sex)

# TO DO: Remove these lines of code if they are not needed.
#
# if (GENE == "CACNA1C") {
#   a <- c(25, 75:91)
#   prepared.pheno[, c(25, 76:92)] <- sapply(prepared.pheno[, a], as.numeric)
# } else if (GENE == "TCF7L2") {
#   a <- 76:92
#   prepared.pheno[, 76:92] <- sapply(prepared.pheno[, a], as.numeric)
# }

# Transform the selected phenotype, if requested.
if (!is.null(transformation)) {
  cat("Applying transformation to ",phenotype,".\n",sep="")
  prepared.pheno[[phenotype]] <- transformation(prepared.pheno[[phenotype]])
}

# Remove specified outlying data points.
if (length(outliers) > 0) {
  to.remove <- match(outliers,prepared.pheno$id)
  prepared.pheno[to.remove,phenotype] <- NA
  cat("Removed ",length(to.remove)," outliers for ",phenotype,".\n",sep="")
}
  
# Fit a linear regression model for phenotype ~ covariates + strain +
# interaction terms.
cat("Fitting linear regression and running ANOVA.\n")
gene.interaction <-
  paste(gene.strain.interaction,"+",gene.sex.interaction,"+",
        strain.sex.interaction,"+",three.way.interaction)
f <- paste(phenotype,"~",paste(covariates,collapse=" + "),
           "+ strain + sex +",gene.interaction)
f <- as.formula(f)
g <- lm(f,prepared.pheno)
out.anova <- anova(g)
print(out.anova)
