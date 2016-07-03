# TO DO: Explain here what this script does.
#
# This script will produce bar graphs from CACNA1C or TCF7L2 phenotype data.

## Script for Single Phenotype Barplotting Based on External Residuals/Raw Data CSV File 
  ## Run the analysis.singlepheno script to obtain the csv used for input into this script
library(plyr)
source("misc.R")
source("read.data.R")
source("data.manip.R")

# *** TEMPORARY ***
phenotype <- "fastglucose"
cohort    <- "tcf7l2"
gene      <- "TCF7L2"

# Load the phenotype data.
raw.pheno <- read.pheno(cohort)

# Compute t-test statistics for assessing strain-specific phenotype effects.
cat("Computing t-test statistics.\n")
strains      <- levels(raw.pheno$strain)
n            <- length(strains)
t.test.stats <- data.frame(wt      = rep(0,n),
                           het     = rep(0,n),
                           diff    = rep(0,n),
                           p.value = rep(0,n),
                           sig     = rep(as.character(NA),n))
rownames(t.test.stats) <- strains
for (i in strains) {
  r   <- subset(raw.pheno,strain == i)
  out <- t.test(y ~ x,data.frame(x = r[[gene]],y = r[[phenotype]]))
  t.test.stats[i,c("wt","het","diff","p.value")] <-
    with(out,
         c(estimate["mean in group WT"],
           estimate["mean in group HET"],
           estimate["mean in group WT"] - estimate["mean in group HET"],
           p.value))
}

# Adds a symbol denoting significance based on the t-test p-value.
x         <- cut(t.test.stats$p.value,c(-1,0.001,0.01,0.05,0.1,Inf))
levels(x) <- c("***","**","*","-"," ")
t.test.stats$sig <- as.character(x)
         
# Sort the strains by difference in phenotype means.
t.test.stats <- t.test.stats[order(t.test.stats$diff,decreasing = TRUE),]
rm(i,n,r,x,out)

# Compute phenotype summary statistics stratified by strain and gene
# (het vs. wild-type).
compute.pheno.stats <- function (x, i) {
  x <- x[[i]]
  n <- sum(!is.na(x))
  return(c(n    = n,
           mean = mean(x,na.rm = TRUE),
           sd   = sd(x,na.rm = TRUE),
           se   = sd(x,na.rm = TRUE)/sqrt(n)))
}
pheno.stats <- ddply(raw.pheno,c(gene,"strain"),compute.pheno.stats,phenotype)
    
stop()

# Retrieve the analysis settings.
analysis         <- model.info[[phenotype]]
transformation   <- analysis$transformation
covariates       <- unique(analysis$covariates)
outliers         <- analysis$outliers
outlier.function <- analysis$outlier.function
  
#Interactions
# GENE.STRAIN.INTERACTION <- paste0(GENE, ":strain")
# GENE.SEX.INTERACTION    <- paste0(GENE, ":sex") 
# STRAIN.SEX.INTERACTION  <- paste0("strain:sex")
# THREE.WAY.INTERACTION   <- paste0(GENE, ":strain:sex")
  

  # Regression and the making of "reduced.pheno"
  ##########################################################################################################
  #Run linear model and get residuals if there are covariates, otherwise use raw phenotype values 
  if (residual) { #if residual = TRUE
    if (is.null(covariates)) {
      resid <- prepared.pheno[[phenotype]] #transformed, not regressed because there are no covariates
    } else {
      f <- paste0(phenotype, " ~ ", paste(covariates, collapse = " + "))
      model  <- lm(f, prepared.pheno, na.action = na.exclude)
      resid <- resid(model)
    }
  } else { #if residual = FALSE
    resid <- prepared.pheno[[phenotype]] #note: keep in mind raw values are saved in a variable called "resid" when no linear regression is run
  }
  
  #Format residuals/raw values data (and export as a text file if needed)
  reduced.pheno <- data.frame(cbind(prepared.pheno$id, prepared.pheno$strain, raw.pheno$group, raw.pheno$sex, raw.pheno[[GENE]], resid))
  names(reduced.pheno) <- c("ID", "strain", "group", "sex", GENE, phenotype)
  reduced.pheno$strain <- raw.pheno$strain
  reduced.pheno$group <- raw.pheno$group
  reduced.pheno$sex <- raw.pheno$sex
  reduced.pheno[[GENE]] <- raw.pheno[[GENE]]
  
  ##########################################################################################################
  #Outliers
  ##########################################################################################################
  #Find mean and standard deviation + lower/upper bounds
  mean <- mean(resid, na.rm = TRUE)
  sd <- sd(resid, na.rm = TRUE)
  
  #Prints index of possible outliers
  mult <- 2.5 #Number of standard deviations away from the mean
  upper <- which(resid > (mean + sd*mult))
  lower <- which(resid < (mean - sd*mult))
  out <- c(upper, lower)
  sum(length(upper), length(lower))
  
  #Prints IDs of possible outliers
  id <- reduced.pheno[out,]
  
  #Scatterplot of residuals with outliers in a different colour
  require(ggplot2)
  graphing.data <- data.frame(cbind(resid, 0))
  names(graphing.data) <- c("Values", "Zero")
  
  left <- min(graphing.data$Values, na.rm=TRUE) - 1.5*sd
  right <- max(graphing.data$Values, na.rm=TRUE) + 1.5*sd
  high <- max(left, right)
  low <- high*-1
  
  g <- ggplot(graphing.data, aes(x = Values, y = Zero, color = "blue")) +
    geom_segment(aes(x = (mean - sd*mult), xend = (mean + sd*mult), y = 0, yend = 0), color = "peachpuff1", lwd=20) + 
    geom_point(shape=16, size=2, colour = "dodgerblue", na.rm = TRUE) +
    scale_x_continuous(limits = c(left, right)) +
    ylab(NULL) +
    xlab("Values") +
    ggtitle(paste0(phenotype, " Values")) +
    theme_bw() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 12, colour = "black"),
          title = element_text(size = 14, face = "bold"),
          axis.text.y = element_blank())
  
  plot(g)
  
  #Removing Outliers 
  out <- analysis$outliers
  prepared.pheno[[phenotype]] <- reduced.pheno[[phenotype]] #This step is necessary if there are covariates and if the data was linearly regressed
  out.rm.pheno <- remove.outliers(prepared.pheno, phenotype, covariates, outlier.function = NULL, out, verbose = TRUE) #not regressed
  
  # Create the gene interaction term to be used in ANOVA
  GENE.INTERACTION <- paste0(GENE, ":strain")

  # Formatted data
  group.variables <- c(GENE, type.of.interaction)
  
  format.data <- function(pheno.data, phenotype, group.variables) {
    reduce.data <- function(x, col) {
      return(c(N    = len(x[[col]], na.rm = TRUE),
               mean = mean(x[[col]], na.rm = TRUE),
               sd   = sd(x[[col]], na.rm = TRUE)))
    }
    
    formatted.data    <- ddply(pheno.data, group.variables, .fun = reduce.data, phenotype)
    formatted.data    <- rename(formatted.data, c("mean" = phenotype))
    formatted.data$se <- formatted.data$sd / sqrt(formatted.data$N)
    
    return(formatted.data)
  }
  formatted.data <- format.data(csv.data, phenotype = phenotype, group.variables = group.variables)
  
  # Data sorted by largest to smallest difference between WT and HET
  xaxis.order <- t.test.data$strain
  
  # Create significance level annotations
  annotations <- annotate.significance(formatted.data, t.test.data, xaxis.order, phenotype,
                                       normalized = normalized, mult.factor = 1.05)
  
  # Produce the bar graph
  produce.graph(phenotype, formatted.data, type.of.data = type.of.data, 
                type.of.interaction = type.of.interaction, xaxis.order = xaxis.order,
                normalized = normalized, annotations = annotations, GENE = GENE)
  
  
