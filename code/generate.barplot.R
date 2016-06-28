# This script will produce bar graphs from CACNA1C or TCF7L2 phenotype data

## Script for Single Phenotype Barplotting Based on External Residuals/Raw Data CSV File 
  ## Run the analysis.singlepheno script to obtain the csv used for input into this script
  
  
source("misc.R")
source("read.data.R")
source("data.manip.R")
# source("Source/transformation.functions.R")
# source("Source/bar.graph.functions.LS.R")
  
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

  #Covariates
  r2 <- sapply (prepared.pheno, FUN = function(x) {
    if (is.numeric(x)) {
      return (cor(x, prepared.pheno[[phenotype]], use = "complete")^2)
    } else{
      return (NA)
    }
  } )
  
  threshold <- 0.02
  print(r2)
  print(r2[r2 >= threshold], na.print = "")
  
  ##########################################################################################################
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
  
  ##########################################################################################################
  # Exporting
  ##########################################################################################################
  
  #Exporting residuals/raw data for barplotting
  if (residual == TRUE && transform == TRUE) {
    if (export.data) {
      reduced.pheno[[phenotype]] <- out.rm.pheno[[phenotype]]
      write.csv(reduced.pheno, paste0(phenotype, ".outrm.txf.residual.csv"))
    }
  } else if (transform == TRUE && residual == FALSE) {
    if (export.data) {
      reduced.pheno[[phenotype]] <- out.rm.pheno[[phenotype]] 
      write.csv(reduced.pheno, paste0(phenotype, ".outrm.txf.raw.csv"))
    }
    
  } else if (residual == FALSE && transform == FALSE) {
    if (export.data) {
      reduced.pheno[[phenotype]] <- out.rm.pheno[[phenotype]] 
      write.csv(reduced.pheno, paste0(phenotype, ".outrm.raw.csv"))
    }
  }
  
  ##########################################################################################################
  #Run GEMMA
  ##########################################################################################################
  if (rungemma == TRUE) {
    # Load packages
    require(data.table)
    require(qtl)
    source("Source/gemma.functions.R")
    
    # Location of GEMMA and other paramters
    gemma.exe <- "/usr/local/src/gemma/bin/gemma"
    seed      <- 1
    
    system("mkdir ~/Desktop/CAC_gemma_output")
    GENE.INTERACTION <- paste0(GENE, ":strain")
    system(paste0("mkdir ~/Desktop/CAC_gemma_output/", phenotype))
    choose <- model.info[[phenotype]]$gemma #Which GEMMAs to run
    
    for (choose.samples in choose) {
      cp0("On ", phenotype, " and ", choose.samples, "\n")
      
      gemmadir <- paste0("~/Desktop/CAC_gemma_output/", phenotype, "/", choose.samples)
      system(paste0("mkdir ", gemmadir))
      resultsfile <- paste0("~/Desktop/CAC_gemma_output/d3.d2totalactivity/",
                            phenotype, ".", choose.samples, ".RData")
      
      #Rename the dataframe
      pheno <- out.rm.pheno
      if (GENE == "CACNA1C") {
        geno.filename <- "data/mda.F1.lg.csv"
      } else if (GENE == "TCF7L2") {
        geno.filename <- "data/mda.F1.csv"
      }
      
      # Remove MA and fix strain name ambiguity
      pheno <- subset(pheno, strain != "MA")
      if (GENE == "CACNA1C") {
        pheno$strain <- as.character.factor(pheno$strain)
        pheno[["strain"]][which(pheno$strain == "RIIS")] <- rep("RII", length(which(pheno$strain == "RIIS")))
        pheno$strain <- factor(pheno$strain)
      }
      pheno        <- pheno[order(pheno$strain),]
      pheno$strain <- droplevels(pheno$strain)
      
      # Only include samples in the analysis for which the phenotype and all
      # covariates are observed.
      cols  <- c(phenotype,covariates)
      rows  <- which(none.missing.row(pheno[cols]))
      pheno <- pheno[rows,]
      
      # Select either: (1) all samples, (2) hets only, or (3) wild-types
      # only. Unless all samples are chosen, it is important to exclude
      # genotype from the set of covariates. If using all samples, then 
      # the WT/HET should be included as a covariate
      if (choose.samples == "het") {
        rows       <- which(pheno[[GENE]] == "HET")
        pheno      <- pheno[rows,]
        covariates <- setdiff(covariates, GENE)
        cat("Selecting",length(rows),"het samples.\n")
      } else if (choose.samples == "wt") {
        rows       <- which(pheno[[GENE]] == "WT")
        pheno      <- pheno[rows,]
        covariates <- setdiff(covariates, GENE)
        cat("Selecting",length(rows),"wild-type samples.\n")  
      } else {
        covariates <- union(covariates, GENE)
        cat("Including all (",nrow(pheno),") samples in analysis.\n",sep="")
      }
      
      # Convert some columns to binary covariates.
      if (GENE == "CACNA1C") {
        levels(pheno$sex)         <- 0:1
        levels(pheno$CACNA1C)     <- 0:1
        levels(pheno$fctimeofday) <- 0:1
      } else if (GENE == "TCF7L2") {
        levels(pheno$sex)         <- 0:1
        levels(pheno$TCF7L2)      <- 0:1
        levels(pheno$FCtimeofday) <- 0:1
      } else {
        stop (cp0("GENE is set to ", GENE, "...it must be either CACNA1C or TCF7L2"))
      }
      
      # LOAD GENOTYPE DATA
      # ------------------
      # Load the genotype data, and sort the rows by strain. I convert the
      # base-pair positions to Megabases (Mb).
      cat("Loading genotype data.\n")
      d    <- read.mda.F1(geno.filename)
      map  <- d$map
      geno <- d$geno
      geno <- geno[order(rownames(geno)),]
      map  <- cbind(data.frame(snp = rownames(map)),map)
      rownames(map) <- NULL
      rm(d)
      
      # Take out 57L for CAC study
      if (GENE == "CACNA1C") {
        which.is.57L <- which(rownames(geno) == "57L")
        geno <- geno[-which.is.57L, ]
      }
      
      # Drop the sex-linked and mitochondrial DNA from the analysis.
      markers <- which(is.element(map$chr,1:19))
      map     <- transform(map[markers,],chr = droplevels(chr))
      geno    <- geno[,markers]
      
      # Also, drop SNPs that are polymorphic in 5 strains or less, because
      # it is unlikely that we will discover phenotype associations with
      # these SNPs.
      markers <- which(pmin(colSums(geno),colSums(1 - geno)) > 5)
      map     <- map[markers,]
      geno    <- geno[,markers]
      
      # Align the phenotypes and genotypes.
      geno <- geno[match(pheno$strain,rownames(geno)),]
      
      # Initialize the random number generator.
      set.seed(seed)
      
      # Run GEMMA
      out <- run.gemma(phenotype,covariates,pheno,geno,map,gemmadir,gemma.exe)
      gwscan.gemma <- out$gwscan
      pve.gemma    <- out$pve
      
      # Save results to file.
      cat("Saving results to file.\n")
      save(list = c("analysis","choose.samples","map","gwscan.gemma","pve.gemma"),
           file = resultsfile)
    }
  }
  
  
  #Which Gene
  if (GENE == "CACNA1C") {
    pheno.filename <- "data/cac.pheno.29jan2015.csv"
    csv.filename <- paste0(phenotype,".",type.of.data,".csv")
    source("Source/cacdefaults.version6.R")
  } else if (GENE == "TCF7L2") {
    pheno.filename <- "data/pheno.12may2014.csv"
    csv.filename <- paste0(phenotype,".",type.of.data,".csv")
    source("Source/tcfdefaults.version0.R")
  }
  
  # Create the gene interaction term to be used in ANOVA
  GENE.INTERACTION <- paste0(GENE, ":strain")
  
  # Read
  read.exported.csv <- function(csv.filename, GENE, verbose = TRUE) {
    if (verbose) cat(paste0("Reading in csv data from: ", csv.filename, ".\n",
                            "Converting genotype column to factor.\n"))
    csv.data <- read.csv(csv.filename, header = TRUE)
    csv.data[, which(colnames(csv.data) == GENE)] <- factor(csv.data[[GENE]], c("WT", "HET"))
    return(csv.data)
  }
  
  csv.data <- read.exported.csv(csv.filename = csv.filename, GENE = GENE, verbose = TRUE)
  
  # T-Test
  t.test.data <- perform.t.tests(csv.data, phenotype, GENE)
  print(t.test.data)
  
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
  
  
  # Exporting as PNG
  if (type.of.data == "outrm.raw") {
    png(filename=paste0(phenotype, ".", GENE, ".outrm.raw.barplot.png"), width = 1500, height = 600)
    produce.graph(phenotype, formatted.data, type.of.data = type.of.data, 
                  type.of.interaction = type.of.interaction, xaxis.order = xaxis.order,
                  normalized = normalized, annotations = annotations, GENE = GENE)
    dev.off()
  } else if (type.of.data == "outrm.txf.raw") {
    png(filename=paste0(phenotype, ".", GENE, ".outrm.txf.raw.barplot.png"), width = 1500, height = 600)
    produce.graph(phenotype, formatted.data, type.of.data = type.of.data, 
                  type.of.interaction = type.of.interaction, xaxis.order = xaxis.order,
                  normalized = normalized, annotations = annotations, GENE = GENE)
    dev.off()
  } else if (type.of.data == "outrm.txf.residual") {
    png(filename=paste0(phenotype, ".", GENE, ".outrm.txf.residual.barplot.png"), width = 1500, height = 600)
    produce.graph(phenotype, formatted.data, type.of.data = type.of.data, 
                  type.of.interaction = type.of.interaction, xaxis.order = xaxis.order,
                  normalized = normalized, annotations = annotations, GENE = GENE)
    dev.off()
  }
  
  
  
