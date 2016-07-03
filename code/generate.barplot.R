# TO DO: Explain here what this script does.
#
# This script will produce bar graphs from CACNA1C or TCF7L2 phenotype data.

## Script for Single Phenotype Barplotting Based on External Residuals/Raw Data CSV File 
  ## Run the analysis.singlepheno script to obtain the csv used for input into this script
library(plyr)
library(ggplot2)
source("misc.R")
source("read.data.R")
source("data.manip.R")

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

# Sort the strains by difference in phenotype means.
rows         <- order(t.test.stats$diff,decreasing = TRUE)
t.test.stats <- t.test.stats[rows,]
strains      <- rownames(t.test.stats)
pheno.stats  <- transform(pheno.stats,
                          strain = factor(strain,strains))
                          
# Create the bar chart summarizing strain-specific phenotype differences.
if (gene == "TCF7L2") {
  legend.title <- expression(italic(Tcf7l2) ~ genotype)
} else if (gene == "CACNA1C") {
  legend.title <- expression(italic(Cacna1c) ~ genotype)
}
  
# Generate the bar chart using ggplot2.
r        <- pheno.stats
names(r) <- c("gene","x","n","y","sd","se")
r        <- r[order(r$gene,r$x),]
r        <- cbind(r,data.frame(sig = rep(t.test.stats$sig,2)))
for (i in strains) {
  rows <- which(r$x == i)
  if (r$y[rows[1]] > r$y[rows[2]]) {
    r[rows[2],"sig"] <- " "
  } else {
    r[rows[1],"sig"] <- " "
  }
}
levels(r$gene) <- c("+/+","+/-")
print(ggplot(r,aes(x,y,fill = gene,label = sig)) +
      geom_bar(stat = "identity",position = position_dodge(width = 0.9),
               color = "black",size = 0.3) +
      geom_errorbar(aes(ymin = y - se, ymax = y + se),
                    position = position_dodge(width = 0.9),width = 0.4,
                    size = 0.3) +
      geom_text(aes(y = y + 0.05*diff(range(r$y)) + se),size = 4) +
      scale_fill_manual(values = c("darkblue","lightgray")) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x      = element_text(size = 8,angle = 45),
            axis.title       = element_text(size = 10),
            legend.title     = element_text(size = 9),
            legend.key.size  = unit(12,"points")) +
      labs(x    = "maternal strain",
           y    = phenotype,
           fill = legend.title))
