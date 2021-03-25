## install devtools, if necessary:
#install.packages("devtools", dep=TRUE)
## install treeWAS from github:
#install_github("caitiecollins/treeWAS", build_vignettes = TRUE)

# load libs
library(devtools)
library(treeWAS)

# running cfml res scripts
# /Library/Frameworks/R.framework/Versions/4.0/Resources/Rscript cfml_results.R cfml_rso

prefix <- snakemake@params[['prefix']]  #"cfml_rso"
#prefix <- "/Users/devseeva/Desktop/work/sm_workflow/snakefiles/outputs/rso_test_cfml/rso_test_cfml"
dat <- read.CFML(prefix=prefix)

meta <- snakemake@input[['meta']]
#meta <- "/Users/devseeva/Desktop/work/sm_workflow/snakefiles/inputs/rso_test_metadata/rso_test_cleaned_data.csv"
output <- snakemake@output[[1]]
#output <- "rso_test_treeWAS_plots_TTTTTT.pdf"

## Required input into treeWAS:
snps <- dat$snps

## Recommended input into treeWAS:
tree <- dat$tree

## Optional input into treeWAS:
n.subs <- dat$n.subs
snps.rec <- dat$snps.rec

## toy phenotype
set.seed(1)
phen <- sample(c(0,1), nrow(snps), replace=TRUE)
phen <- as.factor(phen)
names(phen) <- rownames(snps)

# import custon metadata
# not categorical
# binary or numeric only!
test_phen <- read.table(file = meta, sep = ',', header = TRUE)
test_phen_CD <- test_phen[["collection_date"]]
#View(test_phen_CD)
#test_phen_CD <- as.factor(test_phen_CD)
names(test_phen_CD) <- rownames(snps)
#View (test_phen_CD)

## Exmine toy phenotype
#str(phen)
#table(phen)

###
out <- treeWAS(snps = snps,
               phen = test_phen_CD,
               tree = tree,
               n.subs = n.subs,
               n.snps.sim = ncol(snps)*10,
               chunk.size = ncol(snps),
               test = c("terminal", "simultaneous", "subsequent"),
               snps.reconstruction = "parsimony",
               snps.sim.reconstruction = "parsimony",
               phen.reconstruction = "parsimony",
               phen.type = NULL,
               na.rm = TRUE,
               p.value = 0.01,
               p.value.correct = "bonf",
               p.value.by = "count",
               dist.dna.model = "JC69",
               plot.tree = TRUE,
               plot.manhattan = TRUE,
               plot.null.dist = TRUE,
               plot.dist = TRUE,
               snps.assoc = NULL,
               filename.plot = output,
               seed = 1)
