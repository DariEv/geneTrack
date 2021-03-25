#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")
# install.packages('rhierbaps')

library(rhierbaps)
library(ggtree)

#fasta.file.name <- "/Users/devseeva/Desktop/work/sm_workflow/snakefiles/inputs/panX/rso_test/SNP_whole_matrix_names_corrected.aln"
snp.matrix <- load_fasta(snakemake@input[['msa']])
hb.results <- hierBAPS(snp.matrix, max.depth = 1, n.pops = 20, quiet = TRUE)
head(hb.results$partition.df)

#newick.file.name <- "/Users/devseeva/Desktop/work/pan-genome-visualization/data/rso_test/vis/strain_tree.nwk"
iqtree <- phytools::read.newick(snakemake@input[['tree']])

pdf(snakemake@output[['plot']])
write.csv(x=hb.results$partition.df, file=snakemake@output[['table']])

gg <- ggtree(iqtree) #, layout = "circular")
gg <- gg %<+% hb.results$partition.df
gg <- gg + geom_tippoint(aes(color = factor(`level 1`)))
gg

dev.off()

# visualize the small test dataset on the complete tree 

#newick.file.name <- "/Users/devseeva/Desktop/work/pan-genome-visualization/data/ncbi-rso/vis/strain_tree.nwk"
#iqtree <- phytools::read.newick(newick.file.name)
#gg <- ggtree(iqtree) #, layout = "circular")
#gg <- gg %<+% hb.results$partition.df
#gg <- gg + geom_tippoint(aes(color = factor(`level 3`)))
#gg

## big data set! takes a looooooong time

### saved for previous run 
## write.csv(x=hb.results$partition.df, file="/Users/devseeva/Desktop/ncbi_rso_rhierBAPS.csv")

#fasta.file.name <- "/Users/devseeva/Desktop/SNP_whole_matrix.aln"
#snp.matrix <- load_fasta(fasta.file.name)
#hb.results <- hierBAPS(snp.matrix, max.depth = 3, n.pops = 20, quiet = TRUE)
#head(hb.results$partition.df)

#newick.file.name <- "/Users/devseeva/Desktop/work/pan-genome-visualization/data/ncbi-rso/vis/strain_tree.nwk"
#iqtree <- phytools::read.newick(newick.file.name)

#gg <- ggtree(iqtree)#, layout = "circular")
#gg <- gg %<+% hb.results$partition.df
#gg <- gg + geom_tippoint(aes(color = factor(`level 1`)))
#gg