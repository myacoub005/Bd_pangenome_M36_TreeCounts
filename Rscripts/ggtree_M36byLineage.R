### Percentage of Presence of the M36 loci in each lineage ####

library(dplyr)
library(rvcheck)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggstance)
library(ggtreeExtra)
library(RColorBrewer)

geneCopies <- read.table("M36_loci.by.lineage.tab", header=TRUE, sep="\t", row.names = 1)
geneCopies

tip_metadata <- read.table("M36_metadata.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

tip_metadata

tree <- read.tree("M36Peptidase_BdHp.tree")
tipnames <- tree$tip.label
tipnames
to_drop <- setdiff(tree$tip.label, rownames(geneCopies))
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label

tipnames

p0 <- ggtree(tree, layout="rectangular") +
  geom_tiplab(size=2, color="black")

p0

p1 <- ggtree(straintree, layout="rectangular") +
  geom_tiplab(size=2, color="black")

p1

p <- p0 %<+% tip_metadata + geom_tippoint(aes(color=Species), size=3)

plot(p)

#tip_metadata <- read.table("Lineages.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

#tip_metadata

#p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage), size=2)
#plot(p)

difftable <- setdiff(rownames(geneCopies), straintree$tip.label)
difftable
geneCopiesFilter <- filter(geneCopies,rownames(geneCopies) %in% straintree$tip.label)

geneCopiesFilter
dd = data.frame(id=straintree$tip.label, value=(geneCopiesFilter$CAPE))
dd


#geneCopiesFilter$CAPE <- factor(geneCopiesFilter$CAPE)
#geneCopiesFilter$ASIA <- factor(geneCopiesFilter$ASIA)
#geneCopiesFilter$BRAZIL <- factor(geneCopiesFilter$BRAZIL)
#geneCopiesFilter$GPL <- factor(geneCopiesFilter$GPL)


geneCounts = data.frame(geneCopiesFilter)

geneCounts

geneCopiesCount = as.matrix(geneCopiesFilter)

geneCopiesCount

# Define the number of colors you want
nb.cols <- 21
mycolors <- colorRampPalette(brewer.pal(10, "Set3"))(nb.cols)
# Create a ggplot with 21 colors
# Use scale_fill_manual


phtmap <- gheatmap(p1, geneCopiesCount, width=0.5, offset=0.5, legend_title="M36 Loci", colnames_position="bottom", colnames_angle=90, colnames_offset_y=0, hjust=1, font.size=2) +   
  scale_fill_gradient2("") +
  scale_y_continuous(limits=c(0, 55)) + theme_tree2(legend.position="right")
phtmap


###### Now let's try a heatmap comparing Long Read with Short Read ########

M36 <- read.table("LR_vs_SR_M36s.tab", header=TRUE, sep="\t", row.names = 1)
M36

tipnames <- tree$tip.label
tipnames
to_drop <- setdiff(tree$tip.label, rownames(M36))
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label
tipnames

difftable <- setdiff(rownames(M36), straintree$tip.label)
difftable
M36Filter <- filter(M36,rownames(M36) %in% straintree$tip.label)

M36Filter

geneCounts = data.frame(M36Filter)

geneCounts

geneCopiesCount = as.matrix(M36Filter)

geneCopiesCount
phtmap1 <- gheatmap(p1, geneCopiesCount, width=0.5, offset=0.5, legend_title="M36 Loci", colnames_position="bottom", colnames_angle=90, colnames_offset_y=0, hjust=1, font.size=2) +   
  scale_fill_continuous("") +
  scale_y_continuous(limits=c(0, 50)) + theme_tree2(legend.position="right")
phtmap1
