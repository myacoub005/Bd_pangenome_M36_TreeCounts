######## PAV of M36 Loci from cblaster across Bd Pangenome ##############

library(ggtree)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(reshape)
library(rlang)
library(phytools) # for sims and ASRs
library(ape)
library(geiger)
library(tidyverse)
library(RColorBrewer)
library(ggtree)
library(phyloseq)
library(dplyr)


geneCopies <- read.table("M36_cblast_updated.tab", header=TRUE, sep="\t", row.names=1)
geneCopies

tree <- read.tree("Reroot_tree")
tipnames <- tree$tip.label

to_drop <- setdiff(tree$tip.label, rownames(geneCopies))
to_drop
straintree <- drop.tip(tree, to_drop)
tipnames <- straintree$tip.label


p1 <- ggtree(straintree, layout="circular") +
  geom_tiplab(size=2, color="black")
               
p1

tip_metadata <- read.table("Lineages.tab", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

tip_metadata

p <- p1 %<+% tip_metadata + geom_tippoint(aes(color=Lineage), size=2)
plot(p)

difftable <- setdiff(rownames(geneCopies), straintree$tip.label)
geneCopiesFilter <- filter(geneCopies,rownames(geneCopies) %in% straintree$tip.label)
geneCopiesFilter$BDEG_20379 <- factor(geneCopiesFilter$BDEG_20379)
geneCopiesFilter$BDEG_22774 <- factor(geneCopiesFilter$BDEG_22774)
geneCopiesFilter$BDEG_22775 <- factor(geneCopiesFilter$BDEG_22775)
geneCopiesFilter$X09.OZ_003162 <- factor(geneCopiesFilter$X09.OZ_003162)

geneCopiesFilter$X02.OZ_003076 <- factor(geneCopiesFilter$X02.OZ_003076)
geneCopiesFilter$X11.OZ_002774 <- factor(geneCopiesFilter$X11.OZ_002774)
geneCopiesFilter$BDEG_25898 <- factor(geneCopiesFilter$BDEG_25898)
geneCopiesFilter$X02.OZ_004309 <- factor(geneCopiesFilter$X02.OZ_004309)

geneCopiesFilter$BDEG_27285 <- factor(geneCopiesFilter$BDEG_27285)
geneCopiesFilter$X02.OZ_005874 <- factor(geneCopiesFilter$X02.OZ_005874)
geneCopiesFilter$BDEG_25898 <- factor(geneCopiesFilter$BDEG_25898)
geneCopiesFilter$X11.OZ_004862 <- factor(geneCopiesFilter$X11.OZ_004862)

geneCopiesFilter$BDEG_23209 <- factor(geneCopiesFilter$BDEG_23209)
geneCopiesFilter$X11.OZ_007038 <- factor(geneCopiesFilter$X11.OZ_007038)
geneCopiesFilter$X36.OZ_001712 <- factor(geneCopiesFilter$X36.OZ_001712)
geneCopiesFilter$BDEG_23332 <- factor(geneCopiesFilter$BDEG_23332)
geneCopiesFilter$BDEG_27850 <- factor(geneCopiesFilter$BDEG_27850)
geneCopiesFilter$X03.OZ_003608 <- factor(geneCopiesFilter$X03.OZ_003608)
geneCopiesFilter$X03.OZ_002767 <- factor(geneCopiesFilter$X03.OZ_002767)
geneCopiesFilter$BDEG_27858 <- factor(geneCopiesFilter$BDEG_27858)
geneCopiesFilter$X12.OZ_006060 <- factor(geneCopiesFilter$X12.OZ_006060)
geneCopiesFilter$BDEG_27864 <- factor(geneCopiesFilter$BDEG_27864)
geneCopiesFilter$X38.OZ_006340 <- factor(geneCopiesFilter$X38.OZ_006340)
geneCopiesFilter$BDEG_23888 <- factor(geneCopiesFilter$BDEG_23888)

geneCopiesFilter$X06.OZ_002034 <- factor(geneCopiesFilter$X06.OZ_002034)
geneCopiesFilter$X13.OZ_006875 <- factor(geneCopiesFilter$X13.OZ_006875)
geneCopiesFilter$X14.OZ_007069 <- factor(geneCopiesFilter$X14.OZ_007069)
geneCopiesFilter$BDEG_27923 <- factor(geneCopiesFilter$BDEG_27923)
geneCopiesFilter$BDEG_24855 <- factor(geneCopiesFilter$BDEG_24855)
geneCopiesFilter$X06.OZ_005347 <- factor(geneCopiesFilter$X06.OZ_005347)
geneCopiesFilter$X18.OZ_001716 <- factor(geneCopiesFilter$X18.OZ_001716)
geneCopiesFilter$BDEG_24892 <- factor(geneCopiesFilter$BDEG_24892)
geneCopiesFilter$X0711.1_005485 <- factor(geneCopiesFilter$X0711.1_005485)

geneCopiesFilter$X18.OZ_005506 <- factor(geneCopiesFilter$X18.OZ_005506)
geneCopiesFilter$BDEG_25147 <- factor(geneCopiesFilter$BDEG_25147)
geneCopiesFilter$BDEG_28447 <- factor(geneCopiesFilter$BDEG_28447)
geneCopiesFilter$BDEG_25168 <- factor(geneCopiesFilter$BDEG_25168)
geneCopiesFilter$X22.OZ_002826 <- factor(geneCopiesFilter$X22.OZ_002826)
geneCopiesFilter$BDEG_20405 <- factor(geneCopiesFilter$BDEG_20405)
geneCopiesFilter$X24.OZ_002846 <- factor(geneCopiesFilter$X24.OZ_002846)
geneCopiesFilter$BDEG_25536 <- factor(geneCopiesFilter$BDEG_25536)
geneCopiesFilter$X25.OZ_002430 <- factor(geneCopiesFilter$X25.OZ_002430)
geneCopiesFilter$BDEG_22285 <- factor(geneCopiesFilter$BDEG_22285)

geneCopiesCount = as.matrix(geneCopiesFilter)

geneCopiesCount

# Define the number of colors you want
nb.cols <- 21
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
# Create a ggplot with 21 colors 
# Use scale_fill_manual
pcirc <- gheatmap(p, geneCopiesCount, width=6, offset=0.5, legend_title="M36 Loci", colnames_position="bottom", colnames_angle=90, colnames_offset_y=0, hjust=1, font.size=3) +   
scale_fill_brewer(direction = 1) +
  scale_y_continuous(expand=c(0, 3)) 
pcirc


ggsave("M36_ggtree2.pdf",pcirc, height=15, width=11)
