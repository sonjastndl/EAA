##analyse permutations

library(dplyr)
library(pegas)
library(ggplot2)
library(raster)
library(ggrepel)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(robust)
#BiocManager::install("qvalue")
library(qvalue)
library(ggpubr)
library(geodata)
library('corrr')
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("FactoMineR")
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
#library(dplyr)
library(tidyverse)
library(UpSetR)


folder <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/permutations_geo_control"
permutationFiles <- list.files("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/permutations_geo_control", pattern = "pvalues.csv$",recursive = TRUE,  full.names = TRUE)

permuted_pvalues <- list()

for (file in permutationFiles) {
  data <- read.csv(file)
  number <- gsub("\\D", "", basename(dirname(file)))
  permuted_pvalues[[number]] <- data$x
}

pv <- data.frame(permuted_pvalues)


#meanPvals <- rowMeans(pv)
#pvm <- as.matrix(pv)
#hist(as.numeric(pvm[1,]))

trueRDA <- readRDS("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/RDA_env.rds")
source("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/src/rdadapt.R")
rdadapt_env<-rdadapt(trueRDA, 2)
rownames(pv) <- rownames(rdadapt_env)

SNPs <- colnames(trueRDA$Ybar)

write.csv(rdadapt_env$p.values, paste0(outdir, "/permut", permut,"/rdadapt_pvalues.csv"), row.names = FALSE)

#MeanPermPval<-rowMeans(pv)
#SortMeanPermPval<-sort(MeanPermPval)

##original
#empty<-c()
#SortObs<-sort(rdadapt_env$p.values)
#LEN<-length(SortObs)
#for (i in seq(1,LEN,1)){
#  empty<-rbind(empty,length(SortMeanPermPval[SortMeanPermPval <= SortObs[i]])/i)
#               }
##### naming repeat
#SO <- as.data.frame(rdadapt_env$p.values)
#colnames(SO) <- "pvalues"
#rownames(SO) <- SNPs
#SO_SNPs <- SO[order(SO$pvalues), , drop = FALSE]
#
#empty2 <- c()
#LEN<-nrow(SO_SNPs)
#for (i in seq(1,LEN,1)){
#  empty2<-rbind(empty2,length(SortMeanPermPval[SortMeanPermPval <= SO_SNPs[i,]])/i)
#}
#
#rownames(empty2) <- rownames(SO_SNPs)
#d <- empty2[empty2 <= 0.05,]
#
# Apply strsplit to each element
#split_vec <- strsplit(names(d), "\\.")
#tsv_SNPs <- as.data.frame(do.call(rbind, split_vec))
#write.table(names(d),"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/PermutationOutliers.csv", col.names = FALSE, row.names = FALSE)
#write.table(tsv_SNPs,"/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/PermutationOutliers.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)
#
#### approach 
#
#long_column <- unlist(pv, use.names = FALSE)
#hist(long_column)
#
#hist(rdadapt_env$p.values)
#
#length(long_column[long_column <= sort(rdadapt_env$p.values)[2] ])/200
#

#outliers <- as.data.frame(names(d))

rownames(pv) <- SNPs
non_random_sites <- pv[pv$below_0_05_count <= 5,]
random_sites <- pv[pv$below_0_05_count > 5,]
rownames(rdadapt_env) <- SNPs
outliers_o <- rdadapt_env[order(rdadapt_env$p.values),][1:500,]
outliers <- outliers_o[rownames(outliers_o) %in% rownames(non_random_sites),]

outliers <- data.frame(Loci = rownames(outliers), p.value = outliers$p.values, contig = unlist(lapply(strsplit(rownames(outliers), split = "\\."), function(x) x[1])))




ann="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/annotationdata/annotations.txt"
annotations <- read.table(ann, na.strings = NaN)
annotations <- unique(annotations)

colnames(annotations) <- c("Chr", "Pos", "Gene")
rownames(annotations) <- paste0(annotations$Chr,".", annotations$Pos)





### PLOT permutation threshold outliers

locus_scores <- scores(trueRDA, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = rownames(locus_scores), locus_scores)
TAB_loci$chromosome <- sub("\\..*", "", TAB_loci$names)
TAB_loci$type <- "All Loci"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "Outlier"
#TAB_loci$type[TAB_loci$names%in% rownames(non_random_sites)  ] <- "Outlier"
#TAB_loci$type <- ifelse(TAB_loci$type == "Outlier", 
#                        paste0("Chr ", TAB_loci$chromosome), 
#                        TAB_loci$type)
TAB_loci$type <- factor(TAB_loci$type, levels = c(unique(TAB_loci$type)))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(trueRDA, choices=c(1,2), display="bp")) # pull the biplot scores

TAB_loci$Gene <- ifelse(
  rownames(TAB_loci) %in% rownames(annotations),
  annotations$Gene[match(rownames(TAB_loci), rownames(annotations))],
  NA
)

TAB_var <- as.data.frame(scores(trueRDA, choices=c(1,2), display="bp")) # pull the biplot scores

TAB_loci$color_group <- ifelse(
  TAB_loci$type == "Outlier",
  TAB_loci$chromosome,
  "All Loci"
)

# Ensure 'color_group' is a factor (for consistent color assignment)
TAB_loci$color_group <- factor(TAB_loci$color_group, levels = unique(TAB_loci$color_group))

# Plot
PLOT1 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(
    data = TAB_loci,
    aes(
      x = RDA1 * 20,
      y = RDA2 * 20,
      colour = color_group
    ),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "All Loci" = "gray90",
      "2L" = "#F9A242FF",
      "3R" = "#44ad61",
      "3L" = "#7761d0",
      "2R" = "#4ab09c",
      "X" = "#967dca",
      "Y" = "#588dcc"
    )
  ) +
  
  geom_text_repel(
    data = TAB_loci[TAB_loci$type == "Outlier", ],  # Only label colored genes
    aes(x = RDA1 * 20, y = RDA2 * 20, label = Gene),
    size = 2.5,
    vjust = 1.5,
    hjust = 0.5,
    family = "Times", 
    max.overlaps = 80
  ) +
  
  geom_segment(
    data = TAB_var,
    aes(xend = RDA1, yend = RDA2, x = 0, y = 0),
    colour = "black",
    size = 0.15,
    linetype = 1,
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  
  geom_text(
    data = TAB_var,
    aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = row.names(TAB_var)),
    size = 2.5,
    family = "Times"
  ) +
  
  xlab("RDA 1") + 
  ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times")


ggsave("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/permuted_RDA_SNPs_T56.png", plot = PLOT1, width = 8, height = 6, dpi = 300)


### PLOT top 500 outliers


outliers_o <- data.frame(Loci = rownames(outliers_o), p.value = outliers_o$p.values, contig = unlist(lapply(strsplit(rownames(outliers_o), split = "\\."), function(x) x[1])))

locus_scores <- scores(trueRDA, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = rownames(locus_scores), locus_scores)
TAB_loci$chromosome <- sub("\\..*", "", TAB_loci$names)
TAB_loci$type <- "All Loci"
TAB_loci$type[TAB_loci$names%in%outliers_o$Loci] <- "Outlier"
#TAB_loci$type[TAB_loci$names%in% rownames(non_random_sites)  ] <- "Outlier"
#TAB_loci$type <- ifelse(TAB_loci$type == "Outlier", 
#                        paste0("Chr ", TAB_loci$chromosome), 
#                        TAB_loci$type)
TAB_loci$type <- factor(TAB_loci$type, levels = c(unique(TAB_loci$type)))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(trueRDA, choices=c(1,2), display="bp")) # pull the biplot scores

TAB_loci$Gene <- ifelse(
  rownames(TAB_loci) %in% rownames(annotations),
  annotations$Gene[match(rownames(TAB_loci), rownames(annotations))],
  NA
)

TAB_var <- as.data.frame(scores(trueRDA, choices=c(1,2), display="bp")) # pull the biplot scores

TAB_loci$color_group <- ifelse(
  TAB_loci$type == "Outlier",
  TAB_loci$chromosome,
  "All Loci"
)

# Ensure 'color_group' is a factor (for consistent color assignment)
TAB_loci$color_group <- factor(TAB_loci$color_group, levels = unique(TAB_loci$color_group))

# Plot
PLOT2 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(
    data = TAB_loci,
    aes(
      x = RDA1 * 20,
      y = RDA2 * 20,
      colour = color_group
    ),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "All Loci" = "gray90",
      "2L" = "#F9A242FF",
      "3R" = "#44ad61",
      "3L" = "#7761d0",
      "2R" = "#4ab09c",
      "X" = "#967dca",
      "Y" = "#588dcc"
    )
  ) +
  
  geom_text_repel(
    data = TAB_loci[rownames(TAB_loci) %in% (outliers_o[order(outliers_o$p.value),][1:100,]$Loci),],  # Only label colored genes
    aes(x = RDA1 * 20, y = RDA2 * 20, label = Gene),
    size = 2.5,
    vjust = 1.5,
    hjust = 0.5,
    family = "Times", 
    max.overlaps = 30
  ) +
  
  geom_segment(
    data = TAB_var,
    aes(xend = RDA1, yend = RDA2, x = 0, y = 0),
    colour = "black",
    size = 0.15,
    linetype = 1,
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  
  geom_text(
    data = TAB_var,
    aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = row.names(TAB_var)),
    size = 2.5,
    family = "Times"
  ) +
  
  xlab("RDA 1") + 
  ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times")

ggsave("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/AssociationAnalysis/permuted_RDA_SNPs_T500.png", plot = PLOT2, width = 8, height = 6, dpi = 300)

