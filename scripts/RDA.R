#install.packages('pegas')
args <- commandArgs(TRUE)
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

outdir=args[1]


dir.create(paste0(outdir, "/ordiR2step"))
dir.create(paste0(outdir, "/environment"))
dir.create(paste0(outdir, "/populationStructure"))
dir.create(paste0(outdir, "/partialRDA/plot"), recursive = TRUE)
dir.create(paste0(outdir, "/partialRDA/anova"), recursive = TRUE)
dir.create(paste0(outdir, "/partialRDA/pRDAobjects"), recursive = TRUE)
dir.create(paste0(outdir, "/partialRDA/Rsquared"), recursive = TRUE)
dir.create(paste0(outdir, "/AssociationAnalysis"), recursive = TRUE)

# READING ENVIRONMENTAL DATA

envdata=args[2]
Env <- read.csv(envdata, header=TRUE)
rownames(Env) <- Env$sample
Dates <- as.Date(Env$Date_num)
custom_origin <- as.Date(min(Env$Date_num))
Env$Date_num <- as.numeric(as.Date(Dates) - custom_origin)
## Dates <- as.Date(Env$Date_num + custom_origin)
UncertaintyParameters <- c("SMcVSMU", "SMpVSMU", "AImPP", "SMpDNF", "SMcDNF", "SMaDNF", "SMaPercSatU", "AI340P", "AI354P", "AImHP", "CH4mrP", "CLgCBHP", "SO2gVCP", "O3gOVCP", "NO2gAMFtcpk", "NO2gAMFtcp", "HCOHgFVCP", "COtcP", "CLgCTPP", "CLgCTHP", "CLgCOTP", "CLgCFP", "CLgCBPP")
Env <- Env[!colnames(Env) %in% UncertaintyParameters]

Env_backup <- Env

# READING ALLELE FREQUENCY DATA

af_file="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome.final_DP15.af"
af <- read.table(af_file, header=TRUE, na.strings = NaN)
ann="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/annotationdata/annotations.txt"
annotations <- read.table(ann, na.strings = NaN)
annotations <- unique(annotations)

colnames(annotations) <- c("Chr", "Pos", "Gene")
rownames(annotations) <- paste0(annotations$Chr,".", annotations$Pos)
loci <- paste(af$Chr, af$Pos, sep=".")
af2 <- as.data.frame(t(af %>% dplyr::select(3:ncol(af))))
colnames(af2) <- loci

AllFreq=af2
rownames(AllFreq) <- gsub("\\.", "-", rownames(AllFreq))


# CLEANING ALLELE FREQUENCY DATA

freq_mean <- colMeans(AllFreq)
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]


# INTERSECTING DATA
AllFreq <- AllFreq[rownames(AllFreq) %in% rownames(Env),]
Env <- Env[rownames(Env) %in% rownames(AllFreq),]
Dates <- Env$Date_num
Lat <- Env$lat
Long <- Env$long
Coordinates <- as.data.frame(cbind(Lat, Long))
rownames(Coordinates) <- Env$sample
colnames(Coordinates) <- c("Lat", "Long")
Env <- as.data.frame(Env[,5:ncol(Env)])
Env$Date_num <- Dates

print("Working with Allele Frequencies Means >=0.95 or <=0.05 ")
print(nrow(AllFreq))
print(ncol(AllFreq))

AFsin <- asin(sqrt(AllFreq))

print("SCALING AND CENTERING")
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')


print("FIRST PCA")
#colnames(remains) <- paste(substr(colnames(remains), 1, 4), 1:ncol(remains), sep = "_")
pr_eur <- prcomp(Env, scale. = TRUE)
Env.PCA <- princomp(Env)
RDAPlot <- fviz_pca_var(pr_eur, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

cumvarPlot <- fviz_eig(pr_eur, addlabels = TRUE) + theme_bw()
Env.PCA2 <- pr_eur

#rownames(Env.PCA2$loadings) <- new_labels(rownames(Env.PCA2$loadings))
#RDAPlot <- fviz_pca_var(Env.PCA2, col.var = "black") + theme_bw() 
ggsave(paste0(outdir,"/environment/CumVar.png"), plot = cumvarPlot, width = 8, height = 6, dpi = 300)
ggsave(paste0(outdir,"/environment/EnvRDA.png"), plot = RDAPlot, width = 8, height = 6, dpi = 300)


#plot_list <- list()

summary_pca <- summary(pr_eur)
write.table(capture.output(summary_pca), file = paste0(outdir,"/environment/SummaryPCA0.txt"), sep = "\t", quote = FALSE)


Env <- as.data.frame(Env)
AllFreq <- AFsin

####  ORDI2STEP: VARIABLE SELECTION
###   2.a) RDA0
print("PERFORMING RDA0")
RDA0 <- rda(AllFreq ~ 1,  Env) 

print("PERFORMING RDA_FULL")
###   2.b) RDA FULL
RDAfull <- rda(AllFreq ~ .,  Env) 


print("PERFORMING ORDIR2STEP")
ModFILE <- paste0(outdir, "/ordiR2step/ordiR2step.rds")

if (file.exists(ModFILE)) {
  # Load the file
  mod <- readRDS(file=ModFILE)
  message("File ordiR2step_mod loaded successfully.")
} else {
  # Generate the file
  mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
  saveRDS(mod, file=ModFILE)
  message("File generated and saved successfully.")
}

####   2.c) Stepwise procedure with ordiR2step function
#mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
##saveRDS(mod, file="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDAResearchPlan/ordiR2step/ordiR2step.rds")
#saveRDS(mod, file=paste0(outdir, "/ordiR2step/ordiR2step.rds"))
#
##save mod object 

###best vars and save
best_variables <- all.vars(formula(mod))
best_variables <- setdiff(best_variables, "AllFreq")

xo <- ordiplot(mod, scaling = 1, type = "text")

# Extract site scores (samples)
site_scores <- as.data.frame(scores(mod, display = "sites", scaling = 2))
site_scores$point_name <- rownames(site_scores)
site_scores$group <- substr(site_scores$point_name, 1, 2)

# Extract biplot scores (arrows)
arrow_scores <- as.data.frame(scores(xo, display = "biplot", scaling = 2))
arrow_scores$var_name <- rownames(arrow_scores)

# Normalize arrow lengths to match ordiplot scaling
arrow_scores <- arrow_scores %>%
  mutate(RDA1 = RDA1 * (max(abs(site_scores$RDA1)) / max(abs(RDA1))),
         RDA2 = RDA2 * (max(abs(site_scores$RDA2)) / max(abs(RDA2))))

# Calculate extended plot limits to include all arrows
xlim <- range(c(site_scores$RDA1, arrow_scores$RDA1))
ylim <- range(c(site_scores$RDA2, arrow_scores$RDA2))

# Slightly extend the limits to avoid cutting off arrows
xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
ylim <- ylim + c(-0.1, 0.1) * diff(ylim)

bg2  <- c(
  "#304c62",  # AT
  "#4575b4",  # BY
  "#91bfdb",  # CH
  "#3e5f32",  # DE
  "#abdda4",  # DK
  "#d7191c",  # ES
  "#fdae61",  # FI
  "#4fbbb6",  # FR
  "#9e0142",  # GB
  "#d73027",  # GR
  "#91bfdb",  # HU
  "#4575b4",  # IT
  "#66c2a5",  # NL
  "#f46d43",  # PL
  "gray",  # PT
  "#fee08b",  # RS
  "#ce75ac",  # RU
  "#bd40c0",  # UA
  "#70355e"  #SE
)

OPLOT <- ggplot() +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color = group), size = 3) +
  geom_segment(data = arrow_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 0.4) +
  geom_text(data = arrow_scores, aes(x = RDA1, y = RDA2, label = var_name),
            hjust = 1.5, vjust = -1.5, size = 4, color = "black") +
  coord_fixed(ratio = 1) +  # Ensures equal scaling for x and y axes
  xlim(xlim) + ylim(ylim) +
  scale_color_manual(values = bg2) +  # Apply custom colors
  theme_bw() +
  labs(title = "OrdiR2Step: Variable Selection",
       color = "Group")

#xo <- ordiplot(mod, scaling=1, type="text")
#
png(filename = paste0(outdir,"/ordiR2step/Ordiplot_Colored.png"), width = 800, height = 600)
OPLOT
dev.off()

png(filename = paste0(outdir,"/ordiR2step/Ordiplot.png"), width = 800, height = 600)
ordiplot(mod, scaling=1, type="text")
dev.off()

write.csv(best_variables, file=paste0(outdir,"/ordiR2step/best_variables.csv"))

EnvComp <- Env[best_variables]
#EnvComp <- Env[,colnames(Env) %in% best_variables]
###   3) Generate vector of dates since earliest collection time point

###   4) Combine selected Variables  from OrdiR2step with lat, long, PCs_neutral and time vector

### PREPARING NEUTRAL SNP AND COORDINATES
neutral_SNPS <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome_NeutralSNPS_80.tsv"

AllFreq_neutral <- read.table(neutral_SNPS, header = F)
rownames(AllFreq_neutral) <- paste(AllFreq_neutral$V1, AllFreq_neutral$V2, sep=".")
#AF2_neutral <- as.data.frame(t(AllFreq_neutral %>% dplyr::select(3:ncol(AllFreq_neutral))))
Loci <- rownames(AllFreq_neutral)
AF2_neutral <- AllFreq[,colnames(AllFreq) %in% Loci]
#rownames(AF2_neutral) <- paste(AllFreq_neutral$V1, AllFreq_neutral$V2, sep=".")
#AllFreq=af10k
AllFreq_neutral=AF2_neutral[gsub("\\.", "-", rownames(AF2_neutral)) %in% rownames(Env),]

### PERFORMING PCA NEUTRAL

FilePCANeutral <- paste0(outdir,"/populationStructure/pcaNeutral.rds")

if (file.exists(FilePCANeutral)) {
  # Load the file
  pca_neutral <- readRDS(file=FilePCANeutral)
  message("File pcaNeutral loaded successfully.")
} else {
  # Generate the file
  pca_neutral <- rda(AllFreq_neutral[,-1], scale=T) # PCA in vegan uses the rda() call without any predictors
  saveRDS(pca_neutral, file=FilePCANeutral)
  message("File generated and saved successfully.")
  
}


png(filename = paste0(outdir,"/populationStructure/Screeplot_PopStructure.png"), width = 800, height = 600)
screeplot(pca_neutral, type="barplot")
dev.off()

#here first three are enogh to descirbe neutral population genetic structure??
PCs_neutral <- scores(pca_neutral, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(AllFreq_neutral), PCs_neutral)


print("INTERSECTING COORDINATES AND ENVCOMP AND NEUTRAL PCS")

Variables <- data.frame(Coordinates, EnvComp, PopStruct$PC1, PopStruct$PC2, PopStruct$PC3)
#Variables <- Variables[Variables$sample %in% rownames(AllFreq),]

##only components without lat, long
VariablesN <- Variables[,3:ncol(Variables)]

##included
Factors <- colnames(EnvComp)

### 5) PARTIAL RDA

### FULL MODEL: pRDAfull

# Create the formula string dynamically
formula_string <- paste("AllFreq ~ Lat + Long + PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3 + ", paste(Factors, collapse = " + "))
print(formula_string)
formula <- as.formula(formula_string)

FilepRDAfull <- paste0(outdir,"/partialRDA/pRDAobjects/pRDAfull.rds")
if (file.exists(FilepRDAfull)) {
  # Load the file
  pRDAfull <- readRDS(file=FilepRDAfull)
  message("File 'pRDAfull.rds' loaded successfully.")
} else {
  # Generate the file
  pRDAfull <- rda(formula, Variables)
  saveRDS(pRDAfull, file=FilepRDAfull)
  message("File generated and saved successfully.")
}

#pRDAfull <- rda(formula, Variables)
RsquareAdj(pRDAfull)
R <- as.data.frame(RsquareAdj(pRDAfull))
write.csv(R, file=paste0(outdir,"/partialRDA/Rsquared/Rsquared_full.csv"))
A <- anova(pRDAfull)
write.csv(A, file=paste0(outdir,"/partialRDA/anova/Anova_full.csv"))

png(filename = paste0(outdir,"/partialRDA/plot/pRDA_full.png"), width = 800, height = 600)
plot(pRDAfull)
dev.off()


### SUBMODEL 1: CLIMATE MODEL

formula_string <-paste("AllFreq ~ ", paste(Factors, collapse = " + "), "+ Condition(Lat + Long + PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3)")
formula <- as.formula(formula_string)

FilepRDAclim <- paste0(outdir,"/partialRDA/pRDAobjects/pRDAclim.rds")
if (file.exists(FilepRDAclim)) {
  # Load the file
  pRDAclim <- readRDS(file=FilepRDAclim)
  message("File loaded successfully.")
} else {
  # Generate the file
  pRDAclim <- rda(formula, Variables)
  saveRDS(pRDAclim, file=FilepRDAclim)
  message("File generated and saved successfully.")
}


R <- as.data.frame(RsquareAdj(pRDAclim))
write.csv(R, file=paste0(outdir,"/partialRDA/Rsquared/Rsquared_clim.csv"))
A <- anova(pRDAclim)
write.csv(A, file=paste0(outdir,"/partialRDA/anova/Anova_clim.csv"))

png(filename = paste0(outdir,"/partialRDA/plot/pRDA_climate.png"), width = 800, height = 600)
plot(pRDAclim)
dev.off()



### SUBMODEL 2: POPULATION STRUCTURE MODEL

formula_string <- paste("AllFreq ~ PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3 + Condition(Lat + Long +", paste(Factors, collapse = " + "), ")")
formula_struct <- as.formula(formula_string)

FileRDA_STR <- paste0(outdir,"/partialRDA/pRDAobjects/pRDAstructure.rds")
if (file.exists(FileRDA_STR)) {
  # Load the file
  pRDAstruct <- readRDS(file=FileRDA_STR)
  message("File loaded successfully.")
} else {
  # Generate the file
  pRDAstruct <- rda(formula_struct, Variables)
  saveRDS(pRDAstruct, FileRDA_STR)
  message("File generated and saved successfully.")
}

png(filename = paste0(outdir,"/partialRDA/plot/pRDA_structure.png"), width = 800, height = 600)
plot(pRDAstruct)
dev.off()


R <- as.data.frame(RsquareAdj(pRDAstruct))
write.csv(R, file=paste0(outdir,"/partialRDA/Rsquared/Rsquared_struct.csv"))
A <- anova(pRDAstruct)
write.csv(A, file=paste0(outdir,"/partialRDA/anova/Anova_struct.csv"))



## SUBMODEL 3: GEOGRAPHY  MODEL

print("Making pRDAgeogr")

formula_string <- paste("AllFreq ~ Lat + Long + Condition(", paste(Factors, collapse = " + "), "+ PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3)")
formula_geog <- as.formula(formula_string)


FileRDA_GEO <- paste0(outdir,"/partialRDA/pRDAobjects/pRDAgeo.rds")
if (file.exists(FileRDA_GEO)) {
  # Load the file
  pRDAgeog <- readRDS(file=FileRDA_GEO)
  message("File loaded successfully.")
} else {
  # Generate the file
  pRDAgeog <- rda(formula_geog, Variables)
  saveRDS(pRDAgeog, FileRDA_GEO)
  message("File generated and saved successfully.")
}

png(filename = paste0(outdir,"/partialRDA/plot/pRDA_geography.png"), width = 800, height = 600)
plot(pRDAgeog)
dev.off()

R <- as.data.frame(RsquareAdj(pRDAgeog))
write.csv(R, file=paste0(outdir,"/partialRDA/Rsquared/Rsquared_geog2.csv"))
A <- anova(pRDAgeog)
write.csv(A, file=paste0(outdir,"/partialRDA/anova/Anova_geog2.csv"))

### 5 - TABLE) Inertia 
print("Ready to take stats for table")
### 5.2 - another variance partitioning
vp <- varpart(AllFreq, EnvComp, PopStruct[,2:4], Coordinates)

### 6 - Correlation Plot

png(filename = paste0(outdir,"/partialRDA/plot/CorrelationEnv.png"), width = 800, height = 600)
corrplot::corrplot(cor(Variables))
dev.off()

### 7 - Testing Environmental Associations 



### 7.1) Make an RDA without Geography to invest deeper the climatic factors

formula_string <- paste("AllFreq ~ ", paste(Factors, collapse = " + "), " + Condition(PopStruct.PC1 +PopStruct.PC2 + PopStruct.PC3)")
formula_env <- as.formula(formula_string)

FileRDA_ENV <- paste0(outdir,"/AssociationAnalysis/RDA_env.rds")
if (file.exists(FileRDA_ENV)) {
  # Load the file
  RDA_env <- readRDS(file=FileRDA_ENV)
  message("File loaded successfully.")
} else {
  # Generate the file
  RDA_env <- rda(formula_env, Variables)
  saveRDS(RDA_env, FileRDA_ENV)
  message("File generated and saved successfully.")
}


### 8 - PLOT EnvRDA
#plot(RDA_env, scaling=3)
png(filename = paste0(outdir,"/AssociationAnalysis/RDA_env.png"), width = 800, height = 600)
plot(RDA_env, scaling=3)
dev.off()


CountryCodes <- substr(rownames(Variables),1,2)
Variables$Country <- CountryCodes
unique_levels <- unique(Variables$Country)
levels(CountryCodes) <- unique_levels



wolf.rda <- RDA_env

plot_rda_custom <- function(rda_obj, bg_colors, country_codes) {
  # Basic empty plot
  plot(rda_obj, type = "n", scaling = 3)
  # Add SNPs (species points)
  points(rda_obj, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
  # Add wolves (site points)
  points(rda_obj, display = "sites", pch = 21, cex = 1.3, col = "gray32", bg = bg_colors, scaling = 3)
  # Add predictor arrows with labels
  text(rda_obj, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
  # Add legend
  legend("bottomright", legend = levels(country_codes), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg_colors)
}

##plotting 1 and 3rd axis is also possible (see https://popgen.nescent.org/2018-03-27_RDA_GEA.html)

png(filename = paste0(outdir, "/AssociationAnalysis/RDA_env_CountryColors.png"), width = 800, height = 600)
plot_rda_custom(RDA_env, bg_colors = bg2, country_codes = CountryCodes)
dev.off()


### 9 IDENTIFY CANDIDATE SNPS INVOLVED IN LOCAL ADAPTATION - APPROACH EXCURSE
Variables <- data.frame(Coordinates, EnvComp, PopStruct$PC1, PopStruct$PC2, PopStruct$PC3)

load.rda <- scores(RDA_env, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
#pvalue <- 0.05
#pvalue <- 0.05/length(rdadapt_env$p.values)
#pvalue <- 0.05/ncol(rda_LFMM$Ybar)
zscore <- qnorm(1 - pvalue / 2)



outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],zscore) 
cand2 <- outliers(load.rda[,2],zscore) 
cand3 <- outliers(load.rda[,3],zscore) 

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=ncol(Variables))  # 8 columns for 8 predictors
colnames(foo) <- colnames(Variables)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- AllFreq[,nam]
  foo[i,] <- apply(Variables,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)


length(cand$snp[duplicated(cand$snp)])  # 62 duplicate detections
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
a=ncol(cand)
b=ncol(cand)+1
c=b+1
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,b] <- names(which.max(abs(bar[4:a]))) # gives the variable
  cand[i,c] <- max(abs(bar[4:a]))              # gives the correlation
}

colnames(cand)[b] <- "predictor"
colnames(cand)[c] <- "correlation"

table(cand$predictor)


sel <- cand$snp
env <- cand$predictor
env[env=="AImP"] <- '#c1ec00'
env[env=="AImH"] <- '#c2ec00'
env[env=="Bio_1"] <- '#e0562b' #red-ish
env[env=="Bio_4"] <- '#e31a1c'
env[env=="Bio_5"] <- '#33a02c'
env[env=="Bio_6"] <- '#66df93'
env[env=="Bio_7"] <- '#62171e'
env[env=="Bio_8"] <- '#df4671'
env[env=="Bio_9"] <- '#e1972d'
env[env=="Bio_10"] <- '#d3b074'
env[env=="Bio_11"] <- '#c46752'
env[env=="Bio_17"] <- '#6c6c27'
env[env=="Bio_14"] <- '#a6cee3'
env[env=="Bio_18"] <- '#6a3d9a'
env[env=="Lat"] <- '#ffff33'
env[env=="Long"] <- '#01968e'
env[env=="PasHayFlu_l"] <- '#009ddb'
env[env=="PasHayFlu"] <- '#009ddb'
env[env=="PasHayMet_l"] <- '#9f7590'
env[env=="PO24d_l"] <- '#7a8835'
env[env=="POGly"] <- '#b765a4'
env[env=="PopStruct.PC1"] <- '#1f78b4'
env[env=="PopStruct.PC2"] <- '#f6b639'
env[env=="PopStruct.PC3"] <- '#ff948e'
env[env=="Date_num"] <- '#faf497'
env[env=="PO24d_l"] <- '#a46752'
env[env=="mERA5snowD"] <- '#3a8835'
env[env=="S2L2AB01"] <- '#5482c3'

col.pred <- rownames(wolf.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("\\.",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
#bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99','#b2df8a')
levels(col.pred) <- names(sort(table(cand$predictor)))
col.pred2 <- unique(names(sort(table(col.pred))))

col.pred2 <- col.pred2[-which(col.pred2 == "#f1eef6")]

png(filename = paste0(outdir,"/AssociationAnalysis/OutlierMaxAssoc3.png"), width = 800, height = 600)
plot(wolf.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(wolf.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(wolf.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c(names(sort(table(cand$predictor)))), bty="n", col="gray32", pch=21, cex=1, pt.bg=col.pred2)
dev.off()

############## Adaptive Association via radapt, something does not fit its is > 50k SNPs... :(

source("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/src/rdadapt.R")
rdadapt_env<-rdadapt(RDA_env, 2)

## P-values threshold after Bonferroni correction (st.dev of 3.5 is 0.0005?)
#thres_env <- 0.0027
thres_env <- 0.05/length(rdadapt_env$p.values)

## Identifying the loci that are below the p-value threshold
#outliers <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

#### outliers <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], split = "\\."), function(x) x[1])))
rownames(rdadapt_env) <- colnames(AllFreq)
outliers_o <- rdadapt_env[order(rdadapt_env$p.values),][1:500,]
outliers <- data.frame(Loci = rownames(outliers), p.value = outliers$p.values, contig = unlist(lapply(strsplit(rownames(outliers), split = "\\."), function(x) x[1])))

####intersect outliers####

set_list <- list(Method_Loading = cand$snp,Method_Rdadapt = outliers$Loci)
upset_data <- fromList(set_list)

#include upset library and list code to create the upset_data
u <- upset(upset_data, sets = c("Method_Loading", "Method_Rdadapt"), sets.bar.color = c("#01985a", "#009ddb"),    # Color for set bars
           main.bar.color = c("#ffaf9f", "#01985a", "#009ddb"),                             # Color for intersection bars
           matrix.color = "#ff4500",                               # Color for dots in the matrix
           order.by = "freq",
           mainbar.y.label = "Intersection Size",
           sets.x.label = "Set Size"
)

png(filename = paste0(outdir, "/AssociationAnalysis/OutlierIntersection.png"), width = 800, height = 600)
u
dev.off()


## Top hit outlier per contig
outliers <- outliers[order(outliers$contig, outliers$p.value),]
## List of outlier names
#outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

## Formatting table for ggplot
#outliers %>%
#  group_by(contig) %>%
#  filter(rank(p.value) <= round(n() * 0.1)) %>% # Selecting top 10% outliers per chromosome
#  ungroup() -> top1
#
#top_1p_chromosome_outliers <- as.character(top1$Loci)
#print("Top 1% outliers:")
#print(nrow(top_1perc_outliers_per_chromosome))
#write.csv(top_1p_chromosome_outliers, paste(outdir,"/_RDA_pvals_outliers_top0.01.csv", sep=""),row.names = FALSE)


## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = rownames(locus_scores), locus_scores)
TAB_loci$chromosome <- sub("\\..*", "", TAB_loci$names)
TAB_loci$type <- "All Loci"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "Outlier"
#TAB_loci$type[TAB_loci$names%in%top_1p_chromosome_outliers] <- "Top1"
TAB_loci$type <- ifelse(TAB_loci$type == "Outlier", #top1
                        paste0("Chr ", TAB_loci$chromosome), 
                        TAB_loci$type)
TAB_loci$type <- factor(TAB_loci$type, levels = c(unique(TAB_loci$type)))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

TAB_loci$Gene <- ifelse(
  rownames(TAB_loci) %in% rownames(annotations),
  annotations$Gene[match(rownames(TAB_loci), rownames(annotations))],
  NA
)


PLOT1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color=gray(0.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color=gray(0.80), size=0.6) +
  
  # Conditionally color points by chromosome if type != "neutral"
  geom_point(
    data = TAB_loci,
    aes(
      x = RDA1 * 20,
      y = RDA2 * 20,
      colour = ifelse(type == "All Loci","All Loci", chromosome)
    ),
    size = 2
  ) +
  
  geom_text_repel(
    data = TAB_loci,
    aes(x = RDA1 * 20, y = RDA2 * 20, label = Gene),
    size = 2.5,
    vjust = 1.5,
    hjust = 0.5,
    family = "Times"
  ) +
  
  scale_color_manual(
    values = c("All Loci" = "gray90", "2L"= "#F9A242FF", "3R"="#44ad61", "3L"="#7761d0", "2R"= "#4ab09c", "X"="#967dca", "Y"="#588dcc")
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

ggsave(paste0(outdir,"/AssociationAnalysis/RDA_SNPs2.png"), plot = PLOT1, width = 8, height = 6, dpi = 300)
## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
#dev.off()

## Manhattan plot
Outliers <- rep("Neutral", length(colnames(AllFreq)))
Outliers[colnames(AllFreq)%in%outliers$Loci] <- "All outliers"
#Outliers[colnames(AllFreq)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))

TAB_manhatan <- data.frame(
  Chr = sapply(strsplit(rownames(locus_scores), split = "\\."), `[`, 1),
  Pos = sapply(strsplit(rownames(locus_scores), split = "\\."), `[`, 2),
  pvalues = rdadapt_env$p.values,
  Outliers = Outliers
)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
rownames(TAB_manhatan) <- rownames(TAB_loci)


TAB_manhatan$Gene <- ifelse(rownames(TAB_manhatan) %in% rownames(annotations),
annotations$Gene[match(rownames(TAB_manhatan), rownames(annotations))],
NA
)


PLOT_manhattan <- 
   ggplot(data = TAB_manhatan) +
   geom_point(aes(x = Pos, y = -log10(pvalues), color = Outliers), size = 2) +
   geom_text(data = TAB_manhatan[TAB_manhatan$Outliers != "Neutral",], aes(x = Pos, y = -log10(pvalues), label = Gene),
             size = 3, hjust = 0.5, vjust = -0.5) +
   xlab("Mbp") +
   ylab(paste("-log10(p-value)")) +
   facet_grid(. ~ Chr, scales = "free_x", space = "fixed") +
   geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(0.80), size = 0.6) +
   scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
   theme_bw() +
   theme(
     panel.border = element_rect(color = "black", fill = NA, size = 0.8),
     strip.background = element_rect(color = "black", fill = "white")
   )
ggsave(paste0(outdir,"/AssociationAnalysis/RDA_Manhattan_SNPs.png"), plot = PLOT_manhattan, width = 8, height = 6, dpi = 300)


##LINEARIZE PLOTS
source("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/LinearizeRDAPlot.R")

# Example usage
linearize_rda(
  mod = mod,
  annotations = annotations,
  thres_env = thres_env,
  breakpoint_file_path = "/media/inter/ssteindl/DEST/DEST2_NHM/scripts/inversions_breakpoints_v5v6.txt",
  output_dir = paste0(outdir, "/AssociationAnalysis"))

#############
Linear_Outliers <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/DataAnalysis_after_RDA_1224/results/EuropePass/LinearRegressions/Merged_pvalues.arcsin.csv")
LOC <- Linear_Outliers[,colnames(Linear_Outliers)%in%paste0(bv, ".pval.arcsin")]
LOC$Chr <- Linear_Outliers$Chr
LOC$Pos <- Linear_Outliers$Pos
LOC$Gene <- Linear_Outliers$Gene
rownames(LOC) <- paste0(LOC$Chr, ".", LOC$Pos)
names(LOC) <- gsub(".pval.arcsin","", names(LOC))

#selected_rows <- LOC[apply(LOC[1:(ncol(LOC)-4)], 1, function(row) sum(row < 0.05, na.rm = TRUE) >= 10), ]

#table(selected_rows$Chr)
#2L  2R  3L  3R   X 
#180 185 158 245  34 


####
#LEA_Outliers <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/DataAnalysis_after_RDA_1224/results/EuropePass/LEA/AImH/AImH_LEA_pvals_outliers.csv")
#rownames(LEA_Outliers) <- paste0(LEA_Outliers$Chr, ".", LEA_Outliers$Pos)

LEA_files <- list.files("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/DataAnalysis_after_RDA_1224/results/EuropePass/LEA/", pattern = "outliers.csv$",recursive = TRUE,  full.names = TRUE)
 


### ad

