#attempt to use R for PC and PCA
#first attempt has pH and [H] vaules
setwd("C:/Users/horto/Dropbox/CMUBS Lakes/Data analysis/lake stats")
library(permute)
library(vegan)
library(MASS)
library(Hmisc)
library(ggbiplot)
library(nlme)
library(ggplot2)
library(ggrepel)
library(psych)
library(xlsx)
library(grid)
library(ape)
library(ggvegan)
library(phyloseq)
library(DESeq2)
library(gridExtra)
require(ggmap)
require(tidyverse)

memory.size(max = TRUE)
memory.limit(size = 15000)
memory.limit()

###########
### Map ###
###########

#Here create a text file with all the coords needs a column of lat and long in decimal degrees
Sample_coords<- read.csv("CMUBS_coordinates.csv", header=TRUE, row.names = 1)
Sample_coords <- Sample_coords[1:4,1:2]
map1<-get_googlemap(center=c(-85.57,45.68), zoom=11, size=c(420,640), maptype="satellite")
MapCMUBS<-ggmap(map1)+
  geom_point(data=Sample_coords, aes(x=long, y=lat, shape = rownames(Sample_coords), color = rownames(Sample_coords)), size=2, show.legend = TRUE)+
  scale_shape_manual(values = c(21,22,23,24), aes(fill=sample))+
  scale_color_manual(values = c("white", "white", "white", "white"), aes(fill=sample)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(legend.direction = 'vertical', 
        legend.position = 'right')+
  guides(shape=guide_legend(title="Lake"), color = FALSE)

MapCMUBS

ggsave("CMUBS_Map.pdf", MapCMUBS)


###########
### PCA ###
###########

#chem data
chem <- read.table("C:\\Users\\horto\\Dropbox\\CMUBS Lakes\\Data analysis\\lake stats\\CMUBS dean chem data simple headers reduced.txt",header=TRUE, row.names = 1, sep="\t")

#pearson correlation analysis
#check to see if any variables are autocorrelated
rcorr(as.matrix(chem[,4:10]), type="pearson")

#List (>0.7 and sig at 0.01 level)
#NO3 and SAC340 correlate in non-transformed data (r=.92, p < 0.0000), remove SAC340 
#represented by NO3
#chem <- chem[-c(9)]
chem <- chem[-c(9)]

#apply PCA
#chem.pca <- prcomp(chem[,4:9], center = TRUE, scale. = TRUE)
chem.pca <- prcomp(chem[,4:9], center = TRUE, scale. = TRUE)

#print method
#print(chem.pca)
print(chem.pca)

#plot method
#plot(chem.pca, type = "lines")
plot(chem.pca, type = "lines")

#summary method
#summary(chem.pca)
summary(chem.pca)

##create PCA plot in ggplot

PCAchem <- ggbiplot(chem.pca, obs.scale = 1, var.scale = 1, varname.adjust = 1,
                    varname.size = 4) +
 geom_point(size = 3, aes(fill = factor(chem$Habitat),
                          shape = factor(chem$Lake))) +
 scale_shape_manual(values = c(21,22,23,24), aes(fill=Lake)) +
 scale_fill_manual(values = c("black","white"), aes(fill=Habitat)) +
 geom_text_repel(aes(label = chem$Time)) +
 guides(fill = guide_legend(override.aes = list(shape = c(21,21)))) +
 theme_classic() +
 expand_limits(y=2.5)
PCAchem


#Finding point coordinates on biplot
names(chem.pca)
#chem.pca$x

#reset
#rm(list=ls())

####Alpha Diversity mixed-model ANOVAs testing Lake and S/B
REU.alpha <- read.csv("REU_alpha_singletons.csv")

REU.chao.lm=lm(chao~Lake, REU.alpha)
anova(REU.chao.lm)
REU.chao.lme=lme(fixed=chao~Habitat, random=~1|factor(Lake), data=REU.alpha)
anova(REU.chao.lme)

REU.npshannon.lm=lm(npshannon~Lake,REU.alpha)
anova(REU.npshannon.lm)
REU.npshannon.lme=lme(fixed=npshannon~Habitat, random=~1|factor(Lake), data=REU.alpha)
anova(REU.npshannon.lme)

### look for correlations with environmental variables
corr_meta <- as.matrix(CMUBS_META[,c(4:8,10)])
alphacorr <- rcorr(as.matrix(REU.alpha[,4:5]), corr_meta, type = "spearman")
alphacorr
#no significant correlations between alpha and environmental variables


############################################################
############################################################
### DeSeq, VST, and beta div from these transformed data ###
############################################################
############################################################

#############################################################################
##### Phyloseq, but reducing data further than other diversity analyses #####
#############################################################################

CMUBS_NUT <- read.table("CMUBS dean chem data simple headers reduced.txt", header = TRUE, row.names = 1, sep ="\t")
#CMUBS_OTU <- read.table("REULAKES.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.abund.shared", header = TRUE, row.names = 2, sep="\t")
#CMUBS_OTU <- CMUBS_OTU[1:23,3:20375]
CMUBS_TAX <- read.table("TAX.txt", header = TRUE, row.names = 1, sep = "\t")

#dim(CMUBS_OTU)

#CMUBS_OTU_t <- t(CMUBS_OTU)
#write.csv(CMUBS_OTU_t, "CMUBS_OTU_t.csv")
CMUBS_OTU_t <- read.csv("CMUBS_OTU_t.csv", row.names = 1, header = TRUE)

rownames(CMUBS_OTU_t) <- paste0("OTU", 1:nrow(CMUBS_OTU_t))
rownames(CMUBS_OTU_t)

CMUBS_TAX <- as.matrix(CMUBS_TAX, rownames.force = NA)
rownames(CMUBS_TAX) <- paste0("OTU", 1:nrow(CMUBS_TAX))
rownames(CMUBS_TAX)

CMUBS_OTU = otu_table(CMUBS_OTU_t, taxa_are_rows = TRUE)
CMUBS_TAX = tax_table(CMUBS_TAX)

CMUBS_physeq = phyloseq(CMUBS_OTU,CMUBS_TAX)

CMUBS_META = sample_data(CMUBS_NUT)
rownames(CMUBS_META) <- sample_names(CMUBS_physeq)

CMUBS_META = sample_data(CMUBS_META)

# get all the data in a phyloseq instance, or whatever
CMUBS_ALL = phyloseq(CMUBS_OTU,CMUBS_TAX,CMUBS_META)
CMUBS_ALL

### Prepping data to remove rare or erronious OTUs #####

# We need to de-noise the data by plotting the number of reads on a curve and look for the inflection point

at.least.n.in.m <- function(x, n, m){
  all(x[x>0]>=n)&length(x[x>0])>=m
}
counts<- rep(0,10)
for (i in 1:length(counts)){
  rows.to.keep<- apply(otu_table(CMUBS_ALL, taxa_are_rows = TRUE), 1,at.least.n.in.m, n=i, m=2)
  counts[i]<-sum(rows.to.keep)
}

plot(1:10, counts, xlab= 'Min sequences in 2 samples', ylab= 'Number of taxa remaining')

### Inflection point was 2 ###

# Filter taxa that arent seen more than twice in greater than 20% of the data.
#CMUBS_CUT2<-filter_taxa(CMUBS_ALL, function(x) sum(x > 2) > (0.2*length(x)), TRUE)
#CMUBS_CUT2
write.csv(tax_table(CMUBS_ALL), "CMUBS_seqnames_ALL.csv")
write.csv(otu_table(CMUBS_ALL), "CMUBS_seqOTU_ALL.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_CMUBS_ALL <- phyloseq_to_deseq2(CMUBS_ALL, ~ Lake)

gm_mean_ALL = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_CMUBS_ALL = apply(counts(dds_CMUBS_ALL), 1, gm_mean_ALL)
dds_CMUBS_ALL = estimateSizeFactors(dds_CMUBS_ALL, geoMeans=geoMeans_CMUBS_ALL)
dds_CMUBS_ALL = estimateDispersions(dds_CMUBS_ALL)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_CMUBS_ALL <- varianceStabilizingTransformation(dds_CMUBS_ALL, blind=FALSE)
vstMat_CMUBS_ALL <- assay(vst_CMUBS_ALL)
vstMat_CMUBS_ALL[vstMat_CMUBS_ALL<0]<-0
vst.otu.CMUBS_ALL <- otu_table(vstMat_CMUBS_ALL, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.CMUBS_ALL), "vst.otu.CMUBS_ALL.csv")

###Create a distance matrix 
Lake_Dist_ALL <- phyloseq::distance(vst.otu.CMUBS_ALL, method= 'bray')

##NMDS of all samples
NMDS_lake <- ordinate(vst.otu.CMUBS_ALL, method= 'NMDS', distance = Lake_Dist_ALL, formula = NULL)
NMDS_lake
summary(NMDS_lake)
## plot quick NMDS for comparisons
NMDS_lake_plot <- plot_ordination(CMUBS_ALL, NMDS_lake, shape="Lake", label='Time', color = 'Habitat')
NMDS_lake_plot + geom_point(size=1) + theme_bw()

##################################################################
##### Publication-level NMDS figures and envfit correlations #####
##################################################################

##### Lake NMDS
##### envfit depth
CMUBS_META <- read.table("CMUBS dean chem data simple headers reduced.txt", header = T, row.names = 1, sep="\t")
rcorr(as.matrix(CMUBS_META[,4:10]), type="pearson")
# N03 represents SAC340
ef_lake <- envfit(NMDS_lake, CMUBS_META[,c(4:8,10)], permu=999, na.rm = TRUE)
ef_lake
# NH4 and NO3 not signif. below 0.01, have been removed, redo envfit
ef_lake <- envfit(NMDS_lake, CMUBS_META[,c(4:5,8,10)], permu=999, na.rm = TRUE)
ef_lake
# all values remaining sig < 0.001

vectors_lake<-as.data.frame(ef_lake$vectors$arrows*sqrt(ef_lake$vectors$r))
pvals_lake=ef_lake$vectors$pvals
r_lake=ef_lake$vectors$r
envfit_lake<-cbind(vectors_lake,pvals_lake, r_lake)
envfit_lake
#vectors are too long to fit on plot, transform with /2
envfit_lake2 <- cbind(envfit_lake[,1:2]/2, envfit_lake[,3:4])
envfit_lake2

NMDS_lake_pub = data.frame(MDS1 = NMDS_lake$points[,1], MDS2 = NMDS_lake$points[,2])
NMDS_lake_pub <- cbind(NMDS_lake_pub,CMUBS_META[,1:3])
NMDS_lake_pub
glNMDSplotLake <- ggplot(NMDS_lake_pub, aes(x=NMDS_lake_pub$MDS1, y=NMDS_lake_pub$MDS2)) +
  geom_point(aes(fill = factor(NMDS_lake_pub$Habitat),
                 shape = NMDS_lake_pub$Lake),
             size=5) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme(legend.direction="vertical",
        legend.position="right") +
  theme_bw() +
  scale_shape_manual(values = c(21,22,23,24), aes(fill=Lake),
                     labels=c("BL", "FL", "LG", "LM")) +
  scale_fill_manual(name="Habitat",
                    values = c(1,"white"),
                    guide = FALSE) +
  geom_text(label="Stress = 0.115", x = -0.5, y = -0.3, size = 4) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21)))) +
  coord_equal()
glNMDSplotLake <- glNMDSplotLake + geom_text_repel(aes(label = NMDS_lake_pub$Time))
glNMDSplotLake <- glNMDSplotLake + 
  geom_segment(data = envfit_lake2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=envfit_lake2, aes(x=NMDS1, y=NMDS2, label=rownames(envfit_lake2)),size=4)
glNMDSplotLake

#anosim
Lakevector <- as.vector(CMUBS_META[,2])
Timevector <- as.vector(CMUBS_META[,1])
Habitatvector <- as.vector(CMUBS_META[,3])
anosim <- anosim(Lake_Dist_ALL, grouping = Lakevector, permutations=999)
summary(anosim)
plot(anosim)
#significance with Lakes, nothing too interesting with Habitat and Time

###########
### CCA ###
###########

#OTU table must have no rows with rowSums = 0
otu_tab_cca <- read.csv("vst.otu.CMUBS_ALL.csv", header = TRUE, row.names = 1)
otu_tab_0less <- otu_tab_cca[rowSums(otu_tab_cca)!=0,]
otu_tab_0less_t <- t(otu_tab_0less)
otu_tab_0less_t <- as.matrix(otu_tab_0less_t)
CCA_meta <- as.data.frame(CMUBS_META[,c(4:5,8,10)])
#both otu table and metadata must be 'numeric' and same rowlength
mode(CCA_meta)
mode(otu_tab_0less_t)
class(otu_tab_0less_t)
dim(otu_tab_0less)
dim(CCA_meta)
 
CMUBS_CCA <- cca(otu_tab_0less_t ~ DO + Temp + pH + DOC, data = CCA_meta)
CMUBS_CCA
plot(CMUBS_CCA)

anova_CCA <- anova(CMUBS_CCA, alpha = 0.05, beta = 0.01, step=100, perm.max=9999)
anova_CCA
anova_terms <- anova(CMUBS_CCA, by = "terms", permu = 999, cutoff = 0.001)
anova_terms
anova_axes <- anova(CMUBS_CCA, by = "axis", permu = 999, cutoff = 0.001)
anova_axes
CCApermtest <- permutest(CMUBS_CCA, permutations = 999,
                         model = "reduced")
CCApermtest

#extract vectors
CCAvecs <- as.data.frame(cbind(CMUBS_CCA$CCA$biplot[,1], CMUBS_CCA$CCA$biplot[,2]))

CCApub <- data.frame(CC1 = CMUBS_CCA$CCA$u[,1], CC2 = CMUBS_CCA$CCA$u[,2])
CCApub <- cbind(CCApub,CCA_meta,CMUBS_META[,1:3])
CCAplot <- ggplot(CCApub, aes(x=CCApub$CC1, y=CCApub$CC2)) +
  geom_point(aes(fill = CCApub$Habitat,
                 shape = CCApub$Lake),
             size=5) +
  xlab("CCA1 (31.28%)") +
  ylab("CC2(28.56%)") +
  theme(legend.direction="vertical",
        legend.position="right") +
  theme_bw() +
  scale_shape_manual(values = c(21,22,23,24), aes(fill=Lake),
                     labels=c("BL", "FL", "LG", "LM")) +
  scale_fill_manual(name="Habitat",
                    values = c(1,"white")) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21)))) +
  coord_fixed()
CCAplot <- CCAplot + 
  geom_segment(data = CCAvecs, aes(x=0, xend=V1*2, y=0, yend=V2*2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=CCAvecs, aes(x=V1*2, y=V2*2, label=rownames(CCAvecs)),size=4)
CCAplot <- CCAplot + geom_text_repel(aes(label = CCApub$Time))

CCAplot



### pCCA Without DOC (wasn't significant in CCA)
CMUBS_CCA_DO_xDOC <- cca(otu_tab_0less_t ~ DO + Condition(pH) + Condition(Temp), data = CCA_meta)
CMUBS_CCA_DO_xDOC
autoplot(CMUBS_CCA_DO_xDOC)

anova_CCA_DO_xDOC <- anova(CMUBS_CCA_DO_xDOC, alpha = 0.05, beta = 0.01, step=100, perm.max=9999)
anova_CCA_DO_xDOC
anova_terms_DO_xDOC <- anova(CMUBS_CCA_DO_xDOC, by = "terms", permu = 999, cutoff = 0.001)
anova_terms_DO_xDOC
anova_axes_DO_xDOC <- anova(CMUBS_CCA_DO_xDOC, by = "axis", permu = 999, cutoff = 0.001)
anova_axes_DO_xDOC
CCApermtest_DO_xDOC <- permutest(CMUBS_CCA_DO_xDOC, permutations = 999,
                            model = "reduced")
CCApermtest_DO_xDOC


# DO_xDOC pCCA Plot
CCAvecsDO_xDOC <- as.data.frame(cbind(CMUBS_CCA_DO_xDOC$CCA$biplot[,1],0))

CCApubDO_xDOC <- data.frame(CCA1 = CMUBS_CCA_DO_xDOC$CCA$u[,1], CA1 = CMUBS_CCA_DO_xDOC$CA$u[,1])
CCApubDO_xDOC <- cbind(CCApubDO_xDOC,CCA_meta,CMUBS_META[,1:3])
CCAplotDO_xDOC <- ggplot(CCApubDO_xDOC, aes(x=CCApubDO_xDOC$CCA1, y=CCApubDO_xDOC$CA1)) +
  geom_point(aes(fill = CCApubDO_xDOC$Habitat,
                 shape = CCApubDO_xDOC$Lake),
             size=5) +
  xlab("CCA1 (6.29%)") +
  ylab("CA1(6.79%)") +
  theme_bw() +
  scale_shape_manual(values = c(21,22,23,24), aes(fill=Lake),
                     labels=c("BL", "FL", "LG", "LM")) +
  scale_fill_manual(name="Habitat",
                    values = c(1,"white"),
                    guide = FALSE) +
  theme(legend.direction="vertical",
        legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21)))) +
  coord_fixed()
CCAplotDO_xDOC <- CCAplotDO_xDOC + 
  geom_segment(data = CCAvecsDO_xDOC, aes(x=0, xend=V1*2, y=0, yend=V2*2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",
               linetype = "longdash", inherit.aes = FALSE) +
  geom_text(data=CCAvecsDO_xDOC, aes(x=V1*2, y=V2*2, label= "DO",size=4),
            show.legend = FALSE)
CCAplotDO_xDOC <- CCAplotDO_xDOC + geom_text_repel(aes(label = CCApubDO_xDOC$Time))

CCAplotDO_xDOC

#g_legend<-function(a.gplot){
#  tmp <- ggplot_gtable(ggplot_build(a.gplot))
#  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#  legend <- tmp$grobs[[leg]]
#  return(legend)}

#GLconcLegend<-g_legend(CCAplotDO_xDOC)

#GLChemPlot <- grid.arrange(arrangeGrob(CCAplot + theme(legend.position="none"),
#                                       CCAplotDO_xDOC + theme(legend.position="none"),
#                                       nrow=1),
#                           GLconcLegend, nrow=2,heights=c(10, 1))

############################
##### OTU Correlations #####
############################

# Filter taxa that arent seen more than twice in greater than 20% of the data.
CMUBS_CUT5<-filter_taxa(CMUBS_ALL, function(x) sum(x > 2) > (0.238*length(x)), TRUE)
CMUBS_CUT5
write.csv(tax_table(CMUBS_CUT5), "CMUBS_seqnames.csv")
write.csv(otu_table(CMUBS_CUT5), "CMUBS_seqOTU.csv")

### Deseq2 Normalization
### convert the counts to integer mode and group data to normalize together(e.g. region)
dds_CMUBS <- phyloseq_to_deseq2(CMUBS_CUT5, ~ Lake)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans_CMUBS = apply(counts(dds_CMUBS), 1, gm_mean)
dds_CMUBS = estimateSizeFactors(dds_CMUBS, geoMeans=geoMeans_CMUBS)
dds_CMUBS = estimateDispersions(dds_CMUBS)

### variance stabilizing transformation
#Make it a OTU again for phyloseq (no tax or sample data)
vst_CMUBS <- varianceStabilizingTransformation(dds_CMUBS, blind=FALSE)
vstMat_CMUBS <- assay(vst_CMUBS)
vstMat_CMUBS[vstMat_CMUBS<0]<-0
vst.otu.CMUBS <- otu_table(vstMat_CMUBS, taxa_are_rows=TRUE)
write.csv(otu_table(vst.otu.CMUBS), "vst.otu.CMUBS.csv")

#import otu table from computer to return to non-deseq format
otu_tab_corr <- read.csv("vst.otu.CMUBS.csv", header = TRUE, row.names = 1)
#import environ data if it isn't already
CMUBS_META <- read.table("CMUBS dean chem data simple headers reduced.txt", header = T, row.names = 1, sep="\t")
#recall that NO3 is correlated (p < .05, R > 0.7) to SAC340, and is 
#represented by NO3 from hereon out
corr_meta <- as.matrix(CMUBS_META[,c(4:8,10)])

t_otu_tab_corr <- t(otu_tab_corr)

#Spearman's rank correlations
Spearcorr_CMUBS <- rcorr(t_otu_tab_corr, corr_meta[,1:6], type = "spearman")
Spearcorr_CMUBS
write.csv(Spearcorr_CMUBS$r, "Spearcorr_CMUBS_r.csv")
write.csv(Spearcorr_CMUBS$P, "Spearcorr_CMUBS_P.csv")

#edit spreadsheet in excel and import into R
# included in the below table are retained taxa with p < 0.001 correlation
# to Temp, pH, DO, DOC, or a combination of them. If p > 0.001, the value
# for that R is 0

OTU_corr_dat <- read.csv("corr_fig_dat2.csv", header = TRUE)
heatmap_corr <- ggplot(OTU_corr_dat, aes(x = env, y = OTU)) +
                geom_tile(aes(fill = r),colour = "white") +
                scale_fill_gradient2(low = "steelblue", mid = "white", high = "red")

heatmap_corr

# abundance graph
library(ggplot2)
library(scales)
require(reshape2)
require(plyr)

otu_bar_dat <- read.csv("corr_fig_dat4_noDOC.csv", header = TRUE)

#if one chooses to exclude 'unclassified' bacteria
otu_bar_dat <- otu_bar_dat[-c(78:83),]

otu_bar_dat$Specificity <- factor(otu_bar_dat$Specificity, levels=unique(otu_bar_dat$Specificity))

ggplot(otu_bar_dat, aes(Specificity, direction)) +
  geom_col(size = .25, aes(fill = Phylum, col = corrnum)) + 
  facet_grid(env ~ .) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.title.x=element_blank()) +
  guides(col = FALSE) +
  ylab("Abundance") +
  scale_color_gradient(low = "white", high = "black") +
  geom_segment(aes(y = 0, yend = 0, x = 0, xend = 52))
