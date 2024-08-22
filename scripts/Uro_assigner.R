#####################################
### Assigner code from Erik Sotka ###
#####################################

### do assignment test based on linear discriminant analysis
### analogous to DAPC but for the output from PCAngsd
library(reshape)
library("circlize")
library(RColorBrewer)
library(MASS) ### lda()
library(pophelper)
library(adegenet)
library(vcfR)
######## SNPs ###########

### 1- input the entire PCA
rm(list=ls())

vcf <- vcfR::read.vcfR("Re_analysis_2023/PopGen/Data/final_data.vcf")
pop <- read.table("Re_analysis_2023/PopGen/Data/popmapFiltered_w_coast.txt", header=T, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c("Indiv","STRATA", "Coast")
Coast<-pop[,3]
pop2<-pop[,2]
pop3<-pop[pop$Coast == "ATL",2]

genind <- vcfR::vcfR2genind(vcf)
genind@pop <- pop$STRATA # This line adds the population data. Make sure column header of V2 in pop is named STRATA (or edit here)
genind@pop <- factor(genind@pop, levels = c("DM", "GB", "WH", "CT", "DB", "OY", "BF", "HP",
                                            "FB", "SK", "WP", "HM", "RB", "TO", "CP"))

genl<-vcfR2genlight(vcf)
genl@pop<-pop$STRATA
e <- glPca(genl,center = T, scale = T, nf = 6)

pca1 <- e$scores[,1]; pca2 <- e$scores[,2]; pca3 <- e$scores[,3]; pca4 <- e$scores[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop2,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop2~.,all,subset=(1:183)[Coast=="ATL"])
out <- predict(z,all[(1:183)[Coast=="PAC"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- gsub("-.*","",rownames(out3))
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,10*10),ncol=10)))
tmp <- data.frame(matrix(rep(0,15*5),nrow=5));rownames(tmp) <- c("CP","HM","RB", "TO", "WP");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3<-md3[1:10,1:5]
md3 <- as.matrix(md3)

### 3-make plot
svg('./Re_analysis_2023/Final_plots/assignment-snp.svg',height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- rownames(md3)
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
                    #grid.col=cols.to.use,
                    annotationTrack = "grid", 
                    transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
                    #diffHeight  = -0.04,
                    direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
                    link.arr.length =  0.15)

dev.off()
svg('./Re_analysis_2023/Final_plots/assignment-snps-details_aug2024.svg', width = 10, height = 10)
fig <- chordDiagram(x = md3, directional = 1,
                    direction.type = c("arrows", "diffHeight"),
                    link.arr.type = "big.arrow", diffHeight = -mm_h(2))
dev.off()

#########################################################################################

######################################################
### Re-run with just the outlier SNPs from the RDA ###
######################################################

### 1- input the entire PCA
rm(list=ls())

vcf <- vcfR::read.vcfR("Re_analysis_2023/PopGen/Data/Uro_RDA_SNPs.recode.vcf")
pop <- read.table("Re_analysis_2023/PopGen/Data/popmapFiltered_w_coast.txt", header=T, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c("Indiv","STRATA", "Coast")
Coast<-pop[,3]
pop2<-pop[,2]
pop3<-pop[pop$Coast == "ATL",2]

#genind <- vcfR::vcfR2genind(vcf)
#genind@pop <- pop$STRATA # This line adds the population data. Make sure column header of V2 in pop is named STRATA (or edit here)
#genind@pop <- factor(genind@pop, levels = c("DM", "GB", "WH", "CT", "DB", "OY", "BF", "HP",
#                                            "FB", "SK", "WP", "HM", "RB", "TO", "CP"))

genl<-vcfR2genlight(vcf)
genl@pop<-pop$STRATA
e <- glPca(genl,center = T, scale = T, nf = 6)

pca1 <- e$scores[,1]; pca2 <- e$scores[,2]; pca3 <- e$scores[,3]; pca4 <- e$scores[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop2,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop2~.,all,subset=(1:183)[Coast=="ATL"])
out <- predict(z,all[(1:183)[Coast=="PAC"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- gsub("-.*","",rownames(out3))
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,10*10),ncol=10)))
tmp <- data.frame(matrix(rep(0,15*5),nrow=5));rownames(tmp) <- c("CP","HM","RB", "TO", "WP");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3<-md3[1:10,1:5]
md3 <- as.matrix(md3)

### 3-make plot
svg('./Re_analysis_2023/Final_plots/assignment-snp_RDA_only.svg',height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- rownames(md3)
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
                    #grid.col=cols.to.use,
                    annotationTrack = "grid", 
                    transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
                    #diffHeight  = -0.04,
                    direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
                    link.arr.length =  0.15)

dev.off()
svg('./Re_analysis_2023/Final_plots/assignment-snps-details_RDA_only.svg', width = 10, height = 10)
fig <- chordDiagram(x = md3, directional = 1)
dev.off()

#########################################################################################

#####################################################################
### Re-run with just the outlier SNPs from the RDA (Just min SST) ###
#####################################################################

vcf <- vcfR::read.vcfR("Re_analysis_2023/PopGen/Data/Uro_RDA_MinSST.recode.vcf")
pop <- read.table("Re_analysis_2023/PopGen/Data/popmapFiltered_w_coast.txt", header=T, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c("Indiv","STRATA", "Coast")
Coast<-pop[,3]
pop2<-pop[,2]
pop3<-pop[pop$Coast == "ATL",2]

#genind <- vcfR::vcfR2genind(vcf)
#genind@pop <- pop$STRATA # This line adds the population data. Make sure column header of V2 in pop is named STRATA (or edit here)
#genind@pop <- factor(genind@pop, levels = c("DM", "GB", "WH", "CT", "DB", "OY", "BF", "HP",
#                                            "FB", "SK", "WP", "HM", "RB", "TO", "CP"))

genl<-vcfR2genlight(vcf)
genl@pop<-pop$STRATA
e <- glPca(genl,center = T, scale = T, nf = 6)

pca1 <- e$scores[,1]; pca2 <- e$scores[,2]; pca3 <- e$scores[,3]; pca4 <- e$scores[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop2,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop2~.,all,subset=(1:183)[Coast=="ATL"])
out <- predict(z,all[(1:183)[Coast=="PAC"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- gsub("-.*","",rownames(out3))
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,10*10),ncol=10)))
tmp <- data.frame(matrix(rep(0,15*5),nrow=5));rownames(tmp) <- c("CP","HM","RB", "TO", "WP");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3<-md3[1:10,1:5]
md3 <- as.matrix(md3)

### 3-make plot
svg('./Re_analysis_2023/Final_plots/assignment-snp_RDA_MinSST.svg',height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- rownames(md3)
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
                    #grid.col=cols.to.use,
                    annotationTrack = "grid", 
                    transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
                    #diffHeight  = -0.04,
                    direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
                    link.arr.length =  0.15)

dev.off()
svg('./Re_analysis_2023/Final_plots/assignment-snps-details_RDA_MinSST.svg', width = 10, height = 10)
fig <- chordDiagram(x = md3, directional = 1)
dev.off()

#####################################################################
### Re-run with just the outlier SNPs from the RDA (Just Max SST) ###
#####################################################################

vcf <- vcfR::read.vcfR("Re_analysis_2023/PopGen/Data/Uro_RDA_MaxSST.recode.vcf")
pop <- read.table("Re_analysis_2023/PopGen/Data/popmapFiltered_w_coast.txt", header=T, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c("Indiv","STRATA", "Coast")
Coast<-pop[,3]
pop2<-pop[,2]
pop3<-pop[pop$Coast == "ATL",2]

#genind <- vcfR::vcfR2genind(vcf)
#genind@pop <- pop$STRATA # This line adds the population data. Make sure column header of V2 in pop is named STRATA (or edit here)
#genind@pop <- factor(genind@pop, levels = c("DM", "GB", "WH", "CT", "DB", "OY", "BF", "HP",
#                                            "FB", "SK", "WP", "HM", "RB", "TO", "CP"))

genl<-vcfR2genlight(vcf)
genl@pop<-pop$STRATA
e <- glPca(genl,center = T, scale = T, nf = 6)

pca1 <- e$scores[,1]; pca2 <- e$scores[,2]; pca3 <- e$scores[,3]; pca4 <- e$scores[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop2,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop2~.,all,subset=(1:183)[Coast=="ATL"])
out <- predict(z,all[(1:183)[Coast=="PAC"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- gsub("-.*","",rownames(out3))
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,10*10),ncol=10)))
tmp <- data.frame(matrix(rep(0,15*5),nrow=5));rownames(tmp) <- c("CP","HM","RB", "TO", "WP");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3<-md3[1:10,1:5]
md3 <- as.matrix(md3)

### 3-make plot
svg('./Re_analysis_2023/Final_plots/assignment-snp_RDA_MaxSST.svg',height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- rownames(md3)
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
                    #grid.col=cols.to.use,
                    annotationTrack = "grid", 
                    transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
                    #diffHeight  = -0.04,
                    direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
                    link.arr.length =  0.15)

dev.off()
svg('./Re_analysis_2023/Final_plots/assignment-snps-details_RDA_MaxSST.svg', width = 10, height = 10)
fig <- chordDiagram(x = md3, directional = 1)
dev.off()

######################################################################
### Re-run with just the outlier SNPs from the RDA (Just Mean SST) ###
######################################################################

vcf <- vcfR::read.vcfR("Re_analysis_2023/PopGen/Data/Uro_RDA_MeanSST.recode.vcf")
pop <- read.table("Re_analysis_2023/PopGen/Data/popmapFiltered_w_coast.txt", header=T, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c("Indiv","STRATA", "Coast")
Coast<-pop[,3]
pop2<-pop[,2]
pop3<-pop[pop$Coast == "ATL",2]


genl<-vcfR2genlight(vcf)
genl@pop<-pop$STRATA
e <- glPca(genl,center = T, scale = T, nf = 6)

pca1 <- e$scores[,1]; pca2 <- e$scores[,2]; pca3 <- e$scores[,3]; pca4 <- e$scores[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop2,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop2~.,all,subset=(1:183)[Coast=="ATL"])
out <- predict(z,all[(1:183)[Coast=="PAC"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- gsub("-.*","",rownames(out3))
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,10*10),ncol=10)))
tmp <- data.frame(matrix(rep(0,15*5),nrow=5));rownames(tmp) <- c("CP","HM","RB", "TO", "WP");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3<-md3[1:10,1:5]
md3 <- as.matrix(md3)

### 3-make plot
svg('./Re_analysis_2023/Final_plots/assignment-snp_RDA_MeanSST.svg',height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- rownames(md3)
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
                    #grid.col=cols.to.use,
                    annotationTrack = "grid", 
                    transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
                    #diffHeight  = -0.04,
                    direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
                    link.arr.length =  0.15)

dev.off()
svg('./Re_analysis_2023/Final_plots/assignment-snps-details_RDA_MeanSST.svg', width = 10, height = 10)
fig <- chordDiagram(x = md3, directional = 1)
dev.off()

#######################################################################
### Re-run with just the outlier SNPs from the RDA (Just Range SST) ###
#######################################################################
vcf <- vcfR::read.vcfR("Re_analysis_2023/PopGen/Data/Uro_RDA_RangeSST.recode.vcf")
pop <- read.table("Re_analysis_2023/PopGen/Data/popmapFiltered_w_coast.txt", header=T, sep="\t", stringsAsFactors = TRUE)
colnames(pop)<-c("Indiv","STRATA", "Coast")
Coast<-pop[,3]
pop2<-pop[,2]
pop3<-pop[pop$Coast == "ATL",2]


genl<-vcfR2genlight(vcf)
genl@pop<-pop$STRATA
e <- glPca(genl,center = T, scale = T, nf = 6)

pca1 <- e$scores[,1]; pca2 <- e$scores[,2]; pca3 <- e$scores[,3]; pca4 <- e$scores[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop2,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop2~.,all,subset=(1:183)[Coast=="ATL"])
out <- predict(z,all[(1:183)[Coast=="PAC"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- gsub("-.*","",rownames(out3))
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,10*10),ncol=10)))
tmp <- data.frame(matrix(rep(0,15*5),nrow=5));rownames(tmp) <- c("CP","HM","RB", "TO", "WP");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3<-md3[1:10,1:5]
md3 <- as.matrix(md3)

### 3-make plot
svg('./Re_analysis_2023/Final_plots/assignment-snp_RDA_RangeSST.svg',height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- rownames(md3)
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
                    #grid.col=cols.to.use,
                    annotationTrack = "grid", 
                    transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
                    #diffHeight  = -0.04,
                    direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
                    link.arr.length =  0.15)

dev.off()
svg('./Re_analysis_2023/Final_plots/assignment-snps-details_RDA_RangeSST.svg', width = 10, height = 10)
fig <- chordDiagram(x = md3, directional = 1)
dev.off()
