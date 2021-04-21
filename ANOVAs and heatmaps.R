library(dplyr)
library(FactoMineR) 
library(ggplot2)
library(plotly)
library(fastcluster)
library(pROC)
library(rgl)
library(pracma)
library(funtimes)
library(ClusterR)
library(mclust)
options(stringsAsFactors = FALSE)

#Load data, save trait data seperately
mPFCData <- read.csv("mPFC raw expression.csv")
amygData <- read.csv("Amygdala raw expression.csv")
hippoData <- read.csv("Hippocampus raw expression.csv")

groups <- amygData$Group
hippogroups <- groups[-12]

mPFCData <- mPFCData[,-1]
amygData <- amygData[,-(1:2)]
hippoData <- hippoData[-12,-(1:2)]

subjects <- paste0(groups, c(1:6))
hipposubjects <- subjects[-12]

# heat-maps showing log2 fold-change from control expression in each region, only significant changes shown
library(ComplexHeatmap)
library(circlize)

rownames(mPFCData) <- subjects
rownames(hippoData) <- hipposubjects
rownames(amygData) <- subjects


mPFCvec <- colMeans(mPFCData[13:18,])
amygvec <- colMeans(amygData[13:18,])
hippovec <- colMeans(hippoData[12:17,])

mPFCDataEXPORT <- sweep(mPFCData, 2, mPFCvec, `/`)
amygDataEXPORT <- sweep(amygData, 2, amygvec, `/`)
hippoDataEXPORT <- sweep(hippoData, 2, hippovec, `/`)

write.csv(mPFCDataEXPORT, "mPFC_LFC.csv")
write.csv(amygDataEXPORT, "Amyg_LFC.csv")
write.csv(hippoDataEXPORT, "Hippo_LFC.csv")



mPFCDataNorm <- sweep(mPFCData, 2, mPFCvec, `/`)[-grep("cSPF", rownames(mPFCData)),]
amygDataNorm <- sweep(amygData, 2, amygvec, `/`)[-grep("cSPF", rownames(amygData)),]
hippoDataNorm <- sweep(hippoData, 2, hippovec, `/`)[-grep("cSPF", rownames(hippoData)),]

mPFCDataNorm$groups <- groups[-grep("cSPF", groups)]
amygDataNorm$groups <- groups[-grep("cSPF", groups)]
hippoDataNorm$groups <- hippogroups[-grep("cSPF", hippogroups)]

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


hmap_mPFC <- aggregate(mPFCDataNorm[,c(1:(ncol(mPFCDataNorm)-1))],list(mPFCDataNorm$groups), mean)
rownames(hmap_mPFC) <- hmap_mPFC$Group.1
hmap_mPFC <- t(hmap_mPFC[,-1])
colnames(hmap_mPFC) <- c("Germ Free", "Mono-Colonized", "SPF-Colonized")
hmap_mPFC <- log2(hmap_mPFC)
rownames(hmap_mPFC) <- sapply(tolower(rownames(hmap_mPFC)), simpleCap)

hmap_amyg <- aggregate(amygDataNorm[,c(1:(ncol(amygDataNorm)-1))],list(amygDataNorm$groups), mean)
rownames(hmap_amyg) <- hmap_amyg$Group.1
hmap_amyg <- t(hmap_amyg[,-1])
colnames(hmap_amyg) <- c("Germ Free", "Mono-Colonized", "SPF-Colonized")
hmap_amyg <- log2(hmap_amyg)
rownames(hmap_amyg) <- sapply(tolower(rownames(hmap_amyg)), simpleCap)

hmap_hippo <- aggregate(hippoDataNorm[,c(1:(ncol(hippoDataNorm)-1))],list(hippoDataNorm$groups), mean)
rownames(hmap_hippo) <- hmap_hippo$Group.1
hmap_hippo <- t(hmap_hippo[,-1])
colnames(hmap_hippo) <- c("Germ Free", "Mono-Colonized", "SPF-Colonized")
hmap_hippo <- log2(hmap_hippo)
rownames(hmap_hippo) <- sapply(tolower(rownames(hmap_hippo)), simpleCap)


#Only significant changes - 2 tailed
mPFCData$groups <- groups
amygData$groups <- groups
hippoData$groups <- hippogroups
datalist <- c("mPFCData", "amygData", "hippoData")
library(car)

for (i in 1:length(datalist)){
  levenes <- rep(NA, 21)
  workingdata <- get(datalist[i])
  #levene's test
  for(j in c(1:21)){
    print(colnames(workingdata)[j])
    print(leveneTest(workingdata[,j], group = workingdata[,ncol(workingdata)]))
    workingresults <- leveneTest(workingdata[,j], group = workingdata[,ncol(workingdata)])
    levenes[j] <- workingresults[1,3]
  }
  
  
  ttests <- matrix(ncol=4, nrow=21)
  names(ttests) <- c("Gene", "CON vs GF", "CON vs monoCol", "CON vs SPFcol")
  #t.tests
  for(j in c(1:21)){
    print(colnames(workingdata)[j])
    ttests[j,1] <- colnames(workingdata[j])
    
    
    ttests[j,2] <- t.test(workingdata[grep("cSPF", workingdata[,ncol(workingdata)]),j],
                          workingdata[grep("GF", workingdata[,ncol(workingdata)]),j],
                          alternative = "two.sided", 
                          var.equal = ifelse(levenes[j]<0.05, FALSE, TRUE))$p.value
    
    ttests[j,3] <- t.test(workingdata[grep("cSPF", workingdata[,ncol(workingdata)]),j],
                          workingdata[grep("Mono-colonized", workingdata[,ncol(workingdata)]),j],
                          alternative = "two.sided", 
                          var.equal = ifelse(levenes[j]<0.05, FALSE, TRUE))$p.value
    
    ttests[j,4] <- t.test(workingdata[grep("cSPF", workingdata[,ncol(workingdata)]),j],
                          workingdata[grep("SPF-colonized", workingdata[,ncol(workingdata)]),j],
                          alternative = "two.sided", 
                          var.equal = ifelse(levenes[j]<0.05, FALSE, TRUE))$p.value
  }
  #correcting for multiple comparisons using BH method
  rawpval <- as.vector(c(ttests[,2], ttests[,3], ttests[,4]))
  fdrcorrected <- p.adjust(ttests[,c(2:4)], method = "fdr")
  ttestsFDR <- data.frame("Gene" = ttests[,1], "CON vs GF" = fdrcorrected[1:21], "CON vs monoCol" = fdrcorrected[22:42], "CON vs SPFcol" = fdrcorrected[43:63])
  write.csv(ttests, file=paste(substr(datalist[i], 1, (nchar(datalist[i])-4)), "ttests (uncorrected).csv"), col.names = TRUE, row.names = FALSE)
  write.csv(ttestsFDR, file=paste(substr(datalist[i], 1, (nchar(datalist[i])-4)), "ttests (FDR corrected).csv"), col.names = TRUE, row.names = FALSE)
}



#Only significant changes - 2 tailed, 10% FDR
amygFDR <- read.csv("amyg ttests (FDR corrected).csv")
mPFCFDR <- read.csv("mPFC ttests (FDR corrected).csv")
hippoFDR <- read.csv("hippo ttests (FDR corrected).csv")

hmap_amygFDR <- hmap_amyg
hmap_mPFCFDR <- hmap_mPFC
hmap_hippoFDR <- hmap_hippo

#check that they line up - all good
rownames(hmap_amygFDR) == rownames(hmap_amyg)
rownames(hmap_mPFCFDR) == rownames(hmap_mPFC)
rownames(hmap_hippoFDR) == rownames(hmap_hippo)

for(i in 1:nrow(hmap_amygFDR)){
  for(j in 1:ncol(hmap_amygFDR)){
    #to change FDR later (e.g. 5% or 25% - change here)
    if(amygFDR[i,j+1] >=0.10){hmap_amygFDR[i,j] <- 0}
  }
}

for(i in 1:nrow(hmap_mPFCFDR)){
  for(j in 1:ncol(hmap_mPFCFDR)){
    #to change FDR later (e.g. 5% or 25% - change here)
    if(mPFCFDR[i,j+1] >=0.10){hmap_mPFCFDR[i,j] <- 0}
  }
}

for(i in 1:nrow(hmap_hippoFDR)){
  for(j in 1:ncol(hmap_hippoFDR)){
    #to change FDR later (e.g. 5% or 25% - change here)
    if(hippoFDR[i,j+1] >=0.10){hmap_hippoFDR[i,j] <- 0}
  }
}

#Generate heatmaps, touch-ups etc. done in Inkscape 
Heatmap(hmap_mPFCFDR, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        name="LFC",
        row_names_side = "left",
        row_names_gp = gpar(font=3, cex=1),
        rect_gp = gpar(col = "black", lty = 1, lwd =1),
        column_names_side = "top",
        column_title = "mPFC Significant")

Heatmap(hmap_amygFDR, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col=colorRamp2(c(-2, 0, 2),  c("dodgerblue4", "white", "firebrick3")),
        name="LFC",
        row_names_side = "left",
        row_names_gp = gpar(font=3, cex=1),
        rect_gp = gpar(col = "black", lty = 1, lwd =1),
        column_names_side = "top",
        column_title = "Amygdala Significant")


Heatmap(hmap_hippoFDR, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col=colorRamp2(c(-1, 0, 1),  c("dodgerblue4", "white", "firebrick3")),
        name="LFC",
        row_names_side = "left",
        row_names_gp = gpar(font=3, cex=1),
        rect_gp = gpar(col = "black", lty = 1, lwd =1),
        column_names_side = "top",
        column_title = "Hippocampus Significant")


#ggfacet plot showing boxplots of gene expression in same order as heatmap
library(ggplot2)
library(reshape2)
library(grid)
library(dplyr)

setwd("C:/Users/dwigh/Dropbox/PhD/Side Projects/Vivek - WGCNA")
mPFCdata <- read.csv("mPFC raw expression.csv")
hippodata <- read.csv("Hippocampus raw expression.csv")
amygdata <- read.csv("Amygdala raw expression.csv")


#Reorganizing data for merger
amygdata$ID == mPFCdata$X
mPFCdata$Group <- amygdata$Group


# amygdata <- amygdata[,c(1,3:23)]
# hippodata <- hippodata[,c(1,3:23)]

mPFCdata$ID <- mPFCdata$X
mPFCdata <- mPFCdata[,!colnames(mPFCdata)=="X"]
mPFCdata <- mPFCdata[,c(22,23,1:21)]


colnames(mPFCdata) == colnames(amygdata)
colnames(amygdata) == colnames(hippodata)

mPFCdata$ID == amygdata$ID
amygdata$ID == hippodata$ID

colnames(mPFCdata)[3:23] <- paste0("mPFC_", colnames(mPFCdata)[3:23])
colnames(amygdata)[3:23] <- paste0("Amygdala_", colnames(amygdata)[3:23])
colnames(hippodata)[3:23] <- paste0("Hippocampus_", colnames(hippodata)[3:23])

mergeddata <- merge(mPFCdata, amygdata, by=c("Group", "ID"))
mergeddata <- merge(mergeddata, hippodata, by=c("Group", "ID"))
setwd("C:/Users/dwigh/Dropbox/Gut-Brain Coexpression paper/Manuscript Documents/JPR submission/Response to reviewers/multiple comparisons/")


mergeddata_melt <- melt(mergeddata, id.vars = c("ID", "Group"))
mergeddata_melt$Region <- gsub("_+.*$", "", mergeddata_melt$variable)
mergeddata_melt$Gene <- gsub("^.*_+", "", mergeddata_melt$variable)

mergeddata_melt$Gene <- factor(mergeddata_melt$Gene, levels = unique(mergeddata_melt$Gene))
mergeddata_melt$Region <- factor(mergeddata_melt$Region, levels = c("Hippocampus", "Amygdala", "mPFC"))

mergeddata_melt$variable <- factor(mergeddata_melt$variable, levels=paste0(rep(c("Hippocampus", 'Amygdala', "mPFC"),length(levels(mergeddata_melt$Gene))), "_", rep(levels(mergeddata_melt$Gene), each=3)))

paste0(rep(c("Hippocampus", 'Amygdala', "mPFC"),length(levels(mergeddata_melt$Gene))), "_", rep(levels(mergeddata_melt$Gene), each=3))


ggplot(mergeddata_melt, aes(x=Group, y=value, fill=Group)) +
  geom_boxplot(width=0.9, position=position_dodge(), lwd=0.25, outlier.size = 0.5) +
  facet_wrap(vars(variable), ncol=3, scales = "free")+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing = unit(0.25, "lines"),  panel.grid = element_blank(), axis.text.y = element_text(size = 6, vjust = c(-0.5,0,0.5))) +
  geom_point(data = . %>%    
               group_by(variable) %>% #group by facet variable
               summarise(y.min = pretty(value)[1],
                         y.max = pretty(value)[length(pretty(value))]) %>%
               tidyr::gather(key, value, -variable), 
             aes(x = 1, y = value),
             inherit.aes = FALSE, alpha = 0) +
  scale_y_continuous(breaks = function(x) seq(from = x[1], 
                                              to = x[2], 
                                              length.out = 3), expand = c(0, 0))