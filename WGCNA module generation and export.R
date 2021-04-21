library(WGCNA)
library(dplyr)
options(stringsAsFactors = FALSE)
allowWGCNAthreads()

#Load data, save trait data seperately
expData <- read.csv("Hippocampus raw expression.csv")
rownames(expData) <- expData$ID
traitData <- expData[,c(1:2)]
expData <- expData[,-c(1:2)]

#Check quality of data (missing values, etc)
gsg <- goodSamplesGenes(expData, verbose=3)
gsg$allOK
#Above step returned TRUE, so no further steps required

#Cluster samples to determine any outliers that should be removed
sampleTree <- hclust(dist(expData), method="average")
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

#One (perhaps two) outlier samples, C/0731B is the most dissimilar from the others so it was removed

#Split data into groups to make within-group networks
conData <- expData[13:18,]
GFData <- expData[1:6,]
colMonoData <- expData[7:11,]
colSPFData <- expData[19:24,]

conTraits <- traitData[13:18,]
GFTraits <- traitData[1:6,]
colMonoTraits <- traitData[7:11,]
colSPFTraits <- traitData[19:24,]



#Choosing soft-thresholding power, pick first value which reaches 0.9


###Controls
powers <- c(c(1:10), seq(from=12, to=20, by=2))
sft <- pickSoftThreshold(conData, powerVector=powers, verbose=5, networkType="unsigned")

#Scale-free topology plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence - Controls"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

#Mean connectivity plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity - Controls"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



###Germ-Free
sft <- pickSoftThreshold(GFData, powerVector=powers, verbose=5, networkType="unsigned")

#Scale-free topology plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence - Germ-Free"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

#Mean connectivity plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity - Germ-Free"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




###Mono-Colonized
sft <- pickSoftThreshold(colMonoData, powerVector=powers, verbose=5, networkType="unsigned")

#Scale-free topology plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence - Mono-Colonized"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

#Mean connectivity plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity - Mono-Colonized"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




###SPF-Colonized
sft <- pickSoftThreshold(colSPFData, powerVector=powers, verbose=5, networkType="unsigned")

#Scale-free topology plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence - SPF-Colonized"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

#Mean connectivity plot
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity - SPF-Colonized"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#No network passes criteria, will use powers that yield equivalent scale-free topology (~0):
#figured that permissivity was more important than scale-free here  - so changed thresholds to be the same (5) across all groups

#corType "bicor" is the biweight midcorrelation - correlation based on medians, not means, thus less sensitive to outliers (the equivalent of spearman for our purposes)


#Perform WGCNA, using otherwise normal settings

conNet = blockwiseModules(conData, power = 5,
                          networkType = "unsigned", corType="bicor",
                          TOMType = "unsigned", minModuleSize = 3,
                          reassignThreshold = 0, mergeCutHeight = 0.15,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "con_TOM",
                          verbose = 3, maxBlockSize=30000)

GFNet = blockwiseModules(GFData, power = 5,
                         networkType = "unsigned", corType="bicor",
                         TOMType = "unsigned", minModuleSize = 3,
                         reassignThreshold = 0, mergeCutHeight = 0.15,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "GF_TOM",
                         verbose = 3, maxBlockSize=30000)

colMonoNet = blockwiseModules(colMonoData, power = 5,
                              networkType = "unsigned", corType="bicor",
                              TOMType = "unsigned", minModuleSize = 3,
                              reassignThreshold = 0, mergeCutHeight = 0.15,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "colMono_TOM",
                              verbose = 3, maxBlockSize=30000)

colSPFNet = blockwiseModules(colSPFData, power = 5,
                             networkType = "unsigned", corType="bicor",
                             TOMType = "unsigned", minModuleSize = 3,
                             reassignThreshold = 0, mergeCutHeight = 0.15,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             saveTOMFileBase = "colSPF_TOM",
                             verbose = 3, maxBlockSize=30000)

#Get network plots, # of modules, and genes/module
#Controls
table(conNet$colors)
mergedColors <- labels2colors(conNet$colors)
plotDendroAndColors(conNet$dendrograms[[1]], mergedColors[conNet$blockGenes[[1]]], "Control", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
conmoduleLabels <- conNet$colors
conmoduleColors <- labels2colors(conNet$colors)
conMEs <- conNet$MEs
congeneTree <- conNet$dendrograms[[1]]
save(conmoduleLabels, conmoduleColors, conMEs, congeneTree, file= "con_Network_construction.rData")

#Germ-Free
table(GFNet$colors)
mergedColors <- labels2colors(GFNet$colors)
plotDendroAndColors(GFNet$dendrograms[[1]], mergedColors[GFNet$blockGenes[[1]]], "Germ-Free", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
GFmoduleLabels <- GFNet$colors
GFmoduleColors <- labels2colors(GFNet$colors)
GFMEs <- GFNet$MEs
GFgeneTree <- GFNet$dendrograms[[1]]
save(GFmoduleLabels, GFmoduleColors, GFMEs, GFgeneTree, file= "GF_Network_GFstruction.rData")

#Mono-colonized
table(colMonoNet$colors)
mergedColors <- labels2colors(colMonoNet$colors)
plotDendroAndColors(colMonoNet$dendrograms[[1]], mergedColors[colMonoNet$blockGenes[[1]]], "Mono-Colonized", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
colMonomoduleLabels <- colMonoNet$colors
colMonomoduleColors <- labels2colors(colMonoNet$colors)
colMonoMEs <- colMonoNet$MEs
colMonogeneTree <- colMonoNet$dendrograms[[1]]
save(colMonomoduleLabels, colMonomoduleColors, colMonoMEs, colMonogeneTree, file= "Mono_col_Network_colMonostruction.rData")


#SPF-Colonized
table(colSPFNet$colors)
mergedColors <- labels2colors(colSPFNet$colors)
plotDendroAndColors(colSPFNet$dendrograms[[1]], mergedColors[colSPFNet$blockGenes[[1]]], "SPF-Colonized", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
colSPFmoduleLabels <- colSPFNet$colors
colSPFmoduleColors <- labels2colors(colSPFNet$colors)
colSPFMEs <- colSPFNet$MEs
colSPFgeneTree <- colSPFNet$dendrograms[[1]]
save(colSPFmoduleLabels, colSPFmoduleColors, colSPFMEs, colSPFgeneTree, file= "SPF_col_Network_colSPFstruction.rData")


#Export network constructions to Cytoscape
conModules <- unique(conmoduleColors)
conlist <- c(1:length(conModules))
GFModules <- unique(GFmoduleColors)
GFlist <- c(1:length(GFModules))
colMonoModules <- unique(colMonomoduleColors)
colMonolist <- c(1:length(colMonoModules))
colSPFModules <- unique(colSPFmoduleColors)
colSPFModules <- colSPFModules[-1] #remove 1-member grey module, throwing off export
colSPFlist <- c(1:length(colSPFModules))


#Controls
for (index in conlist){
  TOM <- TOMsimilarityFromExpr(conData, power=5, corType = "bicor", networkType="unsigned" , TOMType="unsigned")
  modules <- conModules[index]
  genes <- names(conData)
  inModule <- is.finite(match(conmoduleColors, modules))
  modGenes <- genes[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modGenes, modGenes)
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = paste("Con_","CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("Con_","CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                  weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames = modGenes, nodeAttr = conmoduleColors[inModule])
}

#Germ-Free
for (index in GFlist){
  TOM <- TOMsimilarityFromExpr(GFData, power=5, corType = "bicor", networkType="unsigned")
  modules <- GFModules[index]
  genes <- names(GFData)
  inModule <- is.finite(match(GFmoduleColors, modules))
  modGenes <- genes[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modGenes, modGenes)
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = paste("GF_","CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("GF_","CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                  weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames = modGenes, nodeAttr = GFmoduleColors[inModule])
}


#Mono-colonized
for (index in colMonolist){
  TOM <- TOMsimilarityFromExpr(colMonoData, power=5, corType = "bicor", networkType="unsigned")
  modules <- colMonoModules[index]
  genes <- names(colMonoData)
  inModule <- is.finite(match(colMonomoduleColors, modules))
  modGenes <- genes[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modGenes, modGenes)
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = paste("colMono_","CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("colMono_","CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                  weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames = modGenes, nodeAttr = colMonomoduleColors[inModule])
}

#SPF-colonized
for (index in colSPFlist){
  (TOM <- TOMsimilarityFromExpr(colSPFData, power=5, corType = "bicor", networkType="unsigned"))
  modules <- colSPFModules[index]
  genes <- names(colSPFData)
  inModule <- is.finite(match(colSPFmoduleColors, modules))
  modGenes <- genes[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modGenes, modGenes)
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = paste("colSPF_","CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("colSPF_","CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                  weighted=TRUE, threshold=0.02, nodeNames=modGenes, altNodeNames = modGenes, nodeAttr = colSPFmoduleColors[inModule])
  
}




###Merging edge files, node files, and adding module membership to node files
#These files take the individual network files for each module and combine into one file for each group
#Allows for (1) consistent visualization parameters and (2) easier viewing of all modules in cytoscape

###########################CONTROLS
#load everything
filestoload <- list.files(pattern="Con")
edges <- filestoload[grep("edges", filestoload)]
nodes <- filestoload[grep("nodes", filestoload)]

alledges <- list(read.table(edges[1], header=TRUE), read.table(edges[2], header=TRUE))
allnodes <- list(read.table(nodes[1], header=TRUE), read.table(nodes[2], header=TRUE))

for (i in 1:length(alledges)){
  name <- edges[i]
  name <- gsub(".*-edges-","",name)
  name <- substr(name, 1, nchar(name)-4)
  alledges[[i]]$module <- name
}

conEdges <- rbind(alledges[[1]], alledges[[2]])
write.table(conEdges, file="Control Network Edges - hippocampus.txt", sep="\t", row.names = FALSE)


#Controls
reducedColors <- cbind(conmoduleColors, conmoduleLabels)
reducedColors <- as.data.frame(reducedColors[!duplicated(reducedColors[,1]),])
reducedColors <- reducedColors[order(match(reducedColors$conmoduleColors, reducedColors$conmoduleLabels)), ]
ColorVec <- reducedColors$conmoduleColors
names(conMEs) <- paste0("ME.",ColorVec)

nSamples <- nrow(conMEs)
conmodnames <- names(conMEs)
congeneModuleMembership <- as.data.frame(cor(conData, conMEs, use="p"))
conMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(congeneModuleMembership), nSamples))
names(congeneModuleMembership) <- paste0("MM.", conmodnames)
names(conMMPvalue) <- paste0("p.MM", conmodnames)



for (i in 1:length(allnodes)){
  name <- nodes[i]
  name <- gsub(".*-nodes-","",name)
  name <- substr(name, 1, nchar(name)-4)
  allnodes[[i]]$module <- name
  targetRowNames <- allnodes[[i]]$nodeName
  targetModule <- allnodes[[i]][1,3]
  MMtoinsert <- congeneModuleMembership[rownames(congeneModuleMembership) %in% targetRowNames, grep(targetModule,names(congeneModuleMembership))]
  allnodes[[i]]$MM <- MMtoinsert
}

conNodes <- rbind(allnodes[[1]], allnodes[[2]])
write.table(conNodes,file="Control Network Nodes - hippocampus.txt", sep="\t", row.names = FALSE)



###########################GF
#load everything
filestoload <- list.files(pattern="GF")
edges <- filestoload[grep("edges", filestoload)]
nodes <- filestoload[grep("nodes", filestoload)]

alledges <- list(read.table(edges[1], header=TRUE), read.table(edges[2], header=TRUE), read.table(edges[3], header=TRUE))
allnodes <- list(read.table(nodes[1], header=TRUE), read.table(nodes[2], header=TRUE), read.table(nodes[3], header=TRUE))

for (i in 1:length(alledges)){
  name <- edges[i]
  name <- gsub(".*-edges-","",name)
  name <- substr(name, 1, nchar(name)-4)
  alledges[[i]]$module <- name
}

GFEdges <- rbind(alledges[[1]], alledges[[2]], alledges[[3]])
write.table(GFEdges, file="GF Network Edges - hippocampus.txt", sep="\t", row.names = FALSE)

reducedColors <- cbind(GFmoduleColors, GFmoduleLabels)
reducedColors <- as.data.frame(reducedColors[!duplicated(reducedColors[,1]),])
reducedColors <- reducedColors[order(match(reducedColors$GFmoduleColors, reducedColors$GFmoduleLabels)), ]
ColorVec <- reducedColors$GFmoduleColors
names(GFMEs) <- paste0("ME.",ColorVec)

nSamples <- nrow(GFMEs)
GFmodnames <- names(GFMEs)
GFgeneModuleMembership <- as.data.frame(cor(GFData, GFMEs, use="p"))
GFMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(GFgeneModuleMembership), nSamples))
names(GFgeneModuleMembership) <- paste0("MM.", GFmodnames)
names(GFMMPvalue) <- paste0("p.MM", GFmodnames)



for (i in 1:length(allnodes)){
  name <- nodes[i]
  name <- gsub(".*-nodes-","",name)
  name <- substr(name, 1, nchar(name)-4)
  allnodes[[i]]$module <- name
  targetRowNames <- allnodes[[i]]$nodeName
  targetModule <- allnodes[[i]][1,3]
  MMtoinsert <- GFgeneModuleMembership[rownames(GFgeneModuleMembership) %in% targetRowNames, grep(targetModule,names(GFgeneModuleMembership))]
  allnodes[[i]]$MM <- MMtoinsert
}

GFNodes <- rbind(allnodes[[1]], allnodes[[2]], allnodes[[3]])
write.table(GFNodes,file="GF Network Nodes - hippocampus.txt", sep="\t", row.names = FALSE)




###########################colMono
#load everything
filestoload <- list.files(pattern="colMono")
edges <- filestoload[grep("edges", filestoload)]
nodes <- filestoload[grep("nodes", filestoload)]

alledges <- list(read.table(edges[1], header=TRUE), read.table(edges[2], header=TRUE), read.table(edges[3], header=TRUE), read.table(edges[4], header=TRUE))
allnodes <- list(read.table(nodes[1], header=TRUE), read.table(nodes[2], header=TRUE), read.table(nodes[3], header=TRUE), read.table(nodes[4], header=TRUE))

for (i in 1:length(alledges)){
  name <- edges[i]
  name <- gsub(".*-edges-","",name)
  name <- substr(name, 1, nchar(name)-4)
  alledges[[i]]$module <- name
}

colMonoEdges <- rbind(alledges[[1]], alledges[[2]], alledges[[3]], alledges[[4]])
write.table(colMonoEdges, file="colMonotrol Network Edges - hippocampus.txt", sep="\t", row.names = FALSE)


#colMonotrols
reducedColors <- cbind(colMonomoduleColors, colMonomoduleLabels)
reducedColors <- as.data.frame(reducedColors[!duplicated(reducedColors[,1]),])
reducedColors <- reducedColors[order(match(reducedColors$colMonomoduleColors, reducedColors$colMonomoduleLabels)), ]
ColorVec <- reducedColors$colMonomoduleColors
names(colMonoMEs) <- paste0("ME.",ColorVec)

nSamples <- nrow(colMonoMEs)
colMonomodnames <- names(colMonoMEs)
colMonogeneModuleMembership <- as.data.frame(cor(colMonoData, colMonoMEs, use="p"))
colMonoMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(colMonogeneModuleMembership), nSamples))
names(colMonogeneModuleMembership) <- paste0("MM.", colMonomodnames)
names(colMonoMMPvalue) <- paste0("p.MM", colMonomodnames)



for (i in 1:length(allnodes)){
  name <- nodes[i]
  name <- gsub(".*-nodes-","",name)
  name <- substr(name, 1, nchar(name)-4)
  allnodes[[i]]$module <- name
  targetRowNames <- allnodes[[i]]$nodeName
  targetModule <- allnodes[[i]][1,3]
  MMtoinsert <- colMonogeneModuleMembership[rownames(colMonogeneModuleMembership) %in% targetRowNames, grep(targetModule,names(colMonogeneModuleMembership))]
  allnodes[[i]]$MM <- MMtoinsert
}

colMonoNodes <- rbind(allnodes[[1]], allnodes[[2]], allnodes[[3]], allnodes[[4]])
write.table(colMonoNodes,file="colMonotrol Network Nodes - hippocampus.txt", sep="\t", row.names = FALSE)



###########################colSPF
#load everything
filestoload <- list.files(pattern="colSPF")
edges <- filestoload[grep("edges", filestoload)]
nodes <- filestoload[grep("nodes", filestoload)]

alledges <- list(read.table(edges[1], header=TRUE), read.table(edges[2], header=TRUE), read.table(edges[3], header=TRUE), read.table(edges[4], header=TRUE))
allnodes <- list(read.table(nodes[1], header=TRUE), read.table(nodes[2], header=TRUE), read.table(nodes[3], header=TRUE), read.table(nodes[4], header=TRUE))

for (i in 1:length(alledges)){
  name <- edges[i]
  name <- gsub(".*-edges-","",name)
  name <- substr(name, 1, nchar(name)-4)
  alledges[[i]]$module <- name
}

colSPFEdges <- rbind(alledges[[1]], alledges[[2]], alledges[[3]], alledges[[4]])
write.table(colSPFEdges, file="colSPFtrol Network Edges - hippocampus.txt", sep="\t", row.names = FALSE)


#colSPFtrols
reducedColors <- cbind(colSPFmoduleColors, colSPFmoduleLabels)
reducedColors <- as.data.frame(reducedColors[!duplicated(reducedColors[,1]),])
reducedColors <- reducedColors[order(match(reducedColors$colSPFmoduleColors, reducedColors$colSPFmoduleLabels)), ]
ColorVec <- reducedColors$colSPFmoduleColors
names(colSPFMEs) <- paste0("ME.",ColorVec)

nSamples <- nrow(colSPFMEs)
colSPFmodnames <- names(colSPFMEs)
colSPFgeneModuleMembership <- as.data.frame(cor(colSPFData, colSPFMEs, use="p"))
colSPFMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(colSPFgeneModuleMembership), nSamples))
names(colSPFgeneModuleMembership) <- paste0("MM.", colSPFmodnames)
names(colSPFMMPvalue) <- paste0("p.MM", colSPFmodnames)



for (i in 1:length(allnodes)){
  name <- nodes[i]
  name <- gsub(".*-nodes-","",name)
  name <- substr(name, 1, nchar(name)-4)
  allnodes[[i]]$module <- name
  targetRowNames <- allnodes[[i]]$nodeName
  targetModule <- allnodes[[i]][1,3]
  MMtoinsert <- colSPFgeneModuleMembership[rownames(colSPFgeneModuleMembership) %in% targetRowNames, grep(targetModule,names(colSPFgeneModuleMembership))]
  allnodes[[i]]$MM <- MMtoinsert
}

colSPFNodes <- rbind(allnodes[[1]], allnodes[[2]], allnodes[[3]], allnodes[[4]])
write.table(colSPFNodes,file="colSPFtrol Network Nodes - hippocampus.txt", sep="\t", row.names = FALSE)



