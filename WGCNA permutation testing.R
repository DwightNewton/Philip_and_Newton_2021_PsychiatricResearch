#Cross tabluation-based preservation measures, analyzed with permutation testing
#Preservation measure is module-wise and asks "how many genes expressed in module A of group X are co-expressed in any modules in group Y?"

library(cluster)
library(scatterplot3d)
library(WGCNA)
library(dplyr)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()


#### Data preparation steps
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

#No obvious outliers here, samples look good, ready for network analysis

#Split data into groups to make within-group networks
conData <- expData[13:18,]
GFData <- expData[1:6,]
colMonoData <- expData[7:12,]
colSPFData <- expData[19:24,]

load("con_Network_construction.rData")
load("GF_Network_GFstruction.rData")
load("Mono_col_Network_colMonostruction.rData")
load("SPF_col_Network_colSPFstruction.rData")


rownames(conData) <- paste0("Con", 1:6)
rownames(GFData) <- paste0("GF", 1:6)
rownames(colMonoData) <- paste0("colMono", 1:6)
rownames(colSPFData) <- paste0("colSPF", 1:6)

conData <- as.matrix(conData)
GFData <- as.matrix(GFData)
colMonoData <- as.matrix(colMonoData)
colSPFData <- as.matrix(colSPFData)

#remove colMonoData outlier (#6)
colMonoData[6,] <- NA

#remove all grey modules
conMEs <- conMEs[,-ncol(conMEs)]
colSPFMEs <- colSPFMEs[,-ncol(colSPFMEs)]




#### Permutation test steps

# Calculate of # of genes from reference modules co-expressed (ref # - comp #)
Reflist <- c("conData")
Complist <- c("conData","GFData", "colMonoData", "colSPFData")
MElist <- c("conMEs", "GFMEs", "colmonoMEs", "colSPFMEs")
Colorlist <- c("conmoduleColors", "GFmoduleColors", "colMonomoduleColors", "colSPFmoduleColors")


for (i in 1:length(Reflist)){
  for(j in 1:length(Complist)){
    #Get of co-expressed genes for each reference module
    Refmods <- unique(get(Colorlist[i]))[unique(get(Colorlist[i])) != "grey"]
    Compmods <- unique(get(Colorlist[j]))[unique(get(Colorlist[j])) != "grey"]
    
    Refgenelist <- list(NA)
    Compgenelist <- list(NA)
    
    #extract gene lists for all modules
    for (k in 1:length(Refmods)){
      Refmodgenes <- get(Colorlist[i]) %in% Refmods[k]
      Refcoex <- colnames(get(Reflist[i]))[Refmodgenes]
      Refgenelist[k+1] <- list(Refcoex)
    }
    for (k in 1:length(Compmods)){
      Compmodgenes <- get(Colorlist[j]) %in% Compmods[k]
      Compcoex <- colnames(get(Complist[j]))[Compmodgenes]
      Compgenelist[k+1] <- list(Compcoex)
    }
    
    #remove NA placeholder
    Refgenelist <- Refgenelist[-1]
    Compgenelist <- Compgenelist[-1]
    #create template to hold preservation results
    maxpres <- matrix(ncol=length(Refgenelist), nrow=length(Compgenelist))
    colnames(maxpres) <- paste0("Ref_", Refmods)
    rownames(maxpres) <- paste0("Comp_", Compmods)
    
    #full-factorial comparison of module preservation
    for (k in 1:length(Refgenelist)){
      for(l in 1:length(Compgenelist)){
        maxpres[l,k] <- sum(Refgenelist[[k]] %in% Compgenelist[[l]])
      }
    }
    
    RealDiff <- rep(NA, length(colnames(maxpres)))
    names(RealDiff) <- colnames(maxpres)
    RealConservation <- rep(NA, length(colnames(maxpres)))
    names(RealConservation) <- colnames(maxpres)
    
    #Calculate maximum conservation for each reference module, and record difference from maximum (ref module size)
    for(k in 1:length(Refgenelist)){
      maxpresresult <- max(maxpres[,k])
      RealDiff[k] <- length(Refgenelist[[k]]) - maxpresresult
      RealConservation[k] <- maxpresresult
    }
    
    #Permutation test : make sure colnames are identical, and in identical order
    RefData <- get(Reflist[i])
    CompData <- get(Complist[j])
    MergedData <- rbind(RefData, CompData)
    
    #Create template for permuted results (each variable has "n" slots (n=# of reference modules), which can hold 1000 permuted values)
    no.perms <- 10000
    FakeDiff <- rep(list(rep(NA, no.perms)), length(Refmods))
    FakeConservation <- rep(list(rep(NA, no.perms)), length(Refmods))
    
    for(p in 1:no.perms){
      #Randomization step where replicates are  sampled for each group
      temp1 <- sample(rownames(MergedData), length(rownames(MergedData)), replace=F)
      RefI <- temp1[1:nrow(RefData)]
      CompI <- temp1[(nrow(RefData)+1):length(temp1)]
      FakeRefData <- MergedData[rownames(MergedData) %in% RefI, ]
      FakeCompData <- MergedData[rownames(MergedData) %in% CompI, ]
      
      #WGCNA step
      
      FakeRefNet <- blockwiseModules(FakeRefData, power = 5,
                                     networkType = "unsigned", corType="bicor",
                                     TOMType = "unsigned", minModuleSize = 3,
                                     reassignThreshold = 0, mergeCutHeight = 0.15,
                                     numericLabels = TRUE, pamRespectsDendro = FALSE,
                                     saveTOMs = FALSE, verbose = 3, maxBlockSize=30000)
      
      FakeCompNet <- blockwiseModules(FakeCompData, power = 5,
                                      networkType = "unsigned", corType="bicor",
                                      TOMType = "unsigned", minModuleSize = 3,
                                      reassignThreshold = 0, mergeCutHeight = 0.15,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      saveTOMs = FALSE, verbose = 3, maxBlockSize=30000)
      
      
      
      FakeRefLabels <- labels2colors(FakeRefNet$colors)
      FakeCompLabels <- labels2colors(FakeCompNet$colors)
      FakeRefmods <- unique(FakeRefLabels)[unique(FakeRefLabels) != "grey"]
      FakeCompmods <- unique(FakeCompLabels)[unique(FakeCompLabels) != "grey"]
      
      FakeRefgenelist <- list(NA)
      FakeCompgenelist <- list(NA)
      
      #extract gene lists for all modules
      for (k in 1:length(FakeRefmods)){
        FakeRefmodgenes <- FakeRefLabels %in% FakeRefmods[k]
        FakeRefcoex <- colnames(FakeRefData)[FakeRefmodgenes]
        FakeRefgenelist[k+1] <- list(FakeRefcoex)
      }
      for (k in 1:length(FakeCompmods)){
        FakeCompmodgenes <- FakeCompLabels %in% FakeCompmods[k]
        FakeCompcoex <- colnames(FakeCompData)[FakeCompmodgenes]
        FakeCompgenelist[k+1] <- list(FakeCompcoex)
      }
      
      #remove NA placeholder
      FakeRefgenelist <- FakeRefgenelist[-1]
      FakeCompgenelist <- FakeCompgenelist[-1]
      
      #create template to hold preservation results
      FakeRefmaxpres <- matrix(ncol=length(Refgenelist), nrow=length(FakeRefgenelist))
      colnames(FakeRefmaxpres) <- paste0("Ref_", Refmods)
      rownames(FakeRefmaxpres) <- paste0("FakeRef_", FakeRefmods)
      #full-factorial Comparison of module preservation in Ref
      for (k in 1:length(Refgenelist)){
        for(l in 1:length(FakeRefgenelist)){
          FakeRefmaxpres[l,k] <- sum(Refgenelist[[k]] %in% FakeRefgenelist[[l]])
        }
      }
      
      #create template to hold preservation results
      FakeCompmaxpres <- matrix(ncol=length(Refgenelist), nrow=length(FakeCompgenelist))
      colnames(FakeCompmaxpres) <- paste0("Ref_", Refmods)
      rownames(FakeCompmaxpres) <- paste0("FakeComp_", FakeCompmods)
      #full-factorial Comparison of module preservation in Comp
      for (k in 1:length(Refgenelist)){
        for(l in 1:length(FakeCompgenelist)){
          FakeCompmaxpres[l,k] <- sum(Refgenelist[[k]] %in% FakeCompgenelist[[l]])
        }
      }
      
      
      #Calculate maximum conservation for each Reference module, and record difference from maximum (Ref preservation amount)
      for(k in 1:length(Refgenelist)){
        FakeRefmaxpresresult <- max(FakeRefmaxpres[,k])
        FakeCompmaxpresresult <- max(FakeCompmaxpres[,k])
        FakeDiff[[k]][p] <- FakeRefmaxpresresult - FakeCompmaxpresresult
        FakeConservation[[k]][p] <- FakeCompmaxpresresult
      }
    }
    DiffSum <- rep(NA, length(Refmods))
    for (k in 1:length(Refmods)){
      DiffSum[k] <- sum(RealDiff[k] <= FakeDiff[[k]])
    }
    
    ConservationSum <-rep(NA, length(Refmods))
    for (k in 1:length(Refmods)){
      ConservationSum[k] <- sum(RealConservation[k] <= FakeConservation[[k]])
    }
    
    p_Diff <- (DiffSum+1)/(no.perms+1)
    p_conservation <- (ConservationSum+1)/(no.perms+1)
    
    p_values <- matrix(ncol=2, nrow=length(Refmods))
    colnames(p_values) <- c("Difference", "Conservation")
    rownames(p_values) <- Refmods
    for(z in 1:length(Refmods)){
      p_values[z,1] <- paste0("p-value for difference of", Refmods[z], "= ", p_Diff[z])
      p_values[z,2] <- paste0("p-value for conservation of", Refmods[z], "= ", p_conservation[z])
    }
    #Save data for later consolidation and analysis
    write.csv(p_values, file=paste0("p-value chart ",Reflist[i], "(Reference) v.s. ", Complist[j], "(comparison).csv"))
    write.csv(maxpres, file=paste0("Preservation chart ",Reflist[i], "(Reference) v.s. ", Complist[j], "(comparison).csv"))
  }
}
