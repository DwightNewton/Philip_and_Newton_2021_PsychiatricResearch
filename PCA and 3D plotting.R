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

#argument "npc" can be changed to change number of components - data is automatically column normalized (Z-noramlized)
#PC1+2 = 55% (hippo), 60% (amyg), 63%(mPFC)
mPFC_PCA <- PCA(mPFCData, ncp=3)
amyg_PCA <- PCA(amygData, ncp=3)
hippo_PCA <- PCA(hippoData, ncp=3) #remove subject 12


mPFC_components <-as.data.frame(mPFC_PCA$ind$coord)
mPFC_components$groups <- groups
amyg_components <- as.data.frame(amyg_PCA$ind$coord)
amyg_components$groups <- groups
hippo_components <- as.data.frame(hippo_PCA$ind$coord)
hippo_components$groups <- hippogroups

subjects <- paste0(groups, c(1:6))
hipposubjects <- subjects[-12]


#plot first 3 components from PCA in 3d scatterplot
library(plotly)
plot_ly(x=hippo_components$Dim.1, y=hippo_components$Dim.2, z=hippo_components$Dim.3, type="scatter3d", mode="markers", color=hippogroups, text= ~paste(rownames(hippo_components)), marker=list(size=12)) %>%layout(
  title="Hippocampus PCA",
  scene=list(
    xaxis = list(backgroundcolor="white",
                 gridcolor="rgb(200,200,200)",
                 showbackground=TRUE,
                 zerolinecolor="black",
                 gridwidth=2,
                 zerolinewidth=2,
                 tickfont=list(size=14),
                 title="1",
                 font=list(size=20)),
    yaxis = list(backgroundcolor="white",
                  gridcolor="rgb(200,200,200)",
                  showbackground=TRUE,
                  zerolinecolor="black",
                  gridwidth=2,
                  zerolinewidth=2,
                  tickfont=list(size=14),
                  title="2",
                  font=list(size=20)),
    zaxis = list(backgroundcolor="white",
                  gridcolor="rgb(200,200,200)",
                  showbackground=TRUE,
                  zerolinecolor="black",
                  gridwidth=2,
                  zerolinewidth=2,
                  tickfont=list(size=14),
                  title="3",
                  font=list(size=20)))
  )

export(hippo_p, file="hippocampus 3d.svg", selenium=RSelenium::rsDriver(browser = "chrome"))
