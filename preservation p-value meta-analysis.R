#Using Stouffer's method of p-value meta analysis to condense multiple module comparisons to a single value, after bonferroni correction

library(metap)
options(stringsAsFactors = FALSE)

NukeAlphabetic <- function(x){
  as.numeric(substr(x, regexpr("=", x) + 2, nchar(x)))
}
setwd("C:/Users/dwigh/Dropbox/PhD/Side Projects/Vivek - WGCNA/Module Preservation analysis/amygdala")
filelist <- list.files(pattern="p-value chart")
###Same for amygdala and hippocampus (check for mPFC)

#set-up each pair of comparisons for p-value meta-analysis
(conVGF <- filelist[c(12)])
(conVcolMono <- filelist[c(9)])
(conVcolSPF <- filelist[c(10)])
(GFVcolMono <- filelist[c(13)])
(GFvcolSPF <- filelist[c(8)])
(colMonoVcolSPF <- filelist[c(5)])


varlist <- c("conVGF", "conVcolMono", "conVcolSPF", "GFVcolMono", "GFvcolSPF", "colMonoVcolSPF")
varnames <- c("Control-GF", "Control-Mono Colonized", "Control-SPF COlonized", "GF-Mono Colonized", "GF-SPF Colonized", "Mono Colonized-SPF Colonized")
workingresults <- as.data.frame(matrix(ncol=2, nrow=1))
names(workingresults) <- c("Contrast", "Stouffer's bon")


for (i in 1:length(varlist)){
  #load each csv, extract p-value vectors
  dir1 <- read.csv(get(varlist[i])[1])
  plist <- as.character(dir1$Difference)
  plist <- sapply(plist, NukeAlphabetic)
  bonlist <- p.adjust(plist, method="bonferroni")
  stoufferbon <- sumz(bonlist)

  resultvector <- c(varnames[i], round(stoufferbon$p,7))
  workingresults <- rbind(workingresults, resultvector)
}
workingresults <- workingresults[-1,]
write.csv(workingresults, file="p-value meta-analysis (single direction).csv", row.names = FALSE)
