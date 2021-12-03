#Loading all the required libraries
BiocManager::install("affy")
BiocManager::install("oligo")
BiocManager::install("Biobase", force =
                       TRUE)
BiocManager::install("GEOquery")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("splitstackshape")
library("GEOquery")
library(oligo)
library(Biobase)
library(affy)
library("splitstackshape")
library("tidyr")
library(dplyr)
library("arrayQualityMetrics")

celFiles <- list.celfiles() #call all the cel files
affyRaw <- read.celfiles(celFiles) # reading those files
eset <- oligo::rma(affyRaw) #normalization
write.exprs(eset, file = "data1.txt") #writing our data in a file (gene info won't be included)
#annotation
mydata <- read.delim ("GSE19697_family.soft", check.names = FALSE)
abc <- data.frame(mydata)
a <- cSplit(abc, "Gene.Symbol", "//")
write.table(a, file = "impl.txt", sep = "\t", row.names = FALSE, quote = FALSE)
data1 <- read.delim("impl.txt", check.names = FALSE)
data2 <- read.delim("data1.txt", check.names = FALSE)
combined <- left_join(data1, data2, by = "ID")
write.csv(combined, "annotated.csv")

library(limma)
pData(eset)
Groups <- c("Tumor", "Tumor", "Tumor", "Tumor", "Tumor", "Normal", "Normal", "Normal", "Normal", "Normal")
design = model.matrix(~factor(Groups))
design
colnames(design) <- c("Tumor", "TumorvsNormal")
fit <- lmFit(eset, design)
fit <- eBayes(fit)
options (digits = 2)
res <- topTable(fit, number = Inf, adjust.method = "none", coef = 1)
write.table(res, "diff_exp.txt", sep = "\t")
