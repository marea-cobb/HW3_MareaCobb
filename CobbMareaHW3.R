#Loads bioconductor packages
# source("http://bioconductor.org/biocLite.R")
library(Biobase)
library(limma)
library(beadarray)
library(pheatmap)
#' Creates the R markdown files.
# library("knitr")
# opts_knit$set(progress = FALSE, verbose = FALSE, message=FALSE)
# spin(hair = "CobbMareaHW3.R", format = "Rmd")
# file.rename("CobbMareaHW3.md", "CobbMareaHW3.Rmd")

#Loads GEO libraries and specific dataset
library(GEOmetadb)
library(GEOquery)

if (!file.exists("GEOmetadb.sqlite")) {
  # Download database only if it's not done already
  getSQLiteFile()
}



#Data was normalized before submission to GEO database. 
#Data is in log form. 
#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1002366
gse <- getGEO("GSE40812", GSEMatrix=TRUE, destdir = "Data/GEO")


# Determines how many data sets are present in the GEO. 
if (length(gse) > 1) idx <- grep("GPL10558", attr(gse, "names")) else idx <- 1
gset <- gse[[idx]]

gset_new <- gset

sanitize_pdata <- function(pd) {
  keepCols <- c("title", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")
  pd <- pd[, keepCols]
  
  pd$title <- sapply(strsplit(as.character(pd$title), "_"), "[", 2)
  pd$characteristics_ch1 <- gsub(".*: ", "", pd$characteristics_ch1)
  pd$characteristics_ch1.1 <- gsub(".*:", "", pd$characteristics_ch1.1)
  pd$characteristics_ch1.2 <- gsub(".*: ", "", pd$characteristics_ch1.2)
  pd <- pd[,c(1,2,3,4)]
  names(pd) <- c("ptid", "HCV", "celltype", "treatment")
  pd
}

pData(gset_new) <- sanitize_pdata(pData(gset_new))

macrophages <- gset_new[, grepl("Monocyte-derived Macrophage", pData(gset_new)$celltype)]
pbmcs <- gset_new[, grepl("PBMC", pData(gset_new)$celltype)]

mm_macrophages <- model.matrix(~treatment+ptid, macrophages) # design matrix
fit_macrophages <- lmFit(macrophages, mm_macrophages) #Fit linear model for each gene given a series of arrays
ebay_macrophages <- eBayes(fit_macrophages) # compute moderated t-statistics, moderated F-statistic, and 29/
colnames(fit_macrophages$coef)

topTreatment_macrophages <- topTable(ebay_macrophages, coef = "treatmentPoly IC H", number = Inf,
                  p.value=0.05, lfc=log2(1.5), sort.by = "p")

selected_macrophages  <- p.adjust(ebay_macrophages$p.value[, 2]) <0.1
esetSel_macrophages <- macrophages [selected_macrophages, ]
heatmap(exprs(esetSel_macrophages))

mm_pbmcs <- model.matrix(~characteristics_ch1.2, pbmcs) # design matrix
fit_pbmcs <- lmFit(pbmcs, mm_pbmcs) #Fit linear model for each gene given a series of arrays
ebay_pbmcs <- eBayes(fit_pbmcs) # compute moderated t-statistics, moderated F-statistic, and 29/

topTreatmentPBMCs <- topTable(ebay_pbmcs, coef = "treatmentPoly IC L", number = Inf,
                         p.value=0.05, lfc=log2(1.5), sort.by = "p")

selected_pbmcs  <- p.adjust(ebay_pbmcs$p.value[, 2]) <0.1
esetSel <- pbmcs [selected_pbmcs, ]
heatmap(topTreatmentPBM)
