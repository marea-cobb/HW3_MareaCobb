#Loads bioconductor packages
# source("http://bioconductor.org/biocLite.R")
library(Biobase)
library(limma)
library(beadarray)
library(pheatmap)

#' Creates the R markdown files.
opts_knit$set(progress = FALSE, verbose = FALSE, message=FALSE)
spin(hair = "CobbMareaHW3.R", format = "Rmd")
file.rename("CobbMareaHW3.md", "CobbMareaHW3.Rmd")

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
  keepCols <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "type", 
                "channel_count", "source_name_ch1", "organism_ch1", "molecule_ch1", "extract_protocol_ch1",
                "label_ch1", "label_protocol_ch1", "taxid_ch1", "hyb_protocol", "scan_protocol", "description",
                "data_processing", "platform_id", "contact_name", "contact_email", "contact_department",
                "contact_institute", "contact_address", "contact_city", "contact_state", "contact_zip/postal_code",
                "contact_country", "supplementary_file", "data_row_count", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")
  pd <- pd[, keepCols]
  pd$x <- NULL
  pd$characteristics_ch1 <- gsub(".*: ", "", pd$characteristics_ch1)
  pd$characteristics_ch1.1 <- gsub(".*:", "", pd$characteristics_ch1.1)
  pd$characteristics_ch1.2 <- gsub(".*: ", "", pd$characteristics_ch1.2)
  pd
}

pData(gset_new) <- sanitize_pdata(pData(gset_new))

macrophages <- gset_new[, grepl("Monocyte-derived Macrophage", pData(gset_new)$characteristics_ch1.1)]
pbmcs <- gset_new[, grepl("PBMC", pData(gset_new)$characteristics_ch1.1)]

# macrophages$characteristics_ch1.2

mm_macrophages <- model.matrix(~characteristics_ch1.2, macrophages) # design matrix
fit_macrophages <- lmFit(macrophages, mm_macrophages) #Fit linear model for each gene given a series of arrays
ebay_macrophages <- eBayes(fit_macrophages) # compute moderated t-statistics, moderated F-statistic, and 29/
colnames(fit_macrophages$coef)

topTreatment_macrophages <- topTable(ebay_macrophages, coef = "characteristics_ch1.2Poly IC H", number = Inf,
                  sort.by = "none")

selected_macrophages  <- p.adjust(ebay_macrophages$p.value[, 2]) <0.1
esetSel_macrophages <- macrophages [selected_macrophages, ]
heatmap(exprs(esetSel_macrophages))


mm_pbmcs <- model.matrix(~characteristics_ch1.2, pbmcs) # design matrix
fit_pbmcs <- lmFit(pbmcs, mm_pbmcs) #Fit linear model for each gene given a series of arrays
ebay_pbmcs <- eBayes(fit_pbmcs) # compute moderated t-statistics, moderated F-statistic, and 29/

topTreatmentPBMCs <- topTable(ebay_pbmcs, coef = "characteristics_ch1.2Poly IC L", number = Inf,
                         sort.by = "none")

selected_pbmcs  <- p.adjust(ebay_pbmcs$p.value[, 2]) <0.1
esetSel <- pbmcs [selected_pbmcs, ]
heatmap(exprs(esetSel))
