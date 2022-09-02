### script for preprocessing and normalizing of human HCC data ###

# load needed libraries
library(ArrayExpress)
library(affy)
library(arrayQualityMetrics)
library(hgu133plus2.db)

# assign variables to needed directories
datadir <- "../data/"
resultsdir <- "../results/"

# loading in pdata_human and assigning file names to rownames
pdata_human <- read.delim(file.path(datadir, "pdata_human.csv"),
                    sep=";",
                    stringsAsFactors = FALSE)
rownames(pdata_human) <- pdata_human$CEL


# write vector with file names for construction of affybatch
file.names <- list.files(datadir_human, pattern=".CEL")

# construct affybatch by reading in the data
abatch <- ReadAffy(filenames = file.path(datadir_human, file.names), 
                   sampleNames = file.names)

# attach pdata to affybatch (metadata is now mapped to abatch)
pData(abatch) <- pdata_human[sampleNames(abatch),]

# save the affybatch object
#saveRDS(abatch,
#        file.path(datadir, "/intermediate",
#                  paste("tumor",
#                        class(abatch),
#                        "raw_human.RDS",
#                        sep="_")))
# read in affybatch
#abatch <- readRDS(file.path(datadir, "/intermediate", "tumor_AffyBatch_raw.RDS"))

################################################################################

# do the QC for the affybatch
QC <- arrayQualityMetrics(abatch,
                          do.logtransform = TRUE,
                          force = TRUE,
                          outdir = file.path(resultsdir, "QC_human.abatch"))

### FILE WILL BE VERY BIG - CARE!!! ###
#saveRDS(QC, file.path(datadir, "intermediate", "QC.abatch.RDS"))
#QC <- readRDS(file.path(datadir, "intermediate", "QC.abatch.RDS"))

outliers <- lapply(QC$modules,
                   function(x){
                     x@outliers@which
                   })

# There are no outliers of note in this data set, especially since the spatial
# distribution looks very good

outliers <- unique(unlist(outliers))
goodSamples <- setdiff(c(1:length(abatch)), outliers)
abatch.outlierremoved <- abatch[,goodSamples]

#saveRDS(abatch.outlierremoved, 
#        file.path(datadir, "/intermediate",
#                  paste("tumor", 
#                        class(abatch.outlierremoved), 
#                        "outlierremoved_human.RDS", 
#                        sep = "_")))

rm(QC)

# normalize the Affybatch with RMA normalization
eset <- rma(abatch,
            verbose = TRUE)

# annotate the probes with the function from Thomas Mohr
fdata <- createAnnotation(annotation(eset))
fData(eset) <- fdata[featureNames(eset),]

#saveRDS(eset,
#        file.path(datadir, "/intermediate",
#                        paste("tumor",
#                              class(eset),
#                              "raw_human.RDS",
#                              sep="_")))

# do QC and remove remaining outliers
QC <- arrayQualityMetrics(eset,
                          do.logtransform = FALSE,
                          force = TRUE,
                          outdir = file.path(resultsdir, "QC_human.eset"))


#saveRDS(QC, file.path(datadir, "intermediate", "QC.eset.RDS"))
#QC <- readRDS(file.path(datadir, "intermediate", "QC.eset.RDS"))

outliers <- lapply(QC$modules,
                   function(x){
                     x@outliers@which
                   })


outliers <- unique(unlist(outliers)) # the remaining outliers do not bother us
goodSamples <- setdiff(c(1:ncol(eset)), outliers)
eset.outlierremoved <- eset[,goodSamples]

# save the eset as it is since the outliers do not bother us
#saveRDS(eset, # no outlierremoved here, see above  
#        file.path(datadir, "/intermediate", 
#                  paste("tumor", 
#                        class(eset), 
#                        "outlierremoved_human.RDS",
#                        sep = "_")))

#rm(eset)
#rm(abatch)

################################################################################

# BiocManager::install("hgu133plus2.db") needed for createAnnotation function to work
# createAnnotation creates a dataframe mapping ProbeIDs to Annotations such as ENTREZID, SYMBOL, etc.
# call: annotation <- createAnnotation(columns, annotation)
# columns: columns of the annotation package
# annotation: annotation package as returned by annotation(eset)
# (c) ScienceConsult - DI Thomas Mohr KG (Thomas Mohr) 2016, licensed under GPL

createAnnotation <- function(annotation, columns = c("ACCNUM",
                                                     "ENTREZID", 
                                                     "ENSEMBL",
                                                     "ENSEMBLTRANS",
                                                     "REFSEQ",
                                                     "SYMBOL", 
                                                     "GENENAME",
                                                     "UNIPROT"),
                             keytype = "PROBEID",
                             keys = NULL){
  package <- paste(annotation, "db", sep = ".")
  require(package, character.only = TRUE)
  DBObject <- get(package)
  names(columns) <- columns
  if (is.null(keys)){
    keys <- keys(DBObject, keytype = keytype)
  }
  coerce <- function(x){
    paste(x, collapse = ", ")
  }
  mapping <- lapply(columns, 
                    function(column, DBObject, keys, keytype, multiVals){
                      map <- mapIds(x = DBObject, 
                                    keys = keys, 
                                    column = column, 
                                    keytype = keytype, 
                                    multiVals = multiVals)
                      return(map)
                    },
                    DBObject = DBObject,
                    keys = keys,
                    keytype = keytype,
                    multiVals = function(x){
                      paste(x, collapse = ", ")
                    })
  mapping <- data.frame(PROBEID = keys,
                        mapping,
                        stringsAsFactors = FALSE)
  return(mapping)
}

################################################################################
