### script for working on limma for the master thesis ###

# load needed libraries
library(limma)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(genefilter) # for featurefilter

# assign variables to needed directories
datadir <- "../data/"
resultsdir <- "../results/"
plotsdir <- "../plots/"

# read in created eset from 01 if needed
eset <- readRDS(file.path(datadir, "intermediate", "tumor_ExpressionSet_outlierremoved_human.RDS"))

### writing design matrix

# change character to factor for model.matrix
eset$INDIVIDUAL <- as.factor(eset$INDIVIDUAL)

# the actual way of generating the design matrix
design <- model.matrix(~ 0 + eset$STAGE + eset$INDIVIDUAL) # Why is number 9 missing???

colnames(design) <- c("adjacent_tissue",
                      "carcinoma",
                      "individual10",
                      "individual11",
                      "individual12",
                      "individual13",
                      "individual14",
                      "individual15",
                      "individual16")
rownames(design) <- eset$SAMPLEID

# contrast fit function with contrast matrix
contrast.matrix <- makeContrasts(carcinoma - adjacent_tissue,
                                 individual10,
                                 individual11,
                                 individual12,
                                 individual13,
                                 individual14,
                                 individual15,
                                 individual16,
                                 levels=design)

################################################################################

### gsva for gene set variation analysis

# downloading entries for species human and saving it to genesets 
# necessary for the term enrichment and the gsva
genesets = msigdbr(species = "human")

# large list containing geneset ids and corresponding entrez gene names
gs <- genesets %>%
  dplyr::select(gs_id, entrez_gene) %>%
  # group by geneset id
  tidyr::nest(gg = -"gs_id") %>%
  # convert into a list of dataframes
  deframe %>%
  # convert the data frame into a vector
  purrr::map(deframe)

# dataframe with feature data from the molsigdb geneset 
gsva.feat <- genesets %>% 
  dplyr::select(gs_id, gs_cat, gs_subcat, gs_name, gs_exact_source) %>% 
  distinct() %>% 
  mutate(rowname = gs_id) %>% 
  column_to_rownames()

# filtering eset to get unique EntrezIDs and removing NAs
# changing feature names to EntrezIDs
eset2 <- featureFilter(eset) 
featureNames(eset2) <- fData(eset2)$ENTREZID
# careful use the newest version of the hgu133.dbr package for the annotation

# gsva
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html
gsva.es <- gsva(eset2, gs, verbose = FALSE)
fData(gsva.es) <- gsva.feat[featureNames(gsva.es),]

# save gsva object
#saveRDS(gsva.es,
#        file.path(datadir, "/intermediate",
#                  paste("tumor",
#                        class(gsva.es),
#                        "gsva.es.human.RDS",
#                        sep = "_")))

#gsva.es <- readRDS(file.path(datadir, "intermediate", "tumor_ExpressionSet_gsva.es.RDS"))

# linear model fit on gsva.es
fit_gsva <- gsva.es %>% 
  lmFit(design) %>% 
  contrasts.fit(contrast.matrix) %>% 
  eBayes()

### save results and visualize them
# save fit as RDS file
#saveRDS(fit_gsva,
#        file.path(datadir, "/intermediate",
#                  paste("tumor",
#                        class(fit_gsva),
#                        "fit_gsva_human.RDS",
#                        sep = "_")))

#results_gsva <- decideTests(fit2_gsva, method = "global", p.value = 0.05, lfc = 2)
#summary(results)

# writing all names of coeffients from the contrast matrix in a new variable
coefs <- colnames(contrast.matrix) %>% 
  set_names(.,.)

# toptable function for all genes of all variables
res.gsva <- coefs %>% 
  purrr::map(function(x) {
    topTable(fit_gsva, coef = x, n = Inf)
  })

# save res as RDS file
#saveRDS(res.gsva,
#        file.path(datadir, "/intermediate",
#                  paste("tumor",
#                        class(res.gsva),
#                        "res_gsva_human.RDS",
#                        sep = "_")))

################################################################################

### differential gene and pathway expression analysis

# linear model fit on eset 2 because of filtered results
fit_eset <- eset2 %>% 
  lmFit(design) %>% 
  contrasts.fit(contrast.matrix) %>% 
  eBayes()

### save results and visualize them
# save fit as RDS file
#saveRDS(fit_eset,
#        file.path(datadir, "/intermediate",
#                  paste("tumor",
#                        class(fit_eset),
#                        "fit_eset_human.RDS",
#                        sep = "_")))

# toptable function for all genes of all variable
res.eset <- coefs %>%
  purrr::map(function(x) {
    topTable(fit_eset, coef = x, n = Inf)
  })

# save res as RDS file
#saveRDS(res.eset,
#        file.path(datadir, "/intermediate",
#                  paste("tumor",
#                        class(res.eset),
#                        "res_eset_human.RDS",
#                        sep = "_")))

################################################################################
