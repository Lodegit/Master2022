### script for working on enrichment analysis for the master thesis ###

# load needed libraries
library(clusterProfiler)
library(msigdbr)
library(tidyverse)
library(limma)
library(writexl)
library(DOSE)
library(enrichplot)
library(EnhancedVolcano)

# assign variables to needed directories
datadir <- "../data/"
resultsdir <- "../results/"
plotsdir <- "../plots/"

# read in created data from 01 and 02 if needed
res.eset <- readRDS(file.path(datadir, "intermediate", "tumor_list_res_eset_human.RDS")) #???
res.gsva <- readRDS(file.path(datadir, "intermediate", "tumor_list_res_gsva_human.RDS")) 

# useful links
# clusterprofiler
# https://support.bioconductor.org/p/102847/
# https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf

# molsigdbr
# https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html

# saving res of eset in a new variable for easier readability
res <- res.eset %>% purrr::map(data.frame)
#res_gsva <- res.gsva %>% purrr::map(data.frame) 

#-------------------------------------------------------------

#Writing the results in .xlsx files 
# turn the relevant data into one dataframe for ease of use
#data <- res$`carcinoma - adjacent_tissue`

# go through all feature annotation columns and replace the strings that are 
# longer than 32k with too long, so it can be written into excel
#data$UNIPROT[nchar(data$UNIPROT[]) > 32000] <- "TOO_LONG"

# writing the enrich results in .xlsx files (needs to be done for each seperatly)
#write_xlsx(data,
#           file.path(resultsdir, 
#                     paste('human_log.xlsx')))

# writing the gsva results in .xlsx files
#write_xlsx(res_gsva$`carcinoma - adjacent_tissue`,
#           file.path(resultsdir, 
#                     paste('human_gsva.xlsx')))


#-------------------------------------------------------------

# creating list of differentially expressed genes for every contrast
degs <- res %>% purrr::map(function(x) {
  x %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::select(SYMBOL) %>%
    deframe()
})

# downloading entries for species human and saving it to genesets 
# necessary for the term enrichment and the gsva
genesets <- msigdbr(species = "human")

# TERM2NAME for enrichment analysis, taking distinct genesets name and id 
TERM2NAME <- genesets %>%
  dplyr::select(gs_id, gs_name) %>%
  distinct()

################################################################################

# Term enrichment analysis with KEGG database everything default
kegg_enrichment <- degs %>%
  purrr::map(function(x) {
    x %>%
      enricher(pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 15,
               maxGSSize = 500,
               TERM2GENE = genesets %>% 
                 dplyr::filter(gs_cat == "C2", gs_subcat == "CP:KEGG") %>%
                 dplyr::select(gs_id, gene_symbol),
               TERM2NAME = TERM2NAME)
  })

# copy and paste for GO:BP/wikipathways/reactome
# Term enrichment analysis with GO:BP database
go_enrichment <- degs %>%
  purrr::map(function(x) {
    x %>%
      enricher(pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 15,
               maxGSSize = 500,
               TERM2GENE = genesets %>% 
                 dplyr::filter(gs_cat == "C5", gs_subcat == "GO:BP") %>%
                 dplyr::select(gs_id, gene_symbol),
               TERM2NAME = TERM2NAME)
  })

# Term enrichment analysis with wikipathways database
wp_enrichment <- degs %>%
  purrr::map(function(x) {
    x %>%
      enricher(pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 15,
               maxGSSize = 500,
               TERM2GENE = genesets %>% 
                 dplyr::filter(gs_cat == "C2", gs_subcat == "CP:WIKIPATHWAYS") %>%
                 dplyr::select(gs_id, gene_symbol),
               TERM2NAME = TERM2NAME)
  })

# Term enrichment analysis with reactome database
reactome_enrichment <- degs %>%
  purrr::map(function(x) {
    x %>%
      enricher(pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 15,
               maxGSSize = 500,
               TERM2GENE = genesets %>% 
                 dplyr::filter(gs_cat == "C2", gs_subcat == "CP:REACTOME") %>%
                 dplyr::select(gs_id, gene_symbol),
               TERM2NAME = TERM2NAME)
  })

################################################################################

### Volcanoplot for differentially expressed Genes with Enhanced volcanoplot

# key-values to differentiate between up- and downregulated genes in terms of logFC and adj. P-value
# adj. p-value cutoff 0.05; logFC cutoff 2, -2 
keyvals <- ifelse(res$`carcinoma - adjacent_tissue`$logFC < -2 & res$`carcinoma - adjacent_tissue`$adj.P.Val < 0.05, 'blue',
           ifelse(res$`carcinoma - adjacent_tissue`$logFC > 2 & res$`carcinoma - adjacent_tissue`$adj.P.Val < 0.05, 'red',
           'grey'))

# assigning labels to the colors 
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'upregulated genes'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'downregulated genes'

# Volcano plot for DEGs
EnhancedVolcano(res$`carcinoma - adjacent_tissue`,
                lab = res$`carcinoma - adjacent_tissue`$SYMBOL,
                x = "logFC",
                y = "adj.P.Val",
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                selectLab = rownames(res)[which(names(keyvals) %in% c('upregulated genes', 'downregulated genes'))],
                pCutoff = 0.05,
                FCcutoff = 2,
                xlim = c(-5.5,5.5),
                ylim = c(0,8),
                title = " ",
                subtitle = " ",
                #legendLabels = c('NS','Log2 FC','Adjusted p-value', 'Adj. p-value & Log2 FC'),
                #legendLabSize = 10,
                colCustom = keyvals,
                border = "full",
                borderWidth = 1,
                borderColour = 'black')

################################################################################

# Volcano plot for GSVA results

# writing gsva results in a new variable as a dataframe
res <- res.gsva %>% purrr::map(data.frame)

# removal of underscores in the gs_name column for better readability
res$`carcinoma - adjacent_tissue`$gs_name <- gsub("_", " ", res$`carcinoma - adjacent_tissue`$gs_name)

# key-values to differentiate between up- and downregulated genes in terms of logFC and adj. P-value
keyvals <- ifelse(res$`carcinoma - adjacent_tissue`$logFC < -0.5 & res$`carcinoma - adjacent_tissue`$adj.P.Val < 0.05, 'blue',
                  ifelse(res$`carcinoma - adjacent_tissue`$logFC > 0.5 & res$`carcinoma - adjacent_tissue`$adj.P.Val < 0.05, 'red',
                         'grey'))

# assigning labels to the colors
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'upregulated pathways'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'downregulated pathways'


# Writing the important pathways into the new variable taken from the xlsx file
carcinoma <- rbind(res$`carcinoma - adjacent_tissue`$gs_name[142],
                   res$`carcinoma - adjacent_tissue`$gs_name[1428],
                   res$`carcinoma - adjacent_tissue`$gs_name[1383],
                   res$`carcinoma - adjacent_tissue`$gs_name[1618],
                   res$`carcinoma - adjacent_tissue`$gs_name[12],
                   res$`carcinoma - adjacent_tissue`$gs_name[1315],
                   res$`carcinoma - adjacent_tissue`$gs_name[165],
                   res$`carcinoma - adjacent_tissue`$gs_name[831],
                   res$`carcinoma - adjacent_tissue`$gs_name[1413])



# Volcano plot for GSVA results
EnhancedVolcano(res$`carcinoma - adjacent_tissue`,
                lab = res$`carcinoma - adjacent_tissue`$gs_name,
                x = "logFC",
                y = "adj.P.Val",
                selectLab = carcinoma,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                title = " ",
                subtitle = " ",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labSize = 3.0,
                labFace = "bold",
                #legendLabels = c('NS','Log2 FC','Adjusted p-value', 'Adj. p-value & Log2 FC'),
                #legendLabSize = 10,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                lengthConnectors = 5.0,
                arrowheads = FALSE,
                colCustom = keyvals,
                border = "full",
                borderWidth = 1,
                borderColour = 'black')


################################################################################

### dotplots with clusterprofiler showing the top ten most enriched pathways with
# count displayed on the x-axis

# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html


dotplot(kegg_enrichment$`carcinoma - adjacent_tissue`, x="count", showCategory=10) + 
  ggtitle("kegg enrichment results")

dotplot(go_enrichment$`carcinoma - adjacent_tissue`, x="count", showCategory=10) + 
  ggtitle("GO enrichment results")

dotplot(wp_enrichment$`carcinoma - adjacent_tissue`, x="count", showCategory=10) + 
  ggtitle("WikiPathways enrichment results")

dotplot(reactome_enrichment$`carcinoma - adjacent_tissue`, x="count", showCategory=10) + 
  ggtitle("Reactome enrichment results")

