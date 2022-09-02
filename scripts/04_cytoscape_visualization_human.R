# Prepare colorized Wikipathways

# load needed libraries
library(limma)
library(tidyverse)
library(rWikiPathways)
library(RCy3)


# assign variables to needed directories
datadir <- "../data/"
resultsdir <- "../results/"
plotsdir <- "../plots/"


## read in limma results

limma_human <- readRDS(file.path(datadir, "intermediate", "tumor_list_res_eset_human.RDS")) %>%
  magrittr::extract(c("carcinoma - adjacent_tissue")) %>%
  purrr::map2(names(.),
              function(x,y){
                x %>% 
                  dplyr::select("SYMBOL", 
                                "logFC",
                                "AveExpr",
                                "t",
                                "P.Value",
                                "adj.P.Val") %>%
                  set_names(c("SYMBOL", paste(c("logFC",
                                                "AveExpr",
                                                "t",
                                                "P.Value",
                                                "adj.P.Val"),
                                              y,
                                              sep = "_")))
              }) %>%
  reduce(left_join, by = "SYMBOL")


################################################################################


## read in gsva results

pathways_human <- readRDS(file.path(datadir, "intermediate", "tumor_list_res_gsva_human.RDS")) %>%
  magrittr::extract(c("carcinoma - adjacent_tissue")) %>%
  purrr::map2(names(.),
              function(x,y){
                x %>% 
                  dplyr::select("gs_exact_source",
                                "logFC",
                                "AveExpr",
                                "t",
                                "P.Value",
                                "adj.P.Val") %>%
                  set_names(c("wpid", paste(c("logFC",
                                              "AveExpr",
                                              "t",
                                              "P.Value",
                                              "adj.P.Val"),
                                            y,
                                            sep = "_")))
              }) %>%
  reduce(left_join, by = "wpid")



################################################################################


# pass all wikipathway ids to vector and merge them in one variable
commands <- paste0('wikipathways import-as-pathway id=',pathways_human$wpid) %>%
  unique()

# get Wikipathway ids from commands
commands <- commands[grep("WP", commands)]
commands <- unique(commands)

commands %>%
  purrr::map(commandsRun)

# draw the networks in cytoscape
networks <- getNetworkList()

limma_human$SYMBOL <- toupper(limma_human$SYMBOL)


networks %>%
  purrr::map(function(x){
    print(x)
    try(loadTableData(limma_human,
                  data.key.column = "SYMBOL",
                  table = "node",
                  table.key.column = "name",
                  namespace = "default",
                  network = x))
  })

networks <- getNetworkList()
names(networks) <- str_replace(networks, " - Homo sapiens", "") %>%
  make.names() %>%
  str_replace_all("\\.", "_")



#----------------------------------------------------------------
### Save the networks of interest from Cytoscape via R console 
# not necessary, look through the networks 
#networks %>%
#  purrr::map2(.,
#              names(.),
#              function(x,y){
#                file_name <- here::here(plotsdir, "PECAM1", paste(y, "png", sep = "."))
#                # setVisualStyle(
#                #   style.name
#                #   = "logFC_CD31", 
#                #   network = x)
#                fitContent(selected.only = FALSE, 
#                           network = x)
#                exportImage(filename = file_name,
#                            type = "PNG",
#                            zoom = 500,
#                            network = x)
#              })

#networks %>%
#  purrr::map2(.,
#              names(.),
#              function(x,y){
#                file_name <- here::here(plotsdir, "ENG", paste(y, "png", sep = "."))
#                setVisualStyle(
#                  style.name
#                  = "logFC_CD105", 
#                  network = x)
#                fitContent(selected.only = FALSE, 
#                           network = x)
#                exportImage(filename = file_name,
#                            type = "PNG",
#                            zoom = 500,
#                            network = x)
#              }) 

# save the networks session to load it later 
#save.image(file = "Network_session_after_loadtabledata_human")
#load(file = "Network_session_human")