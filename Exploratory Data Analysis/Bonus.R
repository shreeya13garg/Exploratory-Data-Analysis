#Pathway analysis
library(enrichR)
dbs <- listEnrichrDbs()
## run enrichment analysis
#run differential expression analysis
if (is.null(dbs)) websiteLive <- FALSE
head(dbs)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)
enriched[["GO_Biological_Process_2015"]]
