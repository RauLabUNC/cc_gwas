# Load required libraries
library(InterMineR)
library(data.table)

# Load in genes list 
genes <- read.csv("results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv")


# Set up MouseMine connection
im <- initInterMine(mine=listMines()["MouseMine"])

template = getTemplates(im)

queryGenePath = getTemplateQuery(
  im = im, 
  name = "_Feature_Phenotype"
)

queryGenePath$where[[4]] <- newConstraint <- list(
  path=c("OntologyAnnotation.subject"),
  op=c("LOOKUP"), 
  value=c(paste(genes$mouse_gene_symbol, collapse=",")), 
  code=c("B")
)

results <- runQuery(im, queryGenePath)

colnames(results)[[which(colnames(results) == "OntologyAnnotation.subject.symbol")]] <- "mouse_gene_symbol"

write.csv(results, "data/processed/joinLoci/relational_tables/mouseGenePhenotypes.csv", row.names = F)
