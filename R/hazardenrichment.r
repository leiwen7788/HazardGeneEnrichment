# R/hazardenrichment.r

#'
#'
#' This function fits a Cox proportional hazards model and extracts
#' the hazard ratios for specified genes.
#'
#' @param geneList The reference gene sets we need to enrich to,for example 'TCGA-BRCA'
#' @param TERM2GENE The gene sets to be analyzed.
#' @return
#' @export
#'
#' @examples
#' library(HazardGeneEnrichment)
#' HazardEnrichment('TCGA-BRCA',TERM2GENE)
HazardEnrichment<-function(geneList,TERM2GENE){
  my_gsea_results <- GSEA(
    geneList = gene_list_all[[geneList]],
    TERM2GENE = TERM2GENE,
    # TERM2NAME = TERM2NAME,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  return(my_gsea_results)
}
