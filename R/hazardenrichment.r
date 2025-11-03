# R/hazardenrichment.r

#'
#'
#' This function fits a Cox proportional hazards model and extracts
#' the hazard ratios for specified genes.
#'
#' @param geneList The reference gene sets we need to enrich to,for example 'TCGA-BRCA'
#' @param TERM2GENE The gene sets to be analyzed.
#' @param pvalueCutoff cut off p-values default 0.05.
#' @return my_gsea_results
#' @export
#'
#' @examples
#' library(HazardGeneEnrichment)
#' HazardEnrichment('TCGA-BRCA',TERM2GENE)
HazardEnrichment<-function(geneList,TERM2GENE,pvalueCutoff=0.05,singlecell_gene){
  geneList<-gene_list_all[[geneList]]
  geneList_gene<-names(geneList)
  geneList_gene<-geneList_gene[geneList_gene%in%singlecell_gene]
  geneList<-geneList[geneList_gene]

  geneList_1<-geneList[geneList>0]
  geneList_2<-geneList[geneList<0]
  geneList_1_max<-abs(max(geneList_1))
  geneList_2_min<-abs(min(geneList_2))
  cutoff_gene<-min(geneList_1_max,geneList_2_min)
  geneList_1<- geneList_1/max(geneList_1)*cutoff_gene
  geneList_2<- -(geneList_2/min(geneList_2))*cutoff_gene
  geneList<-c(geneList_1,geneList_2)
  geneList<-geneList[geneList<=2&geneList>=-2]
  geneList<-sort(geneList, decreasing = TRUE)

  my_gsea_results <- GSEA(
    geneList = geneList,
    TERM2GENE = TERM2GENE,
    # TERM2NAME = TERM2NAME,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  return(my_gsea_results)
}
