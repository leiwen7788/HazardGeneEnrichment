# R/display_geneset.r

#'
#'
#' This function displays the characteristics of the Hazard Ratio gene sets we have curated/included.
#'


#' @return Hazard Ratio gene sets we have curated/included.
#' @export
#'
#' @examples
#' library(HazardGeneEnrichment)
#' display_genesets()
#'
#'
display_genesets<-function(){
  print(names(gene_list_all))
}

