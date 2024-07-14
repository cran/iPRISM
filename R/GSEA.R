#' @title Gene Set Enrichment Analysis (GSEA) Function
#' @name gseafun
#' @description This function performs gene set enrichment analysis using a gene list and a set of pathways.
#'
#' @param genelist A named vector of gene expression values.
#' @param pathlist A list of gene sets (pathways) to test for enrichment.
#' @param nperm Number of permutations for calculating p-values (default is 1000).
#' @param weighted Logical indicating whether to use weighted enrichment scores (default is TRUE).
#'
#' @return A data frame with enrichment scores (ES), p-values, and adjusted p-values.
#'
#' @examples
#' data(path_list, package = "iPRISM")
#' data(genelist_imm, package = "iPRISM")
#' \donttest{
#' res_gsea_imm <- gseafun(genelist = genelist_imm,
#'                         pathlist = path_list[1:2],
#'                         weighted = 1,
#'                         nperm = 1000)
#' print(res_gsea_imm)
#' }
#'
#' @importFrom stats p.adjust
#' @importFrom pbapply pbsapply pblapply
#'
#' @export
gseafun <- function(genelist,
                    pathlist,
                    nperm = 1000,
                    weighted = 1){
  if(weighted == 0){
    pathway_ES <- pbsapply(pathlist, function(s){
      tag.indicator <- sign(match(names(genelist), s, nomatch = 0))
      ES <- ESscore(tag.indicator,correl.vector = genelist)
      return(ES)
    })#Compute enrichment score

    rd_geneid <- pblapply(c(1:nperm),function(r){
      rd_geneid1<-sample(names(genelist),replace = F)
      return(rd_geneid1)
    })#Perturbed gene names
    rd_geneid <- do.call(cbind, rd_geneid)#Randomized gene name matrix

    rd_ES<-pblapply(c(1:nperm),function(r){
      rdES<-sapply(pathlist, function(s){
        tag.indicator <- sign(match(rd_geneid[,r], s, nomatch = 0))
        ES <- ESscore(tag.indicator,correl.vector = genelist)
      })
      return(rdES)
    })
    rd_ES <- do.call(cbind,rd_ES)
  }
  if(weighted == 1){
    pathway_ES <- pbsapply(pathlist, function(s){
      tag.indicator <- sign(match(names(genelist), s, nomatch = 0))
      ES <- ESscore_weighted(tag.indicator,correl.vector = genelist)
      return(ES)
    })#Compute enrichment score

    rd_geneid <- pblapply(c(1:nperm),function(r){
      rd_geneid1<-sample(names(genelist),replace = F)
      return(rd_geneid1)
    })#Perturbed gene names
    rd_geneid <- do.call(cbind, rd_geneid)#Randomized gene name matrix

    rd_ES<-pblapply(c(1:nperm),function(r){
      rdES<-sapply(pathlist, function(s){
        tag.indicator <- sign(match(rd_geneid[,r], s, nomatch = 0))
        ES <- ESscore_weighted(tag.indicator,correl.vector = genelist)
      })
      return(rdES)
    })
    rd_ES <- do.call(cbind,rd_ES)
  }
  p.vals <- NULL
  for(i in 1:length(pathway_ES)){
    if(is.na(pathway_ES[i])==T){
      p.vals[i]<-1
    }else{
      if(pathway_ES[i]>=0){
        p.vals[i]<-sum(rd_ES[i,] >= pathway_ES[i])/nperm
      }else{
        p.vals[i]<-sum(rd_ES[i,] <= pathway_ES[i])/nperm
      }
    }
  }
  gsea_table <- data.frame(ES = pathway_ES,
                           pval = p.vals,
                           padj = p.adjust(p.vals,method = "fdr"))
  return(gsea_table)
}

#' @title Weighted Enrichment Score Calculation
#' @name ESscore_weighted
#' @description
#' Calculates the weighted enrichment score (ES) for a given set of labels and correlation vector.
#'
#' @param labels.list A binary vector indicating membership in a gene set (1 for inclusion, 0 for exclusion).
#' @param correl.vector A vector of correlation values (e.g., gene expression correlations).
#'
#' @return The weighted enrichment score (ES) for the given labels and correlation vector.
#'
#' @importFrom stats p.adjust
#'
ESscore_weighted <- function(labels.list,correl.vector = NULL){
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list==1])
  Nm <-  N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)

  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  return(ES)
}

#' @title Enrichment Score Calculation
#' @name ESscore
#' @description
#' Calculates the enrichment score (ES) for a given set of labels and correlation vector.
#'
#' @param labels.list A binary vector indicating membership in a gene set (1 for inclusion, 0 for exclusion).
#' @param correl.vector A vector of correlation values (e.g., gene expression correlations).
#'
#' @return The enrichment score (ES) for the given labels and correlation vector.
#'
ESscore <- function(labels.list,correl.vector = NULL){
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list==1])
  Nm <-  N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag    <- sum(tag.indicator[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * norm.tag - no.tag.indicator * norm.no.tag)

  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  return(ES)
}
