#' @title Gene Set Enrichment Analysis (GSEA) using Multiplex Networks
#' @name get_gsea_path
#' @description This function performs gene set enrichment analysis (GSEA) based on multiplex network data.
#'
#' @param seed A seed value (optional).
#' @param network A network object (e.g., protein-protein interaction network).
#' @param gamma A parameter for random walk restart (default: 0.7).
#' @param pathlist A predefined list of gene sets (pathways).
#' @param gsea.weight Weight for GSEA (default: 1).
#' @param gsea.nperm Number of permutations for significance testing (default: 1000).
#'
#' @return A GSEA result object.
#'
#' @details The function constructs a multiplex network, performs random walk restart, and calculates gene scores.
#' It then transforms the scores and applies GSEA using the provided gene sets.
#'
#' @examples
#' data(Seeds, package = "iPRISM")
#' data(ppi, package = "iPRISM")
#' data(path_list, package = "iPRISM")
#'
#' \donttest{
#' result <- get_gsea_path(seed = Seeds,
#'                        network = ppi,
#'                        pathlist = path_list[1:2],
#'                        gsea.nperm = 100)
#' print(result)
#' }
#'
#' @importFrom igraph V
#'
#' @keywords gene set enrichment analysis multiplex network pathway
#'
#' @export
get_gsea_path <- function(seed = seed,
                          network = network,
                          gamma = 0.7,
                          pathlist = pathlist,
                          gsea.weight = 1,
                          gsea.nperm = 1000){
  net_MultiplexObject <- create.multiplex(list(net = network))
  AdjMatrix_net <- compute.adjacency.matrix(net_MultiplexObject)
  AdjMatrixNorm_net <- normalize.multiplex.adjacency(AdjMatrix_net)

  RWR_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_net,r = gamma,
                                               net_MultiplexObject,Seeds = intersect(seed,names(V(network))),
                                               DispResults = "Alphabetic")

  walk_result <- RWR_Results$RWRM_Results$Score
  names(walk_result) <- RWR_Results$RWRM_Results$NodeNames

  genelist <- scale(10/(-log10(walk_result)))
  names(genelist) <- names(walk_result)
  genelist <- sort(genelist,decreasing = T)
  res_gsea <- gseafun(genelist = genelist,pathlist = pathlist,weighted = gsea.weight,nperm = gsea.nperm)
  return(res_gsea)
}

