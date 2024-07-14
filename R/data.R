#' @title data_sig
#'
#' @description
#' The `data_sig` represents the sample feature matrix, where rows correspond to samples, and columns correspond to features.
#'
#' @examples
#' library(iPRISM)
#' data(data_sig, package = "iPRISM")
#' plot(data_sig)
#'
"data_sig"

#' @title data.path
#'
#' @description
#' The `data.path` represents the first type of feature matrix used for calculating correlations (in this case, pathway expression levels), where rows correspond to samples, and columns correspond to features.
#'
#' @examples
#' library(iPRISM)
#' data(data.path, package = "iPRISM")
#' plot(data.path)
#'
"data.path"

#' @title data.cell
#'
#' @description
#' The `data.cell` represents the second type of feature matrix used for calculating correlations (in this case, cell abundances), where rows correspond to samples, and columns correspond to features.
#'
#' @examples
#' library(iPRISM)
#' data(data.cell, package = "iPRISM")
#' dim(data.cell)
#'
"data.cell"

#' @title path_list
#'
#' @description
#' The `path_list` contains the gene list associated with pathways.
#'
#' @examples
#' library(iPRISM)
#' data(path_list, package = "iPRISM")
#' length(path_list)
#'
"path_list"

#' @title A protein-protein physical interaction network (PPI network)
#'
#' @description
#' An igraph object containing a protein-protein physical interaction network.
#'
#' @examples
#' library(iPRISM)
#' data(ppi, package = "iPRISM")
#' \donttest{
#' library(igraph)
#' graph <- simplify(ppi)
#' graph_comp <- components(graph)$membership == which.max(components(graph)$csize)
#' graph <- induced_subgraph(graph, V(graph)[graph_comp])
#' plot(graph)
#' }
#'
"ppi"

#' @title Original Class Labels for Samples
#'
#' @description
#' A named vector where each element corresponds to a sample name and represents the original class label.
#'
#' @examples
#' library(iPRISM)
#' data(pred_value, package = "iPRISM")
#' table(pred_value)
#'
"pred_value"

#' @title Original Class Labels for Samples
#'
#' @description
#' A named vector where each element corresponds to a sample name and represents the original class label.
#'
#' @examples
#' library(iPRISM)
#' data(pred_value, package = "iPRISM")
#' table(pred_value)
#'
"pred_value"

#' @title Seed Node Names
#'
#' @description
#' A character vector with seed node names.
#'
#' @examples
#' library(iPRISM)
#' data(Seeds, package = "iPRISM")
#'
"Seeds"

#' @title TME gene list after random walks
#'
#' @description
#' This gene list includes genes from tumor microenvironment (TME). Random Walk with Restart (RWR) is applied to prioritize genes that are relevant to immunotherapy responses.
#'
#' @examples
#' library(iPRISM)
#' data(genelist_cp, package = "iPRISM")
#'
"genelist_cp"

#' @title HLA gene list after random walks
#'
#' @description
#' This gene list includes genes from human leukocyte antigen (HLA). Random Walk with Restart (RWR) is applied to prioritize genes that are relevant to immunotherapy responses.
#'
#' @examples
#' library(iPRISM)
#' data(genelist_hla, package = "iPRISM")
#'
"genelist_hla"

#' @title ICI gene list after random walks
#'
#' @description
#' This gene list includes genes from immune checkpoint inhibitors (ICI). Random Walk with Restart (RWR) is applied to prioritize genes that are relevant to immunotherapy responses.
#'
#' @examples
#' library(iPRISM)
#' data(genelist_imm, package = "iPRISM")
#'
"genelist_imm"

