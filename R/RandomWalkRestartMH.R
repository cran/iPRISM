## Modified from: https://github.com/alberto-valdeolivas/RandomWalkRestartMH

## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Functions to create Multiplex and Multiplex-Heterogeneous objects.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Roxy Documentation comments
#' Create multiplex graphs from individual networks
#'
#' \code{create.multiplex} is a function to create a multiplex network
#' (\code{Multiplex} object) from a list of individual networks defined as
#' \code{igraph} objects. See more details about multiplex networks below.
#' If just one network is provided, a Multiplex object with one layer is
#' therefore created (A monoplex network).
#'
#' @usage create.multiplex(...)
#'
#' @details A multiplex network is a collection of layers (monoplex networks)
#' sharing the same nodes, but in which the edges represent relationships of
#' different nature. At least a list with one element, an igraph object, should
#' be provided.
#'
#' @param LayersList A list containing igraph objects describing monoplex
#' networks in every element. We recommend to give names to the different
#' networks (igraph objects).
#' @param ... Further arguments passed to \code{create.multiplex}
#'
#' @return A Multiplex object. It contains a list of the different graphs
#' integrating the multiplex network, the names and number of its nodes and the
#' number of layers.
#'
#' @seealso \code{\link{isMultiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#'
#' @importFrom igraph is.igraph set_vertex_attr set_edge_attr E vcount
#' @importFrom methods is
#' @noRd
create.multiplex <- function(LayersList, ...){
  if (!is(LayersList, "list")){
    stop("The input object should be a list of graphs.")
  }


  Number_of_Layers <- length(LayersList)
  SeqLayers <- seq(Number_of_Layers)
  Layers_Name <- names(LayersList)

  if (!all(sapply(SeqLayers, function(x) is.igraph(LayersList[[x]])))){
    stop("Not igraph objects")
  }

  Layer_List <- lapply(SeqLayers, function (x) {
    if (is.null(V(LayersList[[x]])$name)){
      LayersList[[x]] <-
        set_vertex_attr(LayersList[[x]],"name",
                        value=seq(1,vcount(LayersList[[x]]),by=1))
    } else {
      LayersList[[x]]
    }
  })

  ## We simplify the layers
  Layer_List <-
    lapply(SeqLayers, function(x) simplify.layers(Layer_List[[x]]))

  ## We set the names of the layers.

  if (is.null(Layers_Name)){
    names(Layer_List) <- paste0("Layer_", SeqLayers)
  } else {
    names(Layer_List) <- Layers_Name
  }

  ## We get a pool of nodes (Nodes in any of the layers.)
  Pool_of_Nodes <-
    sort(unique(unlist(lapply(SeqLayers,
                              function(x) V(Layer_List[[x]])$name))))

  Number_of_Nodes <- length(Pool_of_Nodes)

  Layer_List <-
    lapply(Layer_List, add.missing.nodes,Number_of_Layers,Pool_of_Nodes)

  # We set the attributes of the layer
  counter <- 0
  Layer_List <- lapply(Layer_List, function(x) {
    counter <<- counter + 1;
    set_edge_attr(x,"type",E(x), value = names(Layer_List)[counter])
  })


  MultiplexObject <- c(Layer_List,list(Pool_of_Nodes=Pool_of_Nodes,
                                       Number_of_Nodes_Multiplex=Number_of_Nodes,
                                       Number_of_Layers=Number_of_Layers))

  class(MultiplexObject) <- "Multiplex"

  return(MultiplexObject)
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Functions to compute the matrices and perform the Random Walks.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Roxy Documentation comments
#' Computes the adjacency matrix of a multiplex network
#'
#' \code{compute.adjacency.matrix} is a function to compute the adjacency
#' matrix of a multiplex network provided as a \code{Multiplex} object.
#'
#' @usage compute.adjacency.matrix(x,delta = 0.5)
#'
#' @details The parameter \code{delta} sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#'
#' @param x A \code{Multiplex} object describing a multiplex network generated
#' by the function \code{create.multiplex}.
#' @param delta A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions. It is set by
#' default to 0.5. See more details below.
#'
#' @return A square sparse adjacency matrix created with the \code{Matrix}
#' package.
#'
#' @seealso \code{\link{create.multiplex},\link{normalize.multiplex.adjacency}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' compute.adjacency.matrix(multiObject)
#'
#' @importFrom igraph as_adjacency_matrix is_weighted
#' @importFrom Matrix Diagonal bdiag
#' @importFrom methods as
#' @noRd
compute.adjacency.matrix <- function(x,delta = 0.5)
{
  if (!isMultiplex(x)) {
    stop("Not a Multiplex or Multiplex Heterogeneous object")
  }
  if (delta > 1 || delta <= 0) {
    stop("Delta should be between 0 and 1")
  }

  N <- x$Number_of_Nodes_Multiplex
  L <- x$Number_of_Layers

  ## We impose delta=0 in the monoplex case.
  if (L==1){
    delta = 0
  }

  Layers_Names <- names(x)[seq(L)]

  ## IDEM_MATRIX.
  Idem_Matrix <- Diagonal(N, x = 1)

  counter <- 0
  Layers_List <- lapply(x[Layers_Names],function(x){

    counter <<- counter + 1;
    if (is_weighted(x)){
      Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE,
                                              attr = "weight")
    } else {
      Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE)
    }

    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                       order(colnames(Adjacency_Layer))]
    colnames(Adjacency_Layer) <-
      paste0(colnames(Adjacency_Layer),"_",counter)
    rownames(Adjacency_Layer) <-
      paste0(rownames(Adjacency_Layer),"_",counter)
    Adjacency_Layer
  })

  MyColNames <- unlist(lapply(Layers_List, function (x) unlist(colnames(x))))
  MyRowNames <- unlist(lapply(Layers_List, function (x) unlist(rownames(x))))
  names(MyColNames) <- c()
  names(MyRowNames) <- c()
  SupraAdjacencyMatrix <- (1-delta)*(bdiag(unlist(Layers_List)))
  colnames(SupraAdjacencyMatrix) <-MyColNames
  rownames(SupraAdjacencyMatrix) <-MyRowNames

  offdiag <- (delta/(L-1))*Idem_Matrix

  i <- seq_len(L)
  Position_ini_row <- 1 + (i-1)*N
  Position_end_row <- N + (i-1)*N
  j <- seq_len(L)
  Position_ini_col <- 1 + (j-1)*N
  Position_end_col <- N + (j-1)*N

  for (i in seq_len(L)){
    for (j in seq_len(L)){
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]),
                             (Position_ini_col[j]:Position_end_col[j])] <- offdiag
      }
    }
  }

  SupraAdjacencyMatrix <- as(SupraAdjacencyMatrix, "dgCMatrix")
  return(SupraAdjacencyMatrix)
}


## Roxy Documentation comments
#' Computes column normalization of an adjacency matrix
#'
#' \code{normalize.multiplex.adjacency} is a function to compute the column
#' normalization of a sparse matrix of the package \code{Matrix}.
#'
#' @usage normalize.multiplex.adjacency(x)
#'
#' @param x A \code{Matrix} object describing an adjacency matrix of a network.
#'
#' @return A square sparse column normalized matrix created with the
#' \code{Matrix} package.
#'
#' @seealso \code{\link{compute.adjacency.matrix},
#' \link{Random.Walk.Restart.Multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' normalize.multiplex.adjacency(AdjMatrix)
#'
#' @importFrom Matrix colSums t
#' @importFrom methods is
#' @noRd
normalize.multiplex.adjacency <- function(x)
{
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }

  Adj_Matrix_Norm <- t(t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1,
                                             sparseResult = FALSE)))

  return(Adj_Matrix_Norm)
}

## Roxy Documentation comments
#' Performs Random Walk with Restart on a Multiplex Network
#'
#' \code{Random.Walk.Restart.Multiplex} is a function to perform a Random Walk
#' with Restart on a Multiplex network (on a \code{Multiplex} object). See
#' more details about the algorithm below.
#'
#' @usage Random.Walk.Restart.Multiplex(...)
#'
#' @details Random Walk with Restart simulates an imaginary particle that
#' starts on a seed(s) node(s) and follows randomly the edges of a network. At
#' each step, there is a restart probability, r, meaning that the particle comes
#' back to the seed(s). The extension to multiplex networks allows the particle
#' to explore different monoplex networks (layers). At each step, the particle
#' can also jump to the same node in a different layer.
#'
#'
#' \itemize{
#' \item \code{Seeds}: A vector containing the name of the different seed
#' node(s). It's mandatory to provide at least one seed. The seed(s) node(s)
#' should belong to any of the layers. The length of this vector should be
#' smaller than the total number of nodes in the multiplex network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most common value in this kind of approaches. It means
#' that, at each step, the walker has a 70\% of probability of coming back to
#' one of the seeds.
#' \item \code{tau}: A numeric vector containing the probability of restarting
#' in the nodes of the different layers of the multiplex. In the example below,
#' we define the node 1 as the seed node. However, we can find this node in both
#' layers. Therefore, the walker can restart in any of these seed nodes. It is
#' a way to give different relevance (weight) to the different layers.
#' }
#'
#'
#' @param x An object of the \code{Matrix} package describing a column
#' normalized adjacency matrix of a multiplex network.
#' @param MultiplexObject A \code{Multiplex} object generated by the function
#' \code{create.multiplex} representing a multiplex network.
#' @param Seeds A vector containing the names of the seeds for the Random
#' Walk algorithm. See more details below.
#' @param r A numeric value between 0 and 1. It sets the probability of
#' restarting to a seed node after each step. See more details below.
#' @param tau A vector containing the probability of restart on the seeds
#' of the different layers (layers weights). It must have the same length than
#' the number of layers of the multpiplex network. The sum of its components
#' divided by the number of layers must be 1. See more details below.
#' @param MeanType The user can choose one of the following options:
#' c("Geometric","Arithmetic","Sum"). These options represent the different way
#' to combine the RWR score for the same node in different layers. By default
#' and recommended Geometric (Geometric Mean.). Arithmetic is the arithmetic
#' mean and sum just sum all the scores for the same node across the different
#' layers.
#' @param DispResults The user can choose one of the following options:
#' c("TopScores","Alphabetic"). These options represent the way the RWR results
#' would be presented. By default, and recommended, the nodes would be ordered
#' by score. This option is also required to properly run the
#' \code{create.multiplexNetwork.topResults} and
#' \code{create.multiplexHetNetwork.topResults} functions
#' @param ... Further arguments passed to \code{Random.Walk.Restart.Multiplex}
#'
#' @return A \code{RWRM_Results} object. It contains a sorted ranking of all
#' the nodes of the multiplex network, except the seeds, along with their score.
#' In addition, it contains in a different field the nodes used as seeds.
#'
#' @seealso \code{\link{create.multiplex},
#' \link{compute.adjacency.matrix}, \link{normalize.multiplex.adjacency}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
#' SeedNodes <- c(1)
#' Random.Walk.Restart.Multiplex(AdjMatrixNorm,multiObject,SeedNodes)
#'
#' @importFrom methods is
#' @noRd
Random.Walk.Restart.Multiplex <- function(x, MultiplexObject, Seeds,
                                          r=0.7,tau,MeanType="Geometric", DispResults="TopScores",...){
  ### We control the different values.
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }

  if (!isMultiplex(MultiplexObject)) {
    stop("Not a Multiplex object")
  }

  L <- MultiplexObject$Number_of_Layers
  N <- MultiplexObject$Number_of_Nodes

  Seeds <- as.character(Seeds)
  if (length(Seeds) < 1 | length(Seeds) >= N){
    stop("The length of the vector containing the seed nodes is not
         correct")
  } else {
    if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)){
      stop("Some of the seeds are not nodes of the network")

    }
  }

  if (r >= 1 || r <= 0) {
    stop("Restart partameter should be between 0 and 1")
  }

  if(missing(tau)){
    tau <- rep(1,L)/L
  } else {
    tau <- as.numeric(tau)
    if (sum(tau)/L != 1) {
      stop("The sum of the components of tau divided by the number of
           layers should be 1")
    }
  }

  if(!(MeanType %in% c("Geometric","Arithmetic","Sum"))){
    stop("The type mean should be Geometric, Arithmetic or Sum")
  }

  if(!(DispResults %in% c("TopScores","Alphabetic"))){
    stop("The way to display RWRM results should be TopScores or
         Alphabetic")
  }

  ## We define the threshold and the number maximum of iterations for
  ## the random walker.
  Threeshold <- 1e-10
  NetworkSize <- ncol(x)

  ## We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1

  ## We compute the scores for the different seeds.
  Seeds_Score <- get.seed.scoresMultiplex(Seeds,L,tau)

  ## We define the prox_vector(The vector we will move after the first RWR
  ## iteration. We start from The seed. We have to take in account
  ## that the walker with restart in some of the Seed nodes, depending on
  ## the score we gave in that file).
  prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)

  prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])

  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector

  while(residue >= Threeshold){

    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1;
  }

  NodeNames <- character(length = N)
  Score = numeric(length = N)

  rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
  rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])

  if (MeanType=="Geometric"){
    rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]),L,N)
  } else {
    if (MeanType=="Arithmetic") {
      rank_global$Score <- regular.mean(as.vector(prox_vector[,1]),L,N)
    } else {
      rank_global$Score <- sumValues(as.vector(prox_vector[,1]),L,N)
    }
  }

  if (DispResults=="TopScores"){
    ## We sort the nodes according to their score.
    Global_results <-
      rank_global[with(rank_global, order(-Score, NodeNames)), ]

    ### We remove the seed nodes from the Ranking and we write the results.
    Global_results <-
      Global_results[which(!Global_results$NodeNames %in% Seeds),]
  } else {
    Global_results <- rank_global
  }

  rownames(Global_results) <- c()

  RWRM_ranking <- list(RWRM_Results = Global_results,Seed_Nodes = Seeds)

  class(RWRM_ranking) <- "RWRM_Results"
  return(RWRM_ranking)
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Functions to check the objects used within the package.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Roxy Documentation comments
#' Is this R object a Multiplex object?
#'
#' A Multiplex object is an R object generated as the result of calling
#' the function \code{create.multiplex}
#'+
#' \code{isMultiplex(x)} checks whether an R object is Multiplex.
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a Mutiplex object.
#'
#' @seealso \code{\link{create.multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' isMultiplex(multiObject)
#' isMultiplex(m1)
#'
#' @importFrom methods is
#' @noRd
isMultiplex <- function (x)
{
  is(x,"Multiplex")
}

## Functions to perform Random Walk with Restart on Multiplex Networks.
#' @importFrom igraph as.undirected is_weighted E simplify igraph_opt ecount `E<-`
#' @noRd
simplify.layers <- function(Input_Layer){

  ## Undirected Graphs
  Layer <- as.undirected(Input_Layer, mode = c("collapse"),
                         edge.attr.comb = igraph_opt("edge.attr.comb"))

  ## Unweighted or Weigthed Graphs
  if (is_weighted(Layer)){
    b <- 1
    weigths_layer <- E(Layer)$weight
    if (min(weigths_layer) != max(weigths_layer)){
      a <- min(weigths_layer)/max(weigths_layer)
      range01 <- (b-a)*(weigths_layer-min(weigths_layer))/
        (max(weigths_layer)-min(weigths_layer)) + a
      E(Layer)$weight <- range01
    } else {
      E(Layer)$weight <- rep(1, length(weigths_layer))
    }
  } else {
    E(Layer)$weight <- rep(1, ecount(Layer))
  }

  ## Simple Graphs
  Layer <-
    simplify(Layer,remove.multiple = TRUE,remove.loops = TRUE,
             edge.attr.comb=mean)

  return(Layer)
}

## Add missing nodes in some of the layers.
#' @importFrom igraph add_vertices
#' @noRd
add.missing.nodes <- function (Layers,Nr_Layers,NodeNames) {

  add_vertices(Layers,
               length(NodeNames[which(!NodeNames %in% V(Layers)$name)]),
               name=NodeNames[which(!NodeNames %in%  V(Layers)$name)])
}

## Scores for the seeds of the multiplex network.
get.seed.scoresMultiplex <- function(Seeds,Number_Layers,tau) {

  Nr_Seeds <- length(Seeds)

  Seeds_Seeds_Scores <- rep(tau/Nr_Seeds,Nr_Seeds)
  Seed_Seeds_Layer_Labeled <-
    paste0(rep(Seeds,Number_Layers),sep="_",rep(seq(Number_Layers),
                                                length.out = Nr_Seeds*Number_Layers,each=Nr_Seeds))

  Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
                            Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)

  return(Seeds_Score)
}

geometric.mean <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
  }

  return(FinalScore)
}


regular.mean <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- mean(Scores[seq(from = i, to = N*L, by=N)])
  }

  return(FinalScore)
}

sumValues <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- sum(Scores[seq(from = i, to = N*L, by=N)])
  }

  return(FinalScore)
}
