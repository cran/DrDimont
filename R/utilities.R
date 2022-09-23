graph_metrics <- function(graph, verbose = TRUE, return = FALSE) {
    #' @title [INTERNAL] Analysis of metrics of an iGraph object
    #'
    #' @description [INTERNAL] This helper function prints or returns multiple metrics of arbitrary
    #' iGraph graph object.
    #'
    #' @param graph [igraph] iGraph object to analyze.
    #' @param verbose [bool] If TRUE graph information is printed.
    #' @param return [bool] If TRUE graph information is returned from function.
    #'
    #' @return Named list of metrics including vertex count, edge count, number of components,
    #' size of largest component and the relative frequency of zero degree vertices.
    #' 
    #' @keywords internal
    #' @export

    metrics <- list('n_vertices' = igraph::gorder(graph),
                    'n_edges' = igraph::ecount(graph),
                    'components' = igraph::components(graph),
                    'degree_distribution' = igraph::degree_distribution(graph)
                    )
    ### print graph metrics
    if(verbose == TRUE) {
        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Vertex count: ', metrics$n_vertices)
        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Edge count: ', metrics$n_edges)
        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Number of components: ', metrics$components$no)
        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Size of largest component: ', metrics$components$csize[1])
        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Relative frequency of zero degree vertices: ', metrics$degree_distribution[1])
    }

    if(return == TRUE) {return(metrics)}
}


set_cluster <- function(n_threads) {
    #' @title [INTERNAL] Create and register cluster
    #'
    #' @description [INTERNAL] Helper function to create and register a cluster for parallel
    #' computation of p-value reduction
    #'
    #' @param n_threads [int] Number of nodes in the cluster
    #'
    #' @return No return value, called internally to create cluster
    #' 
    #' @keywords internal
    #' @export
    
    parallel::setDefaultCluster(parallel::makeCluster(n_threads))
}


shutdown_cluster <- function() {
    #' @title [INTERNAL] Shutdown cluster and remove corresponding connections
    #'
    #' @description [INTERNAL] Run this if the pipeline fails during parallel 
    #' computation to clean the state. If a cluster is registered, this functions 
    #' stops it and removes corresponding connections. Ignores errors. Has no effect 
    #' if no cluster is registered.
    #' 
    #' @return No return value, called internally to shutdown cluster
    #' 
    #' @keywords internal
    #' @export
    
    try(parallel::stopCluster(parallel::getDefaultCluster()), silent = TRUE)
    try(closeAllConnections(), silent = TRUE)
}


get_layer <- function(name, layers) {
    #' @title [INTERNAL] Fetch layer by name from layer object
    #' 
    #' @description [INTERNAL] Get a layer by its name from a layer object 
    #' created with \code{\link[DrDimont]{make_layer}}, e.g., \code{\link{layers_example}}.
    #'
    #' @param name The layer to fetch
    #' @param layers A layers object \code{\link{layers_example}}
    #' 
    #' @return Returns the layer along with layer names
    #' 
    #' @keywords internal
    #' @export
    
    layer_names <- sapply(layers, function(l) l[['name']])
    return(layers[[which(layer_names == name)]])
}
