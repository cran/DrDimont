graph_metrics <- function(graph, verbose = TRUE, return = FALSE) {
    #' @title Analysis of metrics of an iGraph object
    #'
    #' @description This helper function prints or returns multiple metrics of arbitrary
    #' iGraph graph object.
    #'
    #' @param graph [igraph] iGraph object to analyze.
    #' @param verbose [bool] If TRUE graph information is printed.
    #' @param return [bool] If TRUE graph information is returned from function.
    #'
    #' @examples
    #' adj_mat <- matrix(rnorm(36), nrow=6)
    #' graph <- igraph::graph_from_adjacency_matrix(adj_mat)
    #' DrDimont::graph_metrics(graph, verbose=TRUE, return=FALSE)
    #'
    #'
    #' @return Named list of metrics including vertex count, edge count, number of components,
    #' size of largest component and the relative frequency of zero degree vertices.
    #' 
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
    #' @export
    
    try(parallel::stopCluster(parallel::getDefaultCluster()), silent = TRUE)
    try(closeAllConnections(), silent = TRUE)
}

install_python_dependencies <- function(package_manager="pip3") {
    #' @title Installs python dependencies needed for interaction score computation
    #'
    #' @description Uses specified pip or conda executable (default: pip3) to 
    #' install all required python modules. When using conda, the currently active 
    #' environment is used. Commands run are `pip install -r requirements` or 
    #' `conda install --file requirements`. Installs the following requirements: 
    #' numpy, tqdm, python-igraph and ray
    #'
    #' @param package_manager [string] The package manager command or path to use (default: pip3)
    #'
    #' @return No return value, called to install python dependencies
    #' 
    #' @export
    
    if (grepl("pip", package_manager, fixed = TRUE)) {
        py_requirements = system.file("requirements_pip.txt", package = "DrDimont")
        system2(package_manager, args = c("install", "-r", py_requirements))
    }
    else if (grepl("conda", package_manager, fixed = TRUE)) {
        py_requirements = system.file("requirements_conda.txt", package = "DrDimont")
        system2(package_manager, args = c("install", "-c", "conda-forge", "--file", py_requirements))
    }
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
    #' @export
    
    layer_names <- sapply(layers, function(l) l[['name']])
    return(layers[[which(layer_names == name)]])
}
