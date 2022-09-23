create_unique_layer_node_ids <- function(identifiersA, identifiersB, layer_name) {
    #' @title [INTERNAL] Assign node IDs to the biological identifiers across a graph layer
    #'
    #' @description [INTERNAL] This function takes two data frames of (biological) identifiers of nodes.
    #' Each data frame corresponds to the identifiers of the components contained in the single-layer
    #' network of a sample group. This function outputs the same data frames, with an added column (`node_id`)
    #' that contains node IDs which can later be used as `name` parameter for an iGraph graph.
    #' Node IDs begin with the defined `prefix` and an underscore. If a molecule is present in both
    #' groups, the node ID will be the same across the whole layer, allowing to easily combine the
    #' graphs of both groups in \code{\link[DrDimont]{generate_differential_score_graph}} to calculate
    #' differential scores of identical nodes in both sample groups.
    #' The function is used by the high-level wrapper \link[DrDimont]{generate_individual_graphs} to
    #' create annotations, which uniquely define nodes across the network layer.
    #'
    #' @param identifiersA,identifiersB [data.frame] Containing the biological identifiers of each
    #' group of the same network layer.
    #' @param layer_name [string] Name of layer that the node ids are created for
    #'
    #' @return Returns an named list. Elements `groupA` and `groupB` contain the input
    #' data frames with an additional column `node_id`. `both` contains all unique node IDs assigned
    #' across the network layer.
    #' 
    #' @keywords internal
    #' @export

    ##### check if only 1 group given
    if (is.null(identifiersB)){

        if (sum(duplicated(identifiersA)) > 0){
            ##### show warning if duplicated node ids given
            message(format(Sys.time(), "[%y-%m-%d %X] "), stringr::str_interp("WARNING: Duplicate node IDs given in layer ${layer_name}, this may cause ERRORS!"))
        }

        assigned_ids <- dplyr::union(identifiersA, identifiersA) %>%
            dplyr::mutate(node_id=paste0(layer_name, "_", dplyr::row_number()), layer=layer_name)

        identifiersA <- suppressMessages(dplyr::left_join(identifiersA, assigned_ids))

        return(list(groupA=identifiersA, groupB=NULL, both=assigned_ids))

    } else {

        if((sum(duplicated(identifiersA)) > 0)|(sum(duplicated(identifiersB)) > 0)){
            ##### show warning if duplicated node ids given
            message(format(Sys.time(), "[%y-%m-%d %X] "), stringr::str_interp("WARNING: Duplicate node IDs given in layer ${layer_name}, this may cause ERRORS!"))
        }
        assigned_ids <- dplyr::union(identifiersA, identifiersB) %>%
            dplyr::mutate(node_id=paste0(layer_name, "_", dplyr::row_number()), layer=layer_name)

        identifiersA <- suppressMessages(dplyr::left_join(identifiersA, assigned_ids))
        identifiersB <- suppressMessages(dplyr::left_join(identifiersB, assigned_ids))

        return(list(groupA=identifiersA, groupB=identifiersB, both=assigned_ids))
    }
}

generate_reduced_graph <- function(adjacency_matrix,
                                   measurement_data,
                                   identifiers,
                                   handling_missing_data='all.obs',
                                   reduction_method='pickHardTreshold',
                                   r_squared_cutoff=0.85,
                                   cut_vector=seq(0.2, 0.8, by=0.01),
                                   mean_number_edges=NULL,
                                   edge_density=NULL,
                                   p_value_adjustment_method='BH',
                                   reduction_alpha=0.05,
                                   n_threads=1,
                                   parallel_chunk_size=10^6,
                                   print_graph_info=TRUE){
    #' @title [INERNAL] Generate a reduced iGraph from adjacency matrices
    #'
    #' @description [INTERNAL] A wrapper functions that calls the functions to generate a network from
    #' correlation data and reduce the network by a given method. Correlation/adjacency matrices are
    #' computed in \code{\link[DrDimont]{compute_correlation_matrices}}. Graph generation uses
    #' \code{\link[igraph]{graph.adjacency}} internally. Methods implemented are
    #' \link[DrDimont]{network_reduction_by_p_value} (reduction by statistical significance of correlation)
    #' and \link[DrDimont]{network_reduction_by_pickHardThreshold} (using WGCNA function
    #' \link[WGCNA]{pickHardThreshold.fromSimilarity} that finds a suitable cutoff value to get a scale-free
    #' network). If no method is given, no reduction will be performed. When using the reduction method
    #' `p_value` the user can specify an alpha significance value and a method for p-value adjustment.
    #' When using the reduction by `pickHardThreshold` a R^2 cutoff and a cut vector can be specified.
    #'
    #' @param adjacency_matrix [matrix] Adjacency matrix of correlations computed using \code{\link[WGCNA]{cor}} in
    #' \code{\link[DrDimont]{compute_correlation_matrices}}
    #' @param measurement_data [data.frame] Data frame containing the respective raw data (e.g. mRNA expression data,
    #' protein abundance, etc.) to the adjacency matrix. Analyzed components (e.g. genes) in rows, samples (e.g. patients)
    #' in columns.
    #' @param identifiers [data.frame] Data frame containing  biological identifiers and the corresponding node ID
    #' created in \code{\link[DrDimont]{compute_correlation_matrices}} via \link[DrDimont]{create_unique_layer_node_ids}.
    #' The column containing node IDs has to be named `node_id`.
    #' @param handling_missing_data ["all.obs"|"pairwise.complete.obs"] Specifying the handling
    #' of missing data during correlation matrix computation. (default: all.obs)
    #' @param reduction_method ["pickHardThreshold"|"p_value"] A character string specifying the method to be used for network
    #' reduction. `p_value` for hard thresholding based on the statistical significance of the
    #' computed correlation. `pickHardThreshold` for a cutoff based on the scale-freeness criterion
    #' (calls \code{\link[WGCNA]{pickHardThreshold}}). (default: pickHardThreshold)
    #' @param r_squared_cutoff [float] A number indicating the desired minimum scale free topology fitting index R^2 for reduction
    #' using \code{\link[WGCNA]{pickHardThreshold}}. (default: 0.85)
    #' @param cut_vector [sequence of float] A vector of hard threshold cuts for which the scale free topology fit indices are to
    #' be calculated during reduction with \code{\link[WGCNA]{pickHardThreshold}}. (default: seq(0.2, 0.8, by = 0.01))
    #' @param mean_number_edges [int] Find a suitable edge weight cutoff employing \code{\link[WGCNA]{pickHardThreshold}} to reduce
    #' the network to at most the specified mean number of edges. Attention: This parameter overwrites the 'r_squared_cutoff' and
    #' 'edge_density' parameters if not set to NULL. (default: NULL)
    #' @param edge_density [float] Find a suitable edge weight cutoff employing \code{\link[WGCNA]{pickHardThreshold}} to reduce the
    #' network to at most the specified edge density. Attention: This parameter overwrites the 'r_squared_cutoff' parameter if not set
    #' to NULL. (default: NULL)
    #' @param p_value_adjustment_method ["holm"|"hochberg"|"hommel"|"bonferroni"|"BH"|"BY"|"fdr"|"none"] String
    #' of the correction method applied to p-values. Passed to \link[stats]{p.adjust}. (default: "BH")
    #' @param reduction_alpha [float] A number indicating the significance value for correlation p-values
    #' during reduction. Not-significant edges are dropped. (default: 0.05)
    #' @param n_threads [int] Number of threads for parallel computation of p-values during p-value reduction.
    #' (default: 1)
    #' @param parallel_chunk_size [int] Number of p-values in smallest work unit when computing in parallel
    #' during network reduction with method `p_value`. (default: 10^6)
    #' @param print_graph_info [bool] Specifying if a summary of the reduced graph should be printed to
    #' the console after network generation. (default: TRUE)
    #' 
    #' @return iGraph graph object of the reduced network.
    #' 
    #' @keywords internal
    #' @export

    # network reduction
    if (reduction_method == 'p_value') {

        message(format(Sys.time(), "[%y-%m-%d %X] "), "Reducing network by p-values...")

        # get number of samples
        message(format(Sys.time(), "[%y-%m-%d %X] "), "Computing sample size...")
        number_of_samples <- sample_size(measurement_data, handling_missing_data)
        message(format(Sys.time(), "[%y-%m-%d %X] "), "Sample size computed. Starting p-value computation...")
        reduced_adjacency_matrix <- network_reduction_by_p_value(adjacency_matrix,
                                                                 number_of_samples,
                                                                 p_value_adjustment_method,
                                                                 reduction_alpha,
                                                                 parallel_chunk_size)
    }
    else if (reduction_method == 'pickHardThreshold') {
        reduced_adjacency_matrix <- network_reduction_by_pickHardThreshold(adjacency_matrix,
                                                                           r_squared_cutoff,
                                                                           cut_vector,
                                                                           mean_number_edges,
                                                                           edge_density)
    }
    else {
        message(format(Sys.time(), "[%y-%m-%d %X] "), 'No network reduction.')
        reduced_adjacency_matrix <- adjacency_matrix
    }

    # add node identifiers
    colnames(reduced_adjacency_matrix) <- identifiers$node_id

    # generate iGraph object
    graph <- igraph::graph.adjacency(adjmatrix=reduced_adjacency_matrix,
                                     weighted=TRUE,
                                     diag=FALSE,
                                     mode='undirected')

    graph <- igraph::delete_edges(graph, which(is.na(igraph::E(graph)$weight)))

    if (print_graph_info == TRUE) {graph_metrics(graph)}

    # remove adjacency matrix from environment
    rm(adjacency_matrix)

    return(graph)
}

sample_size <- function(measurement_data, handling_missing_data) {
    #' @title [INTERNAL] Sample size for correlation computation
    #'
    #' @description [INTERNAL] Depending on how missing data is handled in correlation matrix computation,
    #' the number of samples used is returned. If `all.obs` is specified the number of rows (i.e. samples)
    #' of the original data is returned. If `pairwise.complete.obs` is specified the crossproduct of a
    #' matrix indicating the non-NA values is returned as matrix. This implementation was adopted
    #' from \code{\link[WGCNA]{corAndPvalue}}.
    #'
    #' @param measurement_data [data.frame] Data frame containing the respective raw data (e.g. mRNA expression data,
    #' protein abundance, etc.) to the adjacency matrix. Analyzed components (e.g. genes) in rows, samples (e.g. patients)
    #' in columns.
    #' @param handling_missing_data ["all.obs"|"pairwise.complete.obs"] Specifying the handling
    #' of missing data during correlation matrix computation. (default: all.obs)
    #'
    #' @return For 'all.obs' returns an integer indicating the number of samples in the supplied
    #' matrix (i.e. number of rows). For 'pairwise.complete.obs' returns a matrix in the same size
    #' of the correlation matrix indicating the number of samples for each correlation calculation.
    #' 
    #' @source Method to calculate samples in `pairwise.complete.obs` adopted and improved from
    #' \code{\link[WGCNA]{corAndPvalue}}
    #' 
    #' @keywords internal
    #' @export

    if (handling_missing_data == 'all.obs') { n_samples <- dim(measurement_data)[1]
    }
    else if (handling_missing_data == 'pairwise.complete.obs') {
        finMat <- !is.na(measurement_data)
        n_samples <- Rfast::Crossprod(finMat, finMat)
    }

    return(n_samples)
}
