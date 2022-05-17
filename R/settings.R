drdimont_settings <- function(
                        ### saving
                        saving_path = 'tempdir()',
                        save_data = FALSE,
                        ### network generation
                        correlation_method = "spearman",
                        handling_missing_data = "all.obs",
                        ### network reduction
                        reduction_method = "pickHardThreshold",
                        ### pickHardThreshold
                        r_squared_cutoff = 0.85,
                        cut_vector = seq(0.2, 0.8, by = 0.01),
                        mean_number_edges = NULL,
                        edge_density = NULL,
                        ### p-value
                        p_value_adjustment_method = "BH",
                        reduction_alpha = 0.05,
                        n_threads = 1,
                        parallel_chunk_size = 10^6,
                        print_graph_info = TRUE,
                        ### interaction_score
                        python_executable = "python",
                        conda = FALSE,
                        max_path_length = 3,
                        int_score_mode = "auto",
                        cluster_address = "None",
                        ### drug response score
                        median_drug_response=FALSE,
                        absolute_difference=FALSE,
                        ...){

    #' @title Create global settings variable for DrDimont pipeline
    #'
    #' @description Function that allows creating a global `settings` variable used in the
    #' \code{\link[DrDimont]{run_pipeline}} function and the step-wise DrDimont execution.
    #' Default parameters can be changed within the function call.
    #'
    #' @param correlation_method ["pearson"|"spearman"|"kendall"] Correlation method used
    #' for graph generation. Argument is passed to \code{\link[WGCNA]{cor}}. (default: spearman)
    #' @param handling_missing_data ["all.obs"|"pairwise.complete.obs"] Specifying the handling
    #' of missing data during correlation matrix computation. Argument is passed to \code{\link[WGCNA]{cor}}.
    #' Can be a single character string if the same for all layers, else a named list mapping layer names
    #' to methods, e.g, \code{handling_missing_data=list(mrna="all.obs", protein="pairwise.complete.obs")}.
    #' Layers may be omitted if a method is mapped to `default`, e.g,
    #' handling_missing_data=list(default="pairwise.complete.obs"). (default: all.obs)
    #'
    #' @param print_graph_info [bool] Specifying if a summary of the reduced graph should be printed to
    #' the console after network generation. (default: TRUE)
    #'
    #' @param reduction_method ["pickHardThreshold"|"p_value"]
    #' Reduction method for reducing networks. `p_value` for hard thresholding based on the statistical
    #' significance of the computed correlation. `pickHardThreshold` for a cutoff based on the scale-freeness
    #' criterion (calls \code{\link[WGCNA]{pickHardThreshold}}). Can be a single character string if the same
    #' for all layers, else a named list mapping layer names to methods (see \code{handling_missing_data} setting).
    #' Layers may be omitted if a method is mapped to `default`. (default: pickHardThreshold)
    #'
    #' @param r_squared_cutoff pickHardThreshold setting: [float|named list]
    #' A number indicating the desired minimum scale free topology fitting index R^2 for reduction using
    #' \code{\link[WGCNA]{pickHardThreshold}}.
    #' Can be a single float number if the same for all layers, else a named list mapping layer names to a cutoff
    #' (see \code{handling_missing_data} setting) or a named list in a named list mapping groupA or groupB and layer
    #' names to a cutoff, e.g., \code{r_squared_cutoff=list(groupA=list(mrna=0.85, protein=0.8), groupB=list(mrna=0.9, protein=0.85))}.
    #' Layers/groups may be omitted if a cutoff is mapped to `default`. (default: 0.85)
    #' @param cut_vector pickHardThreshold setting: [sequence of float|named list]
    #' A vector of hard threshold cuts for which the scale free topology fit indices are to be calculated during
    #' reduction with \code{\link[WGCNA]{pickHardThreshold}}.
    #' Can be a single regular sequence if the same for all layers, else a
    #' named list mapping layer names to a cut vector or a named list in a named list mapping groupA or groupB and layer
    #' names to a cut vector (see \code{r_squared_cutoff} setting). Layers/groups may be omitted if a vector is mapped
    #' to `default`. (default: seq(0.2, 0.8, by = 0.01))
    #' @param mean_number_edges pickHardThreshold setting: [int|named list] Find a suitable edge weight cutoff employing
    #' \code{\link[WGCNA]{pickHardThreshold}} to reduce the network to at most the specified mean number of edges.
    #' Can be a single int number if the same for all layers, else a named list mapping layer names to a mean number of edges or
    #' a named list in a named list mapping groupA or groupB and layer names to a cutoff (see \code{r_squared_cutoff} setting).
    #' Attention: This parameter overwrites the 'r_squared_cutoff' and 'edge_density' parameters if not set to NULL. (default: NULL)
    #' @param edge_density pickHardThreshold setting: [float|named list] Find a suitable edge weight cutoff employing
    #' \code{\link[WGCNA]{pickHardThreshold}} to reduce the network to at most the specified edge density.
    #' Can be a single float number if the same for all layers, else a named list mapping layer names to a mean number of edges or
    #' a named list in a named list mapping groupA or groupB and layer names to a cutoff (see \code{r_squared_cutoff} setting).
    #' Attention: This parameter overwrites the 'r_squared_cutoff' parameter if not set to NULL. (default: NULL)
    #'
    #' @param p_value_adjustment_method p_value setting: ["holm"|"hochberg"|"hommel"|"bonferroni"|"BH"|"BY"|"fdr"|"none"] String
    #' of the correction method applied to p-values. Passed to \link[stats]{p.adjust}. (default: "BH")
    #' @param reduction_alpha p_value setting: [float] A number indicating the significance value for correlation p-values
    #' during reduction. Not-significant edges are dropped. (default: 0.05)
    #' @param n_threads p_value setting: [int] Number of threads for parallel computation of p-values during p-value reduction.
    #' (default: 1)
    #' @param parallel_chunk_size p_value setting: [int] Number of p-values in smallest work unit when computing in parallel
    #' during network reduction with method `p_value`. (default: 10^6)
    #'
    #' @param python_executable [string] Path to Python executable used for generating the interaction score graphs.
    #' (default: "python")
    #' @param conda [bool] Specifying if python is installed in a conda environment. Set TRUE if python is installed
    #' with conda. Use \code{python_executable="-n name-of-your-environment python"} (change name-of-your-environment to
    #' your environment) or \code{python_executable="python"} if installed in base environment. (default: FALSE)
    #' @param max_path_length [int] Integer of maximum length of simple paths to include in the
    #' \code{\link[DrDimont]{generate_interaction_score_graphs}} computation. (default: 3)
    #' @param int_score_mode ["auto"|"sequential"|"ray"] Whether to compute interaction
    #' score in parallel using the Ray python library or sequentially. When `auto` it depends on the
    #' graph sizes. (default: "auto")
    #' @param cluster_address [string] Local node IP-address of Ray if executed on a cluster.
    #' On a cluster: Start ray with \code{ray start --head --num-cpus 32} on the console before DrDimont execution.
    #' It should work with "auto", if it does not specifiy IP-address given by the \code{ray start} command. (default: "auto")
    #'
    #' @param median_drug_response [bool] Specifying if the median instead of the mean of a drug's targets differential
    #' scores should be computed (default: FALSE)
    #' @param absolute_difference [bool] Specifying if the absoulte differential scores instead of the actual differential
    #' scores should be used for drug response computation (default: FALSE)
    #'
    #' @param saving_path [string] Path to save intermediate output of DrDimont's functions. Default is current working directory.
    #' @param save_data [bool] Specifying if intermediate data such as correlation_matrices, individual_graphs, etc.
    #' should be saved during \code{\link[DrDimont]{run_pipeline}}. (default: TRUE)
    #' @param ... Supply additional settings.
    #'
    #' @return Named list of settings
    #'
    #' @examples
    #' settings <- drdimont_settings(
    #'                      correlation_method="spearman",
    #'                      handling_missing_data=list(
    #'                                default="pairwise.complete.obs",
    #'                                mrna="all.obs"),
    #'                      reduction_method="pickHardThreshold",
    #'                      max_path_length=3)
    #' 
    #' @export
    
    settings <- c(as.list(environment()), list(...))
    return(settings)

}

get_layer_setting <- function(layer, group, settings, setting_name) {
    #' @title [INTERNAL] Get layer (and group) settings
    #'
    #' @description Returns specified setting for a specific network layer (and group).
    #'
    #' @param layer [list] A network layer created by \code{\link[DrDimont]{make_layer}}
    #' @param group [string] A network group
    #' @param settings [list] Named list of settings created by \code{\link[DrDimont]{drdimont_settings}}
    #' @param setting_name [string] String indicating the setting to return.
    #' @return Setting value(s) for this layer (and group)
    #' 
    #' @export
    
    if (!is.list(settings[[setting_name]])) {return(settings[[setting_name]])}
    else if (!is.null(settings[[setting_name]][[group]][[layer]])) {return(settings[[setting_name]][[group]][[layer]])}
    else if (!is.null(settings[[setting_name]][[layer]])) {return(settings[[setting_name]][[layer]])}
    else if (!is.null(settings[[setting_name]][[group]][['default']])) {return(settings[[setting_name]][[group]][['default']])}
    else if (!is.null(settings[[setting_name]][['default']])){return(settings[[setting_name]][['default']])}
    else {stop(stringr::str_interp("Neither was a setting ${setting_name} given for layer ${layer} nor any default."))}
}
