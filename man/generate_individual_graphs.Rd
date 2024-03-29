% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{generate_individual_graphs}
\alias{generate_individual_graphs}
\title{Builds graphs from specified network layers}
\usage{
generate_individual_graphs(correlation_matrices, layers, settings)
}
\arguments{
\item{correlation_matrices}{[list] List of correlation matrices generated with
code{\link[DrDimont]{compute_correlation_matrices}}}

\item{layers}{[list] Named list with different network layers containing data and
identifiers for both groups (generated from \code{\link[DrDimont]{make_layer}})}

\item{settings}{[list] A named list containing pipeline settings. The settings list has to be
initialized by \code{\link[DrDimont]{drdimont_settings}}. Items in the named list can be
adjusted as desired.}
}
\value{
A nested named list with first-level elements `graphs` and `annotations`. The second
level elements are `groupA` and `groupB` (and `both` at `annotations`). These contain a list of
iGraph objects (`graphs`) and data frames (`annotations`) mapping the graph node IDs to biological
identifiers. The third level elements are layer names given by the user.
}
\description{
Constructs and returns two graphs for each network layer, where nodes
correspond to the rows in the measurement data. Graphs are initially complete and
edges are weighted by correlation values of the measurements across columns. The
number of edges is then reduced by either a threshold on the p-value of the
correlation or a minimum scale-free fit index.
}
\examples{
\dontshow{
WGCNA::disableWGCNAThreads()
}

data(layers_example)
data(correlation_matrices_example)

example_settings <- drdimont_settings(
                        handling_missing_data=list(
                            default="pairwise.complete.obs",
                            mrna="all.obs"),
                        reduction_method="pickHardThreshold",
                        r_squared=list(default=0.65, metabolite=0.1),
                        cut_vector=list(default=seq(0.2, 0.5, 0.01)))

example_individual_graphs <- generate_individual_graphs(
                                 correlation_matrices=correlation_matrices_example,
                                 layers=layers_example, 
                                 settings=example_settings)

graph_metrics(example_individual_graphs$graphs$groupA$mrna)
graph_metrics(example_individual_graphs$graphs$groupB$mrna)

}
