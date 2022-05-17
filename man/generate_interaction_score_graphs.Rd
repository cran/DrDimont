% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{generate_interaction_score_graphs}
\alias{generate_interaction_score_graphs}
\title{Computes interaction score for combined graphs}
\usage{
generate_interaction_score_graphs(graphs, drug_target_edgelists, settings)
}
\arguments{
\item{graphs}{[list] A named list with elements `groupA` and `groupB` containing the combined graphs
of each group as iGraph object (`graphs` from output of \code{\link[DrDimont]{generate_combined_graphs}})}

\item{drug_target_edgelists}{[list] A named list (elements `groupA` and `groupB`). Each element
contains the list of edges adjacent to drug targets as a dataframe (columns `from`, `to` and
`weight`). `edgelists` from output of \code{\link[DrDimont]{determine_drug_targets}}}

\item{settings}{[list] A named list containing pipeline settings. The settings list has to be
initialized by \code{\link[DrDimont]{drdimont_settings}}. Items in the named list can be
adjusted as desired.}
}
\value{
A named list (elements `groupA` and `groupB`). Each element contains an iGraph object
containing the interaction scores as interaction_weight attributes.
}
\description{
Writes the input data (combined graphs for both groups in `gml` format and
lists of edges adjacent to drug targets for both groups) to files and calls a python script
for calculating the score. Output files written by the python script are two graphs in `gml`
format containing the interaction score as additional interaction_weight edge attribute.
These are loaded and returned in a named list.
ATTENTION: Data exchange via files is mandatory and takes a long time for large data. Interaction
score computation is expensive and slow because it involves finding all simple paths up to a
certain length between source and target node of the drug target edges. Don't set `max_path_length`
in settings to a large value and only consider this step if your graphs have up to approximately
2 million edges. Computation is initiated by \code{\link[DrDimont]{calculate_interaction_score}}.
The python script is parallelized using Ray. Use the setting `int_score_mode` to force sequential
or parallel computation. Refer to the Ray documentation if you encounter problems with running
the python script in parallel. DISCLAIMER: Depending on the operating system Python comes
pre-installed or has to be installed manually.
Please pay attention to the version and the executable used (python/python3 or homebrew
python). You can use the `python_executable` setting to specify the command or path.
}
\examples{
data(combined_graphs_example)
data(drug_target_edges_example)

example_settings <- drdimont_settings(
                       save_data=FALSE,
                       python_executable="python")

\donttest{
interaction_score_graphs <- generate_interaction_score_graphs(
                                         graphs=combined_graphs_example$graphs,
                                         drug_target_edgelists=drug_target_edges_example$edgelists,
                                         settings=example_settings)
}

}