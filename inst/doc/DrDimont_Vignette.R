## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(DrDimont)

## ---- echo=TRUE, warning=FALSE, eval=FALSE------------------------------------
#  install_python_dependencies(package_manager="pip")

## ---- echo=TRUE, warning=FALSE, eval=FALSE------------------------------------
#  install_python_dependencies(package_manager="conda")

## ----Load data----------------------------------------------------------------
data("mrna_data")
data("protein_data")
data("phosphosite_data")
data("metabolite_data")
data("metabolite_protein_interactions")
data("drug_gene_interactions")

## -----------------------------------------------------------------------------
# Data inspection
mrna_data$groupA[1:3, 1:5]
protein_data$groupA[1:3, 1:5]
phosphosite_data$groupA[1:3, 1:5]
metabolite_data$groupA[1:3, 1:5]


## ----Create layers------------------------------------------------------------
# Create individual layers
mrna_layer <- make_layer(name="mrna",
                         data_groupA=t(mrna_data$groupA[,-1]),
                         data_groupB=t(mrna_data$groupB[,-1]),
                         identifiers_groupA=data.frame(gene_name=mrna_data$groupA$gene_name),
                         identifiers_groupB=data.frame(gene_name=mrna_data$groupB$gene_name))

protein_layer <- make_layer(name="protein",
                            data_groupA=t(protein_data$groupA[, c(-1,-2)]),
                            data_groupB=t(protein_data$groupB[, c(-1,-2)]),
                            identifiers_groupA=data.frame(gene_name=protein_data$groupA$gene_name, 
                                                         ref_seq=protein_data$groupA$ref_seq),
                            identifiers_groupB=data.frame(gene_name=protein_data$groupB$gene_name, 
                                                         ref_seq=protein_data$groupB$ref_seq))

phosphosite_layer <- make_layer(name="phosphosite",
                                data_groupA=t(phosphosite_data$groupA[, c(-1,-2, -3)]),
                                data_groupB=t(phosphosite_data$groupB[, c(-1,-2, -3)]),
                                identifiers_groupA=data.frame(phosphosite_data$groupA[, 1:3]),
                                identifiers_groupB=data.frame(phosphosite_data$groupB[, 1:3]))

metabolite_layer <- make_layer(name="metabolite",
                               data_groupA=t(metabolite_data$groupA[, c(-1,-2, -3)]),
                               data_groupB=t(metabolite_data$groupB[, c(-1,-2, -3)]),
                               identifiers_groupA=data.frame(metabolite_data$groupA[, 1:3]),
                               identifiers_groupB=data.frame(metabolite_data$groupB[, 1:3]))


## ----Make layers list---------------------------------------------------------
all_layers <- list(mrna_layer, protein_layer, phosphosite_layer, metabolite_layer)

## ----eval=FALSE---------------------------------------------------------------
#  # (i) make inter-layer connection
#  make_connection(from='mrna', to='protein', connect_on='gene_name', weight=1)

## -----------------------------------------------------------------------------
# Data inspection
metabolite_protein_interactions[1:3, ]

## ----eval=FALSE---------------------------------------------------------------
#  # (ii) make inter-layer connection
#  make_connection(from='protein', to='metabolite', connect_on=metabolite_protein_interactions, weight='combined_score')

## ----Inter-layer connections--------------------------------------------------
all_inter_layer_connections = list(
    make_connection(from='mrna', to='protein', connect_on='gene_name', weight=1),
    make_connection(from='protein', to='phosphosite', connect_on='gene_name', weight=1),
    make_connection(from='protein', to='metabolite', 
                    connect_on=metabolite_protein_interactions, weight='combined_score')
)

## -----------------------------------------------------------------------------
# Data inspection
drug_gene_interactions[1:3, ]

## ----Make drug-target interaction---------------------------------------------
all_drug_target_interactions <- make_drug_target(
                                        target_molecules='protein',
                                        interaction_table=drug_gene_interactions,
                                        match_on='gene_name')

## -----------------------------------------------------------------------------
return_errors(check_input(all_layers, all_inter_layer_connections, all_drug_target_interactions))

## ----Settings-----------------------------------------------------------------
example_settings <- drdimont_settings(
                            handling_missing_data = list(
                                default = "pairwise.complete.obs",
                                mrna = "all.obs"
                            ),
                            reduction_method = "pickHardThreshold",
                            r_squared=list(default=0.65, metabolite=0.1),
                            cut_vector=list(default=seq(0.2, 0.65, 0.01)),
                            save_data = FALSE,
                            python_executable = "python3",
                            saving_path = tempdir())
# disable multi-threading for example run; 
# not recommended for actual data processing
WGCNA::disableWGCNAThreads()


## ----Run pipeline, eval=FALSE-------------------------------------------------
#  run_pipeline(layers=all_layers,
#               inter_layer_connections=all_inter_layer_connections,
#               drug_target_interactions=all_drug_target_interactions,
#               settings=example_settings)

## ----Correlation matrices, message=FALSE, results='hide'----------------------
reduced_mrna_layer <- make_layer(name="mrna",
                           data_groupA=t(mrna_data$groupA[1:10,2:11]),
                           data_groupB=t(mrna_data$groupB[1:10,2:11]),
                           identifiers_groupA=data.frame(gene_name=mrna_data$groupA$gene_name[1:10]),
                           identifiers_groupB=data.frame(gene_name=mrna_data$groupB$gene_name[1:10]))

example_correlation_matrices <- compute_correlation_matrices(
                                            layers=list(reduced_mrna_layer), 
                                            settings=example_settings)

## -----------------------------------------------------------------------------
# Data inspection
data("correlation_matrices_example")
correlation_matrices_example$annotations$groupA$protein[1:3, ]

## ----Individual graphs, message=FALSE, results='hide'-------------------------
data("correlation_matrices_example")
example_individual_graphs <- generate_individual_graphs(
                                        correlation_matrices=correlation_matrices_example, 
                                        layers=all_layers, 
                                        settings=example_settings)

## ----Combine graphs-----------------------------------------------------------
example_combined_graphs <- generate_combined_graphs(
                                    graphs=example_individual_graphs[["graphs"]], 
                                    annotations=example_individual_graphs[["annotations"]], 
                                    inter_layer_connections=all_inter_layer_connections, 
                                    settings=example_settings)

## -----------------------------------------------------------------------------
# Data inspection
example_combined_graphs$annotations$both[1:3, ]

## ----Drug targets and their edges---------------------------------------------
example_drug_target_edges <- determine_drug_targets(
                                        graphs=example_combined_graphs[["graphs"]], 
                                        annotations=example_combined_graphs[["annotations"]], 
                                        drug_target_interactions=all_drug_target_interactions, 
                                        settings=example_settings)

## ----Calculate interaction score, eval=FALSE----------------------------------
#  example_interaction_score_graphs <- generate_interaction_score_graphs(
#                                              graphs=example_combined_graphs[["graphs"]],
#                                              drug_target_edgelists=drug_targets[["edgelists"]],
#                                              settings=example_settings)

## ----Calculate differential score---------------------------------------------
data("interaction_score_graphs_example")
example_differential_graph <- generate_differential_score_graph(
                                            interaction_score_graphs=interaction_score_graphs_example, 
                                            settings=example_settings)
#if interaction score graphs were computed use the following:
#example_differential_score_graph <- generate_differential_score_graph(
#                                           interaction_score_graphs=example_interaction_score_graphs, 
#                                           settings=example_settings)

## ----Drug response------------------------------------------------------------
example_drug_response_scores <- compute_drug_response_scores(
                                                     differential_graph=example_differential_graph,
                                                     drug_targets=example_drug_target_edges[["targets"]],
                                                     settings=example_settings)

## ----Result Output------------------------------------------------------------
head(dplyr::filter(example_drug_response_scores, !is.na(drug_response_score)))

