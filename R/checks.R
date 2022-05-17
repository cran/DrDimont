check_layer <- function(layer) {
    #' @title [INTERNAL] Check layer input
    #'
    #' @description [INTERNAL] Checks if the data used to create a network layer is valid and has the right
    #' format
    #'
    #' @param layer [list] Named list of layer to check. Created by \code{\link[DrDimont]{make_layer}}
    #' @return Character string vector containing error messages.
    #' 
    #' @examples
    #' data(protein_data)
    #' protein_layer <- make_layer(
    #'                      name="protein",
    #'                      data_groupA=t(protein_data$groupA[, c(-1,-2)]),
    #'                      data_groupB=t(protein_data$groupB[, c(-1,-2)]),
    #'                      identifiers_groupA=data.frame(gene_name=protein_data$groupA$gene_name,
    #'                                                   ref_seq=protein_data$groupA$ref_seq),
    #'                      identifiers_groupB=data.frame(gene_name=protein_data$groupB$gene_name,
    #'                                                   ref_seq=protein_data$groupB$ref_seq))
    #' return_errors(check_layer(protein_layer))
    #' 
    #' @export
    
    errors <- c()
    for (group in c("groupA", "groupB")) {
    
        ### if group in layer not given skip it
        if (is.null(layer[[group]])){next}
    
        # do length of components and length of identifiers match?
        if (dim(layer[[group]][['data']])[2] != dim(layer[[group]][['identifiers']])[1]) {
            errors <- c(errors, stringr::str_interp("Layer ${layer[['name']]}, ${group}: The number of columns in the supplied data does not equal the number of provided identifiers."))
        }
        # do the identifiers contain an id called 'layer'? layer attribute is assigned in the pipeline
        if ('layer' %in% colnames(layer[[group]][['identifiers']])) {
            errors <- c(errors, stringr::str_interp("Identifiers cannot be named 'layer'. Please choose a different name for the identifier 'layer' in layer ${layer[['name']]}, ${group}."))
        }
        message(format(Sys.time(), "[%y-%m-%d %X] "), stringr::str_interp("Layer ${layer[['name']]}, ${group} contains "), dim(layer[[group]][['data']])[1], " samples and ", dim(layer[[group]][['data']])[2], " genes/proteins/entities.")
    }
    return(errors)
}

check_connection <- function(connection) {
    #' @title [INTERNAL] Check connection
    #'
    #' @description [INTERNAL] Checks if the data given to create an inter-layer connection is valid and has the
    #' right input format
    #'
    #' @param connection [list] Connection to check. Created by \code{\link[DrDimont]{make_connection}}
    #' @return Character string vector containing error messages.
    #' 
    #' @examples
    #' inter_layer_connections = make_connection("mrna", 
    #'                                           "protein", 
    #'                                           connect_on="gene_name")
    #' return_errors(check_connection(inter_layer_connections))
    #' 
    #' @export

    errors <- c()

    # are group arguments 'A', 'B' or 'both'?
    if (as.character(connection$group) != 'A' & as.character(connection$group) != 'B' & connection$group != 'both') {
        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The 'group' argument has to be 'A', 'B' or 'both'. It was set to ${connection$group}"))
    }

    # is argument 'from' a string and length 1?
    if (!is.character(connection$from) | length(connection$from) != 1) {
        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The 'from' argument has to be a character string and not contain more than one element."))
    }

    # is argument 'to' a string and length 1?
    if (!is.character(connection$to) | length(connection$to) > 1) {
        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The 'to' argument has to be a character string and not contain more than one element."))
    }

    if (connection$by == "id") {
        if (!length(connection$connect_on) == 1) {
            errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}):
                               The argument 'connect_on' is a character string and is
                               therefore expected to contain the name of exactly
                               one identifier, but length was not 1.
                               'connect_on' was ${connection$connect_on}."))
        }
        if (!is.numeric(connection$weight)) {
            errors <- c(errors, stringr::str_interp(
            "Connection (${connection$from}-${connection$to}): 
             If 'connect_on' is a character string, 'weight' has to be numeric."))
        }
    }
    else if (connection$by == "table") {

        if (!is.character(connection$weight) || !is.vector(connection$weight) || !length(connection$weight) == 1) {
            errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}):
                             The argument 'weight' is expected to contain the
                             name of exactly one column in the table passed as
                             'connect_on'.
                             'weight' was ${connection$weight}."))
        }

        tryCatch(colnames(connection$connect_on),
                 error = function(e) errors <- c(errors, "Couldn't get column names of 'connect_on'.
                                      Argument seems malformed. Make sure 'connect_on'
                                      is a character string or a table.")
        )

        if (!(connection$weight %in% colnames(connection$connect_on))) {
            errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}):
                             The argument 'weight' is expected to contain the
                             name of exactly one column in the table passed as
                             'connect_on'.
                             'weight' was ${connection$weight}.
                              Column names were ${toString(colnames(connection$connect_on))}"))
        }

        if(!is.numeric(as.matrix(connection$connect_on[ , connection$weight]))) {
            errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The column specified in 'weight' does not contain numeric data. Please provide a numeric column."))
        }
    }

    return(errors)
}

check_sensible_connections <- function(connection, layers) {
    #' @title [INTERNAL] Check connection and layer data
    #'
    #' @description [INTERNAL] Checks if the connection defined in 'connection' makes sense in
    #' context of the defined layers.
    #'
    #' @param connection [list] Connection to check. Created by \code{\link[DrDimont]{make_connection}}
    #' @param layers [list] List of layers to check. Individual layers are created by
    #' \code{\link[DrDimont]{make_layer}} and need to be wrapped in a list.
    #'
    #' @return Character string vector containing error messages.
    #' 
    #' @examples
    #' data(mrna_data)
    #' data(protein_data)
    #' 
    #' mrna_layer <- make_layer(
    #'                     name="mrna",
    #'                     data_groupA=t(mrna_data$groupA[,-1]),
    #'                     data_groupB=t(mrna_data$groupB[,-1]),
    #'                     identifiers_groupA=data.frame(gene_name=mrna_data$groupA$gene_name),
    #'                     identifiers_groupB=data.frame(gene_name=mrna_data$groupB$gene_name))
    #' 
    #' protein_layer <- make_layer(
    #'                      name="protein",
    #'                      data_groupA=t(protein_data$groupA[, c(-1,-2)]),
    #'                      data_groupB=t(protein_data$groupB[, c(-1,-2)]),
    #'                      identifiers_groupA=data.frame(gene_name=protein_data$groupA$gene_name,
    #'                                                   ref_seq=protein_data$groupA$ref_seq),
    #'                      identifiers_groupB=data.frame(gene_name=protein_data$groupB$gene_name,
    #'                                                   ref_seq=protein_data$groupB$ref_seq))
    #' 
    #' inter_layer_connections = make_connection("mrna", 
    #'                                           "protein", 
    #'                                           connect_on="gene_name")
    #' return_errors(check_sensible_connections(inter_layer_connections, 
    #'                                          layers=list(mrna_layer, 
    #'                                                      protein_layer)))
    #'
    #' @export
    
    errors <- c()
    layer_names <- c()
    for (layer in layers) {
        layer_names <- c(layer_names, layer[['name']])
    }

    if (!connection$from %in% layer_names) {
        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The layer given in the 'from' argument cannot be found in the names of the created layers. 'from' was ${connection$from}. Layer names were: ${toString(layer_names)}."))
    }

    if (!connection$to %in% layer_names) {
        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The layer given in the 'to' argument cannot be found in the names of the created layers. 'from' was ${connection$to}. Layer names were: ${toString(layer_names)}."))
    }

    if (connection$to %in% layer_names & connection$from %in% layer_names) {
        if (connection$by == 'id') {
            for (group in c("groupA", "groupB")) {
                for (layer in c(connection$from, connection$to)) {
                    
                    ### if group not given for layer then skip it
                    if (is.null(layers[[which(layer_names == layer)]][[group]])){next}
                
                    identifier_cols <- colnames(layers[[which(layer_names == layer)]][[group]][['identifiers']])
                    if (!connection$connect_on %in% identifier_cols) {
                        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): The identifier specified in 'connect_on' was not found in the identifiers of layer ${layer} of ${group}. Please correct the 'connect_on' argument or the names of the layer identifiers. 'connect_on' was ${connection$connect_on}. Identifiers were: ${toString(identifier_cols)}"))
                    }
                }
            }
        }

        if (connection$by == 'table') {
            for (group in c("groupA", "groupB")) {
                for (layer in c(connection$from, connection$to)){
   
                    ### if group not given for layer then skip it
                    if (is.null(layers[[which(layer_names == layer)]][[group]])){next}
                    
                    identifier_cols <- colnames(layers[[which(layer_names == layer)]][[group]][['identifiers']])
                    if (!any(colnames(connection$connect_on) %in% identifier_cols)) {
                        errors <- c(errors, stringr::str_interp("Connection (${connection$from}-${connection$to}): None of the columnnames supplied in 'connect_on' were found in the identifiers of layer ${layer} of ${group}. Please correct the columnnames of the table passed to 'connect_on' or the names of the layer identifiers. Columnames in 'connect_on' were ${toString(colnames(connection$connect_on))}. Identifiers were: ${toString(identifier_cols)}"))
                    }

                }
            }
        }
    }
    return(errors)
}

check_drug_target <- function(drug_target_interactions) {
    #' @title [INTERNAL] Check drug target interaction data
    #'
    #' @description [INTERNAL] Checks if the data used to define interaction between drugs and
    #' targets is valid and formatted correctly.
    #'
    #' @param drug_target_interactions [list] A named list of the drug interaction data. Created by
    #' \code{\link[DrDimont]{make_drug_target}}
    #'
    #' @return Character string vector containing error messages.
    #' @examples
    #' data(drug_gene_interactions)
    #' drug_target_interactions <- make_drug_target(
    #'                                  target_molecules='protein',
    #'                                  interaction_table=drug_gene_interactions,
    #'                                  match_on='gene_name')
    #' return_errors(check_drug_target(drug_target_interactions))
    #' 
    #' @export

    errors <- c()
    if(!drug_target_interactions$match_on %in% colnames(drug_target_interactions$interaction_table)) {
        errors <- c(errors, stringr::str_interp("Drug-target interaction: The columnname specified in 'match_on' cannot be found in the columnnames of the table supplied in 'interaction_table'. 'match_on' is ${drug_target_interactions$match_on}, columnnames in 'interaction_table' are ${toString(colnames(drug_target_interactions$interaction_table))}"))
    }

    return(errors)
}

check_drug_targets_in_layers <- function(drug_target_interactions, layers) {
    #' @title [INTERNAL] Check drug target and layer data
    #'
    #' @description [INTERNAL] Checks if the parameters supplied in 'drug_target_interactions' makes
    #' sense in the context of the defined layers.
    #' 
    #' @param drug_target_interactions [list] A named list of the drug interaction data. Created by
    #' \code{\link[DrDimont]{make_drug_target}}
    #' @param layers [list] List of layers to check. Individual layers are created by
    #' \code{\link[DrDimont]{make_layer}} and need to be wrapped in a list.
    #'
    #' @return Character string vector containing error messages.
    #' 
    #' @examples
    #' data(layers_example)
    #' data(drug_gene_interactions)
    #' drug_target_interactions <- make_drug_target(
    #'                                  target_molecules='protein',
    #'                                  interaction_table=drug_gene_interactions,
    #'                                  match_on='gene_name')
    #' return_errors(check_drug_targets_in_layers(drug_target_interactions, layers_example))
    #' 
    #' @export
    
    layer_names <- c()
    for (layer in layers) {
        layer_names <- c(layer_names, layer[['name']])
    }

    errors <- c()

    if (!drug_target_interactions$target_molecules %in% layer_names) {
        errors <- c(errors, stringr::str_interp("Drug-target interaction: The defined target was not found in the list of layers. Targets molecules are ${drug_target_interactions$target_molecules}, layers are ${toString(layer_names)}. Please correct target molecule string or the layer names."))
    } else {
        for (group in c("groupA", "groupB")) {
            
            ### if group not given for layer then skip it
            if (is.null(layers[[which(layer_names == drug_target_interactions$target_molecules)]][[group]])){next}
        
            identifier_cols <- colnames(layers[[which(layer_names == drug_target_interactions$target_molecules)]][[group]][['identifiers']])
            if(!drug_target_interactions$match_on %in% identifier_cols) {
                errors <- c(errors, stringr::str_interp("Drug-target interaction: The columnname specified in 'match_on' cannot be found in the columnnames of the identifiers of the target layer. 'match_on' is ${drug_target_interactions$match_on}, columnnames in the identifiers of target layer ${drug_target_interactions$target_molecules} of ${group} are ${toString(identifier_cols)}"))
            }
        }
    }
    return(errors)
}

check_input <- function(layers, inter_layer_connections, drug_target_interactions) {
    #' @title Check pipeline input data for required format
    #'
    #' @description Checks if input data is valid and formatted correctly. This function is a
    #' wrapper for other check functions to be executed as first step of the DrDimont pipeline.
    #'
    #' @param layers [list] List of layers to check. Individual layers were created by
    #' \code{\link[DrDimont]{make_layer}} and need to be wrapped in a list.
    #' @param inter_layer_connections [list] A list containing connections between layers. Each
    #' connection was created by \code{\link[DrDimont]{make_connection}} and wrapped in a list.
    #' @param drug_target_interactions [list] A named list of the drug interaction data. Created by
    #' \code{\link[DrDimont]{make_drug_target}}
    #'
    #' @return Character string vector containing error messages.
    #' 
    #' @export

    errors <- c()
    # check layers
    n <- length(layers)
    for (i in 1:n) {
        errors <- c(errors, check_layer(layers[[i]]))
    }

    # check inter-layer connections
    for (connection in inter_layer_connections) {
        errors <- c(errors, check_connection(connection))
        errors <- c(errors, check_sensible_connections(connection, layers))
    }

    # check drug target
    errors <- c(errors, check_drug_target(drug_target_interactions))
    errors <- c(errors, check_drug_targets_in_layers(drug_target_interactions, layers))

    return(errors)
}

return_errors <- function(errors) {
    #' @title Return detected errors
    #'
    #' @description Throws an error in case errors have been passed to the function. Messages
    #' describing the detected errors are printed.
    #'
    #' @param errors [string] Character string vector containing error messages.
    #' 
    #' @return No return value, writes error messages to console
    #' 
    #' @examples
    #' layer <- DrDimont::layers_example[[2]]
    #' return_errors(check_layer(layer))
    #'
    #' @export

    if (!is.null(errors)) {
        message(format(Sys.time(), "[%y-%m-%d %X] "), stringr::str_interp("\n ----- \n${length(errors)} Error(s) detected:\n ----- \n"))
        message(format(Sys.time(), "[%y-%m-%d %X] "), paste(paste0(seq(1, length(errors), 1), ". ", errors), "\n", collapse = "\n"))
		message(format(Sys.time(), "[%y-%m-%d %X] "),"ERRORs detected! Details are printed above this error message.\n")
        stop("Errors detected! Details are printed above this error message.\n")
    }
}
