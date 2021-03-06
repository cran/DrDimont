% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{check_sensible_connections}
\alias{check_sensible_connections}
\title{[INTERNAL] Check connection and layer data}
\usage{
check_sensible_connections(connection, layers)
}
\arguments{
\item{connection}{[list] Connection to check. Created by \code{\link[DrDimont]{make_connection}}}

\item{layers}{[list] List of layers to check. Individual layers are created by
\code{\link[DrDimont]{make_layer}} and need to be wrapped in a list.}
}
\value{
Character string vector containing error messages.
}
\description{
[INTERNAL] Checks if the connection defined in 'connection' makes sense in
context of the defined layers.
}
\examples{
data(mrna_data)
data(protein_data)

mrna_layer <- make_layer(
                    name="mrna",
                    data_groupA=t(mrna_data$groupA[,-1]),
                    data_groupB=t(mrna_data$groupB[,-1]),
                    identifiers_groupA=data.frame(gene_name=mrna_data$groupA$gene_name),
                    identifiers_groupB=data.frame(gene_name=mrna_data$groupB$gene_name))

protein_layer <- make_layer(
                     name="protein",
                     data_groupA=t(protein_data$groupA[, c(-1,-2)]),
                     data_groupB=t(protein_data$groupB[, c(-1,-2)]),
                     identifiers_groupA=data.frame(gene_name=protein_data$groupA$gene_name,
                                                  ref_seq=protein_data$groupA$ref_seq),
                     identifiers_groupB=data.frame(gene_name=protein_data$groupB$gene_name,
                                                  ref_seq=protein_data$groupB$ref_seq))

inter_layer_connections = make_connection("mrna", 
                                          "protein", 
                                          connect_on="gene_name")
return_errors(check_sensible_connections(inter_layer_connections, 
                                         layers=list(mrna_layer, 
                                                     protein_layer)))

}
