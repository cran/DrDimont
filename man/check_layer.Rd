% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{check_layer}
\alias{check_layer}
\title{[INTERNAL] Check layer input}
\usage{
check_layer(layer)
}
\arguments{
\item{layer}{[list] Named list of layer to check. Created by \code{\link[DrDimont]{make_layer}}}
}
\value{
Character string vector containing error messages.
}
\description{
[INTERNAL] Checks if the data used to create a network layer is valid and has the right
format
}
\examples{
data(protein_data)
protein_layer <- make_layer(
                     name="protein",
                     data_groupA=t(protein_data$groupA[, c(-1,-2)]),
                     data_groupB=t(protein_data$groupB[, c(-1,-2)]),
                     identifiers_groupA=data.frame(gene_name=protein_data$groupA$gene_name,
                                                  ref_seq=protein_data$groupA$ref_seq),
                     identifiers_groupB=data.frame(gene_name=protein_data$groupB$gene_name,
                                                  ref_seq=protein_data$groupB$ref_seq))
return_errors(check_layer(protein_layer))

}
