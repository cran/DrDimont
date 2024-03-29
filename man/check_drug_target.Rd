% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{check_drug_target}
\alias{check_drug_target}
\title{[INTERNAL] Check drug target interaction data}
\usage{
check_drug_target(drug_target_interactions)
}
\arguments{
\item{drug_target_interactions}{[list] A named list of the drug interaction data. Created by
\code{\link[DrDimont]{make_drug_target}}}
}
\value{
Character string vector containing error messages.
}
\description{
[INTERNAL] Checks if the data used to define interaction between drugs and
targets is valid and formatted correctly.
}
\keyword{internal}
