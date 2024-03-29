% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_documentations.R
\docType{data}
\name{metabolite_data}
\alias{metabolite_data}
\title{Metabolomics data}
\format{
\describe{
  \item{groupA}{ER+ data; data.frame: first three columns contain metabolite identifiers
                biochemical_name, metabolon_id and pubchem_id;
                other columns are samples containing the quantified metabolite data per metabolite}
  \item{groupB}{ER- data; data.frame: first three columns contain metabolite identifiers
                biochemical_name, metabolon_id and pubchem_id;
                other columns are samples containing the quantified metabolite data per metabolite}
}
}
\source{
Terunuma, Atsushi et al. “MYC-driven accumulation of 2-hydroxyglutarate is associated
with breast cancer prognosis.”
The Journal of clinical investigation vol. 124,1 (2014): 398-412. doi:10.1172/JCI71180

\url{https://www.metabolon.com}

Pubchem IDs: \url{https://pubchem.ncbi.nlm.nih.gov}

MetaboAnalyst: \url{https://www.metaboanalyst.ca/faces/upload/ConvertView.xhtml}
}
\usage{
metabolite_data
}
\description{
Metabolomics analysis of breast cancer patients data sampled randomly to
generate distributions similar to those reported (e.g., in Terunuma et al. (2014)).
The data is stratified by estrogen receptor (ER) expression status (`groupA` = ER+, `groupB` = ER-). 
The data was reduced to 50 metabolites. 
For each group a data frame is given containing the raw data with the metabolites as rows and the
samples as columns. The first three columns contain the metabolite identifiers
(biochemical_name, metabolon_id and pubchem_id).
}
\keyword{datasets}
