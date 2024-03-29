% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{install_python_dependencies}
\alias{install_python_dependencies}
\title{Installs python dependencies needed for interaction score computation}
\usage{
install_python_dependencies(package_manager = "pip")
}
\arguments{
\item{package_manager}{["pip"|"conda"] Package manager to use (default: pip)}
}
\value{
No return value, called to install python dependencies
}
\description{
Uses pip (default) or conda as specified to 
install all required Python modules. The Python packages are installed 
into a virtual Python or conda environment called 'r-DrDimont'. 
The following requirements are installed: numpy, tqdm, python-igraph and ray.
The environment is created with reticulate.
}
