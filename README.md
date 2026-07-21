# DrDimont: A Pipeline for Drug Response Prediction from Differential Analysis of Multi-omics Networks

## Description
While it has been well established that drugs affect and help patients differently, personalized drug response predictions remain challenging. 
Solutions based on single omics measurements have been proposed, and networks provide means to incorporate molecular interactions into reasoning. 
However, how to integrate the wealth of information contained in multiple omics layers still poses a complex problem. 

We present a novel network analysis pipeline, DrDimont, Drug response prediction from Differential analysis of multi-omics networks. 
It allows for comparative conclusions between two conditions and translates them into differential drug response predictions. 
DrDimont focuses on molecular interactions. It establishes condition-specific networks from correlation within an omics layer that 
are then reduced and combined into heterogeneous, multi-omics molecular networks. A novel semi-local, path-based integration step 
ensures integrative conclusions. Differential predictions are derived from comparing the condition-specific integrated networks. 
DrDimont's predictions are explainable, i.e., molecular differences that are the source of high differential drug scores can be retrieved. 
Our proposed pipeline leverages multi-omics data for differential predictions, e.g. on drug response, and includes prior information on interactions.


## Installation of the package

1. Installing the R package
    - From CRAN: Use `install.packages("DrDimont")` to install the package and the R dependencies (v0.1.3 available)
    - From source: Either clone the repo and use `devtools::install()` within R to install the package, or use `remotes::install_gitlab("PHiort/DrDimont")` without cloning.
2. Installing the python dependencies
    - To use the differential drug response score computation, a Python (>= 3.8) installation is required. Once the DrDimont package is installed, use `DrDimont::install_python_dependencies()` to install the necessary dependencies automatically. You can use the function arguments to customize and use either pip or conda for the installation. If you prefer to install the dependencies manually, check out the requirements files in this repository `inst/requirements_pip.txt`, and `inst/requirements_conda.txt`.


## Exemplary Pipeline Execution
An exemplary pipeline execution with the included data can be found in `doc/DrDimont_Vignette.html` and at https://cran.r-project.org/web/packages/DrDimont/vignettes/DrDimont_Vignette.html.
The supplied case study data uses data published by Krug et al. (2020) (https://www.doi.org/10.1016/j.cell.2020.10.036). The package license applies only to the software and explicitly not to the included data.

## Additional Information

At the moment all functions are exported to make debugging easier. However, many functions are not intended for user-interaction. These functions are marked with `[INTERNAL]` in the function documentation.


The package `DrDimont` is an updated version of the previously published `molnet` package (https://github.com/molnet-org/molnet)