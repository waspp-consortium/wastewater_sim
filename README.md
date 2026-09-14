# wastewater_sim
This repository contains code and data supporting the results presented in the manuscript "The public health value of wastewater surveillance for viruses with pandemic potential: a modelling study". We have also released the package [wastewatchR](https://github.com/mrc-ide/wastewatchR/releases/tag/paper1_version1) which we created to implement our modelling framework that we used to simulate zoonotic spillover, human-to-human transmission, and wastewater/clinical surveillance.

Within this repository is everything you need to recreate the results in our manuscript, with the exception of the SARS-CoV-2 case data shared with us by Dr Joanne Hewitt and Dr Bridget Armstrong which we used to define the relationship between relative viral load shed to wastewater and probability of wastewater detection. 

## Structure of the repository

### /data
- Data extracted during our literature review of viral shedding is found in `data/raw/literature_review_results/`.
- Publicly available published shedding profiles used are found in `data/raw/sc2_shedding_profiles/`.
- The  processed viral load estimates by 1) study and sample type, 2) virus per day, and 3) virus per infection are found in `data/processed/shedding/`.

### /R
- The folder contains support functions specific to the figures and analyses in this manuscript (and so not exported to wastewatchR).

### /scripts
- The scripts that encode the analyses and figure generation are all contained in this fodler.

### /results
- Model fits, estimates and simulations are stored here.

### /figures
- The manuscript figures generated in scripts are stored in here.

## Use of the repository
The code, data and figures can be examined within the repository, or reran from a forked clone of the repository. In order to rerun the majority of the analysis supporting our manuscript you should first fork the repository, create your own local copy and then follow the steps below. 

Please note you will not be able to run scripts `3_estimate_effective_sc2_shedders_NZ.R` or `4_refit_sensitivity_model_sc2_NZ_figS4.R` as these depend on case data from Hewitt et al., however the summary outputs of the refit model are provided to facilitate the rerunning of downstream analyses.

### 1. Install wastewatchR version paper1_version1
The version of the package used to run the analysis in the manuscript can be downloaded using:

```{r}
devtools::install_github("mrc-ide/wastewatchR@paper1_version1")
```

### 2. Run the analysis scripts in order
- The analysis scripts are numbered according to the order they need to be run in to be compatible with their dependencies. However, with the exception of large simulation files, the outputs have already been generated so that you may instead focus on replicating downstream analyses or figure generation.
- The scripts used to generate figures are named `generate_figX.R` where X is the relevant figure number.
  
At the top of each script is a standardised header detailing the purpose of the script, its dependencies and its outputs, e.g.:
```{r}
# This script provides an example (used to generate Figure 3 in the manuscript)
# of how to use wastewatchR to simulate spillover, onward transmission, and
# surveillance for an emerging virus. In this case the virus has a high spillover
# rate and low R0 similar to LASV, with a shedding profile and clinical presentation
# and delays similar to SARS-CoV-2

# It depends on:
### Detection model parameters: results/estimates/ww_detection_params.rds
### Shedding profile: data/processes/shedding/fecal_shedding_profile_sc2.rds
### function "calculate_clinical_detection_prob" in R/simulation_functions

# It produces:
### Simulation to be plotted in generate_fig3:
##### results/simulations/fig3.rds
##### results/simulations/fig3_linelist.rds

# Load libraries
library(wastewatchR)
library(tidyverse)
library(reshape2)
library(igraph)
library(ggraph)
library(cowplot)
library(here)

source(here("R/simulation_functions/estimate_clinical_detection_prob.R"))
#-------------------------------------------------------------------------------------
```

**PLEASE NOTE** that the following scripts will take considerable computational time to run:
- `5_simulations_fig2.R` at **~96 hrs** on a standard laptop
- `7_simulations_fig4AB.R` at **~96 hrs**
- `8_simulations_fig4C_figS2_figS3.R` at **~12 hrs**
  
We've tried to make this repository as complete and reproducible as possible but please contact amy.dighe@imperial.ac.uk when you come across any issues and questions.



