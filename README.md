# Perceived Behavioral Plasticity
This repository includes code to reproduce all analyses and figures in the paper Nielsen, Dablander, Debnath, Emegor, Ghai, Gwozdz, Hahnel, Hofmann, & Bauer (2024). Perceived plasticity of climate-relevant behaviors and policy support among high- and lower-income individuals.

The structure of this repository is as follows:

- **analysis.Rmd** includes the relevant R code
- **analysis.html** is the output of the analysis
- **helpers.R** holds useful R functions
- **Figures/** includes all figures that are included in the manuscript
- **Data/** includes the data
  - **Data/BPdata.RDS** is the main data file
  - **Data/df_preds_loose.csv** are model predictions for the loosely domain matched policies, stored to not rerun the computations
  - **Data/df_sensitivity_matched.csv** are sensitivity analyses results for the policy analysis, stored to not rerun the computations
