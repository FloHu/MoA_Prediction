# MoA_Prediction

## Organisation

### Directories:
`data/`: little data sets picked up from wherever, no big tables, these should stay on the server.
`plots/`: for all plots.
`programmatic_output/`: stuff generated by scripts, don't put stuff there yourself. 
`R/`: all R functions go here (one file per function) so that our notebooks don't get cluttered. 

### Notebooks:
Notebooks are numbered by the order they should be run.


## Some general points:
* Performance measures: contingency tables, precision-recall curves, ROC-AUC, Matthew's correlation coefficient (MCC), Cohen's Kappa
* Do nested cross-validation (starting with 5-fold), splits should be stratified
* Weighting of observations necessary to account for class imbalances? 
* Modes of action 'oxidative_stress', 'pmf', 'protein_qc', and 'unknown' either have too few members or don't make
  sense in training, they are therefore excluded.
* Feature preselection: remove genes w/o significant interactions --> preselect by ANOVA/hierarchical
  clustering/histogram of pairwise correlations with a cutoff
* Florian: random forests
* Léonard: boosted trees 

