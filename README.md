# BayesianHierarchicalModelInNutritionalEpidemiology
This repository is a companion to the manuscript
M. Pittavino, M. Plummer, M. Johansson, E. Riboli, P. Ferrari (2024)
"A Bayesian hierarchical framework to integrate dietary exposure and biomarker measurements into aetiological models."


The file "BayesianHierarchicalModel.BUG" contains the BUG file, with the Bayesian Hierarchical Model described in Section 2.2 and Subsections 2.2.1, 2.2.2 and 2.2.3 and Subsection 2.3.1.

The file "EpicDataResidualMethod" contains the R code, for the data preparation with the Residual Method of the EPIC data, described in Section 2.1.

The file "EpicDataAnalysis" contains the R code, for the implementation and validation of the Bayesian model, presented in Section 2.3.

The file "BayesianHierarchicalModel_DensityMixing" contains the PDF file with the resulting posterior distributions and the mixing of the MCMC simulations chains, discussed in the Section 3.
