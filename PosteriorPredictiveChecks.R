####################################################################
#                                                                  #
# R Code for the Posterior Predictive Checks                       #
#                                                                  #
####################################################################

# Load the packages

library("coda")     # for run.jags()
library("rjags")    # for run.jags()
library("runjags")  # for run.jags()
library("parallel") # for the parallelization of run.jags()
load.module("glm")  # for run.jags
library("bayesplot")
library("ggplot2")
library("abind")
library("cowplot") # for plot_grid


# Load the data

diseaseresidualall_b6foltot <- read.csv(".../OCMdataWCRF.csv")

head( diseaseresidualall_b6foltot)

### Data Preparation for the BUGS Model ###

### Renomination of the data according to the BUGS model

Q_B6  <- diseaseresidualall_b6foltot$ResQEVITB6
Q_FOL <- diseaseresidualall_b6foltot$ResQEFOLATE
R_B6  <- diseaseresidualall_b6foltot$ResREVITB6
R_FOL <- diseaseresidualall_b6foltot$ResREFOLATE
M_B6  <- diseaseresidualall_b6foltot$ResVITB6
M_FOL <- diseaseresidualall_b6foltot$ResFOLATE

N <- length(Q_B6) 

Caseset <- diseaseresidualall_b6foltot$Caseset;

Diet <- cbind( Q_B6,Q_FOL,R_B6,R_FOL,M_B6,M_FOL)

############################################################################################# 
# Data preparation for the Disease Model - Conditional Logistic Regression
#############################################################################################

CaseIndex <- which( diseaseresidualall_b6foltot$Disease=='Incident')
StudyID <- as.factor( diseaseresidualall_b6foltot$STUDY)
casestudy <- StudyID[ CaseIndex]

CasesetIndex <- which( diseaseresidualall_b6foltot$Caseset != 'NA')
CasesetNew <- diseaseresidualall_b6foltot$Caseset[ CasesetIndex]

ContIndex <- which( diseaseresidualall_b6foltot$Disease=='Non-case')
contstudy <- StudyID[ ContIndex]

stopifnot( identical(casestudy,contstudy))  

N_Caseset <- length( CaseIndex)

Bmi <- diseaseresidualall_b6foltot$Bmi_C

Smoking <- as.factor(diseaseresidualall_b6foltot$Smoke_Stat)

## Create unique identifier for each case-control set starting from 1 and going up to
## N_Caseset. The order of the labels does not matter as they are only used for the
## risk set-specific random effects
UniqueRiskSet <- interaction(CasesetNew, StudyID, drop=TRUE)
UniqueRiskSet <- as.numeric(UniqueRiskSet)

Data <- list( "N"=N, "Diet"=Diet, "N_Caseset"=N_Caseset, "CaseIndex"=CaseIndex,
              "ContIndex"=ContIndex, "Study"=StudyID, "Bmi"=Bmi,"Smoking"=Smoking,
              "RiskSet"=UniqueRiskSet);

########################################################################################
# Model Initizialization                                                               
########################################################################################
inits1 <- list("beta_QB6"=0.1, "beta_MB6"=0.2,
               "beta_QFOL"=0.1,"beta_MFOL"=0.1,
               "beta_RFOL"=0.1, "beta_RB6"=0.1,
               "beta_TB6" =c(0,0), "beta_TFOL" =c(0,0),  
               "beta_BMI"=c(0,0), "beta_Smoking" = matrix(nrow=4,ncol=2,c(0,NA,0,0)),  
               .RNG.name="base::Super-Duper", .RNG.seed=200001);

inits2 <- list("beta_QB6"=0.1, "beta_MB6"=0.2,
               "beta_QFOL"=0.1,"beta_MFOL"=0.1,
               "beta_RFOL"=0.1, "beta_RB6"=0.1,
               "beta_TB6" =c(0,0), "beta_TFOL" =c(0,0),
               "beta_BMI"=c(0,0), "beta_Smoking" = matrix(nrow=4,ncol=2,c(0,NA,0,0)),
               .RNG.name="base::Wichmann-Hill", .RNG.seed=200002);

inits3 <- list("beta_QB6"=0.1, "beta_MB6"=0.2,
               "beta_QFOL"=0.1,"beta_MFOL"=0.1, 
               "beta_TB6" =c(0,0), "beta_TFOL" =c(0,0),
               "beta_RFOL"=0.1, "beta_RB6"=0.1,
               "beta_BMI"=c(0,0), "beta_Smoking" = matrix(nrow=4,ncol=2,c(0,NA,0,0)),
               .RNG.name="base::Mersenne-Twister", .RNG.seed=200003);   

########################################################################################       
# Running the BUGS Model for generating the Diet replications for the Posterior Checks                                                             
########################################################################################      
Nreps <- 500 # Number of replicate data sets to generate from each chain
results_bivariatebayesmodel_rep <- run.jags( model="BayesianHierarchicalModel_ForPosteriorPredictiveChecks.BUG", 
                                           monitor=c("Diet.rep"),
                                           data=Data, n.chains=3,
                                           method="parallel",###cl=cl,###
                                           inits=list(inits1,inits2,inits3), modules = "glm",
                                           ### burnin = 400, sample = 100, adapt = 100, thin=1);
                                           burnin = 15000, sample = Nreps, adapt = 10000, thin=100);

### Some data wrangling is required to get the replicate data sets in
### the right format for bayesplot
samples <- results_bivariatebayesmodel_rep$mcmc # Extract samples in coda format
samples <- lapply(samples, array, dim=c(Nreps, N, 6)) # Coerce to array
samples <- abind::abind(samples, along=1) # Collapse samples over all chains

### Extract individual replicate samples
Q_B6_rep <- samples[,,1]
Q_FOL_rep <- samples[,,2]
R_B6_rep <- samples[,,3]
R_FOL_rep <- samples[,,4]
M_B6_rep <- samples[,,5]
M_FOL_rep <- samples[,,6]

## Common parameters for all plots to be passed to ppc_dens_overlay
color_scheme_set("blue")
size <- 1
alpha <- 1 

## Dietary Questionnaire
p1 <- ppc_dens_overlay(Q_B6, Q_B6_rep, size=size, alpha=alpha) +
    ggtitle("Questionnaire vitamin-B6")
p2 <- ppc_dens_overlay(Q_FOL, Q_FOL_rep, size=size, alpha=alpha) +
    ggtitle("Questionnaire folate")
plot_grid(p1, p2)

## 24-hour recall
## It seems ppc_dens_overlay cannot cope with missing data, so we have
## to select the non-missing values
p3 <- ppc_dens_overlay(R_B6[!is.na(R_B6)], R_B6_rep[, !is.na(R_B6)], size=size, alpha=alpha) + ggtitle("24-hour recall vitamin-B6")
p4 <- ppc_dens_overlay(R_FOL[!is.na(R_FOL)], R_FOL_rep[, !is.na(R_FOL)], size=size, alpha=alpha) + ggtitle("24-hour recall folate")
plot_grid(p3, p4)

# Dietary Biomarker
p5 <- ppc_dens_overlay(M_B6, M_B6_rep, size=size, alpha=alpha) + ggtitle("Biomarker vitamin-B6")
p6 <- ppc_dens_overlay(M_FOL, M_FOL_rep, size=size, alpha=alpha) + ggtitle("Biomarker folate")
plot_grid(p5, p6)
