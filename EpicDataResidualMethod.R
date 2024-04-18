#########################################################
#                                                       #
# R code for the Residual Method of the EPIC data for   #
# the Bayesian Hierarchical Model                       #
#                                                       #
#########################################################

# Load the packages

library(lme4)    # for lmer()
library(Matrix)  # for the Matrix decomposition

data <- read.csv("Table5studiesv2.csv")

names(data)

table(data$Study)

data=data[data$Study=="Kidney" |data$Study=="Lung",]

table(data$Batch)
table(data$Study)
table(data$country)

#Log transformation

data$LVITB6  = log(data$VITB6)
data$LFOLATE = log(data$FOLATE)
data$LQVITB6 = log(data$Qe_Vitb6)
data$LQFOL   = log(data$Qe_Fol)
data$LRVITB6 = log(data$RE_VITB6)
data$LRFOL   = log(data$RE_FOL)

data$Study   = as.factor(data$Study)
data$country = as.factor(data$country)
data$Batch   = as.factor(data$Batch)
data$Sex     = as.factor(data$Sex)
data$Age     = as.numeric(data$Age)

summary(data$Alc_Drinker, na.rm = TRUE)
head(data$Alc_Drinker)
table(data$Alc_Drinker)
table(data$FOLATE_Stat)
table(data$Batch)

#Residual method: 6 regression methods to extract the residuals for the data analysis

M1 <- lmer(LVITB6  ~ country + Study + Age + Sex + (1|Study:Batch), data=data,na.action=na.exclude)
M2 <- lmer(LFOLATE ~ country + Study + Age + Sex + (1|Study:Batch), data=data,na.action=na.exclude)
M3 <- lm(LQVITB6   ~ country + Age + Sex, data=data,na.action=na.exclude)
M4 <- lm(LQFOL     ~ country + Age + Sex, data=data,na.action=na.exclude)
M5 <- lm(LRVITB6   ~ Age + Sex, data=data,na.action=na.exclude)
M6 <- lm(LRFOL     ~ Age + Sex, data=data,na.action=na.exclude)

data$RLVITB6  <- resid(M1)
data$RLFOLATE <- resid(M2)
data$RLQVITB6 <- resid(M3)
data$RLQFOL   <- resid(M4)
data$RLRVITB6 <- resid(M5)
data$RLRFOL   <- resid(M6)

#Summary statistics of the data

sd(data$RLQVITB6) 
sd(data$RLRVITB6,na.rm=TRUE) 
sd(data$RLVITB6,na.rm=TRUE)  

sd(data$RLQFOL)   
sd(data$RLRFOL,na.rm=TRUE)   
sd(data$RLFOLATE,na.rm=TRUE) 

summary( data) 
names( data)


