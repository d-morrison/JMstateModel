# Load the packages and the function to estimate joint multi-state models:
library(mstate) # Please use the version 0.2.7
library(JM)
library(nlme)
# source("JMstateModel.R")
# Import two databases which contain longitudinal and multi-state data:
# load("data.RData")
load_all()
# library(JMStateModel)
library(dplyr)
rm(list = ls())

data_long = tibble(data_long)
data_surv = tibble(data_surv)
###############################
#### Longitudinal sub-part ####
###############################
# Fit the longitudinal responses using a linear mixed model:
lmeFit <- lme(fixed = Y ~ (times + I((1 + times)^(-1.2) - 1)) * X,
  data = data_long,
  random = ~ (times + I((1 + times)^(-1.2) - 1)) | id,
  method = "REML",
  control = list(opt = "optim"))
##############################
#### Multi-state sub-part ####
##############################
# Construct the 3*3 matrix of transitions:
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, 3] <- 3
dimnames(tmat) <- list(
  from = c("State_0", "State_1", "State_2"),
  to = c("State_0", "State_1", "State_2"))
tmat
# The transition ’0 -> 1’ is called ’1’,’0 -> 2’ is called ’2’ and
# ’1 -> 2’ is called ’3’.
# Define the covariate in the multi-state sub-part:
covs <- "X"
# The ’msprep()’ function divides the multi-state database in order to have
# one line per transition at risk for each subject, with ’Tstart’ the
# entry time in the current state, and ’Tstop’ the time of transition or
# censorship; ’status’ denotes if the transition has been performed:
data_mstate <- msprep(
  time = c(NA, "t_State_1", "t_State_2"),
  status = c(NA, "State_1", "State_2"),
  data = data_surv,
  trans = tmat,
  keep = covs,
  id = "id")
# ’expand.covs()’ permits to define the set of covariates which impacts
# each transition:
data_mstate <- expand.covs(data_mstate, covs,
  append = TRUE, longnames = FALSE)
# Multi-state model with transition-specific proportional intensities:
coxFit <- coxph(Surv(Tstart, Tstop, status) ~
    X.1 + X.2 + X.3 + strata(trans),
  data = data_mstate,
  method = "breslow",
  x = TRUE,
  model = TRUE)
################################
#### Joint multi-state part ####
################################
# To define the dependency on the slope of the marker, it is necessary
# to specify the derivative of the fixed and random parts in the mixed model,
# and indicate which covariates are kept, :
dForm <- list(fixed = ~ 1 + I((-1.2) * ((1 + times)^(-2.2))) +
    X + I((-1.2) * ((1 + times)^(-2.2))):X,
  indFixed = c(2:3, 5:6),
  random = ~ 1 + I((-1.2) * ((1 + times)^(-2.2))),
  indRandom = 2:3)
# Joint multi-state model with:
# - true current level and true current slope of the marker as dependence function,
# - cubic B-splines with 1 internal knot for each log-baseline intensity,
# - 15 Gauss-Kronrod quadrature points to approximate the integral over time
# (by default),
# - 3 Gauss-Hermite quadrature points in the pseudo-adaptive numerical
# integration to approximate the integral over random effects.
jointFit_1step_GHk3 <-
  JMstateModel(
    lmeObject = lmeFit,
    survObject = coxFit,
    timeVar = "times",
    parameterization = "both",
    method = "spline-PH-aGH",
    interFact = list(
      value = ~strata(trans) - 1,
      slope = ~strata(trans) - 1,
      data = data_mstate),
    derivForm = dForm,
    Mstate = TRUE,
    data.Mstate = data_mstate,
    ID.Mstate = "id",
    verbose = TRUE,
    control = list(GHk = 3, lng.in.kn = 1))
summary(jointFit_1step_GHk3)
# Same joint multi-state model with:
# - 9 Gauss-Hermite quadrature points in the pseudo-adaptive numerical
# integration to approximate the integral over random effects.
jointFit_1step_GHk9 <-
  JMstateModel(lmeObject = lmeFit,
    survObject = coxFit,
    timeVar = "times",
    parameterization = "both",
    method = "spline-PH-aGH",
    interFact = list(value = ~strata(trans) - 1,
      slope = ~strata(trans) - 1,
      data = data_mstate),
    derivForm = dForm,
    Mstate = TRUE,
    data.Mstate = data_mstate,
    ID.Mstate = "id",
    control = list(GHk = 9, lng.in.kn = 1))
summary(jointFit_1step_GHk9)
# To use the multi-step pseudo-adaptive Gauss-Hermite rule, we have to source
# two functions inspired by JM:
# source("modified.log.posterior.b2.R")
# source("modified.ranef.jointModel.R")
# Same joint multi-state model with:
# - 9 and 9 Gauss-Hermite quadrature points in the two-step pseudo-adaptive
# numerical integration to approximate the integral over random effects.
# We can choose the posterior mode (true definition) or the posterior mean
# (faster) of the random effects of the fitted joint model (defined in ’init’)
# to update the quadrature points. Here the mode is used.
jointFit_2step_GHk9_9 <-
  JMstateModel(lmeObject = lmeFit,
    survObject = coxFit,
    timeVar = "times",
    parameterization = "both",
    method = "spline-PH-aGH",
    interFact = list(value = ~strata(trans) - 1,
      slope = ~strata(trans) - 1,
      data = data_mstate),
    derivForm = dForm,
    Mstate = TRUE,
    data.Mstate = data_mstate,
    ID.Mstate = "id",
    control = list(GHk = 9, lng.in.kn = 1),
    init = jointFit_1step_GHk9,
    init.type.ranef = "mode")
summary(jointFit_2step_GHk9_9)
