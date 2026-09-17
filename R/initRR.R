library(TMB)
library(mvtnorm)
library(colorRamps)
library(RColorBrewer)
library(reshape2)

#source("simRR.R")   # Simulate run reconstruction
source(here("R/fitRR.R"))   # Fit run reconstruction
#source("fitSim.R")  # Fit run reconstruction to simmed data
#source("calcRunTiming.R")
source(here("R/tools.R"))   # Background functions
source(here("R/plot.R"))    # Plotting functions

# Compile TMB objective function
compile(here("R/yukonChinookRunRecon.cpp"))
