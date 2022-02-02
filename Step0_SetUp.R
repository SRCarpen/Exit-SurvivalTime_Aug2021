# Setting Up to run Exit Time Analyses
# G.M. Wilkinson, 18 August 2021

#====================================
# NOTE: Make sure you have an updated version of R
# Recommend having R version 4.1.1 (2021-08-10) at minimum
R.version.string

# # If you want to update R on a WINDOWS machine, use the following code:
# if (!require(installr)) install.packages('installr')
# library(installr)
# updateR()

#=====================================
#Check if the required libraries are installed; if not, the library will be installed and then attached. This script will not check for any library updates - do this the first time you run this script. 

if (!require(bvpSolve)) install.packages('bvpSolve')
library('bvpSolve')
if (!require(cubature)) install.packages('cubature')
library('cubature')
if (!require(deSolve)) install.packages('deSolve')
library('deSolve')
if (!require(grDevices)) install.packages('grDevices')
library('grDevices')
if (!require(moments)) install.packages('moments')
library('moments')
if (!require(numDeriv)) install.packages('numDeriv')
library('numDeriv')
if (!require(parallel)) install.packages('parallel')
library('parallel')
if (!require(stats)) install.packages('stats')
library('stats')
if (!require(tseries)) install.packages('tseries')
library('tseries')
