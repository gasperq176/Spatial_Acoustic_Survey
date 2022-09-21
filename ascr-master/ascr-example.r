---
  title: "testing"
author: "Gasper Qian"
date: "5 January 2019"
output: html_document
---
#Fetching the package and fit the model.
library(ascr)
## This will create objects capt, mask, rates, and traps, used below.
library(spatstat)
load("scr-data.RData")
## Fitting the model.
## capt: A list with two components: "bincapt", the 'binary capture histories' and "toa", the times of arrival.
## traps: The microphone coordinates.
## mask: The mask object.
## cue.rates: A vector containing the numbers of times 8 frogs called during a one-minute period.
## survey.length: The length of the survey, in minutes (these data are from only a 25-second recording).
## trace: if TRUE, a bunch of output will be printed relating to optimisation of the likelihood. If this is annoying, set to FALSE.
fit <- fit.ascr(capt = capt, traps = traps, mask = mask,
                cue.rates = rates, survey.length = 25/60, trace = FALSE)
## For parameter estimates.
## Da: Animal density (animals per hectare).
## Dc: Call density (calls per minute per hectare).
## g0, sigma: Detection function parameters.
## sigma.toa: Time-of-arrival error standard deviation (in seconds).
## mu.rates: Expected number of calls per animal per minute.
summary(fit)
## These can also be extracted by parameter names as follows.
coef(fit, c("Da", "sigma"))
## A plot of the estimated detection function.
show.detfn(fit)
## A plot of the estimated location of the second call.
locations(fit, 2)

#Function body for capture histroy and time arrival

pars <- coef(fit, c("Da", "g0", "sigma", "Dc", "sigma.toa", "mu.rates"))
pars["mu.rates"]

sim.scr <- function(pars, traps, buffer){
  Da <- pars["Da"]
  g0 <- pars["g0"]
  sigma <- pars["sigma"]
  Dc <- pars["Dc"]
  sigma.toa <- pars["sigma.toa"]
  mu.rates <- pars["mu.rates"]
  speed <- 343
  surveylen <- 25
  
  ##Generate the region area.
  region.lims <- c(NA,NA,NA,NA)
  region.lims[1]<- min(traps[,1])-buffer
  region.lims[2]<- max(traps[,1])+buffer
  region.lims[3]<- min(traps[,2])-buffer
  region.lims[4]<- max(traps[,2])+buffer
  region.area <- (region.lims[2] - region.lims[1])*
    (region.lims[4] - region.lims[3])/10000
  ## Extracting number of detectors.
  n.traps <- nrow(traps)
  ## Simulating number of animals.
  n.acs <- rpois(1, Da*region.area)
  
  ## Simulating activity centre locations.
  ac.locs.x <- runif(n.acs, region.lims[1], region.lims[2])
  ac.locs.y <- runif(n.acs, region.lims[3], region.lims[4])
  ac.loc <- cbind(ac.locs.x, ac.locs.y)
  ID <- seq.int(nrow(ac.loc))
  ac.locs <- cbind(ac.loc,ID)#location of animals 
  
  ##Simulating call centre locations depend on activity centre locations.
  localrate<- mu.rates*surveylen/60
  dupnum <- rpois(n.acs,localrate)
  ac.locs <- ac.locs[rep(1:nrow(ac.locs),dupnum),]
  
  ## Calculating a matrix of differences between activity centres and detectors.
  dists <- crossdist(ac.locs[, 1], ac.locs[, 2],
                     traps[, 1], traps[, 2])
  ## Calculating detection probabilities.
  det.probs <- g0*exp(-dists^2/(2*sigma^2))
  ## Creating capture histories with ID.
  capts.full <- matrix(rbinom(nrow(ac.locs)*n.traps, 1, det.probs),
                       nrow = nrow(ac.locs), ncol = n.traps)
  capt.full <- cbind(capts.full, ac.locs[,ncol(ac.locs)])
  colnames(capt.full)[ncol(capt.full)] <- "ID"
  ## We only observe capture histories with at least one detection.
  capt$bincapt <- capt.full[apply(capt.full[,1:(ncol(capt.full)-1)], 1, sum) > 0, ]
  
  colfeatures <- c(sprintf("%02d", seq(1,ncol(capt$bincapt)-1)))

  colnames(capt$bincapt)[1:(ncol(capt$bincapt)-1)] <- colfeatures

  ## get the time for each arrival depends on the sd of random error
  timerec <- matrix(data = NA, nrow = nrow(dists), ncol = ncol(dists))
  for(i in 1: nrow(dists)){
    ranstart <- runif(1,0,surveylen)
    for(j in 1: ncol(dists)){
      timerec[i,j] <- ranstart + dists[i,j]/speed + rnorm(1,0,sigma.toa)
    }
  }
  #timerec[timerec > surveylen] <- 0#Clean the data exceeding the period boundry
  timerecs <- cbind(timerec,ac.locs[,ncol(ac.locs)])
  colnames(timerecs)[ncol(timerecs)] <- c("ID")
  timerecs[capt.full[,1:ncol(capt.full-1)] == 0] <- 0
  
  capt$toa <- timerecs[apply(timerecs[,1:(ncol(timerecs)-1)], 1, sum) > 0, ]
  colnames(capt$toa)[1:(ncol(capt$toa)-1)] <- colfeatures

  capt
}

#Simple function test on test data
set.seed(1234)
sim.scr(pars,traps,30)

