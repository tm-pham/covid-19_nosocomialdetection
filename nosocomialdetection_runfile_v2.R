# =========================================================================== #
# Run file for noscomial detection
# =========================================================================== #
library(ggplot2)

setwd("~/PhD/covid-19/nosocomialtransmission/nosocomialdetection/")
source(file="nosocomialdetection_functions_v2.R")

maxday<-50 
cutoff <- 10
# Probability distribution for incubation period
p1<-1.621
p2<-0.418
cum_prob_inc <- plnorm(1:maxday,p1,p2)
prob_inc <- cum_prob_inc-c(0,cum_prob_inc[1:(maxday-1)])

# Probability distribution for LOS (exponential for now)
meanlos <- 7
cum_prob_los <- pexp(1:maxday,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(maxday-1)])
discrete.meanlos<-sum(prob_los*(1:maxday))
prob_los <- cbind(1:maxday, prob_los)

# =========================================================================== #
# Proportion detected in CO-CIN
nosocomial.detection(prob_los, prob_inc, cutoff)


# Probability distribution for delay between symptom onset and getting tested
delay <- read.csv("testondayofsyndrome.csv", header=T, sep=",")
distr_delay <- delay$testondayofsyndrome
prob_delay <- table(distr_delay)/sum(tab_delay)

# =========================================================================== #
# Proportion detected in SUS
nosocomial.detection(prob_los, prob_inc, cutoff, prob_delay)

