# ================================================================= #
# Run file for noscomial detection
# ================================================================= #
setwd("/home/thi.mui.pham/covid-19/nosocomialtransmission/nosocomialdetection/code/")
source(file="nosocomialdetection_functions.R")


maxday<-50 
cutoff <-10 # detection cutoff for definition of nosocomial cases

# Probability distribution for incubation period
p1<-1.621
p2<-0.418
cum_prob_inc <- plnorm(1:maxday,p1,p2)
prob_inc <- cum_prob_inc-c(0,cum_prob_inc[1:(maxday-1)])

# Probality distribution for LOS
meanlos <- 5
cum_prob_los <- pexp(1:maxday,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(maxday-1)])
discrete.meanlos<-sum(prob_los*(1:maxday))

# Mui's function
nosocomial.detection(cutoff,prob_inc,prob_los)
# Ben's function
calc_prob_infection_meets_def_nosocomial(cutoff,prob_inc,prob_los)

# Simulation
sim <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_inc, 
                             cutoff=cutoff)
sim$res
