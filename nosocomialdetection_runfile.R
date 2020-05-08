# ================================================================= #
# Run file for noscomial detection
# ================================================================= #
setwd("/home/thi.mui.pham/covid-19/nosocomialtransmission/code/")
source(file="nosocomialdetection_functions.R")


n<-20 # maximum length of stay
cutoff <-10 # detection cutoff for definition of nosocomial cases

# Exponential LOS distribution
meanlos <- 10
cum_prob_los <- pexp(1:n,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(n-1)])

# Lognormal LOS distribution
# "Data" of COVID-19 patients from GGD Amsterdam
# meanlos <- 9
# sdlos <- 0.9
# cum_prob_los <- plnorm(1:n, meanlog=log(meanlos), sdlog=sdlos)
# prob_los<-cum_prob_los-c(0,cum_prob_los[1:(n-1)])

# Probability distribution for incubation period
prob_inc <- c(0.05,0.1,0.15,0.2,0.1,rep(0.0,15))

# Mui's function
nosocomial.detection(cutoff,prob_inc,prob_los)
# Ben's function
calc_prob_infection_meets_def_nosocomial(cutoff,prob_los,prob_inc)

# Simulation
sim <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_latent, 
                             cutoff=cutoff)
sim$res
