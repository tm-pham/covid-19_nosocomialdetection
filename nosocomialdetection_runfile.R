# ================================================================= #
# Run file for noscomial detection
# ================================================================= #
setwd("/home/thi.mui.pham/covid-19/nosocomialtransmission/nosocomialdetection/code/")
source(file="nosocomialdetection_functions.R")


maxday<-50 
# Probability distribution for incubation period
p1<-1.621
p2<-0.418
cum_prob_inc <- plnorm(1:maxday,p1,p2)
prob_inc <- cum_prob_inc-c(0,cum_prob_inc[1:(maxday-1)])

# Prob. distr. for delay between symptom onset and study enrolment
set.seed(4)
a <- 1.2; b <- 4
delay <- ceiling(rgamma(1000,shape=a,scale=b))
summary(delay)
tab_delay <- table(delay)
values_delay <- as.numeric(names(tab_delay))
prob_delay <- rep(0, max(values_delay))
prob_delay[values_delay] <- tab_delay/sum(tab_delay)

# Probality distribution for LOS
meanlos <- 7
cum_prob_los <- pexp(1:maxday,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(maxday-1)])
discrete.meanlos<-sum(prob_los*(1:maxday))

cutoff <-10 # detection cutoff for definition of nosocomial cases

# Mui's function
nosocomial.detection(cutoff,prob_inc,prob_los)
nosocomial.detection(cutoff,prob_inc,prob_los,prob_delay)
# Ben's function
calc_prob_infection_meets_def_nosocomial(cutoff,prob_inc,prob_los)

# Mui's simulation
sim <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_inc, 
                             delay_distr=prob_delay,
                             cutoff=cutoff)
sim$res

# Ben's simulation
nosocomial.simulation2(N=10000,
                      prob_los=prob_los, 
                      prob_inc=prob_inc, 
                      cutoff=cutoff)




