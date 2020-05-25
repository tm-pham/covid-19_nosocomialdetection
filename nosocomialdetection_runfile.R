# ================================================================= #
# Run file for noscomial detection
# ================================================================= #
setwd("~/PhD/covid-19/nosocomialtransmission/nosocomialdetection/")
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
delay <- floor(rgamma(1000,shape=a,scale=b))
summary(delay)
tab_delay <- table(delay)
values_delay <- as.numeric(names(tab_delay))
prob_delay <- rep(0, max(values_delay))
# First entry in prob_delay corresponds to delay=0
prob_delay[values_delay+1] <- tab_delay/sum(tab_delay)

# First entry in prob_delay corresponds to delay=0
prob_delay <- c(1,rep(0,19))

# Probality distribution for LOS
meanlos <- 5
cum_prob_los <- pexp(1:maxday,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(maxday-1)])
discrete.meanlos<-sum(prob_los*(1:maxday))

cutoff <-10 # detection cutoff for definition of nosocomial cases

# Mui's function
nosocomial.detection(cutoff,prob_inc,prob_los,delay_distr =NULL)
nosocomial.detection(cutoff,prob_inc,prob_los,prob_delay)
# Ben's function
calc_prob_infection_meets_def_nosocomial(cutoff,prob_inc,prob_los)

# Mui's simulation
sim_delay <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_inc, 
                             delay_distr=NULL,
                             cutoff=cutoff)
sim_delay$res

sim_nodelay <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_inc, 
                             delay_distr=prob_delay,
                             cutoff=cutoff)
sim_nodelay$res

# Ben's simulation
nosocomial.simulation2(N=10000,
                      prob_los=prob_los, 
                      prob_inc=prob_inc, 
                      cutoff=cutoff)




