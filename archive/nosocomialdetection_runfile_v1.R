# ================================================================= #
# Run file for noscomial detection
# ================================================================= #
setwd("~/PhD/covid-19/nosocomialtransmission/nosocomialdetection/")
source(file="nosocomialdetection_functions_v1.R")


maxday<-50 
# Probability distribution for incubation period
p1<-1.621
p2<-0.418
cum_prob_inc <- plnorm(1:maxday,p1,p2)
prob_inc <- cum_prob_inc-c(0,cum_prob_inc[1:(maxday-1)])

# Probability distribution for LOS
meanlos <- 7
cum_prob_los <- pexp(1:maxday,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(maxday-1)])
discrete.meanlos<-sum(prob_los*(1:maxday))

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
n <- 100
cum_prob_delay <- pgamma(1:n,shape = 0.811, rate = 0.064)
prob_delay <- cum_prob_delay-c(0,cum_prob_delay[1:(n-1)])
prob_delay

onset_to_discharge_cum_distr <- distr.onset.to.discharge(prob_los,prob_inc)$cum_distr
len_delay <- min(length(prob_delay),length(onset_to_discharge_cum_distr)-1)
# Counterfactual delay distribution (delay from symptom onset to study enrolment)
# that would have occurred if there was no discharge
cf_delay_distr <- prob_delay[1:len_delay]/(1-onset_to_discharge_cum_distr[1:len_delay])
cf_delay_distr <- cf_delay_distr/sum(cf_delay_distr)
cf_delay_distr
sum(cf_delay_distr)

plot(1:(length(cf_delay_distr)),cf_delay_distr[1:(length(cf_delay_distr))],
     xlab="Days",ylab="Probability density",main="Counterfactual delay distribution")

plot(1:(length(prob_delay)),cf_delay_distr[1:(length(prob_delay))],
     xlab="Days",ylab="Probability density",main="Counterfactual delay distribution")



cutoff <-5 # detection cutoff for definition of nosocomial cases

# Mui's function
nosocomial.detection(cutoff,prob_inc,prob_los,delay_distr =NULL)
nosocomial.detection(cutoff,prob_inc,prob_los,prob_delay)
# Ben's function
calc_prob_infection_meets_def_nosocomial(cutoff,prob_inc,prob_los)

# Mui's simulation
sim_nodelay <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_inc, 
                             delay_distr=NULL,
                             cutoff=cutoff)
sim_nodelay$res

sim_delay <- nosocomial.simulation(n_max=10000,
                             los_distr=prob_los, 
                             inc_distr=prob_inc, 
                             delay_distr=prob_delay,
                             cutoff=cutoff)
sim_delay$res

# Ben's simulation
nosocomial.simulation2(N=10000,
                      prob_los=prob_los, 
                      prob_inc=prob_inc, 
                      cutoff=cutoff)




