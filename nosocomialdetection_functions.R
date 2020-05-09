# ================================================================= #
# Nosocomial detection simulation
# ================================================================= #
# Function computes the proportion of all detected nosocomial cases
# INPUT
# n_max = maximum number of nosocomial cases
# los_distr = probability distribution for LOS
# inc_distr = probability distribution for incubation period
# cutoff = detection cutoff for definition of nosocomial case
# 
# DESCRIPTION
# Function generates LOS for n_max patients (according to los_distr)
# Subsequently, infection day is chosen randomly
# Incubation period is drawn from inc_distr
# Compute proportion of cases where symptom onset is before discharge
# and after cutoff = proportion of detected cases
# ================================================================= #
nosocomial.simulation <- function(n_max=1000, 
                                  los_distr, 
                                  inc_distr, 
                                  cutoff=10, 
                                  seed=12345){
  set.seed(seed)
  # Length of stay and discharge times
  t_los <- sample(1:length(los_distr), size=n_max, prob=los_distr, replace=TRUE)
  # Currently: all patients get infected
  # Infection is equally likely on each day in hospital
  t_inf <- sapply(1:length(t_los), function(x) sample(1:(t_los[x]-1),1))
  # Times for symptom onset according to inc_distr
  # Accounts for asymptomatic cases: incubation time = 10000
  inc <- sample(c(1:length(inc_distr),10000), size=length(t_los), prob=c(inc_distr,1-sum(inc_distr)), replace=TRUE)
  t_inc <- t_inf + inc
  # Patients with symptom onset before discharge
  ind_inc_before_discharge <- which(t_inc<=t_los)
  # Patients with symptom onset after cutoff 
  ind_inc_after_cutoff <- which(t_inc>=cutoff)
  # Patients that meet all criteria above
  ind <- intersect(ind_inc_before_discharge, ind_inc_after_cutoff)
  res <- length(ind)/length(t_los)
  return(list(res=res,
              t_los=t_los,
              t_inf=t_inf,
              t_inc=t_inc,
              inc=inc,
              ind=ind))
}

# ================================================================= #
# Separate calculations
# ================================================================= #
# Function computes probability that on a given day, a randomly
# infected patient is detected
# INPUT
# cutoff = detection cutoff for definition of nosocomial cases
# los_distr = probability distribution for LOS
# inc_distr = probability distribution for incubation period
nosocomial.detection <- function(cutoff, inc_distr, los_distr){
  res <- 0
  length_los <- length(los_distr) # max LOS
  mean_los <- sum((1:length_los)*los_distr)
  for(l in cutoff:length_los){ # loop over possible LOS (>= cutoff)
    pl <- los_prob[l]*l/mean_los # proportion of patients with LOS=l
    # pl <- los_distr[l]
    for(t in 1:(l-1)){ # loop over possible infection days
      sum_a <- sum(inc_distr[max((cutoff-t),1):(l-t)]) # possible values for incubation period
      p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
      res <- res + pl*p_inf_on_day_l*sum_a                
    }
  }
  return(res)
}

# ================================================================= #
# Ben's function
# ================================================================= #
calc_prob_infection_meets_def_nosocomial<-function(cutoff, 
                                                   los_vec, 
                                                   inc_vec){
  # input -cutoff for definition of nosocomial infection (delay from admission to onset)
  # los_vec hold probabilities of los of 1 day, 2 days, etc
  # inc_vec holds probabilities incubation  period is 1, 2 day, etc
  # function returns the probabiilty that a randomly infected in-patient in hospital 
  # will be counted as a noscomial case based on definition of X days
  total_prob<-0
  l<-length(los_vec)
  num_days_could_be_infected_for_given_latent_period<-0
  meanlos<- sum((1:l)*los_vec)
  for(i in cutoff:length(los_vec)){ #so loop over possible LoS to be a case
    for(l in 1:(i-1)){  # loop over possible latent period 
      los<-i
      prob_los<-los_vec[i]
      p_i<-prob_los
      pi_i<- p_i*los/ meanlos # prob selecting a patient with given los to infect
      # pi_i <- p_i
      num_days_could_be_infected_for_given_latent_period<- i-cutoff +1
      prob_infected_on_a_given_day_given_infected<-1/los
      prob_meeting_noso_definition<-pi_i*prob_infected_on_a_given_day_given_infected*num_days_could_be_infected_for_given_latent_period*inc_vec[l]
      total_prob<-total_prob+prob_meeting_noso_definition
    }  
  }
  return(total_prob)
}

# ================================================================= #
# Helper functions (not used for now)
# ================================================================= #
# Generate exponential times (for patient arrival times)
exp.times <- function(t_end, rate, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  t <- t_temp <- NULL
  t0 <- t_max <- 0
  while(t_max< t_end){
    temp <- rexp(1, rate=rate)
    t_temp <- c(t_temp, temp)
    t_max <- t_max + temp
    t <-c(t, t_max) 
  }
  if(t[length(t)]>t_end) {
    t <- t[-length(t)]
    t_temp <- t_temp[-length(t_temp)]
  }
  return(list(t=t,t_temp=t_temp))
}


arrival.times<- function(t_end, arrival_distr, seed=12345){
  set.seed(seed)
  t <- t_temp <- NULL
  t0 <- t_max <- 0
  while(t_max< t_end){
    temp <- sample(1:length(arrival_distr), size=1, prob=arrival_distr)
    t_temp <- c(t_temp, temp)
    t_max <- t_max + temp
    t <-c(t, t_max) 
  }
  t <- t[which(t<=t_end)]
  t_temp <- t_temp[1:length(t)]
  return(list(t=t,t_temp=t_temp))
}

