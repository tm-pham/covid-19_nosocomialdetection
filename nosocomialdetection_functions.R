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
                                  delay_distr=NULL,
                                  cutoff=10, 
                                  seed=12345){
  set.seed(seed)
  # Length of stay and discharge times
  t_los <- sample(1:length(los_distr), size=n_max, prob=los_distr, replace=TRUE)
  # Currently: all patients get infected
  # Infection is equally likely on each day in hospital
  t_inf <- sapply(1:length(t_los), function(x) ceiling(runif(1)*t_los[x]))
  # Times for symptom onset according to inc_distr
  # Accounts for asymptomatic cases: incubation time = 10000
  if(sum(inc_distr)==1){
    inc <- sample(1:length(inc_distr), size=length(t_los), prob=inc_distr, replace=TRUE)
  }else inc <- sample(c(1:length(inc_distr),10000), size=length(t_los), prob=c(inc_distr,1-sum(inc_distr)), replace=TRUE)
  t_inc <- t_inf + inc
  # Delay between symptom onset and study enrolment
  if(is.null(delay_distr)){
    delay <- rep(0,length(t_los))
  }else delay <- sample(0:(length(delay_distr)-1),size=length(t_los),prob=delay_distr,replace=TRUE)
  t_detection <- t_inc + delay
  # Patients with symptom onset before discharge
  ind_before_discharge <- which(t_detection<=t_los)
  # Patients with symptom onset after cutoff 
  ind_after_cutoff <- which(t_detection>=cutoff)
  # Patients that meet all criteria above
  ind <- intersect(ind_before_discharge, ind_after_cutoff)
  # Weight patients according to their length-of-stay
  los_weighted <- unlist(sapply(1:length(t_los), function(x) rep(x,t_los[x])))
  ind_res <- which(los_weighted%in%ind)
  res <- length(ind_res)/length(los_weighted)
  return(list(res=res,
              t_los=t_los,
              t_inf=t_inf,
              t_inc=t_inc,
              t_detection=t_detection,
              inc=inc,
              delay=delay,
              ind=ind))
}

# Ben's simulation
nosocomial.simulation2 <- function(N, 
                                   prob_los, 
                                   prob_inc,
                                   cutoff){
  test.data<-data.frame(los=rep(NA,N), time.infected=rep(NA,N), time.symptoms=rep(NA,N), detected=rep(NA,N), proportion_of_infections=rep(NA,N))
  for(i in 1:N){
    test.data$los[i]<-sample(1:maxday, 1, replace=TRUE, prob=prob_los)
    test.data$time.infected[i]<-ceiling(runif(1)*test.data$los[i])
    inc.period<-sample(1:maxday, 1, replace=TRUE, prob=prob_inc)
    test.data$time.symptoms[i]<-test.data$time.infected[i]+inc.period
    test.data$detected[i]<- (test.data$time.symptoms[i] >= cutoff) & ((test.data$time.symptoms[i] <= test.data$los[i]))
    test.data$proportion_of_infections[i] <- test.data$los[i]/discrete.meanlos
    
  }
  detected<-sum(test.data$proportion_of_infections[test.data$detected])
  not.detected<-sum(test.data$proportion_of_infections[!test.data$detected])
  proportion.detected<-detected/(detected+not.detected)
  return(proportion.detected)
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
nosocomial.detection <- function(cutoff, 
                                 inc_distr, 
                                 los_distr, 
                                 delay_distr=NULL){
  res <- 0
  length_los <- length(los_distr) # max LOS
  mean_los <- sum((1:length_los)*los_distr)
  if(is.null(delay_distr)){
    for(l in cutoff:length_los){ # loop over possible LOS (>= cutoff)
      pl <- los_distr[l]*l/mean_los # proportion of patients with LOS=l
      for(t in 1:(l-1)){ # loop over possible infection days
        sum_a <- sum(inc_distr[max((cutoff-t),1):(l-t)]) # possible values for incubation period
        p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
        res <- res + pl*p_inf_on_day_l*sum_a                
      }
    }
  }else{
    for(l in cutoff:length_los){ # loop over possible LOS (>= cutoff)
      pl <- los_distr[l]*l/mean_los # proportion of patients with LOS=l
      for(t in 1:(l-1)){ # loop over possible infection days
        # d = index for delay_distr 
        # maximum delay = length-of-stay - time of infection or maximum delay in delay distr (whatever is smaller)
        # + 1/- 1 due to indices in R starting at 1
        for(d in 1:(min(l-t+1,length(delay_distr)))){
          # Note: Assume that min(incubation period)=1, hence 
          # if cutoff-t-(d-1)=0, then take 1
          sum_a <- sum(inc_distr[max(cutoff-t-(d-1),1):(l-t-(d-1))]) # possible values for incubation period
          p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
          res <- res + pl*p_inf_on_day_l*sum_a*delay_distr[d] 
        }
      }
    }
  }
  return(res)
}

# ================================================================= #
# Ben's function
# ================================================================= #
calc_prob_infection_meets_def_nosocomial <- function(cutoff, 
                                                     inc_vec,
                                                     los_vec){
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
      if(l<cutoff){
        num_days_could_be_infected_for_given_latent_period<- i-cutoff +1
      } else {
        num_days_could_be_infected_for_given_latent_period<- max(0,i-cutoff - (l-cutoff))
      }
      prob_infected_on_a_given_day_given_infected<-1/los
      prob_meeting_noso_definition<-pi_i*prob_infected_on_a_given_day_given_infected*num_days_could_be_infected_for_given_latent_period*inc_vec[l]
      total_prob<-total_prob+prob_meeting_noso_definition
    }
  }
  return(total_prob)
}

# ================================================================= #
# Probability and cumulative distribution of time from symptom 
# onset until discharge
# ================================================================= #
# Note that it does not sum to 1 because it does not return 
# negative values (inc > los)
# Note that prob_distr[1] corresponds to time from symptom onset 
# to discharge = 0
# INPUT
# los_distr = Probability distribution of length-of-stay
# inc_distr = Probability distribution of incubation period
# OUTPUT
# prob_distr = Probability distribution of time of symptom onset to
#              discharge
# cum_distr = Corresponding cumulative distribution
# ================================================================= #
distr.onset.to.discharge <- function(los_distr, inc_distr){
  prob_distr<- NULL
  for(z in 0:length(los_distr)){
    temp <- 0
    for(i in 1:min(length(los_distr)-z, length(inc_distr))){
      temp <- temp + inc_distr[i]*los_distr[z+i]
    }
    prob_distr <- c(prob_distr,temp)
  }
  cum_distr <- sapply(1:length(prob_distr), function(x) sum(prob_distr[1:x]))
  return(list(prob_distr=prob_distr, cum_distr=cum_distr))
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

