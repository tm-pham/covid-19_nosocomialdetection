# Updated functions for nosocomial detection using CO-CIN

# ================================================================= #
# Nosocomial detection, exact computation
# ================================================================= #
# Function computes probability that on a given day, a randomly
# infected patient is detected
# INPUT
# cutoff = detection cutoff for definition of nosocomial cases
# los_distr = probability distribution for LOS
# inc_distr = probability distribution for incubation period
# delay_distr = probability distribution for the delay between symptom onset 
#               and date of test (is empty/null, then no delay is considered)
# OUTPUT
# res = proportion detected
# res_ad = proportion with symptom onset after discharge
# res_bc = proportion with symptom onset before cutoff (1-res-res_ad)
nosocomial.detection <- function(los_distr, 
                                 inc_distr,
                                 cutoff,
                                 delay_distr=NULL){
  if(cutoff<2){
    print("Cutoff has to be at least 2.")
    return(0)
  }
  res_bc <- res <- res_ad <- 0
  length_los <- max(los_distr[,1]) # max LOS
  mean_los <- sum(los_distr[,1]*los_distr[,2])
  if(is.null(delay_distr)){ ##### Calculations for CO-CIN
    if(cutoff<length_los){
      for(l in cutoff:length_los){ # loop over possible LOS (>= cutoff)
        pl <- los_distr[which(los_distr[,1]==l),2]*l/mean_los # proportion of patients with LOS=l
        p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
        for(t in 1:(l-1)){# Proportion detected
          # Symptom onset after cutoff and before discharge 
          sum_res <- sum(inc_distr[max((cutoff-t),1):(l-t)]) 
          res <- res + pl*p_inf_on_day_l*sum_res
        }
      }
    }
    for(l in 2:length_los){# Proportion with symptom onset after discharge
      pl <- los_distr[which(los_distr[,1]==l),2]*l/mean_los # proportion of patients with LOS=l
      p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
      for(t in 1:(l-1)){
        sum_ad <- sum(inc_distr[(l-t+1):length(inc_distr)])
        res_ad <- res_ad + pl*p_inf_on_day_l*sum_ad
      }
      # Patients that get infected on day of discharge
      res_ad <- res_ad + pl*p_inf_on_day_l
    }
    # Proportion with symptom onset before cutoff 
    res_bc <- 1 - res - res_ad 
  }else{ ##### Calculations for SUS
    for(l in cutoff:length_los){ # Proportion
      pl <- los_distr[which(los_distr[,1]==l),2]*l/mean_los # proportion of patients with LOS=l
      p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
      for(t in 1:(l-1)){ # loop over possible infection days
        # d = index for delay_distr   
        # maximum delay = length-of-stay - time of infection or maximum delay in delay distr (whatever is smaller)
        # + 1/- 1 due to indices in R starting at 1
        for(d in 1:(min(l-t,length(delay_distr)))){
          # Note: Assume that min(incubation period)=1, hence 
          # if cutoff-t-(d-1)=0, then take 1
          sum_res <- sum(inc_distr[max(cutoff-t-(d-1),1):(l-t-(d-1))]) # possible values for incubation period
          res <- res + pl*p_inf_on_day_l*sum_res*delay_distr[d] 
        }
      }
    }
    for(l in 2:length_los){# Proportion with symptom onset after discharge
      pl <- los_distr[which(los_distr[,1]==l),2]*l/mean_los # proportion of patients with LOS=l
      p_inf_on_day_l <- 1/l # Assume infection is equally likely on each day
      for(t in 1:(l-1)){
        sum_ad <- sum(inc_distr[(l-t+1):length(inc_distr)])
        res_ad <- res_ad + pl*p_inf_on_day_l*sum_ad
      }
      # Patients that get infected on day of discharge
      res_ad <- res_ad + pl*p_inf_on_day_l
    }
  }
  return(list(res_bc = res_bc,
              res = res,
              res_ad = res_ad))
}


# Just for double-checking
# ================================================================= #
# Nosocomial detection simulation
# ================================================================= #
# Function computes the proportion of all detected nosocomial cases
# INPUT
# n_max = maximum number of nosocomial cases
# los_distr = probability distribution for LOS (can be empirical)
# inc_distr = probability distribution for incubation period
# cutoff = detection cutoff for definition of nosocomial case
# 
# DESCRIPTION
# 1. Function generates LOS for n_max patients (according to los_distr)
# 2. Subsequently, infection day is chosen randomly
# 3. Incubation period is drawn from inc_distr
# 4. Compute proportion of cases where symptom onset is before discharge
#    and after cutoff = proportion of detected cases
# ================================================================= #
nosocomial.simulation <- function(n_max=1000, 
                                  los_distr, 
                                  inc_distr,
                                  cutoff=10, 
                                  seed=12345){
  set.seed(seed)
  # 1. Draw length of stay for each patient (los_distr)
  t_los <- sample(los_distr[,1], size=n_max, prob=los_distr[,2], replace=TRUE)
  # 2. Draw an infection time for each patient
  # Currently: all patients get infected
  # Infection is equally likely on each day in hospital
  t_inf <- sapply(1:length(t_los), function(x) ceiling(runif(1)*t_los[x]))
  # 3. Draw times of symptom onset for each patient (inc_distr)
  # Accounts for asymptomatic cases: incubation time = 10000
  if(sum(inc_distr)==1){
    inc <- sample(1:length(inc_distr), size=length(t_los), prob=inc_distr, replace=TRUE)
  }else inc <- sample(c(1:length(inc_distr),10000), size=length(t_los), prob=c(inc_distr,1-sum(inc_distr)), replace=TRUE)
  t_inc <- t_inf + inc
  # 4. Determine the proportion of patients detected
  # Symptom onset >= cutoff and <= discharge
  # Patients with symptom onset before discharge
  ind_before_discharge <- which(t_inc<=t_los)
  # Patients with symptom onset after cutoff 
  ind_after_cutoff <- which(t_inc>=cutoff)
  # Patients that meet all criteria above
  ind <- intersect(ind_before_discharge, ind_after_cutoff)
  # Weight patients according to their length-of-stay
  los_weighted <- unlist(sapply(1:length(t_los), function(x) rep(x,length(which(t_los==x)))))
  ind_res <- which(los_weighted%in%ind)
  prop_detected <- length(ind_res)/length(los_weighted)
  p_1 <- length(ind)/n_max
  # Patients with symptom onset after discharge
  ind_after_discharge <- intersect(which(t_inc > t_los), which(t_los>1))
  ind_weighted_ad <- which(los_weighted%in%ind_after_discharge)
  prop_after_discharge <- length(ind_weighted_ad)/length(los_weighted)
  p_2 <- length(ind_after_discharge)/n_max
  # Patients with symptom onset before cutoff
  ind_before_cutoff <- intersect(which(t_inc < cutoff), which(t_inc < t_los))
  ind_weighted_c <- which(los_weighted%in%ind_before_cutoff)
  prop_before_cutoff <- length(ind_weighted_c)/length(los_weighted)
  p_3 <- length(ind_before_cutoff)/n_max
  return(list(t_los = t_los, 
              t_inf = t_inf,
              t_inc = t_inc,
              prop_detected = prop_detected, 
              prop_after_discharge = prop_after_discharge,
              prop_before_cutoff = prop_before_cutoff, 
              p_1=p_1, p_2=p_2, p_3=p_3))
}

# Ben's simulation
nosocomial.simulation2 <- function(N, 
                                   prob_los, 
                                   prob_inc,
                                   cutoff){
  test.data<-data.frame(los=rep(NA,N), time.infected=rep(NA,N), time.symptoms=rep(NA,N), detected=rep(NA,N), proportion_of_infections=rep(NA,N))
  for(i in 1:N){
    test.data$los[i]<-sample(1:length(prob_los), 1, replace=TRUE, prob=prob_los)
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








