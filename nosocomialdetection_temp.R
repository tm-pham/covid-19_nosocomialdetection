sim2 <- nosocomial.simulation(n_max=1000, 
                              los_distr=prob_los, 
                              inc_distr=prob_inc, 
                              cutoff=5)
hist(sim2$t_los)
table(sim2$t_los)
sim2$prop_detected


c_vec <- c(1,5,8,10,14,20,25)

prop_det <- NULL; i <- 1
for(cutoff in c_vec){
  prop_det[i] <- nosocomial.simulation(n_max=1000, 
                                       los_distr=prob_los, 
                                       inc_distr=prob_inc, 
                                       cutoff=cutoff)$prop_detected
  i <- i + 1
}

prop_det <- NULL; i <- 1
for(cutoff in c_vec){
  prop_det[i] <- nosocomial.detection(prob_los, prob_inc, cutoff=cutoff)
  i <- i + 1
}


prop_det_data <- as.data.frame(cbind(cutoff=c_vec, prop_det))
prop_det_data

ggplot(data=prop_det_data, aes(x=cutoff, y=prop_det)) + 
  geom_line() + 
  geom_point() + 
  scale_x_continuous(breaks=c_vec, labels=c_vec) + 
  ylab("Proportion detected")