#########################################################################
### Felix Pabon-Rodriguez
### Dissertation R Code
### Bayesian Capture-Recapture Model for Mice/Tick RTV Field Data
#########################################################################


#########################################################################
# Loading libraries
#########################################################################

# Loading libraries
options(repos="https://cran.rstudio.com")
install_load <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

my_packages <- c("readxl", "dplyr", "coda", "rjags", "MASS", "ggplot2",
                 "R2jags", "lattice", "bayesplot", "BayesPostEst", "ggmcmc", 
                 "RCurl", "truncnorm", "kableExtra", "mvtnorm", "rlist", 
                 "extraDistr", "msm", "tmvtnorm", "runjags", "plotrix",
                 "lubridate", "ggpubr", "stringr", "nimble",
                 "igraph", "parallel", "doParallel", "MCMCvis", "tibble")
invisible(install_load(my_packages))


########################################################################
# Reading data
########################################################################

source(file = "Read_RTVField_Data_ExtendedVersion.R")

########################################################################
# Posterior summary
########################################################################

## Function to get other ecological/epid estimates

get_other_params <- function(inc.mice,size.mice,prot.mice,
                            inc.nymphs, sample.nymphs,
                            inc.drags, sample.drags){
  ## Prevalence
  # if     p ~ dunif(0, 1)
  # and    y ~ dbinom(p,n)
  # then,  p|y ~ beta(a+y,b+N-y)
  sim <- 10000
  prev.mice <- matrix(NA,sim,6)
  protected.mice <- matrix(NA,sim,6)
  nip.mice <- matrix(NA,sim,6)
  nip.drag <- matrix(NA,sim,6)
  for(s in 1:6){
    prev.mice[,s] <- rbeta(sim,1+inc.mice[s],
                           1+size.mice[s]-inc.mice[s])
    protected.mice[,s] <- rbeta(sim,1+prot.mice[s],
                                1+size.mice[s]-prot.mice[s])
    nip.mice[,s] <- rbeta(sim,1+inc.nymphs[s],
                          1+sample.nymphs[s]-inc.nymphs[s])
    nip.drag[,s] <- rbeta(sim,1+inc.drags[s],
                          1+sample.drags[s]-inc.drags[s])
  }
  return(list(prev.mice=prev.mice,
              protected.mice=protected.mice,
              nip.mice=nip.mice,
              nip.drag=nip.drag))
}


## 2020

results_file2020 <- readRDS("fit_DP2020_idx386.rds")
results_mcmc2020 <- as.mcmc.list(lapply(1:3, 
                                        function(x){as.mcmc(results_file2020[[x]]$samples)}))
par.names2020 <- colnames(results_mcmc2020[[1]])

MCMCsummary(results_mcmc2020,round = 2)[,c(4,2,3,5)]

MCMCtrace(results_mcmc2020, 
          params = par.names2020,
          ISB = FALSE,
          pdf = TRUE,
          Rhat = TRUE)

cc2020 <- combine.mcmc(results_mcmc2020)
apply(cc2020,2,function(x)mean(x>0))

post2020 <- get_other_params(inc.mice = Final_Incidence_Mice_2020,
                            size.mice = SampleSize2020,
                            prot.mice = Final_Protection_Mice_2020,
                            inc.nymphs = Final_Nymphals_Mice_2020, 
                            sample.nymphs = SampleNymphals2020,
                            inc.drags = Final_Nymphals_Drag_2020, 
                            sample.drags = SampleDrag2020)


## 2021

results_file2021 <- readRDS("fit_DP2021_idx69.rds")
results_mcmc2021 <- as.mcmc.list(lapply(1:3, 
                                        function(x){as.mcmc(results_file2021[[x]]$samples)}))
par.names2021 <- colnames(results_mcmc2021[[1]])

MCMCsummary(results_mcmc2021,round = 2)[,c(4,2,3,5)]

MCMCtrace(results_mcmc2021, 
          params = par.names2021,
          ISB = FALSE,
          pdf = TRUE,
          Rhat = TRUE)

cc2021 <- combine.mcmc(results_mcmc2021)
apply(cc2021,2,function(x)mean(x>0))

post2021 <- get_other_params(inc.mice = Final_Incidence_Mice_2021,
                            size.mice = SampleSize2021,
                            prot.mice = Final_Protection_Mice_2021,
                            inc.nymphs = Final_Nymphals_Mice_2021, 
                            sample.nymphs = SampleNymphals2021,
                            inc.drags = Final_Nymphals_Drag_2021, 
                            sample.drags = SampleDrag2021)

## 2022

results_file2022 <- readRDS("fit_DP2022_idx883.rds")
results_mcmc2022 <- as.mcmc.list(lapply(1:3, 
                                        function(x){as.mcmc(results_file2022[[x]]$samples)}))
par.names2022 <- colnames(results_mcmc2022[[1]])

MCMCsummary(results_mcmc2022,round = 2)[,c(4,2,3,5)]

MCMCtrace(results_mcmc2022, 
          params = par.names2022,
          ISB = FALSE,
          pdf = TRUE,
          Rhat = TRUE)

cc2022 <- combine.mcmc(results_mcmc2022)
apply(cc2022,2,function(x)mean(x>0))

post2022 <- get_other_params(inc.mice = Final_Incidence_Mice_2022,
                            size.mice = SampleSize2022,
                            prot.mice = Final_Protection_Mice_2022,
                            inc.nymphs = Final_Nymphals_Mice_2022, 
                            sample.nymphs = SampleNymphals2022,
                            inc.drags = Final_Nymphals_Drag_2022, 
                            sample.drags = SampleDrag2022)




########################################################################
# Plots
########################################################################

# Site Pairs 
# Control vs Trt
# Orange vs Green
# CPB(1) vs MHB(6) 
# EHH(2) vs HCI(5)
# GFH(3) vs GSV(4)

percent_difference <- function(start,final,decrease=TRUE){
  if(decrease==TRUE){
    return(((start-final)/start))
  }
  return(((final-start)/start))
}


# Plots
df_mice_prev <- data.frame(Pair = rep(c("P1","P2","P3","P3","P2","P1"),3),
                           Group = rep(c("Control","Control","Control",
                                         "Treatment","Treatment","Treatment"),3),
                           Year = c(rep("2020",6),rep("2021",6),rep("2022",6)),
                           Median = c(apply(post2020$prev.mice,2,median),
                                      apply(post2021$prev.mice,2,median),
                                      apply(post2022$prev.mice,2,median)),
                           Lower = c(apply(post2020$prev.mice,2,function(x)quantile(x,prob=0.025)),
                                     apply(post2021$prev.mice,2,function(x)quantile(x,prob=0.025)),
                                     apply(post2022$prev.mice,2,function(x)quantile(x,prob=0.025))),
                           Upper = c(apply(post2020$prev.mice,2,function(x)quantile(x,prob=0.975)),
                                     apply(post2021$prev.mice,2,function(x)quantile(x,prob=0.975)),
                                     apply(post2022$prev.mice,2,function(x)quantile(x,prob=0.975))))


ggplot(df_mice_prev, aes(Year, Median, group = Group, color = Group)) +        
  geom_errorbar(aes(ymin = Lower, ymax = Upper),width=0.3,size=1) +
  geom_point(aes(group = Group),size=2) +
  geom_line(aes(group = Group), size=1) +
  facet_wrap(~Pair) +
  labs(title = "Mice Infection Prevalence \n Control-vs-Treatment Site Pairs",
       y = "Posterior Median (95% Cr-I)", x = "Year") +
  scale_color_manual(breaks = c("Control", "Treatment"),
                     values=c("darkorange", "darkgreen"))



df_nip_mice <- data.frame(Pair = rep(c("P1","P2","P3","P3","P2","P1"),3),
                          Group = rep(c("Control","Control","Control",
                                        "Treatment","Treatment","Treatment"),3),
                           Year = c(rep("2020",6),rep("2021",6),rep("2022",6)),
                           Median = c(apply(post2020$nip.mice,2,median),
                                      apply(post2021$nip.mice,2,median),
                                      apply(post2022$nip.mice,2,median)),
                           Lower = c(apply(post2020$nip.mice,2,function(x)quantile(x,prob=0.025)),
                                     apply(post2021$nip.mice,2,function(x)quantile(x,prob=0.025)),
                                     apply(post2022$nip.mice,2,function(x)quantile(x,prob=0.025))),
                           Upper = c(apply(post2020$nip.mice,2,function(x)quantile(x,prob=0.975)),
                                     apply(post2021$nip.mice,2,function(x)quantile(x,prob=0.975)),
                                     apply(post2022$nip.mice,2,function(x)quantile(x,prob=0.975))))


ggplot(df_nip_mice, aes(Year, Median, group = Group, color = Group)) +        
  geom_errorbar(aes(ymin = Lower, ymax = Upper),width=0.3,size=1) +
  geom_point(aes(group = Group),size=2) +
  geom_line(aes(group = Group), size=1) +
  facet_wrap(~Pair) +
  labs(title = "Nymphal Infection Prevalence (Extracted from Mice) \n Control-vs-Treatment Site Pairs",
       y = "Posterior Median (95% Cr-I)", x = "Year") +
  scale_color_manual(breaks = c("Control", "Treatment"),
                     values=c("darkorange", "darkgreen"))



df_nip_drag <- data.frame(Pair = rep(c("P1","P2","P3","P3","P2","P1"),3),
                          Group = rep(c("Control","Control","Control",
                                        "Treatment","Treatment","Treatment"),3),
                          Year = c(rep("2020",6),rep("2021",6),rep("2022",6)),
                          Median = c(apply(post2020$nip.drag,2,median),
                                     apply(post2021$nip.drag,2,median),
                                     apply(post2022$nip.drag,2,median)),
                          Lower = c(apply(post2020$nip.drag,2,function(x)quantile(x,prob=0.025)),
                                    apply(post2021$nip.drag,2,function(x)quantile(x,prob=0.025)),
                                    apply(post2022$nip.drag,2,function(x)quantile(x,prob=0.025))),
                          Upper = c(apply(post2020$nip.drag,2,function(x)quantile(x,prob=0.975)),
                                    apply(post2021$nip.drag,2,function(x)quantile(x,prob=0.975)),
                                    apply(post2022$nip.drag,2,function(x)quantile(x,prob=0.975))))


ggplot(df_nip_drag, aes(Year, Median, group = Group, color = Group)) +        
  geom_errorbar(aes(ymin = Lower, ymax = Upper),width=0.3,size=1) +
  geom_point(aes(group = Group),size=2) +
  geom_line(aes(group = Group), size=1) +
  facet_wrap(~Pair) +
  labs(title = "Nymphal Infection Prevalence (Drag) \n Control-vs-Treatment Site Pairs",
       y = "Posterior Median (95% Cr-I)", x = "Year") +
  scale_color_manual(breaks = c("Control", "Treatment"),
                     values=c("darkorange", "darkgreen"))



df_prot_mice <- data.frame(Pair = rep(c("P1","P2","P3","P3","P2","P1"),3),
                           Group = rep(c("Control","Control","Control",
                                         "Treatment","Treatment","Treatment"),3),
                          Year = c(rep("2020",6),rep("2021",6),rep("2022",6)),
                          Median = c(apply(post2020$protected.mice,2,median),
                                     apply(post2021$protected.mice,2,median),
                                     apply(post2022$protected.mice,2,median)),
                          Lower = c(apply(post2020$protected.mice,2,function(x)quantile(x,prob=0.025)),
                                    apply(post2021$protected.mice,2,function(x)quantile(x,prob=0.025)),
                                    apply(post2022$protected.mice,2,function(x)quantile(x,prob=0.025))),
                          Upper = c(apply(post2020$protected.mice,2,function(x)quantile(x,prob=0.975)),
                                    apply(post2021$protected.mice,2,function(x)quantile(x,prob=0.975)),
                                    apply(post2022$protected.mice,2,function(x)quantile(x,prob=0.975))))


ggplot(df_prot_mice, aes(Year, Median, group = Group, color = Group)) +        
  geom_errorbar(aes(ymin = Lower, ymax = Upper),width=0.3,size=1) +
  geom_point(aes(group = Group),size=2) +
  geom_line(aes(group = Group), size=1) +
  facet_wrap(~Pair) +
  labs(title = "Proportion of Mice with Protective OspA Levels \n Control-vs-Treatment Site Pairs",
       y = "Posterior Median (95% Cr-I)", x = "Year") +
  scale_color_manual(breaks = c("Control", "Treatment"),
                     values=c("darkorange", "darkgreen"))


########################################################################
# Comparison
########################################################################

# 2020 vs 2021
## Decrease in Mice Infection Prevalence

# Reduction
mean(post2020$prev.mice[,1] - post2021$prev.mice[,1])
mean(post2020$prev.mice[,6] - post2021$prev.mice[,6])

mean(post2020$prev.mice[,2] - post2021$prev.mice[,2])
mean(post2020$prev.mice[,5] - post2021$prev.mice[,5])

mean(post2020$prev.mice[,3] - post2021$prev.mice[,3])
mean(post2020$prev.mice[,4] - post2021$prev.mice[,4])

# Strength
mean(post2020$prev.mice[,1] > post2021$prev.mice[,1])
mean(post2020$prev.mice[,6] > post2021$prev.mice[,6])

mean(post2020$prev.mice[,2] > post2021$prev.mice[,2])
mean(post2020$prev.mice[,5] > post2021$prev.mice[,5])

mean(post2020$prev.mice[,3] > post2021$prev.mice[,3])
mean(post2020$prev.mice[,4] > post2021$prev.mice[,4])


## Decrease in Nymphal Infection Prevalence (Extracted)

# Reduction
mean(post2020$nip.mice[,1] - post2021$nip.mice[,1])
mean(post2020$nip.mice[,6] - post2021$nip.mice[,6])

mean(post2020$nip.mice[,2] - post2021$nip.mice[,2])
mean(post2020$nip.mice[,5] - post2021$nip.mice[,5])

mean(post2020$nip.mice[,3] - post2021$nip.mice[,3])
mean(post2020$nip.mice[,4] - post2021$nip.mice[,4])

# Strength
mean(post2020$nip.mice[,1] > post2021$nip.mice[,1])
mean(post2020$nip.mice[,6] > post2021$nip.mice[,6])

mean(post2020$nip.mice[,2] > post2021$nip.mice[,2])
mean(post2020$nip.mice[,5] > post2021$nip.mice[,5])

mean(post2020$nip.mice[,3] > post2021$nip.mice[,3])
mean(post2020$nip.mice[,4] > post2021$nip.mice[,4])

# Decrease in Nymphal Infection Prevalence (Dragged)

# Reduction
mean(post2020$nip.drag[,1] - post2021$nip.drag[,1])
mean(post2020$nip.drag[,6] - post2021$nip.drag[,6])

mean(post2020$nip.drag[,2] - post2021$nip.drag[,2])
mean(post2020$nip.drag[,5] - post2021$nip.drag[,5])

mean(post2020$nip.drag[,3] - post2021$nip.drag[,3])
mean(post2020$nip.drag[,4] - post2021$nip.drag[,4])

# Strength
mean(post2020$nip.drag[,1] > post2021$nip.drag[,1])
mean(post2020$nip.drag[,6] > post2021$nip.drag[,6])

mean(post2020$nip.drag[,2] > post2021$nip.drag[,2])
mean(post2020$nip.drag[,5] > post2021$nip.drag[,5])

mean(post2020$nip.drag[,3] > post2021$nip.drag[,3])
mean(post2020$nip.drag[,4] > post2021$nip.drag[,4])

# Increase in Proportion of Mice with Protected OspA Levels

# Increase
mean(post2021$protected.mice[,1] - post2020$protected.mice[,1])
mean(post2021$protected.mice[,6] - post2020$protected.mice[,6])

mean(post2021$protected.mice[,2] - post2020$protected.mice[,2])
mean(post2021$protected.mice[,5] - post2020$protected.mice[,5])

mean(post2021$protected.mice[,3] - post2020$protected.mice[,3])
mean(post2021$protected.mice[,4] - post2020$protected.mice[,4])

# Strength
mean(post2020$protected.mice[,1] < post2021$protected.mice[,1])
mean(post2020$protected.mice[,6] < post2021$protected.mice[,6])

mean(post2020$protected.mice[,2] < post2021$protected.mice[,2])
mean(post2020$protected.mice[,5] < post2021$protected.mice[,5])

mean(post2020$protected.mice[,3] < post2021$protected.mice[,3])
mean(post2020$protected.mice[,4] < post2021$protected.mice[,4])


# 2020 vs 2022
## Decrease in Mice Infection Prevalence

# Reduction
mean(post2020$prev.mice[,1] - post2022$prev.mice[,1])
mean(post2020$prev.mice[,6] - post2022$prev.mice[,6])

mean(post2020$prev.mice[,2] - post2022$prev.mice[,2])
mean(post2020$prev.mice[,5] - post2022$prev.mice[,5])

mean(post2020$prev.mice[,3] - post2022$prev.mice[,3])
mean(post2020$prev.mice[,4] - post2022$prev.mice[,4])

# Strength
mean(post2020$prev.mice[,1] > post2022$prev.mice[,1])
mean(post2020$prev.mice[,6] > post2022$prev.mice[,6])

mean(post2020$prev.mice[,2] > post2022$prev.mice[,2])
mean(post2020$prev.mice[,5] > post2022$prev.mice[,5])

mean(post2020$prev.mice[,3] > post2022$prev.mice[,3])
mean(post2020$prev.mice[,4] > post2022$prev.mice[,4])


## Decrease in Nymphal Infection Prevalence (Extracted)

# Reduction
mean(post2020$nip.mice[,1] - post2022$nip.mice[,1])
mean(post2020$nip.mice[,6] - post2022$nip.mice[,6])

mean(post2020$nip.mice[,2] - post2022$nip.mice[,2])
mean(post2020$nip.mice[,5] - post2022$nip.mice[,5])

mean(post2020$nip.mice[,3] - post2022$nip.mice[,3])
mean(post2020$nip.mice[,4] - post2022$nip.mice[,4])

# Strength
mean(post2020$nip.mice[,1] > post2022$nip.mice[,1])
mean(post2020$nip.mice[,6] > post2022$nip.mice[,6])

mean(post2020$nip.mice[,2] > post2022$nip.mice[,2])
mean(post2020$nip.mice[,5] > post2022$nip.mice[,5])

mean(post2020$nip.mice[,3] > post2022$nip.mice[,3])
mean(post2020$nip.mice[,4] > post2022$nip.mice[,4])

# Decrease in Nymphal Infection Prevalence (Dragged)

# Reduction
mean(post2020$nip.drag[,1] - post2022$nip.drag[,1])
mean(post2020$nip.drag[,6] - post2022$nip.drag[,6])

mean(post2020$nip.drag[,2] - post2022$nip.drag[,2])
mean(post2020$nip.drag[,5] - post2022$nip.drag[,5])

mean(post2020$nip.drag[,3] - post2022$nip.drag[,3])
mean(post2020$nip.drag[,4] - post2022$nip.drag[,4])

# Strength
mean(post2020$nip.drag[,1] > post2022$nip.drag[,1])
mean(post2020$nip.drag[,6] > post2022$nip.drag[,6])

mean(post2020$nip.drag[,2] > post2022$nip.drag[,2])
mean(post2020$nip.drag[,5] > post2022$nip.drag[,5])

mean(post2020$nip.drag[,3] > post2022$nip.drag[,3])
mean(post2020$nip.drag[,4] > post2022$nip.drag[,4])

# Increase in Proportion of Mice with Protected OspA Levels 

# Increase
mean(post2022$protected.mice[,1] - post2020$protected.mice[,1])
mean(post2022$protected.mice[,6] - post2020$protected.mice[,6])

mean(post2022$protected.mice[,2] - post2020$protected.mice[,2])
mean(post2022$protected.mice[,5] - post2020$protected.mice[,5])

mean(post2022$protected.mice[,3] - post2020$protected.mice[,3])
mean(post2022$protected.mice[,4] - post2020$protected.mice[,4])

# Strength
mean(post2020$protected.mice[,1] < post2022$protected.mice[,1])
mean(post2020$protected.mice[,6] < post2022$protected.mice[,6])

mean(post2020$protected.mice[,2] < post2022$protected.mice[,2])
mean(post2020$protected.mice[,5] < post2022$protected.mice[,5])

mean(post2020$protected.mice[,3] < post2022$protected.mice[,3])
mean(post2020$protected.mice[,4] < post2022$protected.mice[,4])
