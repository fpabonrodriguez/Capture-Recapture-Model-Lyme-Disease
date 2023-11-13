#########################################################################
### Felix Pabon-Rodriguez
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

# Creating clusters
ncore <- 3    

collected.data <- TRUE 
if(collected.data == TRUE){ # official data (collected)
  #cl <- makeCluster(ncore, outfile = "", type = "FORK")
  cl <- makeCluster(ncore, outfile = "DP2020.log")
  clusterSetRNGStream(cl, iseed = 112301) # official
  registerDoParallel(cl)
  source(file = "Read_RTVField_Data_ExtendedVersion.R")
  
}else{ # permuted/perturbed data (for still blinded team members)
  #cl <- makeCluster(ncore, outfile = "", type = "FORK")
  cl <- makeCluster(ncore, outfile = "DP2020.log")
  clusterSetRNGStream(cl, iseed = 112301) # perturbed
  registerDoParallel(cl)
  source(file = "Read_RTVField_Data_ExtendedVersion_Modified.R")
} 


########################################################################
# Bayesian model
########################################################################

Model_Bayesian <- nimbleCode({
  
  # Priors ====================================================================
  # Infectious Status (Mice)
  # Protective OspA Antibody (Mice)
  # Ticks Count (Mice Skin and Dragged from Site)
  # Infected Ticks (Mice Skin)
  # Infected Ticks (Dragging)
  
  sigma2.site.infected.mice ~ dgamma(1,1)
  sigma2.site.protectiveOspA ~ dgamma(1,1)
  sigma2.site.dragged.ticks ~ dgamma(1,1)
  sigma2.site.infected.dragged.ticks ~ dgamma(1,1)
  for(s in 1:nsites){
    omega.site[s] ~ dnorm(0, sd = sqrt(sigma2.site.infected.mice))
    delta.site[s] ~ dnorm(0, sd = sqrt(sigma2.site.protectiveOspA))
    lambda.site[s] ~ dnorm(0, sd = sqrt(sigma2.site.dragged.ticks))
    theta.site[s] ~ dnorm(0, sd = sqrt(sigma2.site.infected.dragged.ticks))
  }
  sigma2.subj.infected.mice ~ dgamma(1,1)
  sigma2.subj.infticks.mice ~ dgamma(1,1)
  sigma2.subj.protected.mice ~ dgamma(1,1)
  sigma2.subj.infticks.drag ~ dgamma(1,1)
  for(i in 1:N){
    omega.subj[i] ~ dnorm(0, sd = sqrt(sigma2.subj.infected.mice))
    kappa.subj[i] ~ dnorm(0, sd = sqrt(sigma2.subj.infticks.mice))
    gamma.subj[i] ~ dnorm(0, sd = sqrt(sigma2.subj.protected.mice))
    delta.subj[i] ~ dnorm(0, sd = sqrt(sigma2.subj.infticks.drag))
  }
  
  for(i in 1:3){
    delta[i] ~ dbeta(1, 1)
  }
  for(i in 1:2){
    omega[i] ~ dbeta(1, 1)
  }
  sigma2.ExpTicks ~ dgamma(1,1)
  sigma2.ExpDraggedTicks ~ dgamma(1,1)
  rho.ticks ~ dnorm(0, sd = 1/20)
  rho.dragged ~ dnorm(0, sd = 1/4)
  
  # Prior for Missing Values on Covariates or Latent Process ==================
  for(i in 1:N){
    for(k in 1:nOccasions){
      PercentAte[i,k] ~ dbeta(1, 1)
    }
  }
  
  p.OspA_t1 ~ dbeta(1, 1)
  p.InfMice_t1 ~ dbeta(1, 1)
  lam.ExpTicks_t1 ~ dgamma(1,5)
  lam.NymphTicks_t1 ~ dgamma(1,5)
  lam.InfNTicks_t1 ~ dgamma(1,5)
  for(i in 1:N){
    ExpTicks[i,1] ~ dpois(lam.ExpTicks_t1)
    NymphalTicks[i,1] ~ dpois(lam.NymphTicks_t1)
    Infected_Ticks_Mice[i,1] ~ dpois(lam.InfNTicks_t1)
    Infected_Mice[i,1] ~ dbern(p.InfMice_t1)
    ProtectiveOspA[i,1] ~ dbern(p.OspA_t1)
  }
  lam.DragTicks_t1 ~ dgamma(1,5)
  lam.ExpDragTicks_t1 ~ dgamma(1,5)
  lam.DragInfTicks_t1 ~ dgamma(1,5)
  for(s in 1:nsites){
    DraggedTicks[s,1] ~ dpois(lam.DragTicks_t1)
    ExpDraggedTicks[s,1] ~ dpois(lam.ExpDragTicks_t1)
    Infected_Ticks_Drag[s,1] ~ dpois(lam.DragInfTicks_t1)
  }
  
  # MODEL COMPONENTS ==============================================
  
  # Disease Process
  for(i in 1:N){
    for(k in (firstcapture[i]+1):(lastcapture[i])){
      # Infectious Status
      Infected_Mice[i,k] ~ dbern(p.infected.mice[i,k])
      logit(p.infected.mice[i,k]) <- omega.subj[i] + 
        omega.site[Site[i]] +
        omega[1]*Infected_Mice[i,k-1] + 
        omega[2]*ProtectiveOspA[i,k-1]
      
      # Protection Status of OspA Antibody 
      ProtectiveOspA[i,k] ~ dbern(p.protectiveOspA[i,k])
      logit(p.protectiveOspA[i,k]) <- delta.subj[i] + 
        delta.site[Site[i]] + 
        delta[1]*PercentAte[i,k-1] + 
        delta[2]*Infected_Mice[i,k-1] + 
        delta[3]*ProtectiveOspA[i,k-1]
      
      # Number of Nymphal Ticks (Mice Skin)
      NymphalTicks[i,k] ~ dpois(ExpTicks[i,k])
      ExpTicks[i,k] ~ dnorm(mu.NTicks[i,k], sd = sqrt(sigma2.ExpTicks))
      mu.NTicks[i,k] <- gamma.subj[i] + rho.ticks*ExpTicks[i,k-1]
      
      # Infected Nymphal Ticks (Mice Skin)
      Infected_Ticks_Mice[i,k] ~ dbinom(size = NymphalTicks[i,k], 
                                        prob = p.ticks.mice[i,k])
      logit(p.ticks.mice[i,k]) <- kappa.subj[i] 
    }
  }
  
  ### TICKS (Dragged)
  for(s in 1:nsites){
    for(k in 2:nOccasions){
      # Number of Ticks
      DraggedTicks[s,k] ~ dpois(ExpDraggedTicks[s,k])
      ExpDraggedTicks[s,k] ~ dnorm(lambda.site[s] + rho.dragged*ExpDraggedTicks[s,k-1], 
                                   sd = sqrt(sigma2.ExpDraggedTicks))
      
      # Infected Nymphal Ticks 
      Infected_Ticks_Drag[s,k] ~ dbinom(size = DraggedTicks[s,k], 
                                        prob = p.ticks.drag[s,k])
      logit(p.ticks.drag[s,k]) <- theta.site[s]
    }
  }
  
})


########################################################################
# Objects 
########################################################################

parameters <- c("sigma2.site.dragged.ticks", "sigma2.subj.infticks.drag",
                "sigma2.site.infected.mice","sigma2.site.protectiveOspA", 
                "sigma2.site.infected.dragged.ticks", "omega.site",
                "delta.site", "theta.site", "lambda.site", "omega",
                "delta", "sigma2.ExpTicks", "rho.ticks", "rho.dragged",
                "sigma2.subj.infected.mice", "sigma2.subj.infticks.mice", 
                "sigma2.ExpDraggedTicks", "sigma2.subj.protected.mice")

# Use observed data
N.obs <- length(unique(datayear2020$UniqueID))

nimbledata <- list(NymphalTicks = Final_Ticks_Mice2020,
                   DraggedTicks = Final_DraggedTicks2020,
                   Infected_Ticks_Mice = Final_InfTicks_Mice2020,
                   Infected_Ticks_Drag = Final_InfTicks_Drag2020,
                   Infected_Mice = Final_InfectedMice2020,
                   PercentAte = Final_PercentAte2020,
                   ProtectiveOspA = Final_ProtectiveOspA2020)

constants <- list(N = N.obs,
                  nsites = nsites2020,
                  nOccasions = length(full_dates_2020),
                  firstcapture = firstcapture2020[1:N.obs],
                  lastcapture = lastcapture2020[1:N.obs],
                  Site = Site2020[1:N.obs],
                  AgeAdult = AgeAdult2020[1:N.obs],
                  AgeSubAdult = AgeSubAdult2020[1:N.obs],
                  SexMale = SexMale2020[1:N.obs])

inits <- function() {
  
  # Missing Indices for Covariates and Outcomes
  miss.NymphalTicks2020 <- which(is.na(Final_Ticks_Mice2020))
  miss.DraggedTicks_2020 <- which(is.na(Final_DraggedTicks2020))
  miss.InfectedTicks_Mice2020 <- which(is.na(Final_InfTicks_Mice2020))
  miss.InfectedTicks_Drag2020 <- which(is.na(Final_InfTicks_Drag2020))
  miss.Infected2020 <- which(is.na(Final_InfectedMice2020))
  miss.ProtectiveOspA2020 <- which(is.na(Final_ProtectiveOspA2020))
  miss.PercentAte2020 <- which(is.na(Final_PercentAte2020))
  
  # Objects
  init.NymphalTicks2020 <- matrix(NA,nrow(Final_Ticks_Mice2020),ncol(Final_Ticks_Mice2020))
  init.DraggedTicks_2020 <- matrix(NA,nrow(Final_DraggedTicks2020),ncol(Final_DraggedTicks2020))
  init.InfectedTicks_Mice2020 <- matrix(NA,nrow(Final_InfTicks_Mice2020),ncol(Final_InfTicks_Mice2020))
  init.InfectedTicks_Drag2020 <- matrix(NA,nrow(Final_InfTicks_Drag2020),ncol(Final_InfTicks_Drag2020))
  init.Infected2020 <- matrix(NA,nrow(Final_InfectedMice2020),ncol(Final_InfectedMice2020))
  init.ProtectiveOspA2020 <- matrix(NA,nrow(Final_ProtectiveOspA2020),ncol(Final_ProtectiveOspA2020))
  init.PercentAte2020 <- matrix(NA,nrow(Final_PercentAte2020),ncol(Final_PercentAte2020))
  
  # Initial Values for Missing Values
  init.NymphalTicks2020[miss.NymphalTicks2020] <- rep(5,length(miss.NymphalTicks2020))
  init.DraggedTicks_2020[miss.DraggedTicks_2020] <- rep(5,length(miss.DraggedTicks_2020))
  init.InfectedTicks_Mice2020[miss.InfectedTicks_Mice2020] <- rep(0,length(miss.InfectedTicks_Mice2020))
  init.InfectedTicks_Drag2020[miss.InfectedTicks_Drag2020] <- rep(0,length(miss.InfectedTicks_Drag2020))
  init.Infected2020[miss.Infected2020] <- rep(0,length(miss.Infected2020))
  init.ProtectiveOspA2020[miss.ProtectiveOspA2020] <- rep(0,length(miss.ProtectiveOspA2020))
  init.PercentAte2020[miss.PercentAte2020] <- round(runif(length(miss.PercentAte2020),0,1),2)
  
  
  list(sigma2.ExpTicks = runif(1),
       sigma2.ExpDraggedTicks = runif(1),
       sigma2.site.dragged.ticks = runif(1),
       sigma2.site.infected.mice = runif(1),
       sigma2.site.protectiveOspA = runif(1),
       sigma2.site.infected.dragged.ticks = runif(1),
       sigma2.subj.infected.mice = runif(1),
       sigma2.subj.infticks.mice = runif(1),
       sigma2.subj.protected.mice = runif(1),
       sigma2.subj.infticks.drag = runif(1),
       rho.ticks = runif(1,0.04,0.09),
       rho.dragged = runif(1),
       delta = runif(3),
       omega = runif(2),
       kappa.subj = runif(N.obs),
       omega.subj = runif(N.obs),
       gamma.subj = runif(N.obs),
       delta.subj = runif(N.obs),
       omega.site = runif(6), 
       lambda.site = runif(6), 
       delta.site = runif(6), 
       theta.site = runif(6), 
       NymphalTicks = init.NymphalTicks2020,
       DraggedTicks = init.DraggedTicks_2020,
       Infected_Mice = init.Infected2020,
       Infected_Ticks_Mice = init.InfectedTicks_Mice2020,
       Infected_Ticks_Drag = init.InfectedTicks_Drag2020,
       PercentAte = init.PercentAte2020,
       ProtectiveOspA = init.ProtectiveOspA2020,
       ExpTicks = matrix(data = round(runif(nrow(Final_Ticks2020)*ncol(Final_Ticks2020),5,25)),
                         ncol = ncol(Final_Ticks2020), nrow = nrow(Final_Ticks2020)),
       ExpDraggedTicks = matrix(data = round(runif(nrow(Final_DraggedTicks2020)*ncol(Final_DraggedTicks2020),5,25)),
                                ncol = ncol(Final_DraggedTicks2020), nrow = nrow(Final_DraggedTicks2020)),
       p.OspA_t1 = runif(1),
       p.InfMice_t1 = runif(1),
       lam.InfNTicks_t1 = runif(1),  
       lam.ExpTicks_t1 = runif(1),
       lam.NymphTicks_t1 = runif(1),
       lam.DragTicks_t1 = runif(1),
       lam.ExpDragTicks_t1 = runif(1), 
       lam.DragInfTicks_t1 = runif(1))
}

########################################################################
# Running model
########################################################################

time1 <- Sys.time() 
fit <- foreach(x = 1:ncore, .packages = "nimble", .verbose = TRUE) %dopar% {
  nimbleMCMC(code = Model_Bayesian,
             constants = constants,
             data = nimbledata,
             inits = inits(),
             monitors = parameters,
             niter = 300000,  
             nburnin = 100000, 
             nchains = 1,
             thin = 10,
             progressBar = TRUE,
             summary = TRUE,
             samplesAsCodaMCMC = TRUE)
}
time2 <- Sys.time() 
(runtime <- time2-time1)
set_idx <- sample(1:1000,1)
cat(paste0("File name is: fit_DP2020_idx",set_idx,".rds"))
saveRDS(fit, paste0("./fit_DP2020_idx",set_idx,".rds"))
stopCluster(cl)


########################################################################
# Posterior summary
########################################################################

want_summary <- FALSE

if(want_summary == TRUE){
  results_file2020 <- readRDS("fit_DP2020_idx621.rds")
  results_mcmc2020 <- as.mcmc.list(lapply(1:3, 
                                          function(x){as.mcmc(results_file2020[[x]]$samples)}))
  par.names2020 <- colnames(results_mcmc2020[[1]])
  
  MCMCtrace(results_mcmc2020,
            filename = "trace_IDP_2020.pdf", 
            params = par.names2020,
            ISB = FALSE,
            pdf = TRUE,
            Rhat = TRUE,
            iter = (300000-100000)/10)
  
  ## Prevalence
  # if     p ~ dunif(0, 1)
  # and    y ~ dbinom(p,n)
  # then,  p|y ~ beta(a+y,b+N-y)
  sim <- 10000
  prev.mice <- matrix(NA,sim,nsites2020)
  protected.mice <- matrix(NA,sim,nsites2020)
  nip.mice <- matrix(NA,sim,nsites2020)
  nip.drag <- matrix(NA,sim,nsites2020)
  for(s in 1:nsites2020){
    prev.mice[,s] <- rbeta(sim,1+Final_Incidence_Mice_2020[s],
                           1+SampleSize2020[s]-Final_Incidence_Mice_2020[s])
    protected.mice[,s] <- rbeta(sim,1+Final_Protection_Mice_2020[s],
                                1+SampleSize2020[s]-Final_Protection_Mice_2020[s])
    nip.mice[,s] <- rbeta(sim,1+Final_Nymphals_Mice_2020[s],
                          1+SampleNymphals2020[s]-Final_Nymphals_Mice_2020[s])
    nip.drag[,s] <- rbeta(sim,1+Final_Nymphals_Drag_2020[s],
                          1+SampleDrag2020[s]-Final_Nymphals_Drag_2020[s])
  }
  
  outcomes2020 <- list(prev.mice=prev.mice,
                       protected.mice=protected.mice,
                       nip.mice=nip.mice,
                       nip.drag=nip.drag)
  
  saveRDS(object = outcomes2020,file = "disease_outcomes_2020.rds")
  
}















