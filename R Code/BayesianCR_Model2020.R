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

# Creating clusters
ncore <- 3    
cl <- makeCluster(ncore, outfile = "", type = "FORK")
clusterSetRNGStream(cl, iseed = 202303)
registerDoParallel(cl)

source(file = "Read_RTVField_Data_ExtendedVersion.R")

########################################################################
# Bayesian model
########################################################################

Model_Bayesian <- nimbleCode({
  
  # Priors ====================================================================
  
  # Site Membership and Data Augmentation (Mice)
  psi ~ dbeta(1,1)
  sigma2.site.member ~ dgamma(1,1)
  sigma2.site.encounter ~ dgamma(1,1)
  for(s in 1:nsites){
    beta.site[s] ~ dnorm(0, sd = sqrt(sigma2.site.member))
    alpha.site[s] ~ dnorm(0, sd = sqrt(sigma2.site.encounter))
  }
  
  # Subject-specific Covariate Effects
  alpha.adult ~ dnorm(0, sd = 1)  
  alpha.subadult ~ dnorm(0, sd = 1)
  alpha.male ~ dnorm(0, sd = 1)
  alpha.behavior ~ dnorm(0, sd = 1)
  
  # Prior for Missing Values on Covariates or Latent Process ==================
  p.adult ~ dbeta(1, 1)
  p.subadult ~ dbeta(1, 1)
  p.male ~ dbeta(1, 1)
  for(i in 1:M){
    AgeAdult[i] ~ dbern(p.adult)
    AgeSubAdult[i] ~ dbern(p.subadult)
    SexMale[i] ~ dbern(p.male)
  }
  
  
  # MODEL COMPONENTS ==============================================
  
  ### MICE
  # Site Membership
  for(s in 1:nsites){
    siteprob[s] ~ dbeta(1, 1)
  }
  
  # Observation Process
  for(i in 1:M){
    Site[i] ~ dcat(siteprob[1:nsites])
    z[i] ~ dbern(psi)
    site_membership[i] <- Site[i]*z[i]
    for(k in 1:nOccasions){
      logit(p[i,k]) <- alpha.site[Site[i]] + alpha.adult*AgeAdult[i] +
        alpha.subadult*AgeSubAdult[i] + alpha.male*SexMale[i] +
        alpha.behavior*Behavior[i,k]
      Encounter[i,k] ~ dbern(p[i,k] * z[i])  
    }
  }
  
  # Derived Parameters/Quantities
  # Mice Population Size
  for(s in 1:nsites) {N_mice_site[s] <- sum(site_membership[1:M] == s)}
  N_total <- sum(z[1:M])
  
})

########################################################################
# Objects 
########################################################################

parameters <- c("psi", "sigma2.site.member", "sigma2.site.encounter", 
                "alpha.adult", "alpha.site","alpha.subadult", 
                "alpha.male", "alpha.behavior", "p.adult", "p.subadult", 
                "p.male","N_mice_site", "siteprob","N_total")


nimbledata <- list(z = z_ind_2020,
                   Site = Site2020,
                   Encounter = Binary_Encounter2020,
                   AgeAdult = AgeAdult2020,
                   AgeSubAdult = AgeSubAdult2020,
                   SexMale = SexMale2020)

constants <- list(M = nrow(Binary_Encounter2020),
                  nOccasions = length(full_dates_2020),
                  nsites = nsites2020,
                  Behavior = Final_Behavior2020)

inits <- function() {
  
  # Missing Indices for Covariates and Outcomes
  miss.Site2020 <- which(is.na(Site2020))
  miss.z2020 <- which(is.na(z_ind_2020))
  miss.AgeAdult2020 <- which(is.na(AgeAdult2020))
  miss.AgeSubAdult2020 <- which(is.na(AgeSubAdult2020))
  miss.SexMale2020 <- which(is.na(SexMale2020))
  miss.Binary_Encounter2020 <- which(is.na(Binary_Encounter2020))
  
  # Objects
  z_init_2020 <- rep(NA,length(z_ind_2020))
  init.Site2020 <- rep(NA,length(Site2020))
  init.AgeAdult2020 <- rep(NA,length(AgeAdult2020))
  init.AgeSubAdult2020 <- rep(NA,length(AgeSubAdult2020))
  init.SexMale2020 <- rep(NA,length(SexMale2020))
  init.Binary_Encounter2020 <- matrix(NA,nrow(Binary_Encounter2020),
                                      ncol(Binary_Encounter2020))
  
  # Initial Values for Missing Values
  init.Site2020[miss.Site2020] <- sample(1:6,length(miss.Site2020),replace = TRUE)
  z_init_2020[miss.z2020] <- rep(1,length(miss.z2020))
  init.AgeAdult2020[miss.AgeAdult2020] <- rbinom(length(miss.AgeAdult2020),size = 1, prob = 0.5)
  init.AgeSubAdult2020[miss.AgeSubAdult2020] <- 1 - init.AgeAdult2020[miss.AgeAdult2020]
  init.SexMale2020[miss.SexMale2020] <- rbinom(length(miss.SexMale2020),size = 1, prob = 0.5)
  init.Binary_Encounter2020[miss.Binary_Encounter2020] <- rbinom(length(miss.Binary_Encounter2020),
                                                                 size = 1, prob = 0.5)
  
  list(alpha.adult = runif(1), 
       alpha.subadult = runif(1), 
       alpha.male = runif(1),
       alpha.behavior = runif(1), 
       sigma2.site.member = runif(1),
       sigma2.site.encounter = runif(1),
       p.adult = runif(1), 
       p.subadult = runif(1), 
       p.male = runif(1),
       psi = runif(1),
       siteprob = runif(6), 
       alpha.site = runif(6),
       z = z_init_2020,
       Site = init.Site2020,
       Encounter = init.Binary_Encounter2020,
       AgeAdult = init.AgeAdult2020,
       AgeSubAdult = init.AgeSubAdult2020,
       SexMale = init.SexMale2020)
}


########################################################################
# Running model
########################################################################

time1 <- Sys.time() 
fit <- foreach(x = 1:ncore, .packages = "nimble", .verbose = TRUE) %dopar% {
  nimbleMCMC(code = Model_Bayesian,
             constants = constants,
             data = nimbledata,
             monitors = parameters,
             inits = inits(),
             niter = 30000, 
             nburnin = 5000, 
             nchains = 1,
             thin = 5,
             progressBar = TRUE,
             summary = TRUE,
             samplesAsCodaMCMC = TRUE)
}
time2 <- Sys.time() 
(runtime <- time2-time1)
set_idx <- sample(1:1000,1)
cat(paste0("File name is: fit_CR2020_idx",set_idx,".rds"))
saveRDS(fit, paste0("./fit_CR2020_idx",set_idx,".rds"))
stopCluster(cl)


########################################################################
# Posterior summary
########################################################################

want_summary <- FALSE

if(want_summary == TRUE){
  results_file2020 <- readRDS("fit_CR2020_idx523_official.rds")
  results_mcmc2020 <- as.mcmc.list(lapply(1:3, 
                                          function(x){as.mcmc(results_file2020[[x]]$samples)}))
  par.names2020 <- colnames(results_mcmc2020[[1]])
  MCMCsummary(results_mcmc2020,round = 2)[,c(1,2,4,3,5)]
  
  MCMCtrace(results_mcmc2020, 
            params = par.names2020,
            ISB = FALSE,
            pdf = TRUE,
            Rhat = TRUE)
  
  # Ticks size
  mult.factor <- max(c(median(Final_Ticks2020[Final_Ticks2020>0],na.rm = TRUE),
                       median(Final_Ticks2021[Final_Ticks2021>0],na.rm = TRUE),
                       median(Final_Ticks2022[Final_Ticks2022>0],na.rm = TRUE)))
  cc <- as.mcmc.list(lapply(1:3, 
                            function(x){as.mcmc(mult.factor*results_file2020[[x]]$samples)}))
  dft <- MCMCsummary(cc,round = 0)[1:7,c(1,2, 4,3,5)]
  rownames(dft) <- c(paste0("N_ticks_site[",1:6,"]"),"N_total_ticks")
  dft
  
  mcmc_comb_2020 <- combine.mcmc(results_mcmc2020)
  
  # effect on probability of captures
  # behavior effect
  mean(mcmc_comb_2020[,"alpha.behavior"] < 0)

    
}

