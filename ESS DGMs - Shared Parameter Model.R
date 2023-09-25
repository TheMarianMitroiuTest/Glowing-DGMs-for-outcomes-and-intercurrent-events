## License CC-BY-4.0 
## Creative Commons Attribution 4.0 International
## Mitroiu M. et al TBC



## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"


### DGM
#rm(list=ls())
library(gmailr)
library(MASS)
library(nlme)
library(lme4)
library(survival)
library(foreign)
library(tidyverse)
library(janitor)
library(gt)
library(patchwork)



## gmail setup----
# setup to receive e-mails with results of the simulations. Useful to store results, but most importantly to be notified when the simulation is concluded.
# various tutorials can be found to set this up

#google_app <- httr::oauth_app(
#  "renamedapp",
#  key = "xx",
#  secret = "xx"
#)
#
#gm_auth_configure(key = "xx",
#                  secret = "xx")



options(mc.cores = parallel::detectCores()) # M1 with parallel, Old1 without parallel

Scenario <- c("A")

# simulation


## simulation parameters #----

set.seed(2147483629)
n <- 190 # number of patients to be simulated (sample size)
# this is based on a t-test to ensure  90% power at alpha level=0.025 one-sided 

m.iterations <- 500 #416 for the verification of the longitudinal outcomes # number of generated datasets 

# these will be used in the target proportions of intercurrent events used in the function to determine the intercept value in order to obtain the right percentage of intercurrent events
# ranges of probabilities centered around desired percentages of each intercurrent events averaged over all simulated trials
# this is done to increase variability in intercurrent events percentages between trials
p_LoE_sample <-  0.35 #c(0.31, 0.33, 0.34, 0.35, 0.37); mean(p_LoE_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy at trial level
p_AE_Exp_sample <- 0.10 # c(0.055, 0.065, 0.075, 0.085, 0.095)*2; mean(p_AE_Exp_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy in the experimental arm
p_AE_Control_sample <- 0.05 # c(0.05, 0.10, 0.15, 0.20, 0.25)/2; mean(p_AE_Control_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy in the control arm

overlap.adjust <- 1.16 # to adjust for AE overlap with LoE



visits <- as.numeric(c(0, 1, 2, 3, 4, 5, 6))	

delta_spm <- matrix(ncol=1,nrow=m.iterations) # object to store treatment effect estimates at 6 weeks based on MMRM model fitted on each generated dataset
colnames(delta_spm) <-c("TreatmentEffect")
betas_spm <- matrix(ncol=2,nrow=m.iterations) # object to store parameters for the treatment effect at week 6 based on the MMRM model fitted on each generated dataset
colnames(betas_spm) <-c("Treat", "visit42:Treat")

treatmenteffect <- -3.5 #assumed true effect at week 6

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3) # progress bar in percentages relative to the total number of m.iterations

#confint_fit <- matrix(ncol=2,nrow=m.iterations) # object to store the 95% confidence interval bounds for the estimated treatment effect
#colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
#delta_spm_errorz <- matrix(ncol=1,nrow=m.iterations) # standard error of the estimated treatment effect
#colnames(delta_spm_errorz) <- c("SE")

#bias_f <- matrix(ncol=1,nrow=m.iterations) # object to store the bias (estimated treatment effect - true treatment effect)
#colnames(bias_f) <- c("bias_f")

## Randomisation objects
N_Exp  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Exp) <-c("N randomised Exp")

N_Control  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Control) <-c("N randomised Control")

# Intercurrent events objects
# lack of efficacy (LoE)
n_LoE_total <- matrix(ncol=1,nrow=m.iterations) # object to store the number of patients experiencing e.g., treatment discontinuation due to lack of efficacy at trial level
colnames(n_LoE_total) <- c("N LoE Total")

n_LoE_Exp <- matrix(ncol=1,nrow=m.iterations) # object to store the number of patients experiencing e.g., treatment discontinuation due to lack of efficacy in the experimental arm
colnames(n_LoE_Exp) <- c("N LoE Exp")

n_LoE_Control <- matrix(ncol=1,nrow=m.iterations) # object to store the number of patients experiencing e.g., treatment discontinuation due to lack of efficacy in the control arm
colnames(n_LoE_Control) <- c("N LoE Control")


# Adverse events (AE)
n_AE_total <- matrix(ncol=1,nrow=m.iterations) # object to store the number of patients experiencing e.g., treatment discontinuation due to adverse events at trial level
colnames(n_AE_total) <- c("N AE Total")

n_AE_Exp <- matrix(ncol=1,nrow=m.iterations) # object to store the number of patients experiencing e.g., treatment discontinuation due to adverse events in the experimental arm
colnames(n_AE_Exp) <- c("N AE Exp")

n_AE_Control <- matrix(ncol=1,nrow=m.iterations) # object to store the number of patients experiencing e.g., treatment discontinuation due to adverse events in the control arm
colnames(n_AE_Control) <- c("N AE Control")

# AE + LoE Total
n_AE_and_LoE_T <- matrix(ncol=1,nrow=m.iterations) # object to store the total number of patients experiencing e.g., treatment discontinuation due to lack of efficacy or adverse events at trial level
colnames(n_AE_and_LoE_T) <- c("N AE and LoE Total")

# For percentages
# LoE
# below are objects to store percentages data following the above structure
LoE_total_Perc <- matrix(ncol=1,nrow=m.iterations) # 
colnames(LoE_total_Perc) <- c("% LoE Total")

LoE_Exp_Perc <- matrix(ncol=1,nrow=m.iterations) # 
colnames(LoE_Exp_Perc) <- c("% LoE Exp")

LoE_Control_Perc<- matrix(ncol=1,nrow=m.iterations) # 
colnames(LoE_Control_Perc) <- c("% LoE Control")


#AE
AE_total_Perc <-  matrix(ncol=1,nrow=m.iterations) # 
colnames(AE_total_Perc) <- c("% AE Total")

AE_Exp_Perc <- matrix(ncol=1,nrow=m.iterations) # 
colnames(AE_Exp_Perc) <- c("% AE Exp")

AE_Control_Perc <-matrix(ncol=1,nrow=m.iterations) # 
colnames(AE_Control_Perc) <- c("% AE Control")

# AE + LoE percentage
AE_and_LoE_Perc <- matrix(ncol=1,nrow=m.iterations) # 
colnames(AE_and_LoE_Perc) <- c("% AE and LoE Total")

#pb3 <- txtProgressBar(min = 0,  max=length(scaling_factor), style=3) # # progress bar in percentages relative to the total number of scaling factors


# insert model specification from fitted models on the SM super trial



## Suppose we want to simulate binomial data using a logistic regression model of the form logit(Y) ~ b0 + b1 * x
## where x is available and b1 is fixed, but we want to tweak b0 such that the expected proportion of '1's equals
## a prespecified target, the following function can perhaps be a solution (but note that it gives an approximate solution):
tweak.intercept <- function(target, x = x, beta1, K = 10, grid = seq(-200, 200, length.out = 400), show.plot = TRUE){
  res <- rep(NA, length(grid))
  for(i in 1:length(res))    res[i] <- mean(replicate(K, mean(rbinom(length(x), 1, plogis(grid[i] + beta1 * x)))))
  if(show.plot) plot(res ~ grid, ylim = c(0, 1))
  d <- data.frame(intercept = grid, prop1 = res)
  int <- d$intercept[abs(d$prop1 - target) == min(abs(d$prop1 - target))]
  if(show.plot) abline(h = target, col = 2)
  if(show.plot) abline(v = int, col = 2)
  return(int)
}


# Determine intercept for % of LoE at trial level according to the target %
x <- rbinom(10000, 1, 0.5)

c1 <- -0.60629 # coefficient for treatment in the logit model for LoE
# c1 comes from fit_LoE_spm
det.intercept_LoE <- tweak.intercept(target = (p_LoE_sample*overlap.adjust), x = x, beta1 = c1, grid = seq(-5, 5, length.out = 400)) ; det.intercept_LoE
test_LoE <- replicate(10000, mean(rbinom(length(x), 1, plogis(det.intercept_LoE + c1 * x))))
hist(test_LoE)
mean(test_LoE) # at trial level is divided by 2

#describe(d$bi_0)
#hist(d$bi_0)

# Determine intercept for % of AE in experimental arm according to the target %
# should be around the intercept coming from fit_AE_exp_spm
det.intercept_AE_exp <- tweak.intercept(target = (p_AE_Exp_sample*2), x = x, beta1 = 0, grid = seq(-20, 20, length.out = 400)); det.intercept_AE_exp
# p_AE_Exp_sample is multiplied by 2 to adjust for the trial size, not just the arm size (reference is the trial sample size)
test_AE_exp <- replicate(10000, mean(rbinom(length(x), 1, plogis(det.intercept_AE_exp))))
hist(test_AE_exp)
mean(test_AE_exp)# # at trial level is divided by 2


# Determine intercept for % of AE in control arm according to the target %
# fit_AE_control_spm
det.intercept_AE_control <- tweak.intercept(target = (p_AE_Control_sample*2), x = x, beta1 = 0, grid = seq(-20, 20, length.out = 400)); det.intercept_AE_control
# p_AE_Control_sample is multiplied by 2 to adjust for the trial size, not just the arm size
test_AE_control <- replicate(10000, mean(rbinom(length(x), 1, plogis(det.intercept_AE_control))))
hist(test_AE_control)
mean(test_AE_control)




start_time <- Sys.time() # timestamp for the start time of the nested for loop below.
# it was used to have an estimate of time needed for different larger number of trials to be simulated upon scaling up the simulation parameters (e.g., m.iterations)


## Begin for loop----
set.seed(2147483629) # resetting the seed here to get exactly the same longitudinal outcomes in this DGM and in the JM DGM

for(m in 1:m.iterations) {
  
  
  # The specification of the model: random intercept and random slope + some random noise
  # d$MADRS10.true <- (b0 + d$bi_0) + (b1 + d$bi_1) * d$visit +  b2 * d$visit * d$Treat
  # d$MADRS10_collected <- d$MADRS10.true + rnorm(nrow(d), 0, eps.sd)
  
  
  ### Generate longitudinal outcomes----
  #### Generate correlated (random) intercepts and slopes for each individual patient----
  
  # We used as inspiration the trial 003-021 (ref data analysis paper). Here we used a simplified scenario with linear group trajectories.

  
  #### model parameters----
  
  b0 <-  29.374914 #29.79 # intercept
  b1 <-  -0.529752 # -0.55  #slope
  b2 <-  -0.577650 # -0.583 treatment effect 
  # b0, b1 and b2 are taken from a SOURCE TRIAL simulated with SMd DGM with n=8000 patients and re_covm3
  # based on these, errors (eps.sd) are added and the covariance matrix (bi_covm) below is used for the random effects in order to simulate trials.
  
  
  #Fixed effects:  MADRS10 ~ visit + visit:Treat 
  #Value  Std.Error    DF  t-value p-value
  #(Intercept) 29.374914 0.06088202 47998 482.4892       0
  #visit       -0.529752 0.02041562 47998 -25.9484       0
  #visit:Treat -0.577650 0.02898964 47998 -19.9261       0
  
  
  # the bi_covm and eps.sd are taken from a model fitted on a SOURCE TRIAL with n=4000 patients and re_covm3
  bi_means <- c(0, 0)
  bi_covm <- matrix(c(21.140783, 1.307053, 1.307053, 1.044141), nrow = 2)
  
  #bi_covm <- matrix(c(24.611, 0.5869809, 0.5869809, 1.157), nrow = 2) original
  bi <- mvrnorm(n, bi_means, bi_covm)	# generate random effects for n patients, with bi_means and covariance bi_covm
  eps.sd <-3.285955 # 3.247  # residual error
  
  #summary(fit_lme)
  # you should get:
  #>         summary(fit_lme)
  #Linear mixed-effects model fit by REML
  #Data: d 
  #AIC      BIC    logLik
  #327614.9 327677.4 -163800.4
  
  #Random effects:
  #  Formula: ~1 + visit | id
  #Structure: General positive-definite, Log-Cholesky parametrization
  #StdDev   Corr  
  #(Intercept) 4.963853 (Intr)
  #visit       1.142975 0.101 
  #Residual    3.285955       
  
  #Fixed effects:  MADRS10 ~ visit + visit:Treat 
  #Value  Std.Error    DF  t-value p-value
  #(Intercept) 29.374914 0.06088202 47998 482.4892       0
  #visit       -0.529752 0.02041562 47998 -25.9484       0
  #visit:Treat -0.577650 0.02898964 47998 -19.9261       0
  #Correlation: 
  #  (Intr) visit 
  #visit       -0.059       
  #visit:Treat  0.000 -0.702
  
  #Standardized Within-Group Residuals:
  #  Min            Q1           Med            Q3 
  #-3.6560020890 -0.5761900126 -0.0005035796  0.5794871856 
  #Max 
  #3.9895106728 
  
  #Number of Observations: 56000
  #Number of Groups: 8000 
  
  #View(d)
  #mean(d$Baseline)
  #vcov(fit_lme)
  #VarCorr(fit_lme)
  #var(random.effects(fit_lme))
  #            (Intercept)    visit
  #(Intercept)   21.140783 1.307053
  #visit          1.307053 1.044141
  

  d <- data.frame(
    id = rep(1:n, each = length(visits)),
    visit = visits,
    Treat = rep(rbinom(n, 1, 0.5), each = length(visits)),
    bi_0 = rep(bi[,1], each = length(visits)),
    bi_1 = rep(bi[,2], each = length(visits))
  )
  
  # generate the true values of the outcome
  d$MADRS10.true <- (b0 + d$bi_0) + (b1 + d$bi_1) * d$visit +  b2 * d$visit * d$Treat 
  
  # generate the what would be observed values of the outcome = true values + some error 
  d$MADRS10_collected <- d$MADRS10.true + rnorm(nrow(d), 0, eps.sd)
  
  #describe(d$bi_0)
  hist(d$bi_0)
  
  #describe(d$bi_1)
  hist(d$bi_1)
  
  # Plot trajectories
  
  p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id)) 
  p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
    scale_y_continuous(limits = c(-10, 60))
  
  
  # double check with lmer()
  #fit_lmer <- lmer(MADRS10_collected ~ visit + visit:Treat + (1 + visit|id), data = d, REML = T)
  #summary(fit_lmer)
  
  #m<-1
  
  fit_lme <- lme(fixed=MADRS10_collected ~ visit + visit:Treat, 
                random=~1 + visit| id,
                method="REML", 
                data=d)
  
  #std.e <- 0.1749478*length(visits)
  
  #summary(fit_lme)
  betas_spm[m, ] <- fixef(fit_lme)[2:3]
  delta_spm[m, ] <- fixef(fit_lme)[3] * max(visits)
  #-0.788534/0.1696878
  
  #0.1696878*6
  
  #getVarCov(fit_lme, type= "marginal")
  
  #getVarCov(fit_lme,type = c("conditional"))

  #getVarCov(fit_lme,type = c("random.effects"))
  
  #vcov(fit_lme)
  
  #sqrt(9.0522)
  
  #(4*3^2)/(0.1^2)
  #3600 trials vs 502 in the MMRM
  
  
  #m<-1
  # store the number of patients in the objects defined a priori
  Randomised_Exp <- sum(as.numeric(d[,3]))/7 #number of patients in the experimental arm # make sure it is divided by the number of visits (6/7), without or with the baseline included as response
  Randomised_Control <- n-sum(as.numeric(d[,3]))/7  #number of patients in the control arm
  
  N_Exp[m,] <- Randomised_Exp
  N_Control[m,] <- Randomised_Control
  

  
  #####################################################
  # Logit model and Probabilities for LoE at trial level
  # using the intercept value determined above and other model terms from the logit models fitted on the actual trial data
  d$Pr_LoE <- plogis((det.intercept_LoE + d$bi_0/10 + d$bi_1/10) + c1 * d$Treat)
  
  #describe(d$Pr_LoE[d$Treat==1])
  #describe(d$Pr_LoE[d$Treat==0])
  #describe(d$Pr_LoE)
  
  #View(d)
  
  #Use rbinom() directly with each patients characteristics
  
  #d$LoE_yes[d$Treat==1] <- ifelse(Pr_LoE[d$Treat==1]<0.532, 1, 0) * rep(rbinom(length(unique(d$id[d$Treat==1])), 1, 0.7), each = length(visits))
  #d$LoE_yes[d$Treat==0] <- ifelse(Pr_LoE[d$Treat==0]>0.5075, 1, 0) * rep(rbinom(length(unique(d$id[d$Treat==0])), 1, 0.7), each = length(visits))
  
  
  d$LoE_yes <- rep(rbinom(length(unique(d$id)), 1, unique(d$Pr_LoE)), each = length(visits)) 
  
  #describe(d$LoE_yes)
  #reshape to wide

  # check JM model to see how and which RE are in the survival model
  # fit JM model on the SM super trial see if it can be parameterised Weibull and check back here
  
  #unique(d$Pr_LoE)
  
  #describe(d$LoE_yes[d$Treat==0])
  #describe(d$LoE_yes[d$Treat==1])
  
  #describe(d$LoE_yes)
  
  #View(d)
  
  
  #####################################################
  # Logit model and Probabilities for AE in experimental arm
  # using the intercept value determined above and other model terms from the logit models fitted on the actual trial data
  d$Pr_AE[d$Treat==1] <- plogis((det.intercept_AE_exp + d$bi_0[d$Treat==1]/100 +  d$bi_1[d$Treat==1]/100))
  #describe(d$Pr_AE[d$Treat==1])
  
  
  d$AE_yes[d$Treat==1] <- rep(rbinom(length(unique(d$id[d$Treat==1])), 1, unique(d$Pr_AE[d$Treat==1])), each = length(visits)) 
  
  d$AE_yes[d$Treat==1]
  
  #describe(d$AE_yes[d$Treat==1])
  
  
  #####################################################
  # Logit model and Probabilities for AE in control arm
  # using the intercept value determined above and other model terms from the logit models fitted on the actual trial data
  d$Pr_AE[d$Treat==0] <- plogis((det.intercept_AE_control + d$bi_0[d$Treat==0]/100 + d$bi_1[d$Treat==0]/100))
  #describe(d$Pr_AE[d$Treat==0])
  
  
  d$AE_yes[d$Treat==0] <- rep(rbinom(length(unique(d$id[d$Treat==0])), 1, unique(d$Pr_AE[d$Treat==0])), each = length(visits)) 
  
  d$AE_yes[d$Treat==0]
  
  #describe(d$AE_yes[d$Treat==0])

  #describe(d$AE_yes)
  
  #View(d)
  
  class(d$AE_yes)
  #View(d)
  #d <-d[,-6]
  
  #d_w <- d |> spread(visit, MADRS10_collected) # reshape to wide in order to create the CfB variable
  #View(d_w)
  
  class(d$AE_yes)
  class(d$LoE_yes)
  d$AE_yes<- as.numeric(d$AE_yes)
  d$LoE_yes<- as.numeric(d$LoE_yes)


  #d_w$LoE_YES <- ifelse(d_w$AE_yes==1 & d_w$LoE_yes==1, rbinom(1, 1, 0.5), d_w$LoE_yes)
  
  #d_w$AE_YES <- ifelse((d_w$AE_yes==1 & d_w$LoE_yes==1)==T, 1-d_w$LoE_YES, d_w$AE_yes)
  
  
  #d_L <- d_w |> gather(visit, MADRS10_collected, 10:16) # reshape to long format
  #View(d_L)
  
  
  #d_L <- d_L[order(d_L$id, d_L$visit),]; #d # order by subject id and Visit
  
  # determine th adjustment factor
  # set sample size to n=10^5 and simulate 1 trial
  # then determine as below
  #d$compete <- ifelse(d$AE_yes==1 & d$LoE_yes==1, 1, 0)
  
  #sum(ifelse(d$AE_yes==1 & d$LoE_yes==1, 1, 0))/7
  #sum(d$AE_yes)/7
  #(sum(d$compete)/7)/(sum(d$LoE_yes)/7)*100# this is the percentage adjustment that needs to be added to the actual targeted proportion of LoE
  
  d$AE_YES <- d$AE_yes
  
  sum(d$AE_YES)/7
  
  d$LoE_YES <- ifelse(d$AE_YES==0 & d$LoE_yes==1, 1, 0)
  sum(d$LoE_YES)/7
  
  
  #class(d_L$visit)
  #class(d_L$Treat)
  #d_L$Treat <- factor(d_L$Treat)
  #d_L$LoE_yes <- factor(d_L$LoE_yes)
  #d_L$AE_yes <- factor(d_L$AE_yes)
  
  #d$AE_Yes <- ifelse(d[,6]==1, 1, 
  #                        ifelse(d_mis_L[,7]==1, 1, 0))
  
  #d_mis_L$LoE_YES <- ifelse(d_mis_L[,5]==1 & d_mis_L[,10]==0, 1, 0)
  #describe(d_mis_L$LoE_YES)
  
  #View(d_L)
  
  d$Behavior <- ifelse(d$AE_YES==1, "AE",
                       ifelse(d$LoE_YES==1, "LoE", "No IE"))
  
  table(d$Behavior)/7
  
  
  #### assign and save the generated datasets----
  # naming sequence is "SimTrial"_"Method"_"trial sample size"_"iteration number"
  
  assign(paste0("SimTrial_spm", "_", n), d)
  #View(SimTrial_sm_1_5)
  dataset_name.Rdata <- paste0("SimTrial_spm", "_", n, ".Rdata")
  dataset_name <- paste0("SimTrial_spm", "_", n)
  save(dataset_name, file = dataset_name.Rdata)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Intercurrent events descriptives needed for the verification step ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  LoE_Y <- d[ ,c("Treat", "LoE_YES")]
  head(LoE_Y)
  
  LoE_Y_total <- sum(LoE_Y$LoE_YES)/length(visits)
  LoE_Y_Exp <- sum(LoE_Y$LoE_YES[LoE_Y$Treat==1])/length(visits)
  LoE_Y_Control <- sum(LoE_Y$LoE_YES[LoE_Y$Treat==0])/length(visits)
  
  ## LoE ##
  tb_LoE_total <- LoE_Y_total
  tb_LoE_Exp <- LoE_Y_Exp
  tb_LoE_Control <- LoE_Y_Control
  
  n_LoE_total[m, ] <-  tb_LoE_total
  n_LoE_Exp[m, ] <- tb_LoE_Exp
  n_LoE_Control[m, ] <- tb_LoE_Control
  
  AE_Y <- d[ ,c("Treat", "AE_YES")]
  head(AE_Y)
  
  
  AE_Y_total <- sum(AE_Y$AE_YES)/length(visits)
  AE_Y_Exp <- sum(AE_Y$AE_YES[AE_Y$Treat==1])/length(visits)
  AE_Y_Control <- sum(AE_Y$AE_YES[AE_Y$Treat==0])/length(visits)
  
  ## AE ##
  
  tb_AE_total <- AE_Y_total
  tb_AE_Exp <- AE_Y_Exp
  tb_AE_Control <- AE_Y_Control
  
  
  n_AE_total[m, ] <-  tb_AE_total
  n_AE_Exp[m, ] <- tb_AE_Exp
  n_AE_Control[m, ] <- tb_AE_Control
  
  ## LoE ##
  
  LoE_total_Perc[m,] <- round(LoE_Y_total/n*100, digits=2)
  LoE_Exp_Perc[m,] <- round(LoE_Y_Exp/Randomised_Exp*100, digits=2)
  LoE_Control_Perc[m,] <- round(LoE_Y_Control/Randomised_Control*100, digits=2)
  
  #AE
  
  AE_total_Perc[m,] <-  round(AE_Y_total/n*100, digits=2)
  AE_Exp_Perc[m,] <- round(AE_Y_Exp/Randomised_Exp*100, digits=2)
  AE_Control_Perc[m,] <- round(AE_Y_Control/Randomised_Control*100, digits=2)
  
  # Total AE + LoE and percentage relative to the entire study population
  
  
  n_AE_and_LoE_T[m, ] <- LoE_Y_total + AE_Y_total
  AE_and_LoE_Perc[m, ] <- round((LoE_Y_total + AE_Y_total)/n*100, digits=2)

  
  d$Treat <- factor(d$Treat)
  d$visit <- factor(d$visit)
  d$MADRS10_collected <- as.numeric(d$MADRS10_collected) 
  
  d_SPM <- d
  
  # All patients with TRUE trajectory 
  # LoE
  p<- ggplot(data = d_SPM[d_SPM$Behavior=="LoE",], aes(x = visit, y = MADRS10_collected, group = id, color=LoE_YES)) 
  plot_LoE_SPM <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") + facet_wrap(~ Treat)+
    scale_y_continuous(limits = c(-10, 60)) + ggtitle("SPM-LoE pattern") +
    scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_LoE_SPM
  
  # AE
  p<- ggplot(data = d_SPM[d_SPM$Behavior=="AE",], aes(x = visit, y = MADRS10_collected, group = id, color=AE_YES)) 
  plot_AE_SPM <- p + geom_line(size=0.5, color='#F8766D') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
    scale_y_continuous(limits = c(-10, 60)) + ggtitle("SPM-AE pattern")  +
    scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_AE_SPM
  
  #just No IE
  # No IE
  p<- ggplot(data = d_SPM[d_SPM$Behavior=="No IE",], aes(x = visit, y = MADRS10_collected, group = id)) 
  plot_NoIE_SPM <- p + geom_line(size=0.5, color='#619CFF') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="blue") + facet_wrap(~ Treat)+
    scale_y_continuous(limits = c(-10, 60))+ ggtitle("SPM-No IE pattern") +
    scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")) ; plot_NoIE_SPM
  
  
  # All behaviors
  p<- ggplot(data = d_SPM, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
  plot_all_SPM <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="black") + facet_wrap(~ Treat)+
    scale_y_continuous(limits = c(-10, 60))+ ggtitle("SPM-All patterns")+
    scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_all_SPM
  
  #describe(d)
  
  #View(d)
  
  the_plot_SPM <- (plot_all_SPM / plot_LoE_SPM) | (plot_AE_SPM / plot_NoIE_SPM); the_plot_SPM
  

  #View(d_SPM)

  # assign   
  #assign(paste('all_betas_spm'), betas_spm)
  
  #assign(paste('all_delta_spm'), delta_spm)
  
  #setTxtProgressBar(pb3, s)
  setTxtProgressBar(pb1, m)
}
end_time <- Sys.time()

end_time-start_time


# parameters extracted for LMM fitted models on full outcome data
colMeans(betas_spm)
colMeans(delta_spm) ; treatmenteffect


tolerance_margin <- 0.1 
difference_Verification_spm <- abs(treatmenteffect - colMeans(delta_spm))

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_Verification_spm) < tolerance_margin), "Verification SPM *SUCCESSFUL*", "Verification SPM NOT successful :(") 


# compile to report

# Table for the paper ----

table_AE_SPM <- data.frame(
  # descriptives AE  
  n_AE_Control,
  n_AE_Exp); table_AE_SPM

mean(n_AE_Control)
mean(n_AE_Exp)




# descriptives LoE  
table_LoE_SPM <-data.frame(
  n_LoE_Control,
  n_LoE_Exp); table_LoE_SPM

mean(n_LoE_Control)
mean(n_LoE_Exp)


#describe(table_IE_SM)


table_AE_SPM |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "AE") |> 
  rename(N_C_arm=N.AE.Control) |> 
  rename(N_E_arm=N.AE.Exp)

table_LoE_SPM |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "LoE") |> 
  rename(N_C_arm=N.LoE.Control) |> 
  rename(N_E_arm=N.LoE.Exp)


tab_SPM <- tibble(bind_rows(table_AE_SPM |> 
                             as.data.frame() |> 
                             mutate("Intercurrent event" = "AE") |> 
                             rename(N_C_arm=N.AE.Control) |> 
                             rename(N_E_arm=N.AE.Exp), 
                           table_LoE_SPM |> 
                             as.data.frame() |> 
                             mutate("Intercurrent event" = "LoE") |> 
                             rename(N_C_arm=N.LoE.Control) |> 
                             rename(N_E_arm=N.LoE.Exp))); tab_SPM


tab2_SPM <- tab_SPM |> group_by(`Intercurrent event`) |>
  summarise("N" = round(mean(N_C_arm), digits=1), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = round(mean(N_E_arm), digits=1), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = round(mean(N_C_arm + N_E_arm), digits=1),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) |> 
  adorn_totals("row"); tab2_SPM




gt(tab2_SPM) |> 
  tab_header(title = md("Table 8c. Descriptive statistics intercurrent events"), subtitle = md("Shared parameter model DGM")) |>
  tab_source_note(md(paste0("Averaged over", " ", m.iterations,  " ",  "simulated trials.", " ", "Trial sample size = ", " ", n ))) |> 
  tab_spanner(
    label = md("**Control**"),
    columns = c("N", "%")) |> 
  cols_align(
    align = "center",
    columns =  c("N", "%")
  ) |> 
  tab_spanner(
    label = md("**Treatment**"),
    columns = c("N ", "% ")) |> 
  cols_align(
    align = "center",
    columns =  c("N ", "% ")
  ) |> 
  tab_spanner(
    label = md("**Total**"),
    columns = c(" N ", " % ")) |> 
  cols_align(
    align = "center",
    columns =  c(" N ", " % ")
  ) |> 
  data_color(
    columns = c("%", "% ", " % "),
    colors = scales::col_numeric(
      palette = c(
        "#add8e6"),
      domain = NULL)
  ) |> 
  cols_align(
    align = "center",
    columns =  "Intercurrent event"
  ) |> 
  tab_style(
    style = list(
      cell_fill(color = "white"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = "Intercurrent event"
    )
  )



## Session info---- 
sessionInfo()

# build 6 February 2022
#> sessionInfo()
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Monterey 12.1

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#Random number generation:
#  RNG:     Mersenne-Twister 
#Normal:  Inversion 
#Sample:  Rounding 

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] splines   stats     graphics  grDevices utils     datasets 
#[7] methods   base     

#other attached packages:
#  [1] JM_1.4-8        foreign_0.8-82  survival_3.2-13 patchwork_1.1.1
#[5] gt_0.3.1        janitor_2.1.0   lme4_1.1-28     Matrix_1.4-0   
#[9] nlme_3.1-155    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
#[13] purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6   
#[17] ggplot2_3.3.5   tidyverse_1.3.1 MASS_7.3-55     gmailr_1.0.1  
installed.packages()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# determine the number of trials needed to simulate for the verification of the longitudinal outcomes
# Formula from Burton paper
#tolerance_margin <- 0.1 # bias allowed
#std.e <- 1.04 # model-based standard error of the treatment effect estimate from a fitted model on 1 trial

#n.trials_needed <- ceiling(((qnorm(0.975) * std.e)/tolerance_margin)^2) ; n.trials_needed # for the verification 
# 416 trials
# verification of the longitudinal outcomes was successful


