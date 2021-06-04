## License CC-BY-4.0 
## Creative Commons Attribution 4.0 International
## Mitroiu M. et al SMMR


# Simulation study ----
## Scenario A "early separation and treatment effect maintained"----
## Outcomes generated as responses, and analysed with baseline as fixed-effect (as covariate), not as response
## function for DGM, IEGM and analyses methods

## code outline ##----

# admin: load libraries and define objects to store the data
# generate longitudinal outcomes
# ~ fit model and get parameters
# generate intercurrent events
# ~ get descriptives
# assign and save generated datasets
# visualisation of longitudinal outcomes and intercurrent events\

# with differing number of patients per pattern, the overall treatment effect will change obviously.
# Need to finetune the betas for each pattern for all the wanted treatment effects and factoring in the scaling factor.
# create simple function to sample few and visualise them



# to do
# make a comment here to inform the reader/user about the checks.
# rename for each d_orig
# add randomisation objects
# percentages of intercurrent events to be filled in from the SM tables


## load libraries ----
rm(list=ls()) #
# needed for the selection model method
library(gmailr)
library(MASS)
library(tidyverse)
library(nlme)
library(Hmisc)


library(MASS)
library(nlme)
library(survival)
library(foreign)
library(tidyverse)
#install.packages('tinytex')
library(tidyr)
library(haven)
library(Hmisc)
library(nlme)
library(lme4)
library(car)
library(ggplot2)
library(predictmeans)
library(tableone)
library(lattice)
library(mice)
library(data.table)
library(VIM)
library(naniar)
library(patchwork)
library(dplyr)
library(primes)
library(gt)
library(equatiomatic)
library(gmailr)
library(parallel)
#install.packages("clubSandwich",repos="http://cran.r-project.org")
library(gridExtra)
library(metafolio)
library(scales)
library(lavaSearch2)
library(lmerTest)
library(parameters)
library(sandwich)
library(lmtest)
#library(clubSandwich)


## Session info---- 
sessionInfo()
installed.packages()


?# Selection model via marginal model for outcomes-generating model and deterministic rules for generation of intercurrent events

## gmail setup----
# Selection model via marginal model for outcomes-generating model and deterministic rules for generation of intercurrent events
# setup to receive e-mails with results of the simulations. Useful to store results, but most importantly to be notified when the simulation is concluded.
# various tutorials can be found to set this up
google_app <- httr::oauth_app(
  "renamedapp",
  key = "126364165263-nudc2q7h24voutu33a9i6pik9rjou09i.apps.googleusercontent.com",
  secret = "Orgt-B5-eAEplGIbfWmr4Uhy"
)

gm_auth_configure(key = "126364165263-nudc2q7h24voutu33a9i6pik9rjou09i.apps.googleusercontent.com",
                  secret = "Orgt-B5-eAEplGIbfWmr4Uhy")

options(mc.cores = parallel::detectCores()) # M1 with parallel, 2017 macbook pro without parallel setting on
# considerably reduces the time needed to complete the simulations vs parallel setting off
# for other machines it is recommended to test few runs with different iterations and settings to optimise the computation of the simulation

Scenario <- c("A")
# This is scenario A, meaning the target longitudinal trajectories simulated are linear; we took example the Siddiqui paper referenced in the paper.

## simulation parameters #----

set.seed(2147483629) # set seed
#set.seed(2147483399)
m.iterations <- 50# number of generated datasets # number of trials per scaling factor
scaling_factor <-  c(1) # this is used for coding consistency between the four methods.
# In this simulation the scaling factor does not play any role.
# Could be used however to vary difference scenarios,e.g. a range of ratios for the AE:LoE at trial and arm level.

n <- 190# # number of patients to be simulated (sample size)
# this is based on a t-test to ensure  90% power at alpha level=0.025 one-sided

# ranges of probabilities centered around desired percentages of each intercurrent events averaged over all simulated trials
# this is done to increase variability in intercurrent events percentages between trials
prop_LoE <- 0.38 #c(0.36, 0.37, 0.38, 0.39, 0.40); mean(prop_LoE) # c(0.383) 
prop_AE_exp <- 0.04#c(0.020, 0.025, 0.03, 0.035, 0.04); mean(prop_AE_exp) # c(0.0298) 
prop_AE_control <-  0.02#c(0.010, 0.0125, 0.0150, 0.0175, 0.020) ; mean(prop_AE_control) # c(0.0158)
prop_AE <-prop_AE_exp + prop_AE_control

#CFE <- matrix(ncol=4,nrow=length(scaling_factor)*m.iterations)
#colnames(CFE) <-c("N ceiled_floored", "% ceiled_floored", "scaling factor", "simulated trial n")

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline + follow-up measurements
delta <- matrix(ncol=1,nrow=m.iterations) # treatment effect estimate at 6 weeks based on MMRM models fitted on each generated dataset
colnames(delta) <-c("TreatmentEffect")
betas <- matrix(ncol=1,nrow=m.iterations)
colnames(betas) <-c("visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3) # progress bar in percentages relative to the total number of m.iterations

pb3 <- txtProgressBar(min = 0,  max=length(scaling_factor), style=3)

start_time <- Sys.time() # timestamp for the start time of the nested for loop below.
# it was used to have an estimate of time needed for different larger number of trials to be simulated upon scaling up the simulation parameters (e.g., m.iterations)

## Begin for loop----
for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {

    ### percentages of intercurrent events to be filled in from the SM tables ----
    sampled_prop_LoE <- sample(prop_LoE, 1)
    n_LoE_all <- round(sampled_prop_LoE*n * scaling_factor[s], digits = 0) ; n_LoE_all
    
    sampled_prop_AE_exp <- sample(prop_AE_exp, 1)
    n_AE_exp <- round(sampled_prop_AE_exp*n * scaling_factor[s], digits = 0) ; n_AE_exp
    
    sampled_prop_AE_control <- sample(prop_AE_control, 1)
    n_AE_control <- round(sampled_prop_AE_control*n * scaling_factor[s], digits = 0) ; n_AE_control
  
    n_AE_all <- n_AE_exp + n_AE_control ; n_AE_all
    
    n_No_IE <- n - (n_LoE_all + n_AE_exp + n_AE_control); n_No_IE
    prop_No_IE <- 1 - (sampled_prop_LoE + sampled_prop_AE_exp + sampled_prop_AE_control) ; prop_No_IE
  
    
    ### PATTERN OF LoE IN BOTH EXPERIMENTAL AND CONTROL ARMS----
    #### Generate pattern for LoE, both experimental and control arm
    #### Generate correlated errors----
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    
    size_diag <- 0.001
    re_covm_proof <- matrix(
      c(size_diag, 0, 0, 0, 0, 0, 0,
        0, size_diag, 0, 0, 0, 0, 0,
        0, 0, size_diag, 0, 0, 0, 0,
        0, 0, 0, size_diag, 0, 0, 0,
        0, 0, 0, 0, size_diag, 0, 0,
        0, 0, 0, 0, 0, size_diag, 0,
        0, 0, 0, 0, 0, 0, size_diag), nrow = 7)
    
    re_covm_LoE_all <- matrix(
      c(10.029, 18.475, 16.469, 15.586, 12.817, 13.824, 11.079,
        18.475, 49.857, 40.987, 36.110, 26.753, 29.219, 23.940,
        16.469, 40.987, 48.772, 41.651, 35.042, 33.685, 26.541,
        15.586, 36.110, 41.651, 48.975, 37.920, 34.765, 25.684,
        12.817, 26.753, 35.042, 37.920, 44.017, 36.730, 27.583,
        13.824, 29.219, 33.685, 34.765, 36.730, 44.046, 32.626,
        11.079, 23.940, 26.541, 25.684, 27.583, 32.626, 31.158), nrow = 7)
    
    
    re_LoE_all <- mvrnorm(n_LoE_all, re_means, re_covm_LoE_all)	; re_LoE_all
    #View(re_LoE_all)
    
    re_LoE_all <- as.matrix(re_LoE_all)
    colnames(re_LoE_all) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re_LoE_all
    
    d_LoE_all <- data.frame(
      id = rep(1:n_LoE_all, each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n_LoE_all, 1, 0.5), each = length(visits)),
      MADRS10 = rep(NA, n_LoE_all)); d_LoE_all # mean(Treat)
    
    #head(re_LoE_all)
    d_LoE_all <- d_LoE_all[order(d_LoE_all$visit, d_LoE_all$id),]; #d_LoE_all
    #re_LoE_all
    #var(re_LoE_all)
    
    j_LoE_all<-c(re_LoE_all[, 1], re_LoE_all[, 2], re_LoE_all[, 3], re_LoE_all[, 4], re_LoE_all[, 5], re_LoE_all[, 6], re_LoE_all[, 7])
    d_LoE_all$re_LoE_all <- j_LoE_all ; #d_LoE_all
    
    d_LoE_all <- d_LoE_all[order(d_LoE_all$id, d_LoE_all$visit),]; #d_LoE_all
    
    d_LoE_all<-as.matrix(d_LoE_all)
    
    # MMRM parameters to generate the longitudinal outcomes
    # We used as inspiration the trial 003-002 (ref data analysis paper). Here we used a simplified scenario with linear group trajectories.
    
    
    #### model parameters LoE in both arms----
    beta.baseline_LoE_all <- 29.79
    beta_week1_LoE_all <- 1#0.93#-0.117613
    beta_week2_LoE_all <- 1.2#2.3#0.556190
    beta_week3_LoE_all <- 1.2#2.2#0.053362
    beta_week4_LoE_all <- 1.8#3#0.947715 
    beta_week5_LoE_all <- 2.8#4#1.562001
    beta_week6_LoE_all <- 0.95#2#1.450964
    
    beta_v1_treatment_LoE_all <- 0#-0.5#-0.271627
    beta_v2_treatment_LoE_all <- 0.2#-1#-0.024667
    beta_v3_treatment_LoE_all <- -0.5#-1.6#-0.712452
    beta_v4_treatment_LoE_all <- -1#-2.15#-0.816656
    beta_v5_treatment_LoE_all <- -1.6#-2.75#-0.737112
    beta_v6_treatment_LoE_all <- -2.15#0#-2#-1.909406  #-1.318160
    
    treatmenteffect_LoE_all <-  beta_v6_treatment_LoE_all ; treatmenteffect_LoE_all
    
    # following this model:
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6) + 
    #                          Beta_W1_Treat * T + Beta_W2_Treat * T + Beta_W3_Treat * T + Beta_W4_Treat * T + Beta_W5_Treat * T + Beta_W6_Treat * T 
    
    for (i in 1:(n_LoE_all*length(visits))) {
      d_LoE_all[i,4] <- ifelse(d_LoE_all[i, 2]==0, beta.baseline_LoE_all + d_LoE_all[i,5],
                               ifelse(d_LoE_all[i, 2]==7, d_LoE_all[i-1,4] + beta_week1_LoE_all + d_LoE_all[i,5] +  beta_v1_treatment_LoE_all * d_LoE_all[i, 3],
                                      ifelse(d_LoE_all[i, 2]==14, d_LoE_all[i-2,4]+ beta_week2_LoE_all + d_LoE_all[i,5] +  beta_v2_treatment_LoE_all * d_LoE_all[i, 3],
                                             ifelse(d_LoE_all[i, 2]==21, d_LoE_all[i-3,4] + beta_week3_LoE_all + d_LoE_all[i,5] +  beta_v3_treatment_LoE_all * d_LoE_all[i, 3],
                                                    ifelse(d_LoE_all[i, 2]==28, d_LoE_all[i-4,4] + beta_week4_LoE_all + d_LoE_all[i,5] +  beta_v4_treatment_LoE_all * d_LoE_all[i, 3],
                                                           ifelse(d_LoE_all[i, 2]==35, d_LoE_all[i-5,4] + beta_week5_LoE_all + d_LoE_all[i,5] +  beta_v5_treatment_LoE_all * d_LoE_all[i, 3],
                                                                  d_LoE_all[i-6,4] + beta_week6_LoE_all + d_LoE_all[i,5] +  beta_v6_treatment_LoE_all * d_LoE_all[i, 3]))))))
    }
    
    #View(d)
    d_LoE_all <-as.data.frame(d_LoE_all)
    d_LoE_all$visit <-as.factor(d_LoE_all$visit)
    d_LoE_all$Treat <- factor(d_LoE_all$Treat)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MMRM on full outcome data LoE
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # this is the raw dataset used to check the model fit
    d_LoE_all <-d_LoE_all[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d_LoE_all # full outcome data
    
    length(d_LoE_all$id)
    tmp <- sapply(unique(d_LoE_all$id), FUN = function(i) nrow(d_LoE_all[d_LoE_all$id == i,]))
    BaselineMADRS10 <-  rep(d_LoE_all$MADRS10[d_LoE_all$visit == 0], tmp)
    length(BaselineMADRS10)
    d_LoE_all$Baseline <- BaselineMADRS10
    #d_LoE_all<-d_LoE_all[d_LoE_all$visit!=0,]
    
    range(d_LoE_all$MADRS10[d_LoE_all$Treat==1 & d_LoE_all$visit!=0])
    range(d_LoE_all$MADRS10[d_LoE_all$Treat==0 & d_LoE_all$visit!=0])
    
    range(d_LoE_all$MADRS10[d_LoE_all$Treat==1 & d_LoE_all$visit==42])
    range(d_LoE_all$MADRS10[d_LoE_all$Treat==0 & d_LoE_all$visit==42])
    
    #mean(d_LoE_all$MADRS10[d_LoE_all$Treat==1 & d_LoE_all$visit==0])
    #mean(d_LoE_all$MADRS10[d_LoE_all$Treat==0 & d_LoE_all$visit==0])
    
    mean(d_LoE_all$MADRS10[d_LoE_all$Treat==1 & d_LoE_all$visit==42])
    mean(d_LoE_all$MADRS10[d_LoE_all$Treat==0 & d_LoE_all$visit==42])
    
    #### plot of longitudinal outcomes LoE----
    # plot the outcomes to see in big lines how the trajectories look like 
    p<- ggplot(data = d_LoE_all, aes(x = visit, y = MADRS10, group = id)) 
    plot_LoE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("PMMM-LoE pattern"); plot_LoE
    
    #View(d_LoE_all)
    
    d_LoE_all$V0 <- 0
    d_LoE_all$V0[d_LoE_all$visit == 0] <- 1
    
    d_LoE_all$V7 <- 0
    d_LoE_all$V7[d_LoE_all$visit == 7] <- 1
    
    d_LoE_all$V14 <- 0
    d_LoE_all$V14[d_LoE_all$visit == 14] <- 1
    
    d_LoE_all$V21 <- 0
    d_LoE_all$V21[d_LoE_all$visit == 21] <- 1
    
    d_LoE_all$V28 <- 0
    d_LoE_all$V28[d_LoE_all$visit == 28] <- 1
    
    d_LoE_all$V35 <- 0
    d_LoE_all$V35[d_LoE_all$visit == 35] <- 1
    
    d_LoE_all$V42 <- 0
    d_LoE_all$V42[d_LoE_all$visit == 42] <- 1
    #d_LoE_all
    
    # fit a model to check if the estimated parameters are similar/close to the true parameters
    #fit_LoE_all<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
     #                  Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
      #               data=d_LoE_all,
       #              correlation = corSymm(form=~1 | id),
        #             weights = varIdent(form = ~ 1 | visit), 
         #            method="REML")
    
    #summary(fit_LoE_all)
    
    #fit_LoE_all$coefficients[c(13)]
    
    #sum(fit_LoE_all$coefficients[c(13)]); treatmenteffect_LoE_all
    
    # indicator variable
    d_LoE_all$Pattern <- c("LoE_all")
    
    # indicator variable
    d_LoE_all$LoE_Yes <- 1
    #d_LoE_all$AE_exp <- 0
    #d_LoE_all$AE_control <- 0
    d_LoE_all$No_IE <- 0
    d_LoE_all$AE_Yes <- 0
    
    ### PATTERN OF LoE IN BOTH EXPERIMENTAL AND CONTROL ARMS----
    #### Generate pattern for LoE, both experimental and control arm
    #### Generate correlated errors----
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    
    size_diag <- 0.001
    re_covm_proof <- matrix(
      c(size_diag, 0, 0, 0, 0, 0, 0,
        0, size_diag, 0, 0, 0, 0, 0,
        0, 0, size_diag, 0, 0, 0, 0,
        0, 0, 0, size_diag, 0, 0, 0,
        0, 0, 0, 0, size_diag, 0, 0,
        0, 0, 0, 0, 0, size_diag, 0,
        0, 0, 0, 0, 0, 0, size_diag), nrow = 7)
    
    re_covm_AE_all <- matrix(
      c( 6.9032, 10.336,  8.0566,  7.7261,  5.4332,  6.7676,  6.2245,
         10.3360, 27.081, 16.1550, 16.0000, 11.5500, 16.2050, 15.5570,
         8.0566, 16.155, 16.4560, 15.4320, 14.2650, 15.2550, 13.9950,
         7.7261, 16.000, 15.4320, 27.8230, 22.2510, 21.2790, 17.8320,
         5.4332, 11.550, 14.2650, 22.2510, 34.3770, 28.9640, 25.4320,
         6.7676, 16.205, 15.2550, 21.2790, 28.9640, 39.7950, 34.7130,
         6.2245, 15.557, 13.9950, 17.8320, 25.4320, 34.7130, 37.7120), nrow = 7)
    
    re_AE_all <- mvrnorm(n_AE_all, re_means, re_covm_LoE_all)	; re_AE_all
    #View(re_AE_all)
    
    re_AE_all <- as.matrix(re_AE_all)
    colnames(re_AE_all) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re_AE_all
    
    prob_treat_AE <- prop_AE_exp/(prop_AE_exp+prop_AE_control)
    
    d_AE_all <- data.frame(
      id = rep((n_LoE_all+1):(n_LoE_all+n_AE_all), each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n_AE_all, 1, prob_treat_AE), each = length(visits)),
      MADRS10 = rep(NA, n_AE_all)); d_AE_all # mean(Treat)
    
    #describe(d_AE_all$Treat)
    
    #head(re_AE_all)
    d_AE_all <- d_AE_all[order(d_AE_all$visit, d_AE_all$id),]; #d_AE_all
    #re_AE_all
    #var(re_AE_all)
    
    j_AE_all<-c(re_AE_all[, 1], re_AE_all[, 2], re_AE_all[, 3], re_AE_all[, 4], re_AE_all[, 5], re_AE_all[, 6], re_AE_all[, 7])
    d_AE_all$re_AE_all <- j_AE_all ; #d_AE_all
    
    d_AE_all <- d_AE_all[order(d_AE_all$id, d_AE_all$visit),]; #d_AE_all
    
    d_AE_all<-as.matrix(d_AE_all)
    
    #### model parameters AE in both arms----
    beta.baseline_AE_all <- 29.79
    beta_week1_AE_all <- 0.1#0.93#-0.117613
    beta_week2_AE_all <- 0.1#2.3#0.556190
    beta_week3_AE_all <- 0.2#2.2#0.053362
    beta_week4_AE_all <- 0.5#3#0.947715 
    beta_week5_AE_all <- 1#4#1.562001
    beta_week6_AE_all <- 1#2#1.450964
    
    beta_v1_treatment_AE_all <- -7.7#-0.5#-0.271627
    beta_v2_treatment_AE_all <- -12#-1#-0.024667
    beta_v3_treatment_AE_all <- -11#-1.6#-0.712452
    beta_v4_treatment_AE_all <- -12#-2.15#-0.816656
    beta_v5_treatment_AE_all <- -11#-2.75#-0.737112
    beta_v6_treatment_AE_all <- -8#0#-2#-1.909406  #-1.318160
    
    treatmenteffect_AE_all <-  beta_v6_treatment_AE_all ; treatmenteffect_AE_all
    
    # following this model:
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6) + 
    #                          Beta_W1_Treat * T + Beta_W2_Treat * T + Beta_W3_Treat * T + Beta_W4_Treat * T + Beta_W5_Treat * T + Beta_W6_Treat * T 
    
    for (i in 1:(n_AE_all*length(visits))) {
      d_AE_all[i,4] <- ifelse(d_AE_all[i, 2]==0, beta.baseline_AE_all + d_AE_all[i,5],
                               ifelse(d_AE_all[i, 2]==7, d_AE_all[i-1,4] + beta_week1_AE_all + d_AE_all[i,5] +  beta_v1_treatment_AE_all * d_AE_all[i, 3],
                                      ifelse(d_AE_all[i, 2]==14, d_AE_all[i-2,4]+ beta_week2_AE_all + d_AE_all[i,5] +  beta_v2_treatment_AE_all * d_AE_all[i, 3],
                                             ifelse(d_AE_all[i, 2]==21, d_AE_all[i-3,4] + beta_week3_AE_all + d_AE_all[i,5] +  beta_v3_treatment_AE_all * d_AE_all[i, 3],
                                                    ifelse(d_AE_all[i, 2]==28, d_AE_all[i-4,4] + beta_week4_AE_all + d_AE_all[i,5] +  beta_v4_treatment_AE_all * d_AE_all[i, 3],
                                                           ifelse(d_AE_all[i, 2]==35, d_AE_all[i-5,4] + beta_week5_AE_all + d_AE_all[i,5] +  beta_v5_treatment_AE_all * d_AE_all[i, 3],
                                                                  d_AE_all[i-6,4] + beta_week6_AE_all + d_AE_all[i,5] +  beta_v6_treatment_AE_all * d_AE_all[i, 3]))))))
    }
  
    #View(d)
    d_AE_all <-as.data.frame(d_AE_all)
    d_AE_all$visit <-as.factor(d_AE_all$visit)
    d_AE_all$Treat <- factor(d_AE_all$Treat)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MMRM on full outcome data AE
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # this is the raw dataset used to check the model fit
    d_AE_all <-d_AE_all[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # change the numbers of the columns with the names of the variables
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d_AE_all # full outcome data
  
    length(d_AE_all$id)
    tmp <- sapply(unique(d_AE_all$id), FUN = function(i) nrow(d_AE_all[d_AE_all$id == i,]))
    BaselineMADRS10 <-  rep(d_AE_all$MADRS10[d_AE_all$visit == 0], tmp)
    length(BaselineMADRS10)
    d_AE_all$Baseline <- BaselineMADRS10
    #d_AE_all<-d_AE_all[d_AE_all$visit!=0,]
    
    range(d_AE_all$MADRS10[d_AE_all$Treat==1 & d_AE_all$visit!=0])
    range(d_AE_all$MADRS10[d_AE_all$Treat==0 & d_AE_all$visit!=0])
    
    range(d_AE_all$MADRS10[d_AE_all$Treat==1 & d_AE_all$visit==42])
    range(d_AE_all$MADRS10[d_AE_all$Treat==0 & d_AE_all$visit==42])
    
    #mean(d_AE_all$MADRS10[d_AE_all$Treat==1 & d_AE_all$visit==0])
    #mean(d_AE_all$MADRS10[d_AE_all$Treat==0 & d_AE_all$visit==0])
    
    mean(d_AE_all$MADRS10[d_AE_all$Treat==1 & d_AE_all$visit==42])
    mean(d_AE_all$MADRS10[d_AE_all$Treat==0 & d_AE_all$visit==42])
    
    #### plot of longitudinal outcomes AE----
    p<- ggplot(data = d_AE_all, aes(x = visit, y = MADRS10, group = id)) 
    plot_AE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("PMMM-AE pattern") ; plot_AE
    
    #tail(d_LoE_all)
    
    #View(d_AE_all)
    
    d_AE_all$V0 <- 0
    d_AE_all$V0[d_AE_all$visit == 0] <- 1
    
    d_AE_all$V7 <- 0
    d_AE_all$V7[d_AE_all$visit == 7] <- 1
    
    d_AE_all$V14 <- 0
    d_AE_all$V14[d_AE_all$visit == 14] <- 1
    
    d_AE_all$V21 <- 0
    d_AE_all$V21[d_AE_all$visit == 21] <- 1
    
    d_AE_all$V28 <- 0
    d_AE_all$V28[d_AE_all$visit == 28] <- 1
    
    d_AE_all$V35 <- 0
    d_AE_all$V35[d_AE_all$visit == 35] <- 1
    
    d_AE_all$V42 <- 0
    d_AE_all$V42[d_AE_all$visit == 42] <- 1
    #d_AE_all
    
    # fit a model to check if the estimated parameters are similar/close to the true parameters
    #fit_AE_all<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
    #                  Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
    #               data=d_AE_all,
    #              correlation = corSymm(form=~1 | id),
    #             weights = varIdent(form = ~ 1 | visit), 
    #            method="REML")
    
    #summary(fit_AE_all)
    
    #fit_AE_all$coefficients[c(13)]
    
    #sum(fit_AE_all$coefficients[c(13)]); treatmenteffect_AE_all
    
    # indicator variable
    d_AE_all$Pattern <- c("AE_all")
    
    # indicator variable
    d_AE_all$AE_Yes <- 1
    #d_AE_all$AE_exp <- 0
    #d_AE_all$AE_control <- 0
    d_AE_all$No_IE <- 0
    d_AE_all$LoE_Yes <- 0

    ### PATTERN OF No IE IN BOTH EXPERIMENTAL AND CONTROL ARMS----
    #### Generate pattern for No IE, both experimental and control arm
    #### Generate correlated errors----
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    
    re_covm_No_IE <- matrix(
      c(10.040, 18.135, 16.595, 15.814, 13.524, 14.740, 11.582,
        18.135, 48.213, 42.211, 37.687, 30.909, 34.285, 28.580,
        16.595, 42.211, 54.236, 46.775, 44.082, 44.487, 37.226,
        15.814, 37.687, 46.775, 54.825, 46.690, 44.923, 35.807,
        13.524, 30.909, 44.082, 46.690, 56.719, 51.932, 43.362,
        14.740, 34.285, 44.487, 44.923, 51.932, 62.843, 52.026,
        11.582, 28.580, 37.226, 35.807, 43.362, 52.026, 51.272), nrow = 7)
    
    
    re_No_IE <- mvrnorm(n_No_IE, re_means, re_covm_LoE_all)	; re_No_IE
    #View(re)
    
    re_No_IE <- as.matrix(re_No_IE)
    colnames(re_No_IE) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re_No_IE
    
    d_No_IE <- data.frame(
      id = rep((n_LoE_all+n_AE_all+1):(n_LoE_all+n_AE_all+n_No_IE), each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n_No_IE, 1, 0.5), each = length(visits)),
      MADRS10 = rep(NA, n_No_IE)); d_No_IE # mean(Treat)
    
    tail(d_AE_all)
    head(d_No_IE)
    
    d_No_IE <- d_No_IE[order(d_No_IE$visit, d_No_IE$id),]; #d_No_IE
    #re
    
    j_No_IE<-c(re_No_IE[, 1], re_No_IE[, 2], re_No_IE[, 3], re_No_IE[, 4], re_No_IE[, 5], re_No_IE[, 6], re_No_IE[, 7])
    d_No_IE$re_No_IE <- j_No_IE ; #d_No_IE
    
    d_No_IE <- d_No_IE[order(d_No_IE$id, d_No_IE$visit),]; #d_No_IE
    
    d_No_IE<-as.matrix(d_No_IE)
    
    #### model parameters No IE in both arms----
    beta.baseline_No_IE <- 29.79
    beta_week1_No_IE <- -2.5#-2#-2.061548
    beta_week2_No_IE <- -3.5#-3.864898
    beta_week3_No_IE <- -4.4#-4.360199
    beta_week4_No_IE <- -5.75# -5.975185
    beta_week5_No_IE <- -7.4#-7.406465
    beta_week6_No_IE <- -7#-8.161270
    
    beta_v1_treatment_No_IE <- -1.3#-1.5# -0.469287
    beta_v2_treatment_No_IE <- -2#-0.268673
    beta_v3_treatment_No_IE <- -2.5#-0.724476
    beta_v4_treatment_No_IE <- -3#-1.037980
    beta_v5_treatment_No_IE <- -3.5#-1.420660
    beta_v6_treatment_No_IE <- -4.25#-4.25#-4.195854 # -1.904595
  
    treatmenteffect_No_IE <-  beta_v6_treatment_No_IE ; treatmenteffect_No_IE
    
    # following this model:
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6) + 
    #                          Beta_W1_Treat * T + Beta_W2_Treat * T + Beta_W3_Treat * T + Beta_W4_Treat * T + Beta_W5_Treat * T + Beta_W6_Treat * T 
    
    for (i in 1:(n_No_IE*length(visits))) {
      d_No_IE[i,4] <- ifelse(d_No_IE[i, 2]==0, beta.baseline_No_IE + d_No_IE[i,5],
                             ifelse(d_No_IE[i, 2]==7, d_No_IE[i-1,4] + beta_week1_No_IE + d_No_IE[i,5] +  beta_v1_treatment_No_IE * d_No_IE[i, 3],
                                    ifelse(d_No_IE[i, 2]==14, d_No_IE[i-2,4]+ beta_week2_No_IE + d_No_IE[i,5] +  beta_v2_treatment_No_IE * d_No_IE[i, 3],
                                           ifelse(d_No_IE[i, 2]==21, d_No_IE[i-3,4] + beta_week3_No_IE + d_No_IE[i,5] +  beta_v3_treatment_No_IE * d_No_IE[i, 3],
                                                  ifelse(d_No_IE[i, 2]==28, d_No_IE[i-4,4] + beta_week4_No_IE + d_No_IE[i,5] +  beta_v4_treatment_No_IE * d_No_IE[i, 3],
                                                         ifelse(d_No_IE[i, 2]==35, d_No_IE[i-5,4] + beta_week5_No_IE + d_No_IE[i,5] +  beta_v5_treatment_No_IE * d_No_IE[i, 3],
                                                                d_No_IE[i-6,4] + beta_week6_No_IE + d_No_IE[i,5] +  beta_v6_treatment_No_IE * d_No_IE[i, 3]))))))
    }
    
    #View(d)
    d_No_IE <-as.data.frame(d_No_IE)
    d_No_IE$visit <-as.factor(d_No_IE$visit)
    d_No_IE$Treat <- factor(d_No_IE$Treat)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MMRM on full outcome data No IE
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    # this is the raw dataset used to check the model fit
    d_No_IE <-d_No_IE[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d_No_IE # full outcome data
    
    length(d_No_IE$id)
    tmp <- sapply(unique(d_No_IE$id), FUN = function(i) nrow(d_No_IE[d_No_IE$id == i,]))
    BaselineMADRS10 <-  rep(d_No_IE$MADRS10[d_No_IE$visit == 0], tmp)
    length(BaselineMADRS10)
    d_No_IE$Baseline <- BaselineMADRS10
    #d_No_IE
    #d_No_IE<-d_No_IE[d_No_IE$visit!=0,]
    
    range(d_No_IE$MADRS10[d_No_IE$Treat==1 & d_No_IE$visit!=0])
    range(d_No_IE$MADRS10[d_No_IE$Treat==0 & d_No_IE$visit!=0])
    
    range(d_No_IE$MADRS10[d_No_IE$Treat==1 & d_No_IE$visit==42])
    range(d_No_IE$MADRS10[d_No_IE$Treat==0 & d_No_IE$visit==42])
    
    mean(d_No_IE$MADRS10[d_No_IE$Treat==1 & d_No_IE$visit==0])
    mean(d_No_IE$MADRS10[d_No_IE$Treat==0 & d_No_IE$visit==0])
    
    mean(d_No_IE$MADRS10[d_No_IE$Treat==1 & d_No_IE$visit==42])
    mean(d_No_IE$MADRS10[d_No_IE$Treat==0 & d_No_IE$visit==42])
    
    #### plot of longitudinal outcomes No IE----
    # plot the outcomes to see in big lines how the trajectories look like  
    p<- ggplot(data = d_No_IE, aes(x = visit, y = MADRS10, group = id)) 
    plot_NoIE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("PMMM-No IE pattern"); plot_NoIE
    
    #View(d_No_IE)
    
    d_No_IE$V0 <- 0
    d_No_IE$V0[d_No_IE$visit == 0] <- 1
    
    d_No_IE$V7 <- 0
    d_No_IE$V7[d_No_IE$visit == 7] <- 1
    
    d_No_IE$V14 <- 0
    d_No_IE$V14[d_No_IE$visit == 14] <- 1
    
    d_No_IE$V21 <- 0
    d_No_IE$V21[d_No_IE$visit == 21] <- 1
    
    d_No_IE$V28 <- 0
    d_No_IE$V28[d_No_IE$visit == 28] <- 1
    
    d_No_IE$V35 <- 0
    d_No_IE$V35[d_No_IE$visit == 35] <- 1
    
    d_No_IE$V42 <- 0
    d_No_IE$V42[d_No_IE$visit == 42] <- 1
    
    # fit a model to check if the estimated parameters are similar/close to the true parameters
    #fit_No_IE<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
     #                Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
      #             data=d_No_IE,
       #            correlation = corSymm(form=~1 | id),
        #           weights = varIdent(form = ~ 1 | visit), 
         #          method="REML")
    
    
    #summary(fit_No_IE)
    
    #getVarCov(fit_No_IE)
    #vcov(fit_No_IE)
    #re_covm_proof
    
    #fit_No_IE$coefficients[c(13)]
    
    #sum(fit_No_IE$coefficients[c(13)]); treatmenteffect_No_IE
    
    # indicator variable
    d_No_IE$Pattern <- c("No_IE")
    
    # indicator variable
    d_No_IE$LoE_Yes <- 0
    #d_No_IE$AE_exp <- 0
    #d_No_IE$AE_control <- 0
    d_No_IE$No_IE <- 1
    d_No_IE$AE_Yes <- 0
    
    #start with a model with only LoE_all and No_ie and fit the model there.
    
    ### Combine all data from all patterns in one big trial dataset----
    d_pmmm <- rbind(d_LoE_all,
                    d_AE_all,
                    #d_AE_exp,
                    #d_AE_control,
                    d_No_IE)
    
    head(d_LoE_all)
    head(d_AE_all)
    #head(d_AE_exp)
    #head(d_AE_control)
    head(d_No_IE)
    #View(d_pmmm)
    
    describe(d_pmmm)
    head(d_pmmm)
    d_pmmm$Pattern <- as.factor(d_pmmm$Pattern)
    
    d_pmmm$LoE_Yes <- as.factor(d_pmmm$LoE_Yes)
    d_pmmm$AE_Yes <- as.factor(d_pmmm$AE_Yes)
    #d_pmmm$AE_exp <- as.factor(d_pmmm$AE_exp)
    #d_pmmm$AE_control <- as.factor(d_pmmm$AE_control)
    d_pmmm$No_IE <- as.factor(d_pmmm$No_IE)
    levels(d_pmmm$visit)
    
    #### Plot trajectories for all the created patterns----
    # All patients with TRUE trajectory 
    
    ##### LoE_all/any----
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=LoE_Yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    ##### AE_exp----
    #p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=AE_exp)) 
    #p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
     # scale_y_continuous(limits = c(-10, 60))
    
    ##### AE_control----
    #p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=AE_control)) 
    #p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
     # scale_y_continuous(limits = c(-10, 60))
    
    ##### AE_any----
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=AE_Yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    ###### No_IE----
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=No_IE)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    #describe(d_pmmm)
    #view(d_pmmm)
    
    ##### All behaviors----
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=Pattern)) 
    #p + geom_line() + facet_grid(~ Treat) 
    plot_all <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("PMMM - All patterns"); plot_all

    ##### plots for the paper----
    (plot_all / plot_LoE) | (plot_AE / plot_NoIE)
    
    ##### MMRM model on the combined-patterns trial dataset
    fit_pmmm<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                    Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                  data=d_pmmm,
                  correlation = corSymm(form=~1 | id),
                  weights = varIdent(form = ~ 1 | visit), 
                  method="REML")
    
    summary(fit_pmmm)
    
    describe(d_pmmm)
    
    #describe(d_pmmm)
    getVarCov(fit_pmmm, individual = 1)
    
    re_covm_LoE_all
    
    treatmenteffect_pmmm <- sampled_prop_LoE * beta_v6_treatment_LoE_all +
      
    (sampled_prop_AE_exp*(beta_v6_treatment_AE_all + beta_week6_AE_all) - (sampled_prop_AE_control * beta_week6_AE_all)) +
    
     # sampled_prop_AE_exp *  beta_week6_AE_exp +
      
      #sampled_prop_AE_control *  beta_week6_AE_control +
      
    (1-(sampled_prop_LoE + sampled_prop_AE_exp + sampled_prop_AE_control)) *  beta_v6_treatment_No_IE ; treatmenteffect_pmmm
    
    sum(fit_pmmm$coefficients[c(13)]); treatmenteffect_pmmm
    # store parameters from model fit on each dataset
    betas[m, ] <- fit_pmmm$coefficients[c(13)]
    
    delta[m, ] <- sum(fit_pmmm$coefficients[c(13)])
    #mean(delta)
    #bias_f[m, ] <- sum(fit_pmmm$coefficients[c(7,13)]) - treatmenteffect
    
    #delta_error <- sqrt(vcov(fit_pmmm)["Treat1", "Treat1"] + vcov(fit_pmmm)["visit42:Treat1", "visit42:Treat1"] + 2*vcov(fit_pmmm)["Treat1", "visit42:Treat1"]) 
    
    #delta_errorz[m, ] <- delta_error 
    
    #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_error
    #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_error
    
    #### assign and save the generated datasets----
    # naming sequence is "SimTrial"_"Method"_"trial sample size"_"iteration number"_"scaling factor"
    
    assign(paste0("SimTrial_pmmm", "_", n, "_", m, "_", s), d_pmmm)
    #View(SimTrial_pmmm_1_1)
    
    dataset_name.Rdata <- paste0("SimTrial_pmmm", "_", n,"_", m, "_", s, ".Rdata")
    dataset_name <- paste0("SimTrial_pmmm", "_", n, "_", m, "_", s)
    
    save(dataset_name, file = dataset_name.Rdata)
    
    setTxtProgressBar(pb1, m)
  }
  
  
  # parameters extracted for MMRM fitted models on full outcome data
  colMeans(betas)
  colMeans(delta) ; treatmenteffect_pmmm

  # assign   
  assign(paste('all_betas', s, sep="_"), betas)
  
  assign(paste('all_delta', s, sep="_"), delta)
  
  setTxtProgressBar(pb3, s)
}
end_time <- Sys.time()

end_time-start_time

all_betas_1; 
colMeans(all_delta_1); treatmenteffect_pmmm

tolerance_margin <- 0.1 
difference_Verification <- abs(treatmenteffect_pmmm - colMeans(all_delta_1))

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_check) < tolerance_margin), "Verification successful", "Verification NOT successful") 

#hist(treatmenteffect_pmmm - all_betas_1)

hist(all_betas_1)

colMeans(all_betas_1) + 1.96*sd(all_betas_1)/sqrt(n)
colMeans(all_betas_1) - 1.96*sd(all_betas_1)/sqrt(n)

min(all_betas_1)
max(all_betas_1)

describe(SimTrial_pmmm_190_1_1$LoE_Yes[SimTrial_pmmm_190_1_1$Treat==1])
describe(SimTrial_pmmm_190_1_1$LoE_Yes[SimTrial_pmmm_190_1_1$Treat==0])

#### Visualisation of trajectories and patterns

# Plot trajectories
# All patients with TRUE trajectory 

# LoE_all/any
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=LoE_Yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


# AE_exp
#p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=AE_exp)) 
#p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


# AE_control
#p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=AE_control)) 
#p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


# AE_any
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=AE_Yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


# No_IE
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=No_IE)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)

#describe(d_pmmm)
#view(d_pmmm)

# All behaviors
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=Pattern)) 
#p + geom_line() + facet_grid(~ Treat) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)

fit_190<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
               Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
             data=SimTrial_pmmm_190_15_1,
             correlation = corSymm(form=~1 | id),
             weights = varIdent(form = ~ 1 | visit), 
             method="REML")

summary(fit_190)

colMeans(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5))

colMeans(rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5))

cbind(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5),
      rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5),
      rbind(all_delta_errorz_1, all_delta_errorz_2, all_delta_errorz_3, all_delta_errorz_4, all_delta_errorz_5),
      rbind(all_confint_fit_1, all_confint_fit_2, all_confint_fit_3, all_confint_fit_4, all_confint_fit_5),
      rbind(all_N_Exp_1, all_N_Exp_2, all_N_Exp_3, all_N_Exp_4, all_N_Exp_5),
      rbind(all_N_Control_1, all_N_Control_2, all_N_Control_3, all_N_Control_4, all_N_Control_5))
