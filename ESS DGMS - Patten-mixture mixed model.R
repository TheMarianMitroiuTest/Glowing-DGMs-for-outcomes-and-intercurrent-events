## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####
# models for patterns derived from the SM model on 5000 patients, with small covariance matrix


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"



### DGM
rm(list=ls())
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



sessionInfo()
installed.packages()


#View(SimTrial_1_1)


# Selection model via marginal model for outcomes-generating model and deterministic rules for generation of intercurrent events


google_app <- httr::oauth_app(
  "renamedapp",
  key = "126364165263-nudc2q7h24voutu33a9i6pik9rjou09i.apps.googleusercontent.com",
  secret = "Orgt-B5-eAEplGIbfWmr4Uhy"
)

gm_auth_configure(key = "126364165263-nudc2q7h24voutu33a9i6pik9rjou09i.apps.googleusercontent.com",
                  secret = "Orgt-B5-eAEplGIbfWmr4Uhy")



options(mc.cores = parallel::detectCores()) # M1 with parallel, Old1 without parallel

Scenario <- c("A")

# simulation


set.seed(2147483629)
#set.seed(2147483399)
m.iterations <- 1# number of generated datasets # number of trials per scaling factor
scaling_factor <-  c(1) #c(0.20, 0.40, 0.60, 0.80, 1) # to cover a range of IE % from ~5%-45% in total
# total number of simulated trials = m.iterations * length(scaling_factor)
# try with c(0.4, 1.1, 1.8, 2.3, 3)
#c(0.25, 0.5, 1, 2, 2.5) steps for scaling factor
# to multiply for LoE, AE_control and AE_exp ### value 0.5 should yield around 13.5% IEs, 1 ~ %, 2 ~ %, 2.5 ~



n <- 800# number of patients


# percentages/proportions of Intercurrent Events
prop_LoE <- c(0.36, 0.37, 0.38, 0.39, 0.40); mean(prop_LoE) # c(0.383) 
prop_AE_exp <- c(0.020, 0.025, 0.03, 0.035, 0.04); mean(prop_AE_exp) # c(0.0298) 
prop_AE_control <-  c(0.010, 0.0125, 0.0150, 0.0175, 0.020) ; mean(prop_AE_control) # c(0.0158)


#CFE <- matrix(ncol=4,nrow=length(scaling_factor)*m.iterations)
#colnames(CFE) <-c("N ceiled_floored", "% ceiled_floored", "scaling factor", "simulated trial n")

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline + follow-up measurements
delta <- matrix(ncol=1,nrow=m.iterations) # treatment effect estimate at 6 weeks based on MMRM models fitted on each generated dataset
colnames(delta) <-c("TreatmentEffect")
betas <- matrix(ncol=1,nrow=m.iterations)
colnames(betas) <-c("visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3)

pb3 <- txtProgressBar(min = 0,  max=length(scaling_factor), style=3)


s <-1

start_time <- Sys.time()

for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    
    # percentages of intercurrent events To be filled in from the SM tables
    sampled_prop_LoE <- sample(prop_LoE, 1)
    n_LoE_all <- round(sampled_prop_LoE*n * scaling_factor[s], digits = 0) ; n_LoE_all
    
    sampled_prop_AE_exp <- sample(prop_AE_exp, 1)
    n_AE_exp <- round(sampled_prop_AE_exp*n * scaling_factor[s], digits = 0) ; n_AE_exp
    
    sampled_prop_AE_control <- sample(prop_AE_control, 1)
    n_AE_control <- round(sampled_prop_AE_control*n * scaling_factor[s], digits = 0) ; n_AE_control
    
    n_No_IE <- n - (n_LoE_all + n_AE_exp + n_AE_control); n_No_IE
    prop_No_IE <- 1 - (sampled_prop_LoE + sampled_prop_AE_exp + sampled_prop_AE_control) ; prop_No_IE
    
    
    
    
    
    ##### PATTERN OF LoE IN BOTH EXPERIMENTAL AND CONTROL ARMS
    ### Generate pattern for LoE, both experimental and control arm
    ### Generate correlated errors
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
    
    # Scenario A -> pattern of LoE in both experimental and control arms
    beta.baseline_LoE_all <- 29.951049
    beta_week1_LoE_all <- -0.117613
    beta_week2_LoE_all <- 0.556190
    beta_week3_LoE_all <- 0.053362
    beta_week4_LoE_all <- 0.947715 
    beta_week5_LoE_all <- 1.562001
    beta_week6_LoE_all <- 1.450964
    
    beta_v1_treatment_LoE_all <- -0.271627
    beta_v2_treatment_LoE_all <- -0.024667
    beta_v3_treatment_LoE_all <- -0.712452
    beta_v4_treatment_LoE_all <- -0.816656
    beta_v5_treatment_LoE_all <- -0.737112
    beta_v6_treatment_LoE_all <- 0#-1.909406  #-1.318160
    
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
    
    
    ####################################################################################################  
    # MMRM on full outcome data in pattern LoE in both experimental and control arm
    ############################
    
    
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
    
    
    
    p<- ggplot(data = d_LoE_all, aes(x = visit, y = MADRS10, group = id)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-20, 75))
    
    
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
    
    
    
    fit_LoE_all<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                       Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                     data=d_LoE_all,
                     correlation = corSymm(form=~1 | id),
                     weights = varIdent(form = ~ 1 | visit), 
                     method="REML")
    
    
    
    
    summary(fit_LoE_all)
    
    
    fit_LoE_all$coefficients[c(13)]
    
    sum(fit_LoE_all$coefficients[c(13)]); treatmenteffect_LoE_all
    
    # indicator variable
    d_LoE_all$Pattern <- c("LoE_all")
    
    # indicator variable
    d_LoE_all$LoE_Yes <- 1
    d_LoE_all$AE_exp <- 0
    d_LoE_all$AE_control <- 0
    d_LoE_all$No_IE <- 0
    d_LoE_all$AE_Yes <- d_LoE_all$AE_exp + d_LoE_all$AE_control
    
    
    #
    
    
    
    ##### PATTERN OF AEs IN EXPERIMENTAL ARM
    
    
    
    ### Generate pattern for AEs in experimental arm
    ### Generate correlated errors
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    
    re_covm_AE_exp <- matrix(
      c( 6.9032, 10.336,  8.0566,  7.7261,  5.4332,  6.7676,  6.2245,
         10.3360, 27.081, 16.1550, 16.0000, 11.5500, 16.2050, 15.5570,
         8.0566, 16.155, 16.4560, 15.4320, 14.2650, 15.2550, 13.9950,
         7.7261, 16.000, 15.4320, 27.8230, 22.2510, 21.2790, 17.8320,
         5.4332, 11.550, 14.2650, 22.2510, 34.3770, 28.9640, 25.4320,
         6.7676, 16.205, 15.2550, 21.2790, 28.9640, 39.7950, 34.7130,
         6.2245, 15.557, 13.9950, 17.8320, 25.4320, 34.7130, 37.7120), nrow = 7)
    
    
    re_AE_exp <- mvrnorm(n_AE_exp, re_means, re_covm_AE_exp)	; re_AE_exp
    #View(re)
    
    re_AE_exp <- as.matrix(re_AE_exp)
    colnames(re_AE_exp) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re_AE_exp
    
    d_AE_exp <- data.frame(
      id = rep((n_LoE_all+1):(n_LoE_all+n_AE_exp), each = length(visits)),
      visit = visits,
      Treat = rep(1, n_AE_exp*length(visits)),
      MADRS10 = rep(NA, n_AE_exp)); d_AE_exp # mean(Treat)
    
    
    d_AE_exp <- d_AE_exp[order(d_AE_exp$visit, d_AE_exp$id),]; #d_AE_exp
    #re
    
    j_AE_exp<-c(re_AE_exp[, 1], re_AE_exp[, 2], re_AE_exp[, 3], re_AE_exp[, 4], re_AE_exp[, 5], re_AE_exp[, 6], re_AE_exp[, 7])
    d_AE_exp$re_AE_exp <- j_AE_exp ; #d_AE_exp
    
    
    d_AE_exp <- d_AE_exp[order(d_AE_exp$id, d_AE_exp$visit),]; #d_AE_exp
    
    
    d_AE_exp<-as.matrix(d_AE_exp)
    
    # Scenario A -> pattern of AE in  experimental arm
    beta.baseline_AE_exp <- 28.214968
    beta_week1_AE_exp <- -6.477987
    beta_week2_AE_exp <- -11.141631
    beta_week3_AE_exp <- -10.058028
    beta_week4_AE_exp <- -11.484430 
    beta_week5_AE_exp <- -12.509911
    beta_week6_AE_exp <- 0#-12.660152
    
    
    
    # following this model:
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6)
    
    
    
    for (i in 1:(n_AE_exp*length(visits))) {
      d_AE_exp[i,4] <- ifelse(d_AE_exp[i, 2]==0, beta.baseline_AE_exp + d_AE_exp[i,5],
                              ifelse(d_AE_exp[i, 2]==7, d_AE_exp[i-1,4] + beta_week1_AE_exp + d_AE_exp[i,5],
                                     ifelse(d_AE_exp[i, 2]==14, d_AE_exp[i-2,4]+ beta_week2_AE_exp + d_AE_exp[i,5],
                                            ifelse(d_AE_exp[i, 2]==21, d_AE_exp[i-3,4] + beta_week3_AE_exp + d_AE_exp[i,5],
                                                   ifelse(d_AE_exp[i, 2]==28, d_AE_exp[i-4,4] + beta_week4_AE_exp + d_AE_exp[i,5],
                                                          ifelse(d_AE_exp[i, 2]==35, d_AE_exp[i-5,4] + beta_week5_AE_exp + d_AE_exp[i,5],
                                                                 d_AE_exp[i-6,4] + beta_week6_AE_exp + d_AE_exp[i,5]))))))
    }
    
    
    
    #View(d)
    d_AE_exp <-as.data.frame(d_AE_exp)
    d_AE_exp$visit <-as.factor(d_AE_exp$visit)
    d_AE_exp$Treat <- factor(d_AE_exp$Treat)
    
    
    ####################################################################################################  
    # MMRM on full outcome data in pattern AE in experimental arm
    ############################
    
    
    # this is the raw dataset used to check the model fit
    d_AE_exp <-d_AE_exp[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d_AE_exp # full outcome data
    
    
    length(d_AE_exp$id)
    tmp <- sapply(unique(d_AE_exp$id), FUN = function(i) nrow(d_AE_exp[d_AE_exp$id == i,]))
    BaselineMADRS10 <-  rep(d_AE_exp$MADRS10[d_AE_exp$visit == 0], tmp)
    length(BaselineMADRS10)
    d_AE_exp$Baseline <- BaselineMADRS10
    d_AE_exp
    #d_AE_exp<-d_AE_exp[d_AE_exp$visit!=0,]
    
    
    range(d_AE_exp$MADRS10[d_AE_exp$Treat==1 & d_AE_exp$visit!=0])
    range(d_AE_exp$MADRS10[d_AE_exp$Treat==1 & d_AE_exp$visit==42])
    mean(d_AE_exp$MADRS10[d_AE_exp$Treat==1 & d_AE_exp$visit==0])
    mean(d_AE_exp$MADRS10[d_AE_exp$Treat==1 & d_AE_exp$visit==42])
    
    
    
    p<- ggplot(data = d_AE_exp, aes(x = visit, y = MADRS10, group = id)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-20, 75))
    
    
    
    #View(d_AE_exp)
    
    d_AE_exp$V0 <- 0
    d_AE_exp$V0[d_AE_exp$visit == 0] <- 1
    
    d_AE_exp$V7 <- 0
    d_AE_exp$V7[d_AE_exp$visit == 7] <- 1
    
    d_AE_exp$V14 <- 0
    d_AE_exp$V14[d_AE_exp$visit == 14] <- 1
    
    d_AE_exp$V21 <- 0
    d_AE_exp$V21[d_AE_exp$visit == 21] <- 1
    
    d_AE_exp$V28 <- 0
    d_AE_exp$V28[d_AE_exp$visit == 28] <- 1
    
    d_AE_exp$V35 <- 0
    d_AE_exp$V35[d_AE_exp$visit == 35] <- 1
    
    d_AE_exp$V42 <- 0
    d_AE_exp$V42[d_AE_exp$visit == 42] <- 1
    #d_AE_exp
    
    
    
    fit_AE_exp<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42, 
                    data=d_AE_exp,
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | visit), 
                    method="REML")
    
    
    
    
    summary(fit_AE_exp)
    
    # indicator variable
    d_AE_exp$Pattern <- c("AE_exp")
    
    # indicator variable
    
    d_AE_exp$LoE_Yes <- 0
    d_AE_exp$AE_exp <- 1
    d_AE_exp$AE_control <- 0
    d_AE_exp$No_IE <- 0
    d_AE_exp$AE_Yes <- d_AE_exp$AE_exp + d_AE_exp$AE_control
    ###
    
    
    ##### PATTERN OF LoE IN CONTROL ARM
    ### Generate pattern for AE in control arm
    ### Generate correlated errors
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    
    re_covm_AE_control <- matrix(
      c(13.077, 23.179, 17.310, 17.218, 15.430, 16.245, 12.939,
        23.179, 53.052, 33.383, 32.247, 25.094, 30.067, 24.085,
        17.310, 33.383, 30.641, 28.801, 25.875, 28.100, 22.957,
        17.218, 32.247, 28.801, 41.982, 34.827, 35.100, 27.577,
        15.430, 25.094, 25.875, 34.827, 47.057, 47.377, 38.056,
        16.245, 30.067, 28.100, 35.100, 47.377, 59.153, 46.553,
        12.939, 24.085, 22.957, 27.577, 38.056, 46.553, 44.032), nrow = 7)
    
    
    re_AE_control <- mvrnorm(n_AE_control, re_means, re_covm_AE_control)	; re_AE_control
    #View(re_AE_control)
    
    re_AE_control <- matrix(re_AE_control, ncol=7)
    colnames(re_AE_control) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re_AE_control
    
    d_AE_control <- data.frame(
      id = rep((n_LoE_all+n_AE_exp+1):(n_LoE_all+n_AE_exp+n_AE_control), each = length(visits)),
      visit = visits,
      Treat = rep(0, n_AE_control*length(visits)),
      MADRS10 = rep(NA, n_AE_control)); d_AE_control # mean(Treat)
    
    
    
    
    
    d_AE_control <- d_AE_control[order(d_AE_control$visit, d_AE_control$id),]; #d_AE_control
    #re
    
    j_AE_control<-c(re_AE_control[, 1], re_AE_control[, 2], re_AE_control[, 3], re_AE_control[, 4], re_AE_control[, 5], re_AE_control[, 6], re_AE_control[, 7])
    d_AE_control$re_AE_control <- j_AE_control ; #d_AE_control
    
    
    d_AE_control <- d_AE_control[order(d_AE_control$id, d_AE_control$visit),]; #d_AE_control
    
    
    d_AE_control<-as.matrix(d_AE_control)
    
    # Scenario A -> pattern of LoE in both experimental and control arms
    beta.baseline_AE_control <- 31.031848
    beta_week1_AE_control <- 3.056887
    beta_week2_AE_control <- 5.611661
    beta_week3_AE_control <- 2.769695
    beta_week4_AE_control <- 2.991209 
    beta_week5_AE_control <- 3.552138
    beta_week6_AE_control <- 0#0.385184
    
    
    
    
    
    # following this model:
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6) 
    
    
    for (i in 1:(n_AE_control*length(visits))) {
      d_AE_control[i,4] <- ifelse(d_AE_control[i, 2]==0, beta.baseline_AE_control + d_AE_control[i,5],
                                  ifelse(d_AE_control[i, 2]==7, d_AE_control[i-1,4] + beta_week1_AE_control + d_AE_control[i,5],
                                         ifelse(d_AE_control[i, 2]==14, d_AE_control[i-2,4]+ beta_week2_AE_control + d_AE_control[i,5],
                                                ifelse(d_AE_control[i, 2]==21, d_AE_control[i-3,4] + beta_week3_AE_control + d_AE_control[i,5],
                                                       ifelse(d_AE_control[i, 2]==28, d_AE_control[i-4,4] + beta_week4_AE_control + d_AE_control[i,5],
                                                              ifelse(d_AE_control[i, 2]==35, d_AE_control[i-5,4] + beta_week5_AE_control + d_AE_control[i,5],
                                                                     d_AE_control[i-6,4] + beta_week6_AE_control + d_AE_control[i,5]))))))
    }
    
    
    
    #View(d)
    d_AE_control <-as.data.frame(d_AE_control)
    d_AE_control$visit <-as.factor(d_AE_control$visit)
    d_AE_control$Treat <- factor(d_AE_control$Treat)
    
    
    
    
    
    # this is the raw dataset used to check the model fit
    d_AE_control <-d_AE_control[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d_AE_control # full outcome data
    
    
    length(d_AE_control$id)
    tmp <- sapply(unique(d_AE_control$id), FUN = function(i) nrow(d_AE_control[d_AE_control$id == i,]))
    BaselineMADRS10 <-  rep(d_AE_control$MADRS10[d_AE_control$visit == 0], tmp)
    length(BaselineMADRS10)
    d_AE_control$Baseline <- BaselineMADRS10
    d_AE_control
    #d_AE_control<-d_AE_control[d_AE_control$visit!=0,]
    
    
    range(d_AE_control$MADRS10[d_AE_control$Treat==0 & d_AE_control$visit!=0])
    
    
    range(d_AE_control$MADRS10[d_AE_control$Treat==0 & d_AE_control$visit==42])
    
    
    #mean(d_AE_control$MADRS10[d_AE_control$Treat==0 & d_AE_control$visit==0])
    
    
    mean(d_AE_control$MADRS10[d_AE_control$Treat==0 & d_AE_control$visit==42])
    
    
    
    p<- ggplot(data = d_AE_control, aes(x = visit, y = MADRS10, group = id)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 70))
    
    #View(d_AE_control)
    
    d_AE_control$V0 <- 0
    d_AE_control$V0[d_AE_control$visit == 0] <- 1
    
    d_AE_control$V7 <- 0
    d_AE_control$V7[d_AE_control$visit == 7] <- 1
    
    d_AE_control$V14 <- 0
    d_AE_control$V14[d_AE_control$visit == 14] <- 1
    
    d_AE_control$V21 <- 0
    d_AE_control$V21[d_AE_control$visit == 21] <- 1
    
    d_AE_control$V28 <- 0
    d_AE_control$V28[d_AE_control$visit == 28] <- 1
    
    d_AE_control$V35 <- 0
    d_AE_control$V35[d_AE_control$visit == 35] <- 1
    
    d_AE_control$V42 <- 0
    d_AE_control$V42[d_AE_control$visit == 42] <- 1
    #d_AE_control
    
    
    
    
    fit_AE_control<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42, 
                        data=d_AE_control,
                        correlation = corSymm(form=~1 | id),
                        weights = varIdent(form = ~ 1 | visit), 
                        method="REML")
    
    
    
    
    summary(fit_AE_control)
    
    # indicator variable
    d_AE_control$Pattern <- c("AE_control")
    
    # indicator variable
    d_AE_control$LoE_Yes <- 0
    d_AE_control$AE_exp <- 0
    d_AE_control$AE_control <- 1
    d_AE_control$No_IE <- 0
    d_AE_control$AE_Yes <- d_AE_control$AE_exp + d_AE_control$AE_control
    #
    
    ##### PATTERN OF No IE IN BOTH EXPERIMENTAL AND CONTROL ARMS
    ### Generate pattern for LoE, both experimental and control arm
    ### Generate correlated errors
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    
    re_covm_No_IE <- matrix(
      c(10.040, 18.135, 16.595, 15.814, 13.524, 14.740, 11.582,
        18.135, 48.213, 42.211, 37.687, 30.909, 34.285, 28.580,
        16.595, 42.211, 54.236, 46.775, 44.082, 44.487, 37.226,
        15.814, 37.687, 46.775, 54.825, 46.690, 44.923, 35.807,
        13.524, 30.909, 44.082, 46.690, 56.719, 51.932, 43.362,
        14.740, 34.285, 44.487, 44.923, 51.932, 62.843, 52.026,
        11.582, 28.580, 37.226, 35.807, 43.362, 52.026, 51.272), nrow = 7)
    
    
    re_No_IE <- mvrnorm(n_No_IE, re_means, re_covm_No_IE)	; re_No_IE
    #View(re)
    
    re_No_IE <- as.matrix(re_No_IE)
    colnames(re_No_IE) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re_No_IE
    
    d_No_IE <- data.frame(
      id = rep((n_LoE_all+n_AE_exp+n_AE_control+1):(n_LoE_all+n_AE_exp+n_AE_control+n_No_IE), each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n_No_IE, 1, 0.5), each = length(visits)),
      MADRS10 = rep(NA, n_No_IE)); d_No_IE # mean(Treat)
    
    
    
    
    
    d_No_IE <- d_No_IE[order(d_No_IE$visit, d_No_IE$id),]; #d_No_IE
    #re
    
    j_No_IE<-c(re_No_IE[, 1], re_No_IE[, 2], re_No_IE[, 3], re_No_IE[, 4], re_No_IE[, 5], re_No_IE[, 6], re_No_IE[, 7])
    d_No_IE$re_No_IE <- j_No_IE ; #d_No_IE
    
    
    d_No_IE <- d_No_IE[order(d_No_IE$id, d_No_IE$visit),]; #d_No_IE
    
    
    d_No_IE<-as.matrix(d_No_IE)
    
    # Scenario A -> pattern of No IE in both experimental and control arms
    beta.baseline_No_IE <- 29.696133
    beta_week1_No_IE <- -2.061548
    beta_week2_No_IE <- -3.864898
    beta_week3_No_IE <- -4.360199
    beta_week4_No_IE <- -5.975185
    beta_week5_No_IE <- -7.406465
    beta_week6_No_IE <- -8.161270
    
    beta_v1_treatment_No_IE <- -0.469287
    beta_v2_treatment_No_IE <- -0.268673
    beta_v3_treatment_No_IE <- -0.724476
    beta_v4_treatment_No_IE <- -1.037980
    beta_v5_treatment_No_IE <- -1.420660
    beta_v6_treatment_No_IE <- 0#-4.195854 # -1.904595
    
    
    
    
    
    
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
    
    
    ####################################################################################################  
    # MMRM on full outcome data in pattern No IE in both experimental and control arm
    ############################
    
    
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
    
    
    
    p<- ggplot(data = d_No_IE, aes(x = visit, y = MADRS10, group = id)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-20, 70))
    
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
    
    
    fit_No_IE<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                     Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                   data=d_No_IE,
                   correlation = corSymm(form=~1 | id),
                   weights = varIdent(form = ~ 1 | visit), 
                   method="REML")
    
    
    summary(fit_No_IE)
    
    fit_No_IE$coefficients[c(13)]
    
    sum(fit_No_IE$coefficients[c(13)]); treatmenteffect_No_IE
    
    # indicator variable
    d_No_IE$Pattern <- c("No_IE")
    
    
    # indicator variable
    d_No_IE$LoE_Yes <- 0
    d_No_IE$AE_exp <- 0
    d_No_IE$AE_control <- 0
    d_No_IE$No_IE <- 1
    d_No_IE$AE_Yes <- d_No_IE$AE_exp + d_No_IE$AE_control
    
    #
    
    
    
    d_pmmm <- rbind(d_LoE_all,
                    d_AE_exp,
                    d_AE_control,
                    d_No_IE)
    
    
    head(d_LoE_all)
    head(d_AE_exp)
    head(d_AE_control)
    head(d_No_IE)
    #View(d_pmmm)
    
    describe(d_pmmm)
    head(d_pmmm)
    d_pmmm$Pattern <- as.factor(d_pmmm$Pattern)
    
    
    d_pmmm$LoE_Yes <- as.factor(d_pmmm$LoE_Yes)
    d_pmmm$AE_exp <- as.factor(d_pmmm$AE_exp)
    d_pmmm$AE_control <- as.factor(d_pmmm$AE_control)
    d_pmmm$No_IE <- as.factor(d_pmmm$No_IE)
    d_pmmm$AE_Yes <- as.factor(d_pmmm$AE_Yes)
    levels(d_pmmm$visit)
    
    
    
    # Plot trajectories
    # All patients with TRUE trajectory 
    
    
    # LoE_all/any
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=LoE_Yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    # AE_exp
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=AE_exp)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    # AE_control
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=AE_control)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    # AE_any
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=AE_Yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    # No_IE
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=No_IE)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    
    #describe(d_pmmm)
    #view(d_pmmm)
    
    # All behaviors
    p<- ggplot(data = d_pmmm, aes(x = visit, y = MADRS10, group = id, color=Pattern)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    fit_pmmm<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                    Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                  data=d_pmmm,
                  correlation = corCompSymm(form=~1 | id),
                  weights = varIdent(form = ~ 1 | visit), 
                  method="REML")
    
    
    summary(fit_pmmm)
    
    describe(d_pmmm)
    
    t.test(d_pmmm$MADRS10[d_pmmm$Treat==1 & d_pmmm$visit==42], d_pmmm$MADRS10[d_pmmm$Treat==0 & d_pmmm$visit==42])
    
    
    #describe(d_pmmm)
    getVarCov(fit_pmmm, individual = 1)
    
    treatmenteffect_pmmm <- sampled_prop_LoE * beta_v6_treatment_LoE_all +
      
      sampled_prop_AE_exp *  beta_week6_AE_exp +
      
      sampled_prop_AE_control *  beta_week6_AE_control +
      
      (1-(sampled_prop_LoE + sampled_prop_AE_exp + sampled_prop_AE_control)) *  beta_v6_treatment_No_IE
    
    
    sum(fit_pmmm$coefficients[c(13)]); treatmenteffect_pmmm
    
    
    
    treatmenteffect_pmmm <- sampled_prop_LoE * -1.909406 +
      
      sampled_prop_AE_exp * -12.660152 +
      
      sampled_prop_AE_control *  0.385184 +
      
      (1-(sampled_prop_LoE + sampled_prop_AE_exp + sampled_prop_AE_control)) *  -4.195854; treatmenteffect_pmmm
    
    
    
    
    
    # treatment effect in PMMMM based on models fitted on the SM with 5000 patients 
    #-1.318160*0.383+
    #-12.660152*0.0298+
    #0.385184*0.0158+
    #-1.904595*0.5714
    
    
    
    
    
    
    #
    
    
    betas[m, ] <- fit_pmmm$coefficients[c(13)]
    
    
    delta[m, ] <- sum(fit_pmmm$coefficients[c(13)])
    
    #bias_f[m, ] <- sum(fit_pmmm$coefficients[c(7,13)]) - treatmenteffect
    
    #delta_error <- sqrt(vcov(fit_pmmm)["Treat1", "Treat1"] + vcov(fit_pmmm)["visit42:Treat1", "visit42:Treat1"] + 2*vcov(fit_pmmm)["Treat1", "visit42:Treat1"]) 
    
    #delta_errorz[m, ] <- delta_error 
    
    
    #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_error
    #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_error
    
    
    
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




all_betas_1
colMeans(all_delta_1)

describe(SimTrial_pmmm_190_1_1$LoE_Yes[SimTrial_pmmm_190_1_1$Treat==1])
describe(SimTrial_pmmm_190_1_1$LoE_Yes[SimTrial_pmmm_190_1_1$Treat==0])


#### Visualisation of trajectories and patterns

# Plot trajectories
# All patients with TRUE trajectory 


# LoE_all/any
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=LoE_Yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


# AE_exp
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=AE_exp)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


# AE_control
p<- ggplot(data = SimTrial_pmmm_190_1_1, aes(x = visit, y = MADRS10, group = id, color=AE_control)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)


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




# with differing number of patients per pattern, the overall treatment effect will change obviously.
# Need to finetune the betas for each pattern for all the wanted treatment effects and factoring in the scaling factor.



### resume here
## need to save the generated datasets
# create simple function to sample few and visualise them


end_time-start_time

colMeans(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5))

colMeans(rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5))


cbind(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5),
      rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5),
      rbind(all_delta_errorz_1, all_delta_errorz_2, all_delta_errorz_3, all_delta_errorz_4, all_delta_errorz_5),
      rbind(all_confint_fit_1, all_confint_fit_2, all_confint_fit_3, all_confint_fit_4, all_confint_fit_5),
      rbind(all_N_Exp_1, all_N_Exp_2, all_N_Exp_3, all_N_Exp_4, all_N_Exp_5),
      rbind(all_N_Control_1, all_N_Control_2, all_N_Control_3, all_N_Control_4, all_N_Control_5))




