# License CC-BY-4.0 ----
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
# visualisation of longitudinal outcomes and intercurrent events

## load libraries ----
rm(list=ls()) #
# needed for the selection model method
library(gmailr)
library(MASS)#
library(tidyverse)#
library(nlme)#
library(lme4)#
library(Hmisc)
library(janitor)#
library(gt)#
library(patchwork)#


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

n <- 190 # number of patients to be simulated (sample size)
# this is based on a t-test to ensure  90% power at alpha level=0.025 one-sided 

m.iterations <- 382 # 382 is the number of trials needed for the verification of the longitudinal outcomes # number of generated datasets # number of trials per scaling factor
scaling_factor <- c(1) # c(0.5, 1.0, 1.5, 2.0, 2.5) # scaling factor used to vary the percentages of intercurrent events at trial/iteration level
# total number of simulated trials = m.iterations * length(scaling_factor)
# other ranges can be used to ensure variability between simulated trials, as long as they are as envisaged over all simulated trials (e.g., mean percentages)
# and check out the verification step


# ranges of probabilities centered around desired percentages of each intercurrent events averaged over all simulated trials
# this is done to increase variability in intercurrent events percentages between trials
# aim is to obtain 0.35 LoE
# AE exp arm 0.1, and AE control arm 0.5
# the probabilities below should make possible to obtain the desired percentages of intercurrent events

# for 0.35 LoE
p_LoE_sample <- 0.83 #c(0.31, 0.33, 0.34, 0.35, 0.37); mean(p_LoE_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy at trial level
#for 0.1 AE exp arm
p_AE_Exp_sample <- 0.45 #c(0.055, 0.065, 0.075, 0.085, 0.095)*2; mean(p_AE_Exp_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy in the experimental arm
#for 0.05 AE control arm
p_AE_Control_sample <- 0.4 #c(0.05, 0.10, 0.15, 0.20, 0.25)/2; mean(p_AE_Control_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy in the control arm

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline (day 0) + follow-up measurements at day 7, 14, ..., 42
delta <- matrix(ncol=1,nrow=m.iterations) # object to store treatment effect estimates at 6 weeks based on MMRM model fitted on each generated dataset
colnames(delta) <-c("TreatmentEffect")
betas <- matrix(ncol=2,nrow=m.iterations) # object to store parameters for the treatment effect at week 6 based on the MMRM model fitted on each generated dataset
colnames(betas) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3) # progress bar in percentages relative to the total number of m.iterations

confint_fit <- matrix(ncol=2,nrow=m.iterations) # object to store the 95% confidence interval bounds for the estimated treatment effect
colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
delta_errorz <- matrix(ncol=1,nrow=m.iterations) # standard error of the estimated treatment effect
colnames(delta_errorz) <- c("SE")

bias_f <- matrix(ncol=1,nrow=m.iterations) # object to store the bias (estimated treatment effect - true treatment effect)
colnames(bias_f) <- c("bias_f")

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

pb3 <- txtProgressBar(min = 0,  max=length(scaling_factor), style=3) # # progress bar in percentages relative to the total number of scaling factors

    #s<-1
    #m<-1
    
start_time <- Sys.time() # timestamp for the start time of the nested for loop below.
# it was used to have an estimate of time needed for different larger number of trials to be simulated upon scaling up the simulation parameters (e.g., m.iterations)

## Begin for loop----
for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    ### Generate longitudinal outcomes----
    #### Generate correlated residuals----
    re_means <- c(0, 0, 0, 0, 0, 0, 0) # means of residuals at each timepoint
    
    # covariance matrix extracted from trial 003-002 (see reference to the data analysis paper in PST)
    re_covm <- matrix(c(20.2190, 17.149, 14.721, 13.087,  8.4329,  10.854,   4.6417, 
                        17.1490, 48.536, 41.161, 32.151, 24.8400,  30.528,  26.0170,
                        14.7210, 41.161, 72.569, 57.866, 60.2200,  61.974,  54.5400,
                        13.0870, 32.151, 57.866, 74.080, 66.2960,  63.540,  52.1070,
                        8.4329, 24.840, 60.220, 66.296, 97.4730,  90.612,  80.1370,
                        10.8540, 30.528, 61.974, 63.540, 90.6120, 116.410, 102.8300,
                        4.6417, 26.017, 54.540, 52.107, 80.1370, 102.830, 109.5900), nrow = 7)
    
    #Standard Deviations: 4.4965 6.9668 8.5188 8.607 9.8728 10.789 10.468
    
    # covariance matrix with 0 off-diagonal and small variances. This is useful for initial/later checks to see if the simulated data corresponds to target data to be simulated
    #re_covm2 <-matrix(c(0.00001, 0, 0, 0, 0, 0, 0,
     #                   0, 0.00001, 0, 0, 0, 0, 0,
      #                  0, 0, 0.00001, 0, 0, 0, 0,
       #                 0, 0, 0, 0.00001, 0, 0, 0,
        #                0, 0, 0, 0, 0.00001, 0, 0,
         #               0, 0, 0, 0, 0, 0.00001, 0,
          #              0, 0, 0, 0, 0, 0, 0.00001), nrow = 7)
    
    re_covm3 <-re_covm/2 # we scaled the covariance matrix as we observed the trajectories were not similar with the targeted trajectories. This is up to user's decision
    
    re <- mvrnorm(n, re_means, re_covm3)	; re # generate correlated residuals and check them
    #View(re)
    re <- as.matrix(re)
    colnames(re) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re
    
    # put together the dataframe with the id of patients, visit, randomised treatment and MADRS10 values at each visit
    d <- data.frame(
      id = rep(1:n, each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n, 1, 0.5), each = length(visits)),
      MADRS10 = rep(NA, n)); d # mean(Treat)
    
    # data wrangling to add the correlated residuals in the above dataframe such that correlated residuals are for their corresponding patient and visit 
    d <- d[order(d$visit, d$id),]; d
        #re
    
    j<-c(re[, 1], re[, 2], re[, 3], re[, 4], re[, 5], re[, 6], re[, 7])
    d$re <-j; #d
    
        #class(d)
        #View(d)
        #head(re)
    
    d <- d[order(d$id, d$visit),]; #d
    
        #re
    
    d<-as.matrix(d)

    # MMRM parameters to generate the longitudinal outcomes
    # We used as inspiration the trial 003-002 (ref data analysis paper). Here we used a simplified scenario with linear group trajectories.
    
    #### model parameters----
    # for the interpretation of these parameters, please see the table under the selection model method in the body of the paper
    beta.baseline <- 29.79
    beta_week1 <- -1
    beta_week2 <- -1.5
    beta_week3 <- -2
    beta_week4 <- -2.5 
    beta_week5 <- -3
    beta_week6 <- -3.5
    
    beta_v1_treatment <- -1
    beta_v2_treatment <- -1.5
    beta_v3_treatment <- -2
    beta_v4_treatment <- -2.5
    beta_v5_treatment <- -3
    beta_v6_treatment <- -3.5 
    
    treatmenteffect <-  beta_v6_treatment ; treatmenteffect
    
    # The model used to generate the correlated repeated outcomes
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6) + 
    #                          Beta_W1_Treat * T + Beta_W2_Treat * T + Beta_W3_Treat * T + Beta_W4_Treat * T + Beta_W5_Treat * T + Beta_W6_Treat * T 
    
    
    # for loop to generate the correlated repeated outcomes that follows the above model
    for (i in 1:(n*length(visits))) {
      d[i,4] <- ifelse(d[i, 2]==0, beta.baseline + d[i,5],
                       ifelse(d[i, 2]==7, d[i-1,4] + beta_week1 + d[i,5] +  beta_v1_treatment * d[i, 3],
                              ifelse(d[i, 2]==14, d[i-2,4]+ beta_week2 + d[i,5] +  beta_v2_treatment * d[i, 3],
                                     ifelse(d[i, 2]==21, d[i-3,4] + beta_week3 + d[i,5] +  beta_v3_treatment * d[i, 3],
                                            ifelse(d[i, 2]==28, d[i-4,4] + beta_week4 + d[i,5] +  beta_v4_treatment * d[i, 3],
                                                   ifelse(d[i, 2]==35, d[i-5,4] + beta_week5 + d[i,5] +  beta_v5_treatment * d[i, 3],
                                                          d[i-6,4] + beta_week6 + d[i,5] +  beta_v6_treatment * d[i, 3]))))))
    }
    
    # check to see how many values were negative or > 60  
          ## flooring and ceiling
          #d[, 4] <- ifelse(d[, 4] < 0, 0,
          #                ifelse(d[, 4]>60, 60, d[, 4]))
    
    d <-as.data.frame(d)
    d$visit <-as.factor(d$visit)
    #d$Treat <- factor(d$Treat)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MMRM on full outcome data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # this is the raw dataset used to check the model fit
    d <-d[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d # full outcome data
    
    # create a separate column with only the baseline outcomes, if the baseline values will be used as covariate in the model 
    length(d$id)
    tmp <- sapply(unique(d$id), FUN = function(i) nrow(d[d$id == i,]))
    BaselineMADRS10 <-  rep(d$MADRS10[d$visit == 0], tmp)
    length(BaselineMADRS10)
    d$Baseline <- BaselineMADRS10
        #d
    
        #d<-d[d$visit!=0,]
    
        #check ranges of the generated outcomes
        #range(d$MADRS10[d$Treat==1 & d$visit!=0])
        #range(d$MADRS10[d$Treat==0 & d$visit!=0])
    
        #range(d$MADRS10[d$Treat==1 & d$visit==42])
        #range(d$MADRS10[d$Treat==0 & d$visit==42])
    
        #mean(d$MADRS10[d$Treat==1 & d$visit==0])
        #mean(d$MADRS10[d$Treat==0 & d$visit==0])
    
        #mean(d$MADRS10[d$Treat==1 & d$visit==42])
        #mean(d$MADRS10[d$Treat==0 & d$visit==42])
    
        #ceiling_floor <- sum(ifelse(d$MADRS10<0 | d$MADRS10>60, 1, 0))
        #ceiling_floor_perc <- (ceiling_floor/(n*length(visits)))*100
    
        #CFE[s*m,1] <- ceiling_floor
        #CFE[s*m,2] <- ceiling_floor_perc
        #CFE[s*m,3] <- s
        #CFE[s*m,4] <- m
    
    
    
    
    #### plot of longitudinal outcomes ----
    # plot the outcomes to see in big lines how they trajectories look like  
    p<- ggplot(data = d, aes(x = visit, y = MADRS10, group = id)) 
    plot1 <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)) ; plot1
    
    
    # fit a model to check if the estimated parameters are similar/close to the true parameters
    fit <- gls(MADRS10 ~ visit * Treat + Baseline, 
             data=d,
             correlation = corSymm(form=~1 | id),
             #weights = varIdent(form = ~ 1 | visit),
             method="REML")
    
    summary(fit)
    
        #sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"])
        #sqrt(vcov(fit)["Treat", "Treat"])
        #re_covm2
        #fit$coefficients[c(7,13)]
        #sum(fit$coefficients[c(7,13)]); treatmenteffect
        #d$visit <- as.numeric(d$visit)-1
    
        # check the same with another package and models
        #fit_lme <- lme(fixed=MADRS10 ~ visit * Treat + Baseline, 
        #          random=~1 + visit | id,
        #        method="REML", 
        #      correlation = corSymm(form=~1|id),
        #    data=d)
    
        #fit_lmer <- lmer(MADRS10 ~ visit + visit:Treat + (1 |id), data = d, REML = T)
    
        #summary(fit_lmer)
        #summary(fit_lme)
        #model_parameters(fit_lme)
    #### store estimated parameters----
    betas[m, ] <- fit$coefficients[c(8,15)] # store the parameters corresponding to the treatment effect at the end of the trial, at week 6
    
    delta[m, ] <- sum(fit$coefficients[c(8,15)]) # store the treatment effect at the end of the trial, at week 6
    
        #bias_f[m, ] <- sum(fit$coefficients[c(7,13)]) - treatmenteffect
        #delta_error <- sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"]) 
        #delta_errorz[m, ] <- delta_error 
        #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_error
        #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_error
    
    # store the number of patients in the objects defined a priori
    Randomised_Exp <- sum(d[,3])/7 #number of patients in the experimental arm # make sure it is divided by the number of visits (6/7), without or with the baseline included as response
    Randomised_Control <- n-sum(d[,3])/7 #number of patients in the control arm
    #m<-s<-1
    N_Exp[m,] <- Randomised_Exp
    N_Control[m,] <- Randomised_Control

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Intercurrent events generating model (IEGM)----
    # create variable with difference week 6 - baseline
    # if Difference < 5, then missing outcome from visit 3 onwards for LoE
    
    ## take over the original/raw dataset to use for IEGM
    d_mis <-d_orig
    
    d_mis_w <- d_mis %>% spread(visit, MADRS10) # reshape to wide in order to create the CfB variable
        #View(d_mis_w)
    
    colnames(d_mis_w)[3:9] <- c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6"); head(d_mis_w) # name the columns for the wide format dataframe
    
    d_mis_w$CfB <- d_mis_w[,3] - d_mis_w[,9]; d_mis_w # create the CfB (change from baseline) variable. Could also be d_mis_w[,9] - d_mis_w[,3]
    
    d_mis_w$CfW1 <- d_mis_w[,3] - d_mis_w[,4]; d_mis_w # create the CfW1 (change from week 1) variable
        #s<-1
        #View(d_mis_w)
    p_LoE <- sample(p_LoE_sample, 1) * scaling_factor[s] # proportion of e.g., treatment discontinuation due to lack of efficacy multiplied by the scaling factor. 
    # this allows to obtain different smaller/larger percentages of LoE at trial level and cover a broad range of scenarios as would be in reality, and to ensure variability,
    # not all patients that it in the rule will automatically experience the intercurrent event
    
    
    # determine adjustment factors (parameters) in order to obtain desired percentages of intercurrent events
    # simulate a trial with samplesize 10^5, see how many patients fit the IE rules
    # then determine adjustment factors for the Bernoulli process in order to obtain the desired percentages of intercurrent events
    
    
    
    
    d_mis_w$LoE_Yes <-ifelse(d_mis_w$CfB<5, 1, 0)*rbinom(n, 1, p_LoE);  d_mis_w # all patients that fit in the deterministic rule are then adjusted with probability abovementioned
    # to adjust the probabilty of LoE# create the LoE variable
    #35000/sum(ifelse(d_mis_w$CfB<5, 1, 0))
    sum(d_mis_w$LoE_Yes) # check how many patients experienced the intercurrent event
        #View(d_mis_w)  
    
    # same method is applied for e.g., treatment discontinuations due to adverse events at arm level
    # too much efficacy >8 points on MADRS10 at week 2
    
    p_AE_Exp <- sample(p_AE_Exp_sample, 1) * scaling_factor[s]
    
    d_mis_w$AE_Exp_Yes <-ifelse(d_mis_w$Treat==1 & d_mis_w$CfW1>3, 1, 0)*rbinom(n, 1, p_AE_Exp);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    #sum(ifelse(d_mis_w$Treat==1 & d_mis_w$CfW1>3, 1, 0))
    
    sum(d_mis_w$AE_Exp_Yes)
    
    p_AE_Control <- sample(p_AE_Control_sample, 1) * scaling_factor[s]
    
    d_mis_w$AE_Control_Yes <-ifelse(d_mis_w$Treat==0 & d_mis_w$CfW1< (-2), 1, 0)*rbinom(n, 1, p_AE_Control);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    #sum(ifelse(d_mis_w$Treat==0 & d_mis_w$CfW1< (-2), 1, 0))
    sum(d_mis_w$AE_Control_Yes)

    sum(d_mis_w$LoE_Yes)
    
    describe(ifelse(d_mis_w$Treat==0 & d_mis_w$CfW1< (-2), 1, 0))
    
    #
    
    d_mis_L <- d_mis_w %>% gather(Visit, MADRS10, Baseline:Week6) # reshape to long format
    
    d_mis_L <- d_mis_L[order(d_mis_L$id, d_mis_L$Visit),]; #d # order by subject id and Visit
  
        # not needed for the purpose of this paper. this is needed for the next paper with application of this DGMs paper. This can be skipped as well as the step below involving missing data
        # create the NAs by LoE variable
        #d_mis_L[,10] <-ifelse(d_mis_L[,5]==1 & d_mis_L[,8]=="Week4", "NA", 
        #                     ifelse(d_mis_L[,5]==1 & d_mis_L[,8]=="Week5", "NA",
        #                           ifelse(d_mis_L[,5]==1 & d_mis_L[,8]=="Week6", "NA", d_mis_L[,9])))
    
        #View(d_mis_L)
    
        #colnames(d_mis_L)
        # [1] "id"             "Treat"          "CfB"           
        # [4] "CfW2"           "LoE_Yes"        "AE_Exp_Yes"    
        # [7] "AE_Control_Yes" "Visit"          "MADRS10"       
        #[10] "V10"  
    
    
    #View(d_mis_L)
    

    
    d_mis_L$AE_Yes <- ifelse(d_mis_L$AE_Exp_Yes==1, 1, 
                             ifelse(d_mis_L$AE_Control_Yes==1, 1, 0))
    
        # create the NAs by AE variable
        #d_mis_L[,10] <-ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week3", "NA",
        #                     ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week4", "NA", 
        #                           ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week5", "NA",
        #                                 ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week6", "NA", d_mis_L[,10]))))
    
        #View(d_mis_L)
    
    # determine competing AE and LoE on the same patients
    #d_mis_L$compete <- ifelse(d_mis_L$AE_Yes==1 & d_mis_L$LoE_Yes==1, 1, 0)
    
    #sum(ifelse(d_mis_L$AE_Yes==1 & d_mis_L$LoE_Yes==1, 1, 0))/7
    
    #sum(as.numeric(d_mis_L[,c("AE_Yes")])-1)/7
    
    #(sum(d_mis_L$compete)/7)/(sum(as.numeric(d_mis_L$LoE_Yes)-1)/7)*100# this is the percentage adjustment that needs to be added to achieve targeted proportion of LoE
    
    
    d_mis_L$LoE_YES <- ifelse(d_mis_L$LoE_Yes==1 & d_mis_L$AE_Yes==0, 1, 0)
    describe(d_mis_L$LoE_YES)
    
    d_mis_L <- d_mis_L[order(d_mis_L$id),]; #d # order by subject id and Visit
        #colnames(d_mis_L)[10] <- c("MADRS10_mis"); d_mis_L # name the column. See above comment about missing data. This is not needed here, but in the companion paper
    
    # check classes of variable and reclassify as needed
    class(d_mis_L$id)
    class(d_mis_L$Treat)
    d_mis_L$Treat<-as.factor(d_mis_L$Treat)
    class(d_mis_L$Visit)
    d_mis_L$LoE_Yes <- as.factor(d_mis_L$LoE_Yes)
    d_mis_L$LoE_YES <- as.factor(d_mis_L$LoE_YES)
    d_mis_L$AE_Yes <- as.factor(d_mis_L$AE_Yes)
    
    # create the Behavior column/variable
    d_mis_L$Behavior <- ifelse(d_mis_L[,10]==1, "AE",
                               ifelse(d_mis_L[,11]==1, "LoE", "No IE"))
    
    d_mis_L$Behavior <-factor(d_mis_L$Behavior, levels = c("AE", "LoE", "No IE"))
    
    class(d_mis_L$Visit)
    d_mis_L$Visit <- as.factor(d_mis_L$Visit)
    rownames(d_mis_L) <-NULL 
    
        # check range of values
        #range(d_mis_L$MADRS10[d_mis_L$Treat==1])
        #range(d_mis_L$MADRS10[d_mis_L$Treat==0])
    
    #### assign and save the generated datasets----
    # naming sequence is "SimTrial"_"Method"_"trial sample size"_"iteration number"_"scaling factor"

    assign(paste0("SimTrial_sm", "_", n,"_", m, "_", s), d_mis_L)
    #View(SimTrial_sm_1_5)
    dataset_name.Rdata <- paste0("SimTrial_sm", "_", n,"_", m, "_", s, ".Rdata")
    dataset_name <- paste0("SimTrial_sm", "_", n,"_", m, "_", s)
    save(dataset_name, file = dataset_name.Rdata)
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Intercurrent events descriptives needed for the verification step ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    LoE_Y <- d_mis_L[,c("Treat", "LoE_YES")]
    LoE_Y$LoE_YES <- as.numeric(LoE_Y$LoE_YES)-1
    
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

    AE_Y <- d_mis_L[,c("Treat", "AE_Yes")]
    AE_Y$AE_Yes <- as.numeric(AE_Y$AE_Yes)-1
    
    AE_Y_total <- sum(AE_Y$AE_Yes)/length(visits)
    AE_Y_Exp <- sum(AE_Y$AE_Yes[AE_Y$Treat==1])/length(visits)
    AE_Y_Control <- sum(AE_Y$AE_Yes[AE_Y$Treat==0])/length(visits)
    
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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Plot trajectories for all and all intercurrent events----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # All patients with their trajectory 
  
    # LoE #
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=LoE_YES)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    # AE #
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=AE_Yes)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
        #describe(d_mis_L$Behavior)
  
    # All behaviors #
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=Behavior)) 
    #p + geom_line() + facet_grid(~ Treat) 
    plot_all_SM <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="black") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SM-All patterns") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_all_SM
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##### Plots by subsets/patterns of patients that experienced a certain intercurrent event----
    
    d_mis_L_LoE <- d_mis_L[d_mis_L$LoE_YES==1,] # subset only patients that experienced LoE
    d_mis_L_AE <- d_mis_L[d_mis_L$AE_Yes==1,] # subset only patients that experienced AE
    d_mis_L_NoIE <- d_mis_L[d_mis_L$LoE_YES==0 & d_mis_L$AE_Yes==0,] # subset only patients that did not experience any IE
    
    # just (subset) LoE patients with their trajectory #
    p<- ggplot(data = d_mis_L_LoE, aes(x = Visit, y = MADRS10, group = id))
    plot_LoE_SM <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SM-LoE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_LoE_SM

    
    # just AE patients with their trajectory #
    p<- ggplot(data = d_mis_L_AE, aes(x = Visit, y = MADRS10, group = id))
    plot_AE_SM <- p + geom_line(size=0.5, color='#F8766D') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SM-AE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_AE_SM
    
    # just No IE patients with their trajectory # 
    p<- ggplot(data = d_mis_L_NoIE, aes(x = Visit, y = MADRS10, group = id))
    plot_NoIE_SM <- p + geom_line(size=0.5, color='#619CFF') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="blue") + facet_wrap(~ Treat) +
     scale_y_continuous(limits = c(-10, 60)) + ggtitle("SM-No IEs pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_NoIE_SM
    
    
    
    
    setTxtProgressBar(pb1, m)
  }

  # parameters extracted for MMRM fitted models on full outcome data
  colMeans(betas) # average of treatment effect parameters estimated from the model
  colMeans(delta) ; treatmenteffect # average treatment effect estimated from the model
  
  # assign and save all parameters----
  # pertaining to the treatment effect, the estimated treatment effect, the standard errors, 95% CI and intercurrent event descriptives
  assign(paste('all_betas', s, sep="_"), betas)
  
  assign(paste('all_delta', s, sep="_"), delta)
  
  #assign(paste('all_delta_errorz', s, sep="_"), delta_errorz)
  
  #assign(paste('all_confint_fit', s, sep="_"), confint_fit)
  
  #assign(paste('all_N_Exp', s, sep="_"), N_Exp)
  
  #assign(paste('all_N_Control', s, sep="_"), N_Control)
  

  setTxtProgressBar(pb3, s)
}

# end for loop----
end_time <- Sys.time() # timestamp for end of simulation
end_time-start_time # total time to complete the simulation



betas; 
colMeans(delta); treatmenteffect

tolerance_margin <- 0.1 
difference_Verification <- abs(treatmenteffect - colMeans(delta))

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_Verification) < tolerance_margin), "Verification SUCCESSFUL", "Verification NOT successful") 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot relevant graphs for the paper----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

the_plot_SM <- (plot_all_SM / plot_LoE_SM) | (plot_AE_SM / plot_NoIE_SM); the_plot_SM






# Table for the paper ----

table_AE_SMd <- data.frame(
# descriptives AE  
n_AE_Control,
n_AE_Exp); table_AE_SM

mean(n_AE_Control)
mean(n_AE_Exp)




# descriptives LoE  
table_LoE_SMd <-data.frame(
n_LoE_Control,
n_LoE_Exp); table_LoE_SMd

mean(n_LoE_Control)
mean(n_LoE_Exp)


#describe(table_IE_SM)


table_AE_SMd %>% 
  as.data.frame() %>% 
  mutate("Intercurrent event" = "AE") %>% 
  rename(N_C_arm=N.AE.Control) %>% 
  rename(N_E_arm=N.AE.Exp)

table_LoE_SMd %>% 
  as.data.frame() %>% 
  mutate("Intercurrent event" = "LoE") %>% 
  rename(N_C_arm=N.LoE.Control) %>% 
  rename(N_E_arm=N.LoE.Exp)


tab_SMd <- tibble(bind_rows(table_AE_SMd %>% 
            as.data.frame() %>% 
            mutate("Intercurrent event" = "AE") %>% 
            rename(N_C_arm=N.AE.Control) %>% 
            rename(N_E_arm=N.AE.Exp), 
          table_LoE_SMd %>% 
            as.data.frame() %>% 
            mutate("Intercurrent event" = "LoE") %>% 
            rename(N_C_arm=N.LoE.Control) %>% 
            rename(N_E_arm=N.LoE.Exp))); tab_SMd

  

tab2_SMd<- tab_SMd %>% group_by(`Intercurrent event`) %>%
  summarise("N" = round(mean(N_C_arm), digits=1), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = round(mean(N_E_arm), digits=1), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = round(mean(N_C_arm + N_E_arm), digits=1),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) %>% 
  adorn_totals("row"); tab2_SMd



  
gt(tab2_SMd) %>% 
  tab_header(title = md("Table 4. Descriptive statistics intercurrent events"), subtitle = md("Selection model DGM - deterministic rule")) %>%
  tab_source_note(md(paste0("Averaged over", " ", m.iterations,  " ",  "simulated trials.", " ", "Trial sample size = ", " ", n ))) %>% 
tab_spanner(
  label = md("**Control**"),
  columns = c("N", "%")) %>% 
  cols_align(
    align = "center",
    columns =  c("N", "%")
  ) %>% 
  tab_spanner(
    label = md("**Treatment**"),
    columns = c("N ", "% ")) %>% 
  cols_align(
    align = "center",
    columns =  c("N ", "% ")
  ) %>% 
  tab_spanner(
    label = md("**Total**"),
    columns = c(" N ", " % ")) %>% 
  cols_align(
    align = "center",
    columns =  c(" N ", " % ")
  ) %>% 
  data_color(
    columns = c("%", "% ", " % "),
    colors = scales::col_numeric(
      palette = c(
      "light blue"),
      domain = NULL)
  ) %>% 
  cols_align(
    align = "center",
    columns =  "Intercurrent event"
  ) %>% 
  tab_style(
    style = list(
      cell_fill(color = "white"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = "Intercurrent event"
    )
)







# determine the number of trials needed to simulate for the verification of the longitudinal outcomes
# Formula from Burton paper
#tolerance_margin <- 0.1 # bias allowed
#std.e <- 0.997 # model-based standard error of the treatment effect estimate from a fitted model on 1 trial

#n.trials_needed <- ceiling(((qnorm(0.975) * std.e)/tolerance_margin)^2) ; n.trials_needed # for the verification 
# 382 trials
# verification of the longitudinal outcomes was successful


## Session info---- 
sessionInfo()
installed.packages()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The code used to extract the MMRM (longitudinal) models for each pattern in preparation for PMMMM method----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Source trial with 1000 patients per arm
# Pattern for LoE at trial level (as per rule above)
# Pattern for AE in experimental arm (as per rule above)
# Pattern for AE in control arm (as per rule above)
# Pattern for No intercurrent event at trial level (remaining patients)

# 50% and then use scaling factors to get to 5% the entire range of IE percentages

# Use simulated trial with highest percentage of IE. Any other trial or parameters can be used. We chose it for operational reasons and model-fitting optimisation
# e.g., to have more patients in each pattern such that models fit

    #View(SimTrial_sm_2000_1_5)
SimTrial_sm_2000_1_5 <-SimTrial_sm_2000_1_1
describe(SimTrial_sm_2000_1_5)

# prepare the dataset for reparameterisation to not have any treatment coefficient at baseline
# as per the formulated model (see above)
SimTrial_sm_2000_1_5$V0 <- 0
SimTrial_sm_2000_1_5$V0[SimTrial_sm_2000_1_5$Visit == "Baseline"] <- 1

SimTrial_sm_2000_1_5$V7 <- 0
SimTrial_sm_2000_1_5$V7[SimTrial_sm_2000_1_5$Visit == "Week1"] <- 1

SimTrial_sm_2000_1_5$V14 <- 0
SimTrial_sm_2000_1_5$V14[SimTrial_sm_2000_1_5$Visit == "Week2"] <- 1

SimTrial_sm_2000_1_5$V21 <- 0
SimTrial_sm_2000_1_5$V21[SimTrial_sm_2000_1_5$Visit == "Week3"] <- 1

SimTrial_sm_2000_1_5$V28 <- 0
SimTrial_sm_2000_1_5$V28[SimTrial_sm_2000_1_5$Visit == "Week4"] <- 1

SimTrial_sm_2000_1_5$V35 <- 0
SimTrial_sm_2000_1_5$V35[SimTrial_sm_2000_1_5$Visit == "Week5"] <- 1

SimTrial_sm_2000_1_5$V42 <- 0
SimTrial_sm_2000_1_5$V42[SimTrial_sm_2000_1_5$Visit == "Week6"] <- 1
SimTrial_sm_2000_1_5

    #head(SimTrial_sm_2000_1_5)
    #class(SimTrial_sm_2000_1_5$Treat)

fit<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
           Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
         data= SimTrial_sm_2000_1_5,
         correlation = corSymm(form=~1 | id),
         #weights = varIdent(form = ~ 1 | Visit), 
         method="REML")

    #View(SimTrial_sm_2000_1_5)

summary(fit)
    #getVarCov(fit, individual = 1)
    #describe(SimTrial_sm_2000_1_5)


## Pattern for LoE at trial level----
fit_LoE_trial<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                     Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
          data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$LoE_YES==1,] ,
          correlation = corSymm(form=~1 | id),
          #weights = varIdent(form = ~ 1 | Visit),
          method="REML")

    #View(SimTrial_sm_2000_1_5)


summary(fit_LoE_trial)
getVarCov(fit_LoE_trial, individual = 8)

## Pattern for AE in experimental arm----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fit_AE_exp<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42, 
                    data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Exp_Yes==1,] ,
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | Visit),
                    method="REML")

summary(fit_AE_exp)
getVarCov(fit_AE_exp, individual = '2')

    #View(SimTrial_1_5[SimTrial_1_5$AE_Exp_Yes==1,])
    #corMatrix(fit_AE_exp$modelStruct$corStruct)[[1]]  

## Pattern for AE in control arm----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fit_AE_control<-gls(MADRS10 ~ MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                      Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                    data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Control_Yes==1,] ,
                correlation = corSymm(form=~1 | id),
                weights = varIdent(form = ~ 1 | Visit),
                method="REML")

summary(fit_AE_control)
getVarCov(fit_AE_control, individual = '28')

    #describe(SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Control_Yes==1,])
    #View(SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Control_Yes==1,])

## Pattern for no intercurrent events at trial level----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fit_no_IE<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                 Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
               data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Behavior== "No IE", ] ,
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | Visit),
                    method="REML")

summary(fit_no_IE)
    #getVarCov(fit_no_IE, individual = 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  The code used to extract the logit models for LoE at trial level, AE in experimental arm and AE in control arm----
## obtain the logistic regression models for LoE at trial level, AE in experimental arm and AE in control arm ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class(SimTrial_sm_2000_1_5$AE_Control_Yes)
SimTrial_sm_2000_1_5$AE_Control_Yes <- factor(SimTrial_sm_2000_1_5$AE_Control_Yes)

class(SimTrial_sm_2000_1_5$AE_Exp_Yes)
SimTrial_sm_2000_1_5$AE_Exp_Yes <- factor(SimTrial_sm_2000_1_5$AE_Exp_Yes)

trial_AE_exp <- SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Treat==1,]
trial_AE_control <- SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Treat==0,]


## LoE at trial level----
#logit(Pr(LoE))= fi_LoE * (Y_i0 - Y_i6)
logit_LoE <- glm(LoE_Yes ~ CfB,
                 data = SimTrial_sm_2000_1_5,
                 family = "binomial")

summary(logit_LoE)
#View(SimTrial_sm_2000_1_5)

predicted_LoE <- predict(logit_LoE, type="response", newdata=SimTrial_sm_2000_1_5)

SimTrial_sm_2000_1_5$predicted_LoE <- predicted_LoE
    #View(SimTrial_sm_2000_1_5[,c(5, 20)])


probz_predicted_LoE_Yes <- SimTrial_sm_2000_1_5$predicted_LoE[SimTrial_sm_2000_1_5$LoE_Yes==1]
up_boundary_prob_LoE <- max(probz_predicted_LoE_Yes)  ; up_boundary_prob_LoE
low_boundary_prob_LoE <- min(probz_predicted_LoE_Yes) ; low_boundary_prob_LoE

## AE experimental arm----
#logit(Pr(AE_exp))= fi_AE_exp * (Y_i0 - Y_i2)
logit_AE_exp <- glm(AE_Exp_Yes ~ CfW1,
                 data = trial_AE_exp,
                 family = "binomial")

    #summary(logit_AE_exp)

predicted_AE_exp <- predict(logit_AE_exp, type="response", newdata = trial_AE_exp)

trial_AE_exp$predicted_AE_exp <- predicted_AE_exp

    #View(trial_AE_exp[,c(6, 20)])
    #View(trial_AE_exp)

probz_predicted_AE_exp_Yes <- trial_AE_exp$predicted_AE_exp[trial_AE_exp$AE_Exp_Yes==1]
up_boundary_prob_AE_exp <- max(probz_predicted_AE_exp_Yes)  ; up_boundary_prob_AE_exp
low_boundary_prob_AE_exp <- min(probz_predicted_AE_exp_Yes)  ;  low_boundary_prob_AE_exp

## AE  control arm----
#logit(Pr(AE_control))= fi_AE_control * (Y_i0 - Y_i2)
logit_AE_control <- glm(AE_Control_Yes ~ CfW1,
                    data = trial_AE_control,
                    family = "binomial")

summary(logit_AE_control)

predicted_AE_control <- predict(logit_AE_control, type="response", newdata = trial_AE_control)

trial_AE_control$predicted_AE_control <- predicted_AE_control
    #View(trial_AE_control[,c(7, 20)])

probz_predicted_AE_control_Yes <- trial_AE_control$predicted_AE_control[trial_AE_control$AE_Control_Yes==1]
up_boundary_prob_AE_control <- max(probz_predicted_AE_control_Yes)  ; up_boundary_prob_AE_control
low_boundary_prob_AE_control <- min(probz_predicted_AE_control_Yes)  ; low_boundary_prob_AE_control

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The code for SPM method. On the same trial generated with the selection model (above), we fitted a mixed-effects model with random intercept and random slope # 
# fit LMM model on the simulated trial
# fit the glm to get the intercurrent events models with random effects in the linear predictor

SimTrial_sm_2000_1_5 <- SimTrial_sm_2000_1_1
class(SimTrial_sm_2000_1_5$Visit)
SimTrial_sm_2000_1_5$Visit <- as.numeric(SimTrial_sm_2000_1_5$Visit)-1

SimTrial_sm_2000_1_5$MADRS10 <- ifelse(SimTrial_sm_2000_1_5$Visit>2 & SimTrial_sm_2000_1_5$Behavior== "AE", NA,
                                       ifelse(SimTrial_sm_2000_1_5$Visit>3 & SimTrial_sm_2000_1_5$Behavior== "LoE", NA, SimTrial_sm_2000_1_5$MADRS10))



fit_lmer <- lmer(MADRS10 ~ Visit + Visit : Treat + (1 + Visit|id), data = SimTrial_sm_2000_1_5, REML = T)
summary(fit_lmer)

fit_lme <- lme(fixed = MADRS10 ~ Visit + Visit : Treat, 
              random = ~ 1 + Visit| id,
              method ="REML", 
               data = SimTrial_sm_2000_1_5,
              na.action = na.omit)
#View(SimTrial_sm_2000_1_5)

summary(fit_lme)

vcov(fit_lme)
VarCorr(fit_lme)
random.effects(fit_lme)[,1]
var(random.effects(fit_lme))

#View(SimTrial_sm_2000_1_5)
# Fit the Survival model: time to LoE
# but first, bring the dataset in the right format, namely all LoE at week 3 and all at week 2.
SimTrial_sm_2000_1_5_surv <- SimTrial_sm_2000_1_5
#View(SimTrial_sm_2000_1_5_surv)
SimTrial_sm_2000_1_5_surv$LoE_yes_surv <- ifelse(SimTrial_sm_2000_1_5_surv$LoE_YES==1 & SimTrial_sm_2000_1_5_surv$Visit==3, 1, 0)
#View(SimTrial_sm_2000_1_5_surv)
SimTrial_sm_2000_1_5_surv$AE_yes_surv <- ifelse(SimTrial_sm_2000_1_5_surv$AE_Yes==1 & SimTrial_sm_2000_1_5_surv$Visit==2, 1, 0)
SimTrial_sm_2000_1_5_surv_cox_LoE <- SimTrial_sm_2000_1_5_surv[SimTrial_sm_2000_1_5_surv$Visit==3, ]
SimTrial_sm_2000_1_5_surv_cox_LoE$Visit <- ifelse(SimTrial_sm_2000_1_5_surv_cox_LoE$LoE_yes_surv==1, 3, 6)
#View(SimTrial_sm_2000_1_5_surv_cox_LoE)


SimTrial_sm_2000_1_5_surv_cox_AE <- SimTrial_sm_2000_1_5_surv[SimTrial_sm_2000_1_5_surv$Visit==2, ]
SimTrial_sm_2000_1_5_surv_cox_AE$Visit <- ifelse(SimTrial_sm_2000_1_5_surv_cox_AE$AE_yes_surv==1, 2, 6)
#View(SimTrial_sm_2000_1_5_surv_cox_AE)

# now fit the survival models on the source trials (own source trial)
cox_fit_LoE <- coxph(Surv(Visit, LoE_yes_surv) ~ Treat, data = SimTrial_sm_2000_1_5_surv_cox_LoE, x = TRUE)
summary(cox_fit_LoE)



SimTrial_sm_2000_1_5_surv_cox_AE_exp <-SimTrial_sm_2000_1_5_surv_cox_AE[SimTrial_sm_2000_1_5_surv_cox_AE$Treat==1,]
SimTrial_sm_2000_1_5_surv_cox_AE_ctrl <-SimTrial_sm_2000_1_5_surv_cox_AE[SimTrial_sm_2000_1_5_surv_cox_AE$Treat==0,]


cox_fit_AE <- coxph(Surv(Visit, AE_yes_surv) ~ Treat, data = SimTrial_sm_2000_1_5_surv_cox_AE, x = TRUE)
summary(cox_fit_AE)




jointFit.LoE <- jointModel(fit_lme, cox_fit_LoE, timeVar = "Visit", 
                           method = "weibull-AFT-aGH")


summary(jointFit.LoE)




jointFit.AE <- jointModel(fit_lme, cox_fit_AE, timeVar = "Visit", 
                           method = "weibull-AFT-aGH")


summary(jointFit.AE)



    # logit(Pr(IE_ij)) = (beta_0 + b0_ij) [aka baseline MADRS] + beta_1 * Treat
    # baseline MADRS10 is made of the general intercept + each individual random intercept

d_re <- d_mis_w
    #View(d_re)


# add back the Baseline and CfB and CfW2 in the model specifications

fit_LoE_spm <- glm(LoE_Yes ~ Treat,
                   data = d_re, 
                   family = "binomial")

summary(fit_LoE_spm)
#Call:
#glm(formula = LoE_Yes ~ Treat, family = "binomial", data = d_re)

#Deviance Residuals: 
#  Min      1Q  Median      3Q     Max  
#-1.200  -0.956  -0.956   1.155   1.416  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.05233    0.06346   0.825     0.41    
#Treat       -0.59830    0.09114  -6.565 5.22e-11 ***

fit_AE_exp_spm <- glm(AE_Exp_Yes ~ 1,
                   data = d_re[d_re$Treat==1,], 
                   family = "binomial")

summary(fit_AE_exp_spm)

#Call:
#glm(formula = AE_Exp_Yes ~ 1, family = "binomial", data = d_re[d_re$Treat == 
#                                                                 1, ])

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-0.4527  -0.4527  -0.4527  -0.4527   2.1581  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -2.2263     0.1063  -20.94   <2e-16 ***


fit_AE_control_spm <- glm(AE_Control_Yes ~ 1,
                      data = d_re[d_re$Treat==0,], 
                      family = "binomial")

summary(fit_AE_control_spm)
#Call:
#glm(formula = AE_Control_Yes ~ 1, family = "binomial", data = d_re[d_re$Treat == 
#                                                                     0, ])

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-0.2517  -0.2517  -0.2517  -0.2517   2.6335  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -3.4361     0.1825  -18.83   <2e-16 ***

# These coefficients were used in the shared-parameter model (SPM) method

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


######### The code for SPM and glmm for logit models for LoE  at trial level, AE in experimental arm and AE in control arm

#View(SimTrial_sm_2000_1_5)
# estimate the model LoE
fit_glmm_LoE <- glmer(LoE_YES ~ Visit + Treat +  
                        (1 | id), data = SimTrial_sm_2000_1_5, family = binomial, control = glmerControl(optimizer = "bobyqa"),
                      nAGQ = 10)

summary(fit_glmm_LoE)




# estimate the model AE exp arm
fit_glmm_AE_exp <- glmer(AE_Exp_Yes ~ 1 +  
                        (1 | id), data = SimTrial_sm_2000_1_5[d$Treat==1,], family = binomial, control = glmerControl(optimizer = "bobyqa"),
                      nAGQ = 10)

summary(fit_glmm_AE_exp)


# estimate the model AE exp arm
fit_glmm_AE_control <- glmer(AE_Control_Yes ~ 1 +  
                           (1 | id), data = SimTrial_sm_2000_1_5[d$Treat==0,], family = binomial, control = glmerControl(optimizer = "bobyqa"),
                         nAGQ = 10)

summary(fit_glmm_AE_control)

# difficult to use due to overdispersion



#### the code for Survival models fitted on the large trial










# Selection model via marginal model----
# for outcomes-generating model and FULLY STOCHASTIC models for generation of intercurrent events
# for notation and description, please see above in the deterministic approach. The difference hereonwards is that we used logistic models to model the probability of intercurrent events, vs the deterministic rules as applied in the above approach

 # most of the code below is very similar with the code above, in terms of scope and operational reasons, as well as meaning and interpretation
 # code for longitudinal outcomes generation is as above
google_app <- httr::oauth_app(
  "appwhatever",
  key = "126364165263-nudc2q7h24voutu33a9i6pik9rjou09i.apps.googleusercontent.com",
  secret = "4ytmmz0zaj0709XBajk5Xyfb"
)
gm_auth_configure(key = "126364165263-nudc2q7h24voutu33a9i6pik9rjou09i.apps.googleusercontent.com",
                  secret = "4ytmmz0zaj0709XBajk5Xyfb")

options(mc.cores = parallel::detectCores()) # M1 with parallel, Old1 without parallel

Scenario <- c("A")

# simulation
set.seed(2147483629)
#set.seed(2147483399)
m.iterations <- 382# 382number of generated datasets # number of trials per scaling factor
scaling_factor <-  c(1)#c(0.5, 1.0, 1.5, 2.0, 2.5)
# total number of simulated trials = m.iterations * length(scaling_factor)
# try with c(0.4, 1.1, 1.8, 2.3, 3)
#c(0.25, 0.5, 1, 2, 2.5) steps for scaling factor
# to multiply for LoE, AE_control and AE_exp ### value 0.5 should yield around 13.5% IEs, 1 ~ %, 2 ~ %, 2.5 ~

n <- 190# number of patients at trial level to be randomised 1:1, experimental and control arm

#CFE <- matrix(ncol=4,nrow=length(scaling_factor)*m.iterations)
#colnames(CFE) <-c("N ceiled_floored", "% ceiled_floored", "scaling factor", "simulated trial n")

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline + follow-up measurements
delta <- matrix(ncol=1,nrow=m.iterations) # treatment effect estimate at 6 weeks based on MMRM models fitted on each generated dataset
colnames(delta) <-c("TreatmentEffect")
betas <- matrix(ncol=2,nrow=m.iterations)
colnames(betas) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3)

#confint_fit <- matrix(ncol=2,nrow=m.iterations) 
#colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
#delta_errorz <- matrix(ncol=1,nrow=m.iterations)
#colnames(delta_errorz) <- c("SE")

#bias_f <- matrix(ncol=1,nrow=m.iterations)
#colnames(bias_f) <- c("bias_f")

## randomisation objects, allocation
N_Exp  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Exp) <-c("N randomised Exp")

N_Control  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Control) <-c("N randomised Control")


# Intercurrent events objects
# LoE
n_LoE_total <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_LoE_total) <- c("N LoE Total")

n_LoE_Exp <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_LoE_Exp) <- c("N LoE Exp")

n_LoE_Control <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_LoE_Control) <- c("N LoE Control")


# AE
n_AE_total <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_AE_total) <- c("N AE Total")

n_AE_Exp <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_AE_Exp) <- c("N AE Exp")

n_AE_Control <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_AE_Control) <- c("N AE Control")

# AE + LoE Total

n_AE_and_LoE_T <- matrix(ncol=1,nrow=m.iterations) # 
colnames(n_AE_and_LoE_T) <- c("N AE and LoE Total")


# for percentages
# LoE
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

pb3 <- txtProgressBar(min = 0,  max=length(scaling_factor), style=3)

start_time <- Sys.time()

## Begin for loop----
for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    ### Generate longitudinal outcomes----
    ### Generate random effects
    re_means <- c(0, 0, 0, 0, 0, 0, 0)
    re_covm <- matrix(c(20.2190, 17.149, 14.721, 13.087,  8.4329,  10.854,   4.6417,
                        17.1490, 48.536, 41.161, 32.151, 24.8400,  30.528,  26.0170,
                        14.7210, 41.161, 72.569, 57.866, 60.2200,  61.974,  54.5400,
                        13.0870, 32.151, 57.866, 74.080, 66.2960,  63.540,  52.1070,
                        8.4329, 24.840, 60.220, 66.296, 97.4730,  90.612,  80.1370,
                        10.8540, 30.528, 61.974, 63.540, 90.6120, 116.410, 102.8300,
                        4.6417, 26.017, 54.540, 52.107, 80.1370, 102.830, 109.5900), nrow = 7)
    
       #Standard Deviations: 4.4965 6.9668 8.5188 8.607 9.8728 10.789 10.468
    
    re_covm2 <-matrix(c(0.00001, 0, 0, 0, 0, 0, 0,
                        0, 0.00001, 0, 0, 0, 0, 0,
                        0, 0, 0.00001, 0, 0, 0, 0,
                        0, 0, 0, 0.00001, 0, 0, 0,
                        0, 0, 0, 0, 0.00001, 0, 0,
                        0, 0, 0, 0, 0, 0.00001, 0,
                        0, 0, 0, 0, 0, 0, 0.00001), nrow = 7)
    
    re_covm3 <-re_covm/2
    
    re <- mvrnorm(n, re_means, re_covm3)	; re
      #View(re)
    
    re <- as.matrix(re)
    colnames(re) <-c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6") ; re
    
    d <- data.frame(
      id = rep(1:n, each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n, 1, 0.5), each = length(visits)),
      MADRS10 = rep(NA, n)); d # mean(Treat)
    
    d <- d[order(d$visit, d$id),]; #d
       #re
    
    j<-c(re[, 1], re[, 2], re[, 3], re[, 4], re[, 5], re[, 6], re[, 7])
    d$re <-j; #d
    
        #class(d)
        #View(d)
        #head(re)
    
    d <- d[order(d$id, d$visit),]; #d
    
        #re
    
    d<-as.matrix(d)
    
    # Scenario A
    beta.baseline <- 29.79
    beta_week1 <- -1
    beta_week2 <- -1.5
    beta_week3 <- -2
    beta_week4 <- -2.5 
    beta_week5 <- -3
    beta_week6 <- -3.5
    
    beta_v1_treatment <- -1
    beta_v2_treatment <- -1.5
    beta_v3_treatment <- -2
    beta_v4_treatment <- -2.5
    beta_v5_treatment <- -3
    
    beta_v6_treatment <- -3.5 # 
    
    treatmenteffect <-  beta_v6_treatment ; treatmenteffect
    
    # Y_ij = (Beta_0 + bi0) + (BetaWeek1 + bi1) + (BetaWeek2 + bi2) + (BetaWeek3 + bi3) + (BetaWeek4 + bi4) + (BetaWeek5 + bi5) + (BetaWeek6 + bi6) + 
    #                          Beta_W1_Treat * T + Beta_W2_Treat * T + Beta_W3_Treat * T + Beta_W4_Treat * T + Beta_W5_Treat * T + Beta_W6_Treat * T 
    
    
    for (i in 1:(n*length(visits))) {
      d[i,4] <- ifelse(d[i, 2]==0, beta.baseline + d[i,5],
                       ifelse(d[i, 2]==7, d[i-1,4] + beta_week1 + d[i,5] +  beta_v1_treatment * d[i, 3],
                              ifelse(d[i, 2]==14, d[i-2,4]+ beta_week2 + d[i,5] +  beta_v2_treatment * d[i, 3],
                                     ifelse(d[i, 2]==21, d[i-3,4] + beta_week3 + d[i,5] +  beta_v3_treatment * d[i, 3],
                                            ifelse(d[i, 2]==28, d[i-4,4] + beta_week4 + d[i,5] +  beta_v4_treatment * d[i, 3],
                                                   ifelse(d[i, 2]==35, d[i-5,4] + beta_week5 + d[i,5] +  beta_v5_treatment * d[i, 3],
                                                          d[i-6,4] + beta_week6 + d[i,5] +  beta_v6_treatment * d[i, 3]))))))
    }
    
    
       ## flooring and ceiling
       #d[, 4] <- ifelse(d[, 4] < 0, 0,
       #                ifelse(d[, 4]>60, 60, d[, 4]))
    
      #View(d)
    d <-as.data.frame(d)
    d$visit <-as.factor(d$visit)
    #d$Treat <- factor(d$Treat)

    # MMRM on full outcome data
        # this is the raw dataset used to check the model fit
    d <-d[,-5] # remove re (residuals) column from the dataset, they have been added to betas
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d # full outcome data
    
    length(d$id)
    tmp <- sapply(unique(d$id), FUN = function(i) nrow(d[d$id == i,]))
    BaselineMADRS10 <-  rep(d$MADRS10[d$visit == 0], tmp)
    length(BaselineMADRS10)
    d$Baseline <- BaselineMADRS10
    d
    d<-d[d$visit!=0,]
    
    range(d$MADRS10[d$Treat==1 & d$visit!=0])
    range(d$MADRS10[d$Treat==0 & d$visit!=0])
    
    range(d$MADRS10[d$Treat==1 & d$visit==42])
    range(d$MADRS10[d$Treat==0 & d$visit==42])
    
    mean(d$MADRS10[d$Treat==1 & d$visit==0])
    mean(d$MADRS10[d$Treat==0 & d$visit==0])
    
    mean(d$MADRS10[d$Treat==1 & d$visit==42])
    mean(d$MADRS10[d$Treat==0 & d$visit==42])
    
    
        #ceiling_floor <- sum(ifelse(d$MADRS10<0 | d$MADRS10>60, 1, 0))
        #ceiling_floor_perc <- (ceiling_floor/(n*length(visits)))*100
    
        #CFE[s*m,1] <- ceiling_floor
        #CFE[s*m,2] <- ceiling_floor_perc
        #CFE[s*m,3] <- s
        #CFE[s*m,4] <- m
    
    ### Plot trajectories
    p<- ggplot(data = d, aes(x = visit, y = MADRS10, group = id)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60))
  
    fit<-gls(MADRS10 ~ visit * Treat + Baseline, 
             data=d,
             correlation = corSymm(form=~1 | id),
             weights = varIdent(form = ~ 1 | visit), 
             method="REML")
    
    summary(fit)
    
    fit$coefficients[c(7,13)]

    sum(fit$coefficients[c(7,13)]); treatmenteffect
    
        #fit_lme <- lme(fixed=MADRS10 ~ visit * Treat + Baseline, 
        #           random=~1 | id,
        #          method="REML", 
        #         correlation = corSymm(form=~1|id),
        #        na.action=na.omit,
        #       data=d)
    
        #summary(fit_lme)
        #model_parameters(fit_lme)
    #m<-1
    ### Store estimated parameters
    betas[m, ] <- fit$coefficients[c(7,13)]
    
    delta[m, ] <- sum(fit$coefficients[c(7,13)])
    
    #bias_f[m, ] <- sum(fit$coefficients[c(7,13)]) - treatmenteffect
    
    #delta_error <- sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"]) 
    
    #delta_errorz[m, ] <- delta_error 
    
    
    #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_error
    #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_error
    
    Randomised_Exp <- sum(d[,3])/6 #number of patients in the experimental arm
    Randomised_Control <- n-sum(d[,3])/6 #number of patients in the control arm
    
    N_Exp[m,] <- Randomised_Exp
    N_Control[m,] <- Randomised_Control

    
    ## IEGM (intercurrent events generating model)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Generate Intercurrent events----
    
    # Use propensity models derived from original trial or from the dataset above 
    # One model for LoE for all patients (as derived above)
    # One model for AE in experimental arm (as derived above)
    # One model for AE in control arm (as derived above)
    
    ## take over the original/raw dataset to use for IEGM/MDGM
    d_mis <-d_orig
    
    d_mis_w <- d_mis %>% spread(visit, MADRS10) # reshape to wide in order to create the CfB variable
         #View(d_mis_w)
    
    colnames(d_mis_w)[3:9] <- c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6"); head(d_mis_w)
    
    d_mis_w$CfB <- d_mis_w[,"Baseline",] - d_mis_w[,"Week6"]; d_mis_w # create the CfB variable
    
    #View(d_mis_w)
    
    d_mis_w$CfW1 <- d_mis_w[,"Baseline"] - d_mis_w[,"Week1"]; d_mis_w # create the CfW1 variable
    
    d_mis_L <- d_mis_w %>% gather(Visit, MADRS10, Baseline:Week6) # reshape to long format
    
    d_mis_L <- d_mis_L[order(d_mis_L$id, d_mis_L$Visit),]; #d # order by subject id and Visit
    summary(logit_LoE)
    
    d_mis_L$predicted_LoE <- predict(logit_LoE, type="response", newdata=d_mis_L) # predicted occurrence of intercurrent events based on estimated probabilities from the logit models for e.g., treatment discontinuation due to lack of efficacy at trial level
        d#describe(predict(logit_LoE, type="response", newdata=d_mis_L)) # check 
    #View(d_mis_L)
    
    d_mis_L$LoE_yes <- rep(rbinom(length(unique(d_mis_L$id)), 1, unique(d_mis_L$predicted_LoE)), each = length(visits)) 
    d_mis_L$CfB
    #View(d_mis_L)
    
    trial_AE_X <- d_mis_L[d_mis_L$Treat==1,]
    
    trial_AE_X$predicted_AE <- predict(logit_AE_exp, type="response", newdata = trial_AE_X) #  predicted occurrence of intercurrent events based on estimated probabilities from the logit models for e.g., treatment discontinuation due to adverse events in the experimental arm
        #describe(predict(logit_LoE, type="response", newdata=d_mis_L))
    
    trial_AE_X$AE_yes <- rep(rbinom(length(unique(trial_AE_X$id)), 1, unique(trial_AE_X$predicted_AE)), each = length(visits)) 
    
    trial_AE_X$CfW1
    
    trial_AE_C <- d_mis_L[d_mis_L$Treat==0,]
    
    trial_AE_C$predicted_AE <- predict(logit_AE_control, type="response", newdata = trial_AE_C) # predicted occurrence of intercurrent events based on estimated probabilities from the logit models for e.g., treatment discontinuation due to adverse events in the control arm

    ### need to change here to not use a cutoff anymore, but directly rbinom with means equal to the estimated probabilities from the logit models
    trial_AE_C$AE_yes <- rep(rbinom(length(unique(trial_AE_C$id)), 1, unique(trial_AE_C$predicted_AE)), each = length(visits)) 
    
    trial_AE_C$CfW1
    
    
    # AE and LoE in experimental arm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #View(trial_AE_X)
    trial_AE_X$AE_YES <- ifelse(trial_AE_X[,"AE_yes"]==1, 1, 0)
    trial_AE_X$LoE_YES <- ifelse(trial_AE_X[,"AE_yes"]==0 & trial_AE_X[,"LoE_yes"]==1, 1, 0)
    
    #View(trial_AE_C)
    # AE and LoE in control arm #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~#
    trial_AE_C$AE_YES <- ifelse(trial_AE_C[,"AE_yes"]==1, 1, 0)
    trial_AE_C$LoE_YES <- ifelse(trial_AE_C[,"AE_yes"]==0 & trial_AE_C[,"LoE_yes"]==1, 1, 0)

  
    d_mis_LL <- rbind(trial_AE_X, trial_AE_C)
        #View(d_mis_LL)
    
    d_mis_LL
    
  
    d_mis_LL <- d_mis_LL[order(d_mis_LL$id),]; #d # order by subject id and Visit
       #colnames(d_mis_L)[10] <- c("MADRS10_mis"); d_mis_L # name the column
    
       # check classes of variable and transform in order 
    class(d_mis_LL$id)
    class(d_mis_LL$Treat)
    d_mis_LL$Treat <- as.factor(d_mis_LL$Treat)
    class(d_mis_LL$Visit)
    d_mis_LL$LoE_YES <- as.factor(d_mis_LL$LoE_YES)
    d_mis_LL$AE_YES <- as.factor(d_mis_LL$AE_YES)
    
    d_mis_LL$Behavior <- ifelse(d_mis_LL[,11]==1, "AE",
                               ifelse(d_mis_LL[,12]==1, "LoE", "No IE"))
    
    d_mis_LL$Behavior <-factor(d_mis_LL$Behavior, levels = c("AE", "LoE", "No IE"))
    
    class(d_mis_LL$Visit)
    d_mis_LL$Visit <- as.factor(d_mis_LL$Visit)
    rownames(d_mis_LL) <-NULL 
    
        # check range of values
        #range(d_mis_LL$MADRS10[d_mis_LL$Treat==1], na.rm = T)
        #range(d_mis_LL$MADRS10[d_mis_LL$Treat==0], na.rm = T)
    
    describe(d_mis_LL$AE_YES[d_mis_LL$Treat==0])
    describe(d_mis_LL$AE_YES[d_mis_LL$Treat==1])
    describe(d_mis_LL$AE_YES)
    
    describe(d_mis_LL$Behavior[d_mis_LL$Treat==0])
    describe(d_mis_LL$Behavior[d_mis_LL$Treat==1])
    describe(d_mis_LL$Behavior)
    
    describe(d_mis_LL$LoE_YES[d_mis_LL$Treat==0])
    describe(d_mis_LL$LoE_YES[d_mis_LL$Treat==1])
    describe(d_mis_LL$LoE_YES)
        #View(d_mis_LL)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Intercurrent events descriptives needed for the verification step ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    describe(d_mis_LL)
    LoE_Y <- d_mis_LL[ ,c("Treat", "LoE_YES")]
    head(LoE_Y)
    
    LoE_Y_total <- sum(as.numeric(LoE_Y$LoE_YES)-1)/length(visits)
    LoE_Y_Exp <- sum(as.numeric(LoE_Y$LoE_YES[LoE_Y$Treat==1])-1)/length(visits)
    LoE_Y_Control <- sum(as.numeric(LoE_Y$LoE_YES[LoE_Y$Treat==0])-1)/length(visits)
    
    ## LoE ##
    tb_LoE_total <- LoE_Y_total
    tb_LoE_Exp <- LoE_Y_Exp
    tb_LoE_Control <- LoE_Y_Control
    
    n_LoE_total[m, ] <-  tb_LoE_total
    n_LoE_Exp[m, ] <- tb_LoE_Exp
    n_LoE_Control[m, ] <- tb_LoE_Control
    
    AE_Y <- d_mis_LL[ ,c("Treat", "AE_YES")]
    head(AE_Y)
    
    
    AE_Y_total <- sum(as.numeric(AE_Y$AE_YES)-1)/length(visits)
    AE_Y_Exp <- sum(as.numeric(AE_Y$AE_YES[AE_Y$Treat==1])-1)/length(visits)
    AE_Y_Control <- sum(as.numeric(AE_Y$AE_YES[AE_Y$Treat==0])-1)/length(visits)
    
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
    
    
    
    
    ## Plot trajectories for all and each intercurrent event----
    # All patients with their trajectory 
    # LoE
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=LoE_YES)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    # AE
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=AE_YES)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    # All behaviors
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=Behavior)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  #  head(d_mis_LL)
    d_mis_LL_LoE <- d_mis_LL[d_mis_LL$LoE_YES==1,] # subset only patients that experienced LoE
    d_mis_LL_AE <- d_mis_LL[d_mis_LL$AE_YES==1,] # subset only patients that experienced AE
    d_mis_LL_NoIE <- d_mis_LL[d_mis_LL$LoE_YES==0 & d_mis_LL$AE_YES==0,] # subset only patients that did not experience any IE
    
        # View(d_mis_L)
    
    # All patients with their trajectory with different colours by LoE (Y/N)
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=LoE_YES))
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark red") + facet_wrap(~ Treat)
    
    # just LoE patients with their trajectory
    p<- ggplot(data = d_mis_LL_LoE, aes(x = Visit, y = MADRS10, group = id))
    #p + geom_line() + facet_grid(~ Treat) 
    plot_LoE_SMs <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") + facet_wrap(~ Treat)+ 
      ggtitle("SMs-LoE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_LoE_SMs
    
    #View(d_mis_L)
    
    
    # All patients with their trajectory with different colours by AE (Y/N)
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=AE_YES))
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark red") + facet_wrap(~ Treat)
    
    # just AE patients with their trajectory
    p<- ggplot(data = d_mis_LL_AE, aes(x = Visit, y = MADRS10, group = id))
    #p + geom_line() + facet_grid(~ Treat) 
    plot_AE_SMs <- p + geom_line(size=0.5, color='#F8766D') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + 
      facet_wrap(~ Treat)+ ggtitle("SPM-AE pattern")  +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_AE_SMs
    
    
    
    #just No IE
    # No IE
    p<- ggplot(data = d_mis_LL[d_mis_LL$Behavior=="No IE",], aes(x = Visit, y = MADRS10, group = id)) 
    plot_NoIE_SMs <- p + geom_line(size=0.5, color='#619CFF') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="blue") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("SMs-No IE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")) ; plot_NoIE_SMs
    
    
    # All behaviors
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=Behavior)) 
    plot_all_SMs <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="black") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("SMs-All patterns")+
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_all_SMs
    
    

    setTxtProgressBar(pb1, m)
  }
  

  # parameters extracted for MMRM fitted models on full outcome data
  colMeans(betas)
  colMeans(delta) ; treatmenteffect
  
  
  ### assign and save parameters---- 
  assign(paste('all_betas', s, sep="_"), betas)
  assign(paste('all_delta', s, sep="_"), delta)
  assign(paste('all_delta_errorz', s, sep="_"), delta_errorz)
  assign(paste('all_confint_fit', s, sep="_"), confint_fit)
  assign(paste('all_N_Exp', s, sep="_"), N_Exp)
  assign(paste('all_N_Control', s, sep="_"), N_Control)
  
  setTxtProgressBar(pb3, s)
}

## End for loop----
end_time <- Sys.time()
end_time-start_time

colMeans(all_betas_1)


tolerance_margin <- 0.1 

difference_Verification <- abs(treatmenteffect - colMeans(all_delta_1))

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_Verification) < tolerance_margin), "Verification SUCCESSFUL", "Verification NOT successful") 



min(all_betas_1)
max(all_betas_1)



# Table for the paper ----



table_AE_SMs <- data.frame(
  # descriptives AE  
  n_AE_Control,
  n_AE_Exp); table_AE_SMs

mean(n_AE_Control)
mean(n_AE_Exp)



# descriptives LoE  
table_LoE_SMs <-data.frame(
  n_LoE_Control,
  n_LoE_Exp); table_LoE_SMs

mean(n_LoE_Control)
mean(n_LoE_Exp)


#describe(table_IE_SMs)


table_AE_SMs %>% 
  as.data.frame() %>% 
  mutate("Intercurrent event" = "AE") %>% 
  rename(N_C_arm=N.AE.Control) %>% 
  rename(N_E_arm=N.AE.Exp)

table_LoE_SMs %>% 
  as.data.frame() %>% 
  mutate("Intercurrent event" = "LoE") %>% 
  rename(N_C_arm=N.LoE.Control) %>% 
  rename(N_E_arm=N.LoE.Exp)


tab_SMs <- tibble(bind_rows(table_AE_SMs %>% 
                               as.data.frame() %>% 
                               mutate("Intercurrent event" = "AE") %>% 
                               rename(N_C_arm=N.AE.Control) %>% 
                               rename(N_E_arm=N.AE.Exp), 
                             table_LoE_SMs %>% 
                               as.data.frame() %>% 
                               mutate("Intercurrent event" = "LoE") %>% 
                               rename(N_C_arm=N.LoE.Control) %>% 
                               rename(N_E_arm=N.LoE.Exp))); tab_SMs



tab2_SMs <- tab_SMs %>% group_by(`Intercurrent event`) %>%
  summarise("N" = round(mean(N_C_arm), digits=1), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = round(mean(N_E_arm), digits=1), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = round(mean(N_C_arm + N_E_arm), digits=1),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) %>% 
  adorn_totals("row"); tab2_SMs




gt(tab2_SMs) %>% 
  tab_header(title = md("Table 4. Descriptive statistics intercurrent events"), subtitle = md("Selection model DGM - stochastic")) %>%
  tab_source_note(md(paste0("Averaged over", " ", m.iterations*length(scaling_factor),  " ",  "simulated trials.", " ", "Trial sample size = ", " ", n ))) %>% 
  tab_spanner(
    label = md("**Control**"),
    columns = c("N", "%")) %>% 
  cols_align(
    align = "center",
    columns =  c("N", "%")
  ) %>% 
  tab_spanner(
    label = md("**Treatment**"),
    columns = c("N ", "% ")) %>% 
  cols_align(
    align = "center",
    columns =  c("N ", "% ")
  ) %>% 
  tab_spanner(
    label = md("**Total**"),
    columns = c(" N ", " % ")) %>% 
  cols_align(
    align = "center",
    columns =  c(" N ", " % ")
  ) %>% 
  data_color(
    columns = c("%", "% ", " % "),
    colors = scales::col_numeric(
      palette = c(
        "light blue"),
      domain = NULL)
  ) %>% 
  cols_align(
    align = "center",
    columns =  "Intercurrent event"
  ) %>% 
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
installed.packages()

