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


#install.packages("remotes")
#library(remotes)
#install_version("MASS", "7.3.54")

#install.packages("devtools")
#library(devtools)

#devtools::install_version("tidyverse", version = "1.3.0") 
#devtools::install_version("nlme", version = "3.1.50") 





#install_version("lme4", "1.1.27")
#install_version("Hmisc", "4.5.0")
#install_version("janitor", "2.0.1")
#install_version("gt", "0.3.0")
#install_version("patchwork", "1.1.0")


## load libraries ----
#rm(list=ls()) #
library(gmailr)
library(MASS)#
library(tidyverse)#
library(nlme)#
library(lme4)#
#library(Hmisc)
library(janitor)#
library(gt)#
library(patchwork)#
library(JM)

## gmail setup----
# Selection model via marginal model for outcomes-generating model and deterministic/stochastic rules for generation of intercurrent events
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


#set.seed(2147483399)

n <- 190#8000 # number of patients to be simulated (sample size)
# this is based on a t-test to ensure  90% power at alpha level=0.025 one-sided 

m.iterations <- 500 # 382 is the number of trials needed for the verification of the longitudinal outcomes # number of generated datasets # number of trials per scaling factor
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
delta_SMd <- matrix(ncol=1,nrow=m.iterations) # object to store treatment effect estimates at 6 weeks based on MMRM model fitted on each generated dataset
colnames(delta_SMd) <-c("TreatmentEffect")
betas_SMd <- matrix(ncol=2,nrow=m.iterations) # object to store parameters for the treatment effect at week 6 based on the MMRM model fitted on each generated dataset
colnames(betas_SMd) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3) # progress bar in percentages relative to the total number of m.iterations

confint_fit <- matrix(ncol=2,nrow=m.iterations) # object to store the 95% confidence interval bounds for the estimated treatment effect
colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
delta_SMd_errorz <- matrix(ncol=1,nrow=m.iterations) # standard error of the estimated treatment effect
colnames(delta_SMd_errorz) <- c("SE")

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
set.seed(2147483629) # set seed
## Begin for loop----
for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    ### Generate longitudinal outcomes----
    #### Generate correlated residuals----
    re_means <- c(0, 0, 0, 0, 0, 0, 0) # means of residuals at each timepoint
    
    # covariance matrix extracted from trial 003-002 (see reference to the data analysis paper in PST)
    # this dataset is not available to share. A MMRM model was fitted and this is the corresponding covariance matrix
    re_covm <- matrix(c(20.2190, 17.149, 14.721, 13.087,  8.4329,  10.854,   4.6417, 
                        17.1490, 48.536, 41.161, 32.151, 24.8400,  30.528,  26.0170,
                        14.7210, 41.161, 72.569, 57.866, 60.2200,  61.974,  54.5400,
                        13.0870, 32.151, 57.866, 74.080, 66.2960,  63.540,  52.1070,
                        8.4329, 24.840, 60.220, 66.296, 97.4730,  90.612,  80.1370,
                        10.8540, 30.528, 61.974, 63.540, 90.6120, 116.410, 102.8300,
                        4.6417, 26.017, 54.540, 52.107, 80.1370, 102.830, 109.5900), nrow = 7)
    
    #Standard Deviations: 4.4965 6.9668 8.5188 8.607 9.8728 10.789 10.468
    
    # covariance matrix with 0 off-diagonal and small variances. This is useful for initial/later checks to see if the simulated data corresponds to target data to be simulated
    size_diagonal <- 0.0000001
    re_covm2 <-matrix(c(size_diagonal, 0, 0, 0, 0, 0, 0,
                        0, size_diagonal, 0, 0, 0, 0, 0,
                        0, 0, size_diagonal, 0, 0, 0, 0,
                        0, 0, 0, size_diagonal, 0, 0, 0,
                        0, 0, 0, 0, size_diagonal, 0, 0,
                        0, 0, 0, 0, 0, size_diagonal, 0,
                        0, 0, 0, 0, 0, 0, size_diagonal), nrow = 7)
    
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
    d <-d[,-5] # remove re (residuals) column from the dataset, they have been added to betas_SMd
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d # full outcome data
    
    # create a separate column with only the baseline outcomes, if the baseline values will be used as covariate in the model 
    length(d$id)
    tmp <- sapply(unique(d$id), FUN = function(i) nrow(d[d$id == i,]))
    BaselineMADRS10 <-  rep(d$MADRS10[d$visit == 0], tmp)
    length(BaselineMADRS10)
    d$Baseline <- BaselineMADRS10
        #View(d)
    
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
    fit <- gls(MADRS10 ~ visit*Treat, 
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~#
    #### The code for SPM DGM
    #~~~~~~~~~~~~~~~~~~~~~~~~~#
    # On a trial generated with the selection model (in this for loop with 8000 patients and 1 iteration), we fitted a mixed-effects model with random intercept and random slope # 
    # fit LMM model on the simulated trial (n=8000, re_covm3), to get good estimates for the outcomes model
    # fit the glm models to get the intercurrent events models. On these, random effects will be added in the linear predictor
    
        d$visit <- as.numeric(d$visit)-1
    
        fit_lme <- lme(fixed=MADRS10 ~ visit + visit:Treat, 
                  random=~1 + visit | id,
                 method="REML", 
              #correlation = corSymm(form=~1|id),
            data=d)
        summary(fit_lme)
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
        vcov(fit_lme)
        VarCorr(fit_lme)
        var(random.effects(fit_lme))
        #            (Intercept)    visit
        #(Intercept)   21.140783 1.307053
        #visit          1.307053 1.044141
        
        
        #fit_lmer <- lmer(MADRS10 ~visit + visit:Treat + (1 + visit |id), data = d, REML = T)
    
        #summary(fit_lmer)
      
        #model_parameters(fit_lme)
    #### store estimated parameters----
        summary(fit)
    betas_SMd[m, ] <- fit$coefficients[c(8,14)] # store the parameters corresponding to the treatment effect at the end of the trial, at week 6
    
    delta_SMd[m, ] <- sum(fit$coefficients[c(8,14)]) # store the treatment effect at the end of the trial, at week 6
    
        #bias_f[m, ] <- sum(fit$coefficients[c(7,13)]) - treatmenteffect
        #delta_SMd_error <- sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"]) 
        #delta_SMd_errorz[m, ] <- delta_SMd_error 
        #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_SMd_error
        #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_SMd_error
    
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
    
    d_mis_w <- d_mis |> spread(visit, MADRS10) # reshape to wide in order to create the CfB variable
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
    
   #describe(ifelse(d_mis_w$Treat==0 & d_mis_w$CfW1< (-2), 1, 0))
    
    #
    
    d_mis_L <- d_mis_w |> gather(Visit, MADRS10, Baseline:Week6) # reshape to long format
    
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
    #describe(d_mis_L$LoE_YES)
    
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
    plot_all_SMd <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="black") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SMd-All patterns") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_all_SMd
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##### Plots by subsets/patterns of patients that experienced a certain intercurrent event----
    
    d_mis_L_LoE <- d_mis_L[d_mis_L$LoE_YES==1,] # subset only patients that experienced LoE
    d_mis_L_AE <- d_mis_L[d_mis_L$AE_Yes==1,] # subset only patients that experienced AE
    d_mis_L_NoIE <- d_mis_L[d_mis_L$LoE_YES==0 & d_mis_L$AE_Yes==0,] # subset only patients that did not experience any IE
    
    # just (subset) LoE patients with their trajectory #
    p<- ggplot(data = d_mis_L_LoE, aes(x = Visit, y = MADRS10, group = id))
    plot_LoE_SMd <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SMd-LoE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_LoE_SMd

    
    # just AE patients with their trajectory #
    p<- ggplot(data = d_mis_L_AE, aes(x = Visit, y = MADRS10, group = id))
    plot_AE_SMd <- p + geom_line(size=0.5, color='#F8766D') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SMd-AE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_AE_SMd
    
    # just No IE patients with their trajectory # 
    p<- ggplot(data = d_mis_L_NoIE, aes(x = Visit, y = MADRS10, group = id))
    plot_NoIE_SMd <- p + geom_line(size=0.5, color='#619CFF') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="blue") + facet_wrap(~ Treat) +
     scale_y_continuous(limits = c(-10, 60)) + ggtitle("SMd-No IE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_NoIE_SMd
    
    
    
    
    setTxtProgressBar(pb1, m)
  }

  # parameters extracted for MMRM fitted models on full outcome data
  colMeans(betas_SMd) # average of treatment effect parameters estimated from the model
  colMeans(delta_SMd) ; treatmenteffect # average treatment effect estimated from the model
  
  # assign and save all parameters----
  # pertaining to the treatment effect, the estimated treatment effect, the standard errors, 95% CI and intercurrent event descriptives
  assign(paste('all_betas_SMd', s, sep="_"), betas_SMd)
  
  assign(paste('all_delta_SMd', s, sep="_"), delta_SMd)
  
  #assign(paste('all_delta_SMd_errorz', s, sep="_"), delta_SMd_errorz)
  
  #assign(paste('all_confint_fit', s, sep="_"), confint_fit)
  
  #assign(paste('all_N_Exp', s, sep="_"), N_Exp)
  
  #assign(paste('all_N_Control', s, sep="_"), N_Control)
  

  setTxtProgressBar(pb3, s)
}

# end for loop----
end_time <- Sys.time() # timestamp for end of simulation
end_time-start_time # total time to complete the simulation



betas_SMd; 
colMeans(delta_SMd); treatmenteffect# change the name of delta_SMd_SMd

tolerance_margin <- 0.1 
difference_Verification_SMd <- abs(treatmenteffect - colMeans(delta_SMd))# check parameterisation gls()
#rename each to avoid duplication of verification conclusion due to the same name.

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_Verification_SMd) < tolerance_margin), "Verification SMd *SUCCESSFUL*", "Verification SMd NOT successful :(") 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot relevant graphs for the paper----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

the_plot_SMd <- (plot_all_SMd / plot_LoE_SMd) | (plot_AE_SMd / plot_NoIE_SMd); the_plot_SMd






# Table for the paper ----

table_AE_SMd <- data.frame(
# descriptives AE  
n_AE_Control,
n_AE_Exp); table_AE_SMd

mean(n_AE_Control)
mean(n_AE_Exp)




# descriptives LoE  
table_LoE_SMd <-data.frame(
n_LoE_Control,
n_LoE_Exp); table_LoE_SMd

mean(n_LoE_Control)
mean(n_LoE_Exp)


#describe(table_IE_SM)


table_AE_SMd |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "AE") |> 
  rename(N_C_arm=N.AE.Control) |> 
  rename(N_E_arm=N.AE.Exp)

table_LoE_SMd |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "LoE") |> 
  rename(N_C_arm=N.LoE.Control) |> 
  rename(N_E_arm=N.LoE.Exp)


tab_SMd <- tibble(bind_rows(table_AE_SMd |> 
            as.data.frame() |> 
            mutate("Intercurrent event" = "AE") |> 
            rename(N_C_arm=N.AE.Control) |> 
            rename(N_E_arm=N.AE.Exp), 
          table_LoE_SMd |> 
            as.data.frame() |> 
            mutate("Intercurrent event" = "LoE") |> 
            rename(N_C_arm=N.LoE.Control) |> 
            rename(N_E_arm=N.LoE.Exp))); tab_SMd

  

tab2_SMd<- tab_SMd |> group_by(`Intercurrent event`) |>
  summarise("N" = round(mean(N_C_arm), digits=1), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = round(mean(N_E_arm), digits=1), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = round(mean(N_C_arm + N_E_arm), digits=1),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) |> 
  adorn_totals("row"); tab2_SMd



  
gt(tab2_SMd) |> 
  tab_header(title = md("Table 8a. Descriptive statistics intercurrent events"), subtitle = md("Selection model DGM - deterministic rule")) |>
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

# 5 February 2022
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
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets 
#[6] methods   base     

#other attached packages:
#  [1] patchwork_1.1.1 gt_0.3.1        janitor_2.1.0  
#[4] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
#[7] purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
#[10] tibble_3.1.6    ggplot2_3.3.5   tidyverse_1.3.1
#[13] foreign_0.8-82  survival_3.2-13 lme4_1.1-28    
#[16] Matrix_1.4-0    nlme_3.1-155    MASS_7.3-55    
#[19] gmailr_1.0.1   

#installed.packages()






# Determine the number of trials needed to simulate for the verification of the longitudinal outcomes
# Formula from Burton paper
#tolerance_margin <- 0.1 # bias allowed
#std.e <- 0.997 # model-based standard error of the treatment effect estimate from a fitted model on 1 trial

#n.trials_needed <- ceiling(((qnorm(0.975) * std.e)/tolerance_margin)^2) ; n.trials_needed # for the verification 
# 382 trials
# verification of the longitudinal outcomes was successful


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The code used to extract the MMRM (longitudinal) models for each pattern in preparation for PMMMM method----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Source trial with 1000 patients per arm (total sample size = 2000 patients, re_covm3)
# Pattern for LoE at trial level (as per rule above)
# Pattern for AE in experimental and control arms (as per rules above in SMd)
# Pattern for No intercurrent event at trial level (remaining patients)


    #View(SimTrial_sm_2000_1_5)
SimTrial_sm_2000_1_5 <-SimTrial_sm_2000_1_1 #
#This is 2000 for most of the models, for SPM DGM for the linear mixed effects model n=8000 (see above in the for loop)
# to avoid duplicating the code, which is already lengthy
#describe(SimTrial_sm_2000_1_5)

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
# you should get
#> summary(fit_LoE_trial)
#Generalized least squares fit by REML
#Model: MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 + Treat:V7 + Treat:V14 +      Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42 
#Data: SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$LoE_YES == 1, ] 
#AIC      BIC    logLik
#26597.48 26824.18 -13263.74

#Correlation Structure: General
#Formula: ~1 | id 
#Parameter estimate(s):
#  Correlation: 
#  1     2     3     4     5     6    
#2 0.669                              
#3 0.499 0.765                        
#4 0.463 0.642 0.772                  
#5 0.440 0.497 0.696 0.772            
#6 0.540 0.533 0.632 0.665 0.785      
#7 0.641 0.541 0.586 0.574 0.688 0.850

#Coefficients:
#  Value Std.Error   t-value p-value
#(Intercept) 29.877763 0.2226789 134.17420  0.0000
#V7          -1.266730 0.2260979  -5.60257  0.0000
#V14         -0.387887 0.2731455  -1.42008  0.1557
#V21         -0.835799 0.2816776  -2.96722  0.0030
#V28          0.507517 0.2870108   1.76829  0.0771
#V35          0.879907 0.2630121   3.34550  0.0008
#V42          0.952685 0.2349301   4.05518  0.0001
#V7:Treat1    2.190204 0.3375204   6.48910  0.0000
#V14:Treat1   1.863464 0.3936713   4.73355  0.0000
#V21:Treat1   0.924152 0.4026746   2.29503  0.0218
#V28:Treat1   0.139989 0.4080712   0.34305  0.7316
#V35:Treat1   0.685633 0.3824412   1.79278  0.0731
#V42:Treat1  -0.337289 0.3487655  -0.96710  0.3335

#Correlation: 
#  (Intr) V7     V14    V21    V28    V35   
#V7         -0.326                                   
#V14        -0.408  0.711                            
#V21        -0.424  0.569  0.756                     
#V28        -0.435  0.401  0.680  0.767              
#V35        -0.390  0.366  0.577  0.629  0.770       
#V42        -0.340  0.285  0.482  0.491  0.647  0.807
#V7:Treat1   0.000 -0.599 -0.387 -0.289 -0.173 -0.160
#V14:Treat1  0.000 -0.401 -0.578 -0.404 -0.349 -0.290
#V21:Treat1  0.000 -0.302 -0.407 -0.573 -0.407 -0.324
#V28:Treat1  0.000 -0.182 -0.354 -0.410 -0.570 -0.423
#V35:Treat1  0.000 -0.165 -0.287 -0.319 -0.413 -0.583
#V42:Treat1  0.000 -0.118 -0.231 -0.234 -0.336 -0.455
#V42    V7:Tr1 V14:T1 V21:T1 V28:T1 V35:T1
#V7                                                  
#V14                                                 
#V21                                                 
#V28                                                 
#V35                                                 
#V42                                                 
#V7:Treat1  -0.117                                   
#V14:Treat1 -0.238  0.669                            
#V21:Treat1 -0.243  0.503  0.704                     
#V28:Treat1 -0.351  0.304  0.612  0.714              
#V35:Treat1 -0.464  0.275  0.497  0.556  0.725       
#V42:Treat1 -0.596  0.196  0.399  0.407  0.589  0.779

#Standardized residuals:
#  Min          Q1         Med          Q3 
#-3.83413438 -0.66126098 -0.04652564  0.58405776 
#Max 
#5.58338731 

#Residual standard error: 5.840813 
#Degrees of freedom: 4816 total; 4803 residual



## Pattern for AE in experimental arm----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fit_AE_all<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                  Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                    data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Yes==1,],
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | Visit),
                    method="REML")


summary(fit_AE_all)
# you should get:
#> summary(fit_AE_all)
#Generalized least squares fit by REML
#Model: MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 + Treat:V7 + Treat:V14 +      Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42 
#Data: SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Yes == 1, ] 
#AIC      BIC    logLik
#11498.99 11731.46 -5708.495

#Correlation Structure: General
#Formula: ~1 | id 
#Parameter estimate(s):
#  Correlation: 
#  1     2     3     4     5     6    
#2 0.864                              
#3 0.679 0.759                        
#4 0.594 0.641 0.833                  
#5 0.467 0.493 0.768 0.829            
#6 0.467 0.498 0.730 0.754 0.872      
#7 0.358 0.415 0.675 0.665 0.801 0.923
#Variance function:
#  Structure: Different standard deviations per stratum
#Formula: ~1 | Visit 
#Parameter estimates:
#  Baseline    Week1    Week2    Week3    Week4    Week5 
#1.000000 1.527025 1.972734 2.114375 2.381852 2.612523 
#Week6 
#2.408459 

#Coefficients:
#  Value Std.Error   t-value p-value
#(Intercept) 29.580079 0.1825747 162.01628  0.0000
#V7           4.028240 0.2472695  16.29089  0.0000
#V14          2.739713 0.4570672   5.99411  0.0000
#V21          1.577882 0.5340134   2.95476  0.0032
#V28         -0.384598 0.6589913  -0.58362  0.5595
#V35         -0.307635 0.7234364  -0.42524  0.6707
#V42         -0.352707 0.7034991  -0.50136  0.6162
#V7:Treat1   -9.843380 0.2959763 -33.25733  0.0000
#V14:Treat1  -8.976663 0.5577959 -16.09310  0.0000
#V21:Treat1  -8.198202 0.6552633 -12.51131  0.0000
#V28:Treat1  -6.771080 0.8113319  -8.34564  0.0000
#V35:Treat1  -7.991637 0.8897234  -8.98216  0.0000
#V42:Treat1  -8.795151 0.8659999 -10.15606  0.0000

#Correlation: 
#  (Intr) V7     V14    V21    V28    V35    V42   
#V7          0.236                                          
#V14         0.136  0.482                                   
#V21         0.087  0.326  0.730                            
#V28         0.031  0.203  0.693  0.774                     
#V35         0.056  0.219  0.637  0.672  0.837              
#V42        -0.035  0.208  0.619  0.597  0.766  0.912       
#V7:Treat1   0.000 -0.789 -0.376 -0.255 -0.163 -0.172 -0.181
#V14:Treat1  0.000 -0.368 -0.804 -0.588 -0.564 -0.516 -0.511
#V21:Treat1  0.000 -0.249 -0.585 -0.809 -0.629 -0.543 -0.489
#V28:Treat1  0.000 -0.159 -0.559 -0.627 -0.811 -0.678 -0.623
#V35:Treat1  0.000 -0.168 -0.512 -0.542 -0.679 -0.811 -0.743
#V42:Treat1  0.000 -0.176 -0.506 -0.487 -0.623 -0.742 -0.811
#V7:Tr1 V14:T1 V21:T1 V28:T1 V35:T1
#V7                                           
#V14                                          
#V21                                          
#V28                                          
#V35                                          
#V42                                          
#V7:Treat1                                    
#V14:Treat1  0.467                            
#V21:Treat1  0.315  0.728                     
#V28:Treat1  0.201  0.695  0.775              
#V35:Treat1  0.212  0.637  0.670  0.837       
#V42:Treat1  0.223  0.630  0.602  0.768  0.916

#Standardized residuals:
#  Min          Q1         Med          Q3         Max 
#-3.70652558 -0.69297245  0.00234884  0.67972462  3.52416401 

#Residual standard error: 3.204173 
#Degrees of freedom: 2156 total; 2143 residual

getVarCov(fit_AE_all, individual = '6')
#you should get:
#> getVarCov(fit_AE_all, individual = '6')
#Marginal variance covariance matrix
#[,1]   [,2]   [,3]   [,4]   [,5]   [,6]    [,7]
#[1,] 10.2670 13.548 13.753 12.890 11.415 12.529  8.8644
#[2,] 13.5480 23.940 23.483 21.243 18.407 20.406 15.6570
#[3,] 13.7530 23.483 39.955 35.670 37.071 38.650 32.9240
#[4,] 12.8900 21.243 35.670 45.898 42.849 42.779 34.7890
#[5,] 11.4150 18.407 37.071 42.849 58.245 55.731 47.1980
#[6,] 12.5290 20.406 38.650 42.779 55.731 70.073 59.6430
#[7,]  8.8644 15.657 32.924 34.789 47.198 59.643 59.5540
#Standard Deviations: 3.2042 4.8929 6.321 6.7748 7.6319 8.371 7.7171 




## Pattern for no intercurrent events at trial level----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fit_no_IE<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                 Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
               data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Behavior== "No IE", ] ,
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | Visit),
                    method="REML")

summary(fit_no_IE)
#you should get:
#> summary(fit_no_IE)
#Generalized least squares fit by REML
#Model: MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 + Treat:V7 + Treat:V14 +      Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42 
#Data: SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Behavior == "No IE",      ] 
#AIC      BIC    logLik
#37858.53 38139.61 -18888.26

#Correlation Structure: General
#Formula: ~1 | id 
#Parameter estimate(s):
#  Correlation: 
#  1     2     3     4     5     6    
#2 0.822                              
#3 0.707 0.820                        
#4 0.684 0.728 0.846                  
#5 0.559 0.582 0.784 0.826            
#6 0.565 0.598 0.750 0.756 0.877      
#7 0.498 0.555 0.705 0.678 0.812 0.927
#Variance function:
#  Structure: Different standard deviations per stratum
#Formula: ~1 | Visit 
#Parameter estimates:
#  Baseline    Week1    Week2    Week3    Week4    Week5    Week6 
#1.000000 2.180655 2.312579 2.300759 2.388043 2.577183 2.335500 

#Coefficients:
#  Value Std.Error   t-value p-value
#(Intercept) 29.514812 0.0983291 300.16353  0.0000
#V7          -2.529964 0.1923158 -13.15525  0.0000
#V14         -4.052008 0.2398648 -16.89288  0.0000
#V21         -4.431931 0.2442224 -18.14711  0.0000
#V28         -6.018855 0.2823669 -21.31572  0.0000
#V35         -7.601607 0.3044194 -24.97084  0.0000
#V42         -8.457351 0.2871898 -29.44865  0.0000
#V7:Treat1    0.397323 0.2442874   1.62646  0.1039
#V14:Treat1  -0.309205 0.3217999  -0.96086  0.3367
#V21:Treat1  -0.902596 0.3301530  -2.73387  0.0063
#V28:Treat1  -1.183286 0.3896808  -3.03655  0.0024
#V35:Treat1  -1.083183 0.4183702  -2.58905  0.0096
#V42:Treat1  -1.424697 0.3984256  -3.57582  0.0004

#Correlation: 
#  (Intr) V7     V14    V21    V28    V35    V42    V7:Tr1
#V7          0.405                                                 
#V14         0.260  0.628                                          
#V21         0.231  0.449  0.721                                   
#V28         0.117  0.282  0.666  0.736                            
#V35         0.148  0.316  0.612  0.625  0.824                     
#V42         0.056  0.291  0.570  0.531  0.742  0.899              
#V7:Treat1   0.000 -0.658 -0.411 -0.280 -0.185 -0.202 -0.211       
#V14:Treat1  0.000 -0.389 -0.695 -0.492 -0.474 -0.428 -0.414  0.592
#V21:Treat1  0.000 -0.263 -0.489 -0.700 -0.524 -0.437 -0.383  0.399
#V28:Treat1  0.000 -0.170 -0.460 -0.514 -0.715 -0.584 -0.533  0.259
#V35:Treat1  0.000 -0.186 -0.418 -0.430 -0.587 -0.712 -0.648  0.283
#V42:Treat1  0.000 -0.193 -0.400 -0.373 -0.530 -0.642 -0.719  0.294
#V14:T1 V21:T1 V28:T1 V35:T1
#V7                                    
#V14                                   
#V21                                   
#V28                                   
#V35                                   
#V42                                   
#V7:Treat1                             
#V14:Treat1                            
#V21:Treat1  0.703                     
#V28:Treat1  0.663  0.734              
#V35:Treat1  0.601  0.614  0.821       
#V42:Treat1  0.576  0.533  0.741  0.902

#Standardized residuals:
#  Min          Q1         Med          Q3         Max 
#-4.10416049 -0.65565895 -0.02936365  0.63449169  4.18037828 

#Residual standard error: 3.115652 
#Degrees of freedom: 7028 total; 7015 residual

getVarCov(fit_no_IE, individual = 3)
#you should get:
#> getVarCov(fit_no_IE, individual = 3)
#Marginal variance covariance matrix
#[,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]
#[1,]  9.7073 17.403 15.873 15.283 12.956 14.140 11.301
#[2,] 17.4030 46.161 40.123 35.473 29.397 32.609 27.429
#[3,] 15.8730 40.123 51.915 43.718 42.015 43.401 36.981
#[4,] 15.2830 35.473 43.718 51.385 44.059 43.510 35.372
#[5,] 12.9560 29.397 42.015 44.059 55.358 52.421 43.936
#[6,] 14.1400 32.609 43.401 43.510 52.421 64.475 54.154
#[7,] 11.3010 27.429 36.981 35.372 43.936 54.154 52.949
#Standard Deviations: 3.1157 6.7942 7.2052 7.1684 7.4403 8.0296 7.2766 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOR SM Stochastic implementation
#  The code used to extract the logit models for LoE at trial level, AE in experimental arm and AE in control arm----
## obtain the logistic regression models for LoE at trial level, AE in experimental arm and AE in control arm ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# n= 2000 patients, re_covm3

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
#you should get
#> summary(logit_LoE)

#Call:
#  glm(formula = LoE_Yes ~ CfB, family = "binomial", data = SimTrial_sm_2000_1_5)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-3.9859  -0.5398  -0.1679   0.6016   1.4843  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#  (Intercept)  1.077405   0.033011   32.64   <2e-16 ***
#  CfB         -0.355617   0.006021  -59.07   <2e-16 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 18839  on 13999  degrees of freedom
#Residual deviance: 10473  on 13998  degrees of freedom
#AIC: 10477

#Number of Fisher Scoring iterations: 6


#View(SimTrial_sm_2000_1_5)




predicted_LoE <- predict(logit_LoE, type="response", newdata=SimTrial_sm_2000_1_5)

SimTrial_sm_2000_1_5$predicted_LoE <- predicted_LoE
    #View(SimTrial_sm_2000_1_5[,c(5, 20)])


#probz_predicted_LoE_Yes <- SimTrial_sm_2000_1_5$predicted_LoE[SimTrial_sm_2000_1_5$LoE_Yes==1]
#up_boundary_prob_LoE <- max(probz_predicted_LoE_Yes)  ; up_boundary_prob_LoE
#low_boundary_prob_LoE <- min(probz_predicted_LoE_Yes) ; low_boundary_prob_LoE

## AE experimental arm----
#logit(Pr(AE_exp))= fi_AE_exp * (Y_i0 - Y_i2)
logit_AE_exp <- glm(AE_Exp_Yes ~ CfW1,
                 data = trial_AE_exp,
                 family = "binomial")

summary(logit_AE_exp)
#you should get
#> summary(logit_AE_exp)

#Call:
#  glm(formula = AE_Exp_Yes ~ CfW1, family = "binomial", data = trial_AE_exp)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.9455  -0.6155  -0.3926  -0.1830   1.8673  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -2.370965   0.053292  -44.49   <2e-16 ***
#  CfW1         0.272188   0.008634   31.53   <2e-16 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 7060.3  on 6992  degrees of freedom
#Residual deviance: 5654.8  on 6991  degrees of freedom
#AIC: 5658.8

#Number of Fisher Scoring iterations: 5


#predicted_AE_exp <- predict(logit_AE_exp, type="response", newdata = trial_AE_exp)

#trial_AE_exp$predicted_AE_exp <- predicted_AE_exp

    #View(trial_AE_exp[,c(6, 20)])
    #View(trial_AE_exp)

#probz_predicted_AE_exp_Yes <- trial_AE_exp$predicted_AE_exp[trial_AE_exp$AE_Exp_Yes==1]
#up_boundary_prob_AE_exp <- max(probz_predicted_AE_exp_Yes)  ; up_boundary_prob_AE_exp
#low_boundary_prob_AE_exp <- min(probz_predicted_AE_exp_Yes)  ;  low_boundary_prob_AE_exp

## AE  control arm----
#logit(Pr(AE_control))= fi_AE_control * (Y_i0 - Y_i2)
logit_AE_control <- glm(AE_Control_Yes ~ CfW1,
                    data = trial_AE_control,
                    family = "binomial")

summary(logit_AE_control)
#you should get:
#> summary(logit_AE_control)

#Call:
#  glm(formula = AE_Control_Yes ~ CfW1, family = "binomial", data = trial_AE_control)

#Deviance Residuals: 
#  Min        1Q    Median        3Q       Max  
#-2.40987  -0.38049  -0.21449  -0.09241   2.04156  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#  (Intercept) -2.78108    0.06472  -42.97   <2e-16 ***
#  CfW1        -0.39944    0.01364  -29.28   <2e-16 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 4704.6  on 7006  degrees of freedom
#Residual deviance: 3233.2  on 7005  degrees of freedom
#AIC: 3237.2

#Number of Fisher Scoring iterations: 6



#predicted_AE_control <- predict(logit_AE_control, type="response", newdata = trial_AE_control)

#trial_AE_control$predicted_AE_control <- predicted_AE_control
    #View(trial_AE_control[,c(7, 20)])

#probz_predicted_AE_control_Yes <- trial_AE_control$predicted_AE_control[trial_AE_control$AE_Control_Yes==1]
#up_boundary_prob_AE_control <- max(probz_predicted_AE_control_Yes)  ; up_boundary_prob_AE_control
#low_boundary_prob_AE_control <- min(probz_predicted_AE_control_Yes)  ; low_boundary_prob_AE_control

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Logit models for SPM DGM     ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# logit(Pr(IE_ij)) = (beta_0 + b0_ij) [aka baseline MADRS] + beta_1 * Treat
# baseline MADRS10 is made of the general intercept + each individual random intercept

## these logistic regression models are run to inform the parameters (e.g., parameter for Treatment) used in the SPM DGM for intercurrent events
# fit_LoE_spm is used in order to determine the treatment parameter and to have an idea for the general intercept
# dataset with n=2000 and re_covm3

d_re <- d_mis_w
    #View(d_re)

fit_LoE_spm <- glm(LoE_Yes ~ Treat,
                   data = d_re, 
                   family = "binomial")

summary(fit_LoE_spm)
#you should get:
#summary(fit_LoE_spm)

#Call:
#  glm(formula = LoE_Yes ~ Treat, family = "binomial", data = d_re)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.1294  -1.1294  -0.8905   1.2262   1.4945  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -0.11401    0.06332  -1.801   0.0718 .  
#Treat       -0.60629    0.09249  -6.555 5.57e-11 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 2691.2  on 1999  degrees of freedom
#Residual deviance: 2647.7  on 1998  degrees of freedom
#AIC: 2651.7

#Number of Fisher Scoring iterations: 4

fit_AE_exp_spm <- glm(AE_Exp_Yes ~ 1,
                   data = d_re[d_re$Treat==1,], 
                   family = "binomial")

summary(fit_AE_exp_spm)
# you should get:
#> summary(fit_AE_exp_spm)
#Call:
#  glm(formula = AE_Exp_Yes ~ 1, family = "binomial", data = d_re[d_re$Treat == 
    #                                                               1, ])

#Deviance Residuals: 
#  Min      1Q  Median      3Q     Max  
#-0.674  -0.674  -0.674  -0.674   1.785  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -1.36639    0.07863  -17.38   <2e-16 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 1008.6  on 998  degrees of freedom
#Residual deviance: 1008.6  on 998  degrees of freedom
#AIC: 1010.6

#Number of Fisher Scoring iterations: 4


fit_AE_control_spm <- glm(AE_Control_Yes ~ 1,
                      data = d_re[d_re$Treat==0,], 
                      family = "binomial")

summary(fit_AE_control_spm)
#length(d_re$id) # = 2000
#you should get:
#> summary(fit_AE_control_spm)
#Call:
#  glm(formula = AE_Control_Yes ~ 1, family = "binomial", data = d_re[d_re$Treat == 
#                                                                       0, ])

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-0.4708  -0.4708  -0.4708  -0.4708   2.1236  

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -2.1440     0.1031  -20.79   <2e-16 ***
#  ---
#  Signif. codes:  
#  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 672.09  on 1000  degrees of freedom
#Residual deviance: 672.09  on 1000  degrees of freedom
#AIC: 674.09

#Number of Fisher Scoring iterations: 4
# These coefficients were used in the logit models in the SPM DGM.

# the parameters from these 3 models were used against the tweak.intercept function in the SPM DGM code, they are very close
# we did this in order to check and make sure the intercepts are close to the model-derived parameter estimates
# tweak.intercept can be used for any percentage of intercurrent events. The logit models for SPM were fitted on the SMd source trial with specific percentages of intercurrent events.
# tweak.intercep function is versatile

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### the code for Survival models fitted on the large trial


#View(SimTrial_sm_2000_1_5)
# Fit the Survival model: time to LoE
# but first, bring the dataset in the right format, namely all LoE at week 3 and all AE at week 2.
#3 n= 8000, Source trial via SMd, re_covm3
# reformat the dataset
# fit lme models corresponding to each intercurrent event.
SimTrial_sm_2000_1_5_surv <- SimTrial_sm_8000_1_1
class(SimTrial_sm_2000_1_5_surv$Visit)
SimTrial_sm_2000_1_5_surv$Visit <- as.numeric(SimTrial_sm_2000_1_5_surv$Visit)-1
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



# Now fit the survival models on the source trials (own source trial)
##### LoE joint model ##

cox_fit_LoE <- coxph(Surv(Visit, LoE_yes_surv) ~ Treat, data = SimTrial_sm_2000_1_5_surv_cox_LoE, x = TRUE)
summary(cox_fit_LoE)

fit_lme_LoE <- lme(fixed=MADRS10 ~ Visit * Treat, 
               random=~1 | id,
               method="REML", 
               data=SimTrial_sm_2000_1_5_surv_cox_LoE)

summary(fit_lme_LoE)

jointFit.LoE <- jointModel(fit_lme_LoE, cox_fit_LoE, timeVar = "Visit", 
                           method = "weibull-AFT-aGH")

summary(jointFit.LoE)
# you should get:
#> summary(jointFit.LoE)

#Call:
#  jointModel(lmeObject = fit_lme_LoE, survObject = cox_fit_LoE, 
#             timeVar = "Visit", method = "weibull-AFT-aGH")

#Data Descriptives:
#  Longitudinal Process		Event Process
#Number of Observations: 8000	Number of Events: 2850 (35.6%)
#Number of Groups: 8000

#Joint Model Summary:
#  Longitudinal Process: Linear mixed-effects model
#Event Process: Weibull accelerated failure time model
#Parameterization: Time-dependent 

#log.Lik     AIC      BIC
#-37166.6 74353.2 74423.07

#Variance Components:
#  StdDev
#(Intercept) 7.054264
#Residual    2.495667

#Coefficients:
#  Longitudinal Process
#Value Std.Err z-value p-value
#(Intercept)  28.6515  0.7659 37.4110 <0.0001
#Visit        -0.1646  0.1608 -1.0237  0.3060
#Treat1        2.8411  0.5716  4.9704 <0.0001
#Visit:Treat1 -0.9399  0.1151 -8.1647 <0.0001

#Event Process
#Value Std.Err z-value p-value
#(Intercept)  2.8789  0.0971 29.6500 <0.0001
#Treat1       0.3111  0.0239 13.0051 <0.0001
#Assoct      -0.0279  0.0035 -7.9562 <0.0001
#log(shape)   0.5663  0.0170 33.3322 <0.0001

#Scale: 1.7618 

#Integration:
#  method: (pseudo) adaptive Gauss-Hermite
#quadrature points: 3 

#Optimization:
#  Convergence: 0 







##### AE joint model ##


cox_fit_AE <- coxph(Surv(Visit, AE_yes_surv) ~ Treat, data = SimTrial_sm_2000_1_5_surv_cox_AE, x = TRUE)
summary(cox_fit_AE)

fit_lme_AE <- lme(fixed=MADRS10 ~ Visit * Treat, 
                   random=~1 | id,
                   method="REML", 
                   data=SimTrial_sm_2000_1_5_surv_cox_AE)

summary(fit_lme_AE)


jointFit.AE <- jointModel(fit_lme_AE, cox_fit_AE, timeVar = "Visit", 
                          method = "weibull-AFT-aGH")

summary(jointFit.AE)
# you should get:
#> summary(jointFit.AE)

#Call:
#  jointModel(lmeObject = fit_lme_AE, survObject = cox_fit_AE, timeVar = "Visit", 
#             method = "weibull-AFT-aGH")

#Data Descriptives:
#  Longitudinal Process		Event Process
#Number of Observations: 8000	Number of Events: 1171 (14.6%)
#Number of Groups: 8000

#Joint Model Summary:
#  Longitudinal Process: Linear mixed-effects model
#Event Process: Weibull accelerated failure time model
#Parameterization: Time-dependent 

#log.Lik      AIC      BIC
#-32675.62 65371.24 65441.11

#Variance Components:
#  StdDev
#(Intercept) 6.898532
#Residual    2.592738

#Coefficients:
#  Longitudinal Process
#Value Std.Err  z-value p-value
#(Intercept)   39.2867  0.6016  65.3010 <0.0001
#Visit         -1.9818  0.1060 -18.6942 <0.0001
#Treat1       -20.9535  0.6742 -31.0790 <0.0001
#Visit:Treat1   3.5870  0.1201  29.8554 <0.0001

#Event Process
#Value Std.Err z-value p-value
#(Intercept)  3.7550  0.2074 18.1051 <0.0001
#Treat1      -0.5394  0.0881 -6.1195 <0.0001
#Assoct       0.0076  0.0057  1.3465  0.1782
#log(shape)  -0.0285  0.0283 -1.0054  0.3147

#Scale: 0.9719 

#Integration:
#  method: (pseudo) adaptive Gauss-Hermite
#quadrature points: 3 

#Optimization:
# Convergence: 0 





# Selection model via marginal model----
# for outcomes-generating model and FULLY STOCHASTIC models for generation of intercurrent events
# for notation and description, please see above in the deterministic approach. The difference hereonwards is that we used logistic models to model the probability of intercurrent events, vs the deterministic rules as applied in the above approach

 # most of the code below is very similar with the code above, in terms of scope and operational reasons, as well as meaning and interpretation
 # code for longitudinal outcomes generation is as above

#rm(list=ls()) #
# needed for the selection model method
library(gmailr)
library(MASS)#
library(tidyverse)#
library(nlme)#
library(lme4)#
#library(Hmisc)
library(janitor)#
library(gt)#
library(patchwork)#

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

#set.seed(2147483399)
n <- 190# number of patients at trial level to be randomised 1:1, experimental and control arm. 190 in total
m.iterations <- 500# 382number of generated datasets needed for verification of longitudinal outcomes# number of trials per scaling factor
scaling_factor <-  c(1)#c(0.5, 1.0, 1.5, 2.0, 2.5)
# total number of simulated trials = m.iterations * length(scaling_factor)
# try with c(0.4, 1.1, 1.8, 2.3, 3)
#c(0.25, 0.5, 1, 2, 2.5) steps for scaling factor
# to multiply for LoE, AE_control and AE_exp ### value 0.5 should yield around 13.5% IEs, 1 ~ %, 2 ~ %, 2.5 ~

#CFE <- matrix(ncol=4,nrow=length(scaling_factor)*m.iterations)
#colnames(CFE) <-c("N ceiled_floored", "% ceiled_floored", "scaling factor", "simulated trial n")

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline + follow-up measurements
delta_SMs <- matrix(ncol=1,nrow=m.iterations) # treatment effect estimate at 6 weeks based on MMRM models fitted on each generated dataset
colnames(delta_SMs) <-c("TreatmentEffect")
betas_SMs <- matrix(ncol=2,nrow=m.iterations)
colnames(betas_SMs) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3)

#confint_fit <- matrix(ncol=2,nrow=m.iterations) 
#colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
#delta_SMs_errorz <- matrix(ncol=1,nrow=m.iterations)
#colnames(delta_SMs_errorz) <- c("SE")

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

start_time <- Sys.time() # timestamp for the start time of the nested for loop below.
# it was used to have an estimate of time needed for different larger number of trials to be simulated upon scaling up the simulation parameters (e.g., m.iterations)
set.seed(2147483629) # set seed
## Begin for loop----
for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    ### Generate longitudinal outcomes----
    #### Generate correlated residuals----
    re_means <- c(0, 0, 0, 0, 0, 0, 0) # means of residuals at each timepoint
    
    # covariance matrix extracted from trial 003-002 (see reference to the data analysis paper in PST)
    # this dataset is not available to share. A MMRM model was fitted and this is the corresponding covariance matrix
    re_covm <- matrix(c(20.2190, 17.149, 14.721, 13.087,  8.4329,  10.854,   4.6417, 
                        17.1490, 48.536, 41.161, 32.151, 24.8400,  30.528,  26.0170,
                        14.7210, 41.161, 72.569, 57.866, 60.2200,  61.974,  54.5400,
                        13.0870, 32.151, 57.866, 74.080, 66.2960,  63.540,  52.1070,
                        8.4329, 24.840, 60.220, 66.296, 97.4730,  90.612,  80.1370,
                        10.8540, 30.528, 61.974, 63.540, 90.6120, 116.410, 102.8300,
                        4.6417, 26.017, 54.540, 52.107, 80.1370, 102.830, 109.5900), nrow = 7)
    
    #Standard Deviations: 4.4965 6.9668 8.5188 8.607 9.8728 10.789 10.468
    
    # covariance matrix with 0 off-diagonal and small variances. This is useful for initial/later checks to see if the simulated data corresponds to target data to be simulated
    size_diagonal <- 0.0000001
    re_covm2 <-matrix(c(size_diagonal, 0, 0, 0, 0, 0, 0,
                        0, size_diagonal, 0, 0, 0, 0, 0,
                        0, 0, size_diagonal, 0, 0, 0, 0,
                        0, 0, 0, size_diagonal, 0, 0, 0,
                        0, 0, 0, 0, size_diagonal, 0, 0,
                        0, 0, 0, 0, 0, size_diagonal, 0,
                        0, 0, 0, 0, 0, 0, size_diagonal), nrow = 7)
    
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
    d <-d[,-5] # remove re (residuals) column from the dataset, they have been added to betas_SMd
    
    # assign this to another object to make sure each time for each analysis the dataset used is the same
    d_orig<-d # full outcome data
    
    # create a separate column with only the baseline outcomes, if the baseline values will be used as covariate in the model 
    length(d$id)
    tmp <- sapply(unique(d$id), FUN = function(i) nrow(d[d$id == i,]))
    BaselineMADRS10 <-  rep(d$MADRS10[d$visit == 0], tmp)
    length(BaselineMADRS10)
    d$Baseline <- BaselineMADRS10
    #View(d)
    
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
    fit <- gls(MADRS10 ~ visit*Treat, 
               data=d,
               correlation = corSymm(form=~1 | id),
               #weights = varIdent(form = ~ 1 | visit),
               method="REML")
    
    summary(fit)
    
    fit$coefficients[c(8,14)]

    sum(fit$coefficients[c(8,14)]); treatmenteffect
    
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
    betas_SMs[m, ] <- fit$coefficients[c(8,14)]
    
    delta_SMs[m, ] <- sum(fit$coefficients[c(8,14)])
    
    #bias_f[m, ] <- sum(fit$coefficients[c(7,13)]) - treatmenteffect
    
    #delta_SMs_error <- sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"]) 
    
    #delta_SMs_errorz[m, ] <- delta_SMs_error 
    
    
    #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_SMs_error
    #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_SMs_error
    
    Randomised_Exp <- sum(d[,3])/7 #number of patients in the experimental arm
    Randomised_Control <- n-sum(d[,3])/7 #number of patients in the control arm
    #View(d)
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
    
    d_mis_w <- d_mis |> spread(visit, MADRS10) # reshape to wide in order to create the CfB variable
         #View(d_mis_w)
    
    colnames(d_mis_w)[3:9] <- c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6"); head(d_mis_w)
    
    d_mis_w$CfB <- d_mis_w[,"Baseline",] - d_mis_w[,"Week6"]; d_mis_w # create the CfB variable
    
    #View(d_mis_w)
    
    d_mis_w$CfW1 <- d_mis_w[,"Baseline"] - d_mis_w[,"Week1"]; d_mis_w # create the CfW1 variable
    
    d_mis_L <- d_mis_w |> gather(Visit, MADRS10, Baseline:Week6) # reshape to long format
    
    d_mis_L <- d_mis_L[order(d_mis_L$id, d_mis_L$Visit),]; #d # order by subject id and Visit
    #summary(logit_LoE)
    
    d_mis_L$predicted_LoE <- plogis((1.077405 + d_mis_L$CfB*-0.355617))
    #summary(logit_LoE)
    #d_mis_L$predicted_LoE <- predict(logit_LoE, type="response", newdata=d_mis_L) # predicted occurrence of intercurrent events based on estimated probabilities from the logit models for e.g., treatment discontinuation due to lack of efficacy at trial level
        #d#describe(predict(logit_LoE, type="response", newdata=d_mis_L)) # check 
    #View(d_mis_L)
    
    d_mis_L$LoE_yes <- rep(rbinom(length(unique(d_mis_L$id)), 1, unique(d_mis_L$predicted_LoE)), each = length(visits)) 
    #d_mis_L$CfB
    #View(d_mis_L)
    
    trial_AE_X <- d_mis_L[d_mis_L$Treat==1,]
    trial_AE_X$predicted_AE <- plogis(-2.370965 + trial_AE_X$CfW1*0.272188)
    #summary(logit_AE_exp)
    #trial_AE_X$predicted_AE <- predict(logit_AE_exp, type="response", newdata = trial_AE_X) #  predicted occurrence of intercurrent events based on estimated probabilities from the logit models for e.g., treatment discontinuation due to adverse events in the experimental arm
        #describe(predict(logit_LoE, type="response", newdata=d_mis_L))
    
    trial_AE_X$AE_yes <- rep(rbinom(length(unique(trial_AE_X$id)), 1, unique(trial_AE_X$predicted_AE)), each = length(visits)) 
    
    trial_AE_X$CfW1
    
    trial_AE_C <- d_mis_L[d_mis_L$Treat==0,]
    
    trial_AE_C$predicted_AE <- plogis(-2.78108 + trial_AE_C$CfW1*-0.39944)
    #summary(logit_AE_control)
    #trial_AE_C$predicted_AE <- predict(logit_AE_control, type="response", newdata = trial_AE_C) # predicted occurrence of intercurrent events based on estimated probabilities from the logit models for e.g., treatment discontinuation due to adverse events in the control arm

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
    
    # Check for 
    #View(d_mis_LL)
    #d_mis_LL$compete <- ifelse(d_mis_LL$AE_yes==1 & d_mis_LL$LoE_yes==1, 1, 0)
    #sum(as.numeric(d_mis_LL$AE_YES)-1)/length(visits)/n*100
    #sum(as.numeric(d_mis_LL$LoE_YES)-1)/length(visits)/n*100
    #describe(d_mis_LL$Behavior)
    
    
    d_mis_LL$Behavior <- ifelse(d_mis_LL[,"AE_YES"]==1, "AE",
                               ifelse(d_mis_LL[,"LoE_YES"]==1, "LoE", "No IE"))
    
    d_mis_LL$Behavior <-factor(d_mis_LL$Behavior, levels = c("AE", "LoE", "No IE"))
    
    class(d_mis_LL$Visit)
    d_mis_LL$Visit <- as.factor(d_mis_LL$Visit)
    rownames(d_mis_LL) <-NULL 
    
        # check range of values
        #range(d_mis_LL$MADRS10[d_mis_LL$Treat==1], na.rm = T)
        #range(d_mis_LL$MADRS10[d_mis_LL$Treat==0], na.rm = T)
    
    #describe(d_mis_LL$AE_YES[d_mis_LL$Treat==0])
    #describe(d_mis_LL$AE_YES[d_mis_LL$Treat==1])
    #describe(d_mis_LL$AE_YES)
    
    #describe(d_mis_LL$Behavior[d_mis_LL$Treat==0])
    #describe(d_mis_LL$Behavior[d_mis_LL$Treat==1])
    #describe(d_mis_LL$Behavior)
    
    #describe(d_mis_LL$LoE_YES[d_mis_LL$Treat==0])
    #describe(d_mis_LL$LoE_YES[d_mis_LL$Treat==1])
    #describe(d_mis_LL$LoE_YES)
        #View(d_mis_LL)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Intercurrent events descriptives needed for the verification step ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #describe(d_mis_LL)
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
    
    ## AE ##

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
    plot_LoE_SMs <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") +
      facet_wrap(~ Treat)+ 
      ggtitle("SMs-LoE pattern") +
      scale_y_continuous(limits = c(-10, 60))+
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
      facet_wrap(~ Treat)+ ggtitle("SMs-AE pattern")  +
      scale_y_continuous(limits = c(-10, 60))+
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
    
    
    (plot_all_SMs + plot_AE_SMs)/(plot_LoE_SMs+plot_NoIE_SMs)
    

    setTxtProgressBar(pb1, m)
  }
  

  # parameters extracted for MMRM fitted models on full outcome data
  colMeans(betas_SMs)
  colMeans(delta_SMs) ; treatmenteffect
  
  
  ### assign and save parameters---- 
  assign(paste('all_betas_SMs', s, sep="_"), betas_SMs)
  assign(paste('all_delta_SMs', s, sep="_"), delta_SMs)
  #assign(paste('delta_SMs_errorz', s, sep="_"), delta_SMs_errorz)
  #assign(paste('all_confint_fit', s, sep="_"), confint_fit)
  assign(paste('all_N_Exp', s, sep="_"), N_Exp)
  assign(paste('all_N_Control', s, sep="_"), N_Control)
  
  setTxtProgressBar(pb3, s)
}

## End for loop----
end_time <- Sys.time()
end_time-start_time

colMeans(all_betas_SMs_1)
#> colMeans(all_betas_SMs_1)
#Treat visit42:Treat 
#-1.637695     -1.783089 

tolerance_margin <- 0.1 

difference_Verification_SMs <- abs(treatmenteffect - colMeans(all_delta_SMs_1))

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_Verification_SMs) < tolerance_margin), "Verification SMs *SUCCESSFUL*", "Verification SMs NOT successful :(") 



min(all_betas_SMs_1)
max(all_betas_SMs_1)



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


table_AE_SMs |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "AE") |> 
  rename(N_C_arm=N.AE.Control) |> 
  rename(N_E_arm=N.AE.Exp)

table_LoE_SMs |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "LoE") |> 
  rename(N_C_arm=N.LoE.Control) |> 
  rename(N_E_arm=N.LoE.Exp)


tab_SMs <- tibble(bind_rows(table_AE_SMs |> 
                               as.data.frame() |> 
                               mutate("Intercurrent event" = "AE") |> 
                               rename(N_C_arm=N.AE.Control) |> 
                               rename(N_E_arm=N.AE.Exp), 
                             table_LoE_SMs |> 
                               as.data.frame() |> 
                               mutate("Intercurrent event" = "LoE") |> 
                               rename(N_C_arm=N.LoE.Control) |> 
                               rename(N_E_arm=N.LoE.Exp))); tab_SMs



tab2_SMs <- tab_SMs |> group_by(`Intercurrent event`) |>
  summarise("N" = round(mean(N_C_arm), digits=1), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = round(mean(N_E_arm), digits=1), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = round(mean(N_C_arm + N_E_arm), digits=1),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) |> 
  adorn_totals("row"); tab2_SMs




gt(tab2_SMs) |> 
  tab_header(title = md("Table 8e. Descriptive statistics intercurrent events"), subtitle = md("Selection model DGM - stochastic")) |>
  tab_source_note(md(paste0("Averaged over", " ", m.iterations*length(scaling_factor),  " ",  "simulated trials.", " ", "Trial sample size = ", " ", n ))) |> 
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
# 7 February 2022
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
#  [1] splines   stats     graphics  grDevices utils    
#[6] datasets  methods   base     

#other attached packages:
#  [1] JM_1.4-8        foreign_0.8-82  survival_3.2-13
#[4] patchwork_1.1.1 gt_0.3.1        janitor_2.1.0  
#[7] lme4_1.1-28     Matrix_1.4-0    nlme_3.1-155   
#[10] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
#[13] purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
#[16] tibble_3.1.6    ggplot2_3.3.5   tidyverse_1.3.1
#[19] MASS_7.3-55     gmailr_1.0.1   


installed.packages()

