## License CC-BY-4.0 
## Creative Commons Attribution 4.0 International
## Mitroiu M. et al SMMR





## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"



### DGM
#rm(list=ls())
library(MASS)
library(nlme)
library(survival)
library(foreign)
library(tidyverse)
library(janitor)

#install.packages('tinytex')
library(tidyr)
library(haven)
library(Hmisc)
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





n <- 190# number of patients

#0.11*sqrt(24.611*1.157)
    
#0.12*sqrt(24.893*1.049)
    # resume here, integrate this code in this structure


## simulation parameters #----

set.seed(2147483629)
    
    
m.iterations <- 10 # number of generated datasets # number of trials per scaling factor
scaling_factor <-  c(0.5, 1.0, 1.5, 2.0, 2.5) # scaling factor used to vary the percentages of intercurrent events at trial/iteration level
    # total number of simulated trials = m.iterations * length(scaling_factor)
    # other ranges can be used to ensure variability between simulated trials, as long as they are as envisaged over all simulated trials (e.g., mean percentages)
    # and check out the verification step
    
# these will be used in the target proportions of intercurrent events used in the function to determine the intercept value in order to obtain the right percentage of intercurrent events
# ranges of probabilities centered around desired percentages of each intercurrent events averaged over all simulated trials
# this is done to increase variability in intercurrent events percentages between trials
p_LoE_sample <-c(0.31, 0.33, 0.34, 0.35, 0.37); mean(p_LoE_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy at trial level
p_AE_Exp_sample <- c(0.055, 0.065, 0.075, 0.085, 0.095)*2; mean(p_AE_Exp_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy in the experimental arm
p_AE_Control_sample <- c(0.05, 0.10, 0.15, 0.20, 0.25)/2; mean(p_AE_Control_sample) # proportion of e.g.,  treatment discontinuation due to lack of efficacy in the control arm


    
    visits <- as.numeric(c(0, 1, 2, 3, 4, 5, 6))	

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
    
    
    start_time <- Sys.time() # timestamp for the start time of the nested for loop below.
    # it was used to have an estimate of time needed for different larger number of trials to be simulated upon scaling up the simulation parameters (e.g., m.iterations)
    
    ## Begin for loop----
    
    for (s in 1:length(scaling_factor)) {
      for(m in 1:m.iterations) {
    
        
        # The specification of the model: random intercept and random slope + some random noise
        # d$MADRS10.true <- (b0 + d$bi_0) + (b1 + d$bi_1) * d$visit +  b2 * d$visit * d$Treat
        # d$MADRS10_collected <- d$MADRS10.true + rnorm(nrow(d), 0, eps.sd)
        
        
        ### Generate longitudinal outcomes----
        #### Generate correlated (random) intercepts and slopes for each individual patient----
        
        # We used as inspiration the trial 003-002 (ref data analysis paper). Here we used a simplified scenario with linear group trajectories.
        
        #### model parameters----
        
        b0 <- 29.79
        b1 <- -0.55
        b2 <- -0.583
        bi_means <- c(0, 0)
        bi_covm <- matrix(c(24.611, 0.5869809, 0.5869809, 1.157), nrow = 2)
        bi <- mvrnorm(n, bi_means, bi_covm)	
        eps.sd <-3.247  
        
        c1 <- 0.1
        
    
    d <- data.frame(
      id = rep(1:n, each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n, 1, 0.5), each = length(visits)),
      bi_0 = rep(bi[,1], each = length(visits)),
      bi_1 = rep(bi[,2], each = length(visits))
    )
    
    
    
    # insert model specification from fitted models on the SM super trial
    # insert Rutger's code and adjust it
    
  
    
    #####################################################
    # Logit model and Probabilities for LoE at trial level
    logit_Pr_LoE <- (b0 + d$bi_1)/1000 +  c1 * d$Treat
    Pr_LoE <- 1/(1+exp(-logit_Pr_LoE))
    
    describe(Pr_LoE[d$Treat==1])
    describe(Pr_LoE[d$Treat==0])
    
    d$Pr_LoE <- Pr_LoE
    
    d$LoE_yes[d$Treat==1] <- ifelse(Pr_LoE[d$Treat==1]<0.532, 1, 0) * rep(rbinom(length(unique(d$id[d$Treat==1])), 1, 0.7), each = length(visits))
    d$LoE_yes[d$Treat==0] <- ifelse(Pr_LoE[d$Treat==0]>0.5075, 1, 0) * rep(rbinom(length(unique(d$id[d$Treat==0])), 1, 0.7), each = length(visits))
  
    #####################################################
    # Logit model and Probabilities for AE in experimental arm
    logit_Pr_AE_exp <- (d$bi_1[d$Treat==1])*5  # for the subset of experimental arm patients
    
    
    
    
    #View(cbind(logit_Pr_AE_exp, Pr_AE_exp, exp(-logit_Pr_AE_exp)))
    
    #View(d)
    
    describe(d$bi_1[d$Treat==1])
    hist(d$bi_1[d$Treat==1])
    
    #1/(1+exp(-(-0.45*5)))
    
    
    Pr_AE_exp <- 1/(1+exp(-logit_Pr_AE_exp))
    
    describe(Pr_AE_exp)
    hist(Pr_AE_exp)
    d$Pr_AE_exp[d$Treat==1] <- Pr_AE_exp
    
    
    d$AE_yes[d$Treat==1] <- ifelse(Pr_AE_exp>0.9944527, 1, 0)
    
    
    
    
    
    #####################################################
    # Logit model and Probabilities for AE in control arm
    logit_Pr_AE_control <- (b1 + d$bi_1[d$Treat==0])*2 # for the subset of control arm patients
    
    Pr_AE_control <- 1/(1+exp(-logit_Pr_AE_control))
    
    describe(Pr_AE_control)
    hist(Pr_AE_control)
    d$Pr_AE_control[d$Treat==0] <- Pr_AE_control
    
    
    d$AE_yes[d$Treat==0] <- ifelse(Pr_AE_control>0.55, 1, 0) # add rbinom probabilities
    

    
    #View(d)
    d
    
    class(d$visit)
    class(d$Treat)
    d$Treat <- factor(d$Treat)
    d$LoE_yes <- factor(d$LoE_yes)
    d$AE_yes <- factor(d$AE_yes)
    
    
    d$AE_YES <- ifelse(d$AE_yes==1, 1, 0)
    d$LoE_YES <- ifelse(d$AE_yes==0 & d$LoE_yes==1, 1, 0)
    
    
    
    d$Behavior <- ifelse(d[,13]==1, "AE",
                                ifelse(d[,14]==1, "LoE", "No IE"))
    
    
    
    
    p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60))
    
    
    
    # Plot trajectories
    
    # just LoE
    p<- ggplot(data = d[d$LoE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
    plot_LoE <-  p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("SPM-LoE pattern")
    
    
    #just AE
    # AE
    p<- ggplot(data = d[d$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
    plot_AE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("SPM-AE pattern")
    
    

    
    
    d_SPM <- d

    
    # All patients with TRUE trajectory 
    # LoE
    p<- ggplot(data = d_SPM[d_SPM$Behavior=="LoE",], aes(x = factor(visit), y = MADRS10_collected, group = id, color=LoE_YES)) 
    plot_LoE_SPM <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SPM-LoE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_LoE_SPM
    
    # AE
    p<- ggplot(data = d_SPM[d_SPM$Behavior=="AE",], aes(x = factor(visit), y = MADRS10_collected, group = id, color=AE_YES)) 
    plot_AE_SPM <- p + geom_line(size=0.5, color='#F8766D') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60)) + ggtitle("SPM-AE pattern")  +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_AE_SPM
    
    #just No IE
    # No IE
    p<- ggplot(data = d_SPM[d_SPM$Behavior=="No IE",], aes(x = factor(visit), y = MADRS10_collected, group = id)) 
    plot_NoIE_SPM <- p + geom_line(size=0.5, color='#619CFF') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="blue") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("SPM-No IE pattern") +
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")) ; plot_NoIE_SPM
    
    
    # All behaviors
    p<- ggplot(data = d_SPM, aes(x = factor(visit), y = MADRS10_collected, group = id, color=Behavior)) 
    plot_all_SPM <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="black") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))+ ggtitle("SPM-All patterns")+
      scale_x_discrete(labels=c("0", "1", "2", "3", "4", "5", "6")); plot_all_SPM
    
#describe(d)

#View(d)
    
the_plot_SPM <- (plot_all_SPM / plot_LoE_SPM) | (plot_AE_SPM / plot_NoIE_SPM); the_plot_SPM
    
#plot_LoE_SPM + plot_LoE_JM
    
    fit_lmer <- lmer(MADRS10_collected ~ visit + visit:Treat + (1 + visit|id), data = d, REML = T)
    
    summary(fit_lmer)
    

    
    
    fit_lme <- lme(fixed=MADRS10_collected ~ visit + visit:Treat, 
                   random=~1 + visit| id,
                   method="REML", 
                   data=d)
    
    summary(fit_lme)
  
    -0.788534/0.1696878
    
    
    0.1696878*6
    
    getVarCov(fit_lme, type= "marginal")
    
    
    getVarCov(fit_lme,type = c("conditional"))
    
    
    getVarCov(fit_lme,type = c("random.effects"))
    
    #vcov(fit_lme)
    
    
    #sqrt(9.0522)
    
    
    
    #(4*3^2)/(0.1^2)
    #3600 trials vs 502 in the MMRM
    # need to derive the marginal covariance structure and from there to take the variance?
    
    
    #t(matrix(c(0, 1, 1, 0), nrow=2))
    
    
    
# visualise
# get descriptive statistics

    
    

    
    
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
    
    

    p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 70))
    
    
    
    
    fit<-gls(MADRS10_collected ~ visit * Treat + Baseline, 
             data=d,
             correlation = corSymm(form=~1 | id),
             #weights = varIdent(form = ~ 1 | visit), 
             na.action=na.exclude, 
             method="REML")
    
    
    
    
    summary(fit)
    
    
    fit$coefficients[c(7,13)]
    
    sum(fit$coefficients[c(7,13)]); treatmenteffect
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # LMM on full outcome data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    #fit_lme <- lme(fixed=MADRS10 ~ visit * Treat + Baseline, 
    #           random=~1 | id,
    #          method="REML", 
    #         correlation = corSymm(form=~1|id),
    #        na.action=na.omit,
    #       data=d)
    
    #summary(fit_lme)
    #model_parameters(fit_lme)
    
    #refine 
    
    betas[m, ] <- fit$coefficients[c(7,13)]
    
    delta[m, ] <- sum(fit$coefficients[c(7,13)])
    
    bias_f[m, ] <- sum(fit$coefficients[c(7,13)]) - treatmenteffect
    
    delta_error <- sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"]) 
    
    delta_errorz[m, ] <- delta_error 
    
    
    confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_error
    confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_error
    
    
    Randomised_Exp <- sum(d[,3])/6 #number of patients in the experimental arm
    Randomised_Control <- n-sum(d[,3])/6 #number of patients in the control arm
    
    N_Exp[m,] <- Randomised_Exp
    N_Control[m,] <- Randomised_Control
    
    
    
    
    
    
    ####################################################################################################  
    ####################################################################################################
    
    # IEGM
    # Generate missing outcome dataset
    # create variable with difference week 6 - baseline
    # if Difference < 5, then missing outcome from visit 3 onwards for LoE
    
    ## take over the original/raw dataset to use for IEGM/MDGM
    d_mis <-d_orig
    
    d_mis_w <- d_mis %>% spread(visit, MADRS10) # reshape to wide in order to create the CfB variable
    #View(d_mis_w)
    
    colnames(d_mis_w)[3:9] <- c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6"); head(d_mis_w)
    
    d_mis_w$CfB <- d_mis_w[,3] - d_mis_w[,9]; d_mis_w # create the CfB variable
    
    d_mis_w$CfW2 <- d_mis_w[,3] - d_mis_w[,5]; d_mis_w # create the CfW2 variable
    
    #View(d_mis_w)
    
    p_LoE <-sample(p_LoE_sample, 1) * scaling_factor[s]
    
    d_mis_w$LoE_Yes <-ifelse(d_mis_w$CfB<5, 1, 0)*rbinom(n, 1, p_LoE);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    #sum(d_mis_w$LoE_Yes)
    #View(d_mis_w)
    
    # integrate for AE
    # too much efficacy >8 points on MADRS10 at week 2
    
    p_AE_Exp <- sample(p_AE_Exp_sample, 1) * scaling_factor[s]
    
    d_mis_w$AE_Exp_Yes <-ifelse(d_mis_w$Treat==1 & d_mis_w$CfW2>8, 1, 0)*rbinom(n, 1, p_AE_Exp);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    sum(d_mis_w$AE_Exp_Yes)
    
    p_AE_Control <- sample(p_AE_Control_sample, 1) * scaling_factor[s]
    
    d_mis_w$AE_Control_Yes <-ifelse(d_mis_w$Treat==0 & d_mis_w$CfW2< (-2), 1, 0)*rbinom(n, 1, p_AE_Control);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    sum(d_mis_w$AE_Control_Yes)
    
    #
    
    d_mis_L <- d_mis_w %>% gather(Visit, MADRS10, Baseline:Week6) # reshape to long format
    
    d_mis_L <- d_mis_L[order(d_mis_L$id, d_mis_L$Visit),]; #d # order by subject id and Visit
    
    
    
    #weeks <-c(2, 3, 4, 5)
    #weeks_IE <- sample(weeks, sum(d_mis_w$LoE_Yes), replace = TRUE)
    # resume here
    
    
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
    
    
    d_mis_L$AE_Yes <- ifelse(d_mis_L[,6]==1, 1, 
                             ifelse(d_mis_L[,7]==1, 1, 0))
    
    
    
    
    # create the NAs by AE variable
    #d_mis_L[,10] <-ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week3", "NA",
    #                     ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week4", "NA", 
    #                           ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week5", "NA",
    #                                 ifelse(d_mis_L[,11]==1 & d_mis_L[,8]=="Week6", "NA", d_mis_L[,10]))))
    
    #View(d_mis_L)
    
    d_mis_L$LoE_YES <- ifelse(d_mis_L[,5]==1 & d_mis_L[,10]==0, 1, 0)
    
    
    
    #View(d_mis_L)
    
    d_mis_L <- d_mis_L[order(d_mis_L$id),]; #d # order by subject id and Visit
    #colnames(d_mis_L)[10] <- c("MADRS10_mis"); d_mis_L # name the column
    
    
    # check classes of variable and transform in order 
    class(d_mis_L$id)
    class(d_mis_L$Treat)
    d_mis_L$Treat<-as.factor(d_mis_L$Treat)
    class(d_mis_L$Visit)
    d_mis_L$LoE_Yes <- as.factor(d_mis_L$LoE_Yes)
    d_mis_L$LoE_YES <- as.factor(d_mis_L$LoE_YES)
    d_mis_L$AE_Yes <- as.factor(d_mis_L$AE_Yes)
    
    d_mis_L$Behavior <- ifelse(d_mis_L[,10]==1, "AE",
                               ifelse(d_mis_L[,11]==1, "LoE", "No IE"))
    
    
    d_mis_L$Behavior <-factor(d_mis_L$Behavior, levels = c("AE", "LoE", "No IE"))
    
    class(d_mis_L$Visit)
    d_mis_L$Visit <- as.factor(d_mis_L$Visit)
    rownames(d_mis_L) <-NULL 
    
    # check range of values
    range(d_mis_L$MADRS10[d_mis_L$Treat==1], na.rm = T)
    range(d_mis_L$MADRS10[d_mis_L$Treat==0], na.rm = T)
    
    
    
    
    
    #do with colsums instead of table to avoid Error subscript out of bounds when there is no LoE generated
    #View(d_mis_L)
    
    LoE_Y <- d_mis_L[,c(2, 11)]
    LoE_Y$LoE_YES <- as.numeric(LoE_Y$LoE_YES)-1
    
    
    LoE_Y_total <- sum(LoE_Y$LoE_YES)/length(visits)
    LoE_Y_Exp <- sum(LoE_Y$LoE_YES[LoE_Y$Treat==1])/length(visits)
    LoE_Y_Control <- sum(LoE_Y$LoE_YES[LoE_Y$Treat==0])/length(visits)
    
    
    ## LoE
    tb_LoE_total <- LoE_Y_total
    tb_LoE_Exp <- LoE_Y_Exp
    tb_LoE_Control <- LoE_Y_Control
    
    
    n_LoE_total[m, ] <-  tb_LoE_total
    
    n_LoE_Exp[m, ] <- tb_LoE_Exp
    
    n_LoE_Control[m, ] <- tb_LoE_Control
    
    ###
    
    
    #do with colsums instead of table to avoid Error subscript out of bounds when there is no AE generated
    
    
    AE_Y <- d_mis_L[,c(2, 10)]
    AE_Y$AE_Yes <- as.numeric(AE_Y$AE_Yes)-1
    
    
    AE_Y_total <- sum(AE_Y$AE_Yes)/length(visits)
    AE_Y_Exp <- sum(AE_Y$AE_Yes[AE_Y$Treat==1])/length(visits)
    AE_Y_Control <- sum(AE_Y$AE_Yes[AE_Y$Treat==0])/length(visits)
    
    
    
    ## AE
    tb_AE_total <- AE_Y_total
    tb_AE_Exp <- AE_Y_Exp
    tb_AE_Control <- AE_Y_Control
    
    
    n_AE_total[m, ] <-  tb_AE_total
    
    
    n_AE_Exp[m, ] <- tb_AE_Exp
    
    
    n_AE_Control[m, ] <- tb_AE_Control
    
    
    
    #LoE
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
    
    
    
    # Plot trajectories
    # All patients with TRUE trajectory 
    # LoE
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=LoE_YES)) 
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    
    # AE
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=AE_Yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    # All behaviors
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=Behavior)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    
    
    
    
    
    d_mis_L_LoE <- d_mis_L[d_mis_L$LoE_YES==1,] # subset only patients that experienced LoE
    d_mis_L_AE <- d_mis_L[d_mis_L$AE_Yes==1,] # subset only patients that experienced AE
    
    
    # All patients with true trajectory with different colours by LoE (Y/N)
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=LoE_YES))
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark red") + facet_wrap(~ Treat)
    
    # just LoE patients with true trajectory
    p<- ggplot(data = d_mis_L_LoE, aes(x = Visit, y = MADRS10, group = id))
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    #View(d_mis_L)
    
    
    # All patients with true trajectory with different colours by AE (Y/N)
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=AE_Yes))
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark red") + facet_wrap(~ Treat)
    
    # just AE patients with true trajectory
    p<- ggplot(data = d_mis_L_AE, aes(x = Visit, y = MADRS10, group = id))
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    
    
    
    setTxtProgressBar(pb1, m)
  }
  
  
  
  
  
  
  # parameters extracted for MMRM fitted models on full outcome data
  colMeans(betas)
  colMeans(delta) ; treatmenteffect
  
  
  
  
  # assign   
  assign(paste('all_betas', s, sep="_"), betas)
  
  assign(paste('all_delta', s, sep="_"), delta)
  
  assign(paste('all_delta_errorz', s, sep="_"), delta_errorz)
  
  assign(paste('all_confint_fit', s, sep="_"), confint_fit)
  
  assign(paste('all_N_Exp', s, sep="_"), N_Exp)
  
  assign(paste('all_N_Control', s, sep="_"), N_Control)
  
  
  
  ### plot bias
  
  
  
  setTxtProgressBar(pb3, s)
}


end_time <- Sys.time()

end_time-start_time

colMeans(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5))

colMeans(rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5))





# compile to report
cbind(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5),
      rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5),
      rbind(all_delta_errorz_1, all_delta_errorz_2, all_delta_errorz_3, all_delta_errorz_4, all_delta_errorz_5),
      rbind(all_confint_fit_1, all_confint_fit_2, all_confint_fit_3, all_confint_fit_4, all_confint_fit_5),
      rbind(all_N_Exp_1, all_N_Exp_2, all_N_Exp_3, all_N_Exp_4, all_N_Exp_5),
      rbind(all_N_Control_1, all_N_Control_2, all_N_Control_3, all_N_Control_4, all_N_Control_5))





# Table for the paper ----

table_AE_SM <- data.frame(
  # descriptives AE  
  n_AE_Control,
  n_AE_Exp); table_AE_SM

mean(n_AE_Control)
mean(n_AE_Exp)




# descriptives LoE  
table_LoE_SM <-data.frame(
  n_LoE_Control,
  n_LoE_Exp); table_LoE_SM

mean(n_LoE_Control)
mean(n_LoE_Exp)


#describe(table_IE_SM)


table_AE_SM %>% 
  as.data.frame() %>% 
  mutate("Intercurrent event" = "AE") %>% 
  rename(N_C_arm=N.AE.Control) %>% 
  rename(N_E_arm=N.AE.Exp)

table_LoE_SM %>% 
  as.data.frame() %>% 
  mutate("Intercurrent event" = "LoE") %>% 
  rename(N_C_arm=N.LoE.Control) %>% 
  rename(N_E_arm=N.LoE.Exp)


tab_SM <- tibble(bind_rows(table_AE_SM %>% 
                             as.data.frame() %>% 
                             mutate("Intercurrent event" = "AE") %>% 
                             rename(N_C_arm=N.AE.Control) %>% 
                             rename(N_E_arm=N.AE.Exp), 
                           table_LoE_SM %>% 
                             as.data.frame() %>% 
                             mutate("Intercurrent event" = "LoE") %>% 
                             rename(N_C_arm=N.LoE.Control) %>% 
                             rename(N_E_arm=N.LoE.Exp))); tab_SM



tab2_SM <- tab_SM %>% group_by(`Intercurrent event`) %>%
  summarise("N" = mean(N_C_arm), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = mean(N_E_arm), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = mean(N_C_arm + N_E_arm),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) %>% 
  adorn_totals("row"); tab2_SM




gt(tab2_SM) %>% 
  tab_header(title = md("Table 4. Descriptive statistics intercurrent events"), subtitle = md("Selection model DGM")) %>%
  tab_source_note(md("Averaged over n simulated trials")) %>% 
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



