## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"


### DGM
rm(list=ls())
library(MASS)
library(nlme)
library(survival)
library(foreign)
library(tidyverse)
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
library(optimx)
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

# simulation parameters

set.seed(2147483629) # set seed
#set.seed(2147483399)
m.iterations <- 1 # number of generated datasets # number of trials per scaling factor
scaling_factor <-  c(0.5, 1.0, 1.5, 2.0, 2.5) # scaling factor used to vary the percentages of intercurrent events at trial/iteration level
# total number of simulated trials = m.iterations * length(scaling_factor)
# try with c(0.4, 1.1, 1.8, 2.3, 3)
#c(0.25, 0.5, 1, 2, 2.5) steps for scaling factor
# to multiply for LoE, AE_control and AE_exp ### value 0.5 should yield around 13.5% IEs, 1 ~ %, 2 ~ %, 2.5 ~



n <- 2000# number of patients

CFE <- matrix(ncol=4,nrow=length(scaling_factor)*m.iterations)
colnames(CFE) <-c("N ceiled_floored", "% ceiled_floored", "scaling factor", "simulated trial n")


p_LoE_sample <-c(0.31, 0.33, 0.34, 0.35, 0.37); mean(p_LoE_sample)
p_AE_Exp_sample <- c(0.055, 0.065, 0.075, 0.085, 0.095)*2; mean(p_AE_Exp_sample)
p_AE_Control_sample <- c(0.05, 0.10, 0.15, 0.20, 0.25)/2; mean(p_AE_Control_sample)

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline + follow-up measurements
delta <- matrix(ncol=1,nrow=m.iterations) # treatment effect estimate at 6 weeks based on MMRM models fitted on each generated dataset
colnames(delta) <-c("TreatmentEffect")
betas <- matrix(ncol=2,nrow=m.iterations)
colnames(betas) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3)

confint_fit <- matrix(ncol=2,nrow=m.iterations) 
colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
delta_errorz <- matrix(ncol=1,nrow=m.iterations)
colnames(delta_errorz) <- c("SE")

bias_f <- matrix(ncol=1,nrow=m.iterations)
colnames(bias_f) <- c("bias_f")

## randomisation objects, allocation
N_Exp  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Exp) <-c("N randomised Exp")

N_Control  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Control) <-c("N randomised Control")


# IE objects
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

for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    
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
    
    re
    
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
    beta_v6_treatment <- -3.5 # this is roughly 1/2 of variance. If 1:1 then 16 to 16, and if 0.5 to 1, then 64 to 64.
    
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
    
    
    #sqrt(var(d$MADRS10[d$visit==42]))
    
    ## flooring and ceiling
    #d[, 4] <- ifelse(d[, 4] < 0, 0,
    #                ifelse(d[, 4]>60, 60, d[, 4]))
    
    
    #View(d)
    d <-as.data.frame(d)
    d$visit <-as.factor(d$visit)
    #d$Treat <- factor(d$Treat)
    
    
    ####################################################################################################  
    # MMRM on full outcome data
    ############################
    
    
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
    #d<-d[d$visit!=0,]
    
    
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
    
    
 
    p<- ggplot(data = d, aes(x = visit, y = MADRS10, group = id)) 
    plot1 <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)) ; plot1
    
    
    
    
    fit<-gls(MADRS10 ~ visit * Treat + Baseline, 
             data=d,
             correlation = corSymm(form=~1 | id),
             #weights = varIdent(form = ~ 1 | visit),
             method="REML")
    

    

    summary(fit)
    
    
    #sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"])
    
    #sqrt(vcov(fit)["Treat", "Treat"])
    
    #sqrt(getVarCov(fit)[1,1]/n)
    
    
    
    re_covm2
    #fit$coefficients[c(7,13)]
    
    #sum(fit$coefficients[c(7,13)]); treatmenteffect
    
    d$visit <- as.numeric(d$visit)-1
    
    #fit_lme <- lme(fixed=MADRS10 ~ visit * Treat + Baseline, 
     #          random=~1 + visit | id,
      #        method="REML", 
       #      correlation = corSymm(form=~1|id),
        #    data=d)
    
    
    
    
    #fit_lmer <- lmer(MADRS10 ~ visit + visit:Treat + (1 |id), data = d, REML = T)
    
    
    #summary(fit_lmer)
    
    #summary(fit_lme)
    
    #model_parameters(fit_lme)
    
    #betas[m, ] <- fit$coefficients[c(7,13)]
    
    #delta[m, ] <- sum(fit$coefficients[c(7,13)])
    
    #bias_f[m, ] <- sum(fit$coefficients[c(7,13)]) - treatmenteffect
    
    #delta_error <- sqrt(vcov(fit)["Treat", "Treat"] + vcov(fit)["visit42:Treat", "visit42:Treat"] + 2*vcov(fit)["Treat", "visit42:Treat"]) 
    
    #delta_errorz[m, ] <- delta_error 
    
    
    #confint_fit[m,1] <- sum(fit$coefficients[c(7,13)])-qnorm(0.975)*delta_error
    #confint_fit[m,2] <- sum(fit$coefficients[c(7,13)])+qnorm(0.975)*delta_error
    
    
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
    #s<-1
    p_LoE <-sample(p_LoE_sample, 1) * scaling_factor[s]
    
    d_mis_w$LoE_Yes <-ifelse(d_mis_w$CfB<5, 1, 0)*rbinom(n, 1, p_LoE);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    sum(d_mis_w$LoE_Yes)
    #View(d_mis_w)
    
    # integrate for AE
    # too much efficacy >8 points on MADRS10 at week 2
    
    p_AE_Exp <- sample(p_AE_Exp_sample, 1) * scaling_factor[s]
    
    d_mis_w$AE_Exp_Yes <-ifelse(d_mis_w$Treat==1 & d_mis_w$CfW2>8, 1, 0)*rbinom(n, 1, p_AE_Exp);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    sum(d_mis_w$AE_Exp_Yes)
    
    p_AE_Control <- sample(p_AE_Control_sample, 1) * scaling_factor[s]
    
    d_mis_w$AE_Control_Yes <-ifelse(d_mis_w$Treat==0 & d_mis_w$CfW2< (-2), 1, 0)*rbinom(n, 1, p_AE_Control);  d_mis_w # # to adjust the probabilty of LoE# create the LoE variable
    sum(d_mis_w$AE_Control_Yes)
    sum(d_mis_w$LoE_Yes)
    
    
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
    describe(d_mis_L$LoE_YES)
    
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
    
    
    #View(d_mis_L)
    
    assign(paste0("SimTrial_sm", "_", n,"_", m, "_", s), d_mis_L)
    #View(SimTrial_sm_1_5)
    
    
    dataset_name.Rdata <- paste0("SimTrial_sm", "_", n,"_", m, "_", s, ".Rdata")
    dataset_name <- paste0("SimTrial_sm", "_", n,"_", m, "_", s)
  
    save(dataset_name, file = dataset_name.Rdata)
  
    
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
    
    #View(d_mis_L)
    
    
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
    #p + geom_line() + facet_grid(~ Treat) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)
    
    describe(d_mis_L$Behavior)
  
    
    # All behaviors
    p<- ggplot(data = d_mis_L, aes(x = Visit, y = MADRS10, group = id, color=Behavior)) 
    #p + geom_line() + facet_grid(~ Treat) 
    plot_all <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)); plot_all
    
    

   
    
    d_mis_L_LoE <- d_mis_L[d_mis_L$LoE_YES==1,] # subset only patients that experienced LoE
    d_mis_L_AE <- d_mis_L[d_mis_L$AE_Yes==1,] # subset only patients that experienced AE
    d_mis_L_NoIE <- d_mis_L[d_mis_L$LoE_YES==0 & d_mis_L$AE_Yes==0,] # subset only patients that did not experience any IE
    

    
    # just LoE patients with true trajectory
    p<- ggplot(data = d_mis_L_LoE, aes(x = Visit, y = MADRS10, group = id))
    plot_LoE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
      scale_y_continuous(limits = c(-10, 60)); plot_LoE
    
    #View(d_mis_L)
    
    

    # just AE patients with true trajectory
    p<- ggplot(data = d_mis_L_AE, aes(x = Visit, y = MADRS10, group = id))
    plot_AE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60)) ; plot_AE
    
    
    # just No IE patients with true trajectory
    p<- ggplot(data = d_mis_L_NoIE, aes(x = Visit, y = MADRS10, group = id))
   plot_NoIE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
     scale_y_continuous(limits = c(-10, 60)); plot_NoIE
    
    
    
    
    
    
    
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
  
## continue to add here the other statistics
  
  
  LoE_Control_Perc
  
  all_delta_errorz_1

  
  ### plot bias
  

  
  setTxtProgressBar(pb3, s)
}





(plot1 / plot_all) /

(plot_LoE / plot_AE / plot_NoIE)



plot_all / plot_LoE / plot_AE / plot_NoIE


(plot_all / plot_all) | (plot_all/plot_all)


end_time <- Sys.time()

end_time-start_time

colMeans(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5))

colMeans(rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5))


cbind(rbind(all_betas_1,all_betas_2, all_betas_3, all_betas_4,all_betas_5),
      rbind(all_delta_1,all_delta_2, all_delta_3, all_delta_4,all_delta_5),
      rbind(all_delta_errorz_1, all_delta_errorz_2, all_delta_errorz_3, all_delta_errorz_4, all_delta_errorz_5),
      rbind(all_confint_fit_1, all_confint_fit_2, all_confint_fit_3, all_confint_fit_4, all_confint_fit_5),
      rbind(all_N_Exp_1, all_N_Exp_2, all_N_Exp_3, all_N_Exp_4, all_N_Exp_5),
      rbind(all_N_Control_1, all_N_Control_2, all_N_Control_3, all_N_Control_4, all_N_Control_5))




# find an efficient way to store all descriptive statistics, betas, deltas, percentages for each intercurrent event



############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################



#### code to extract the models for each pattern in preparation for PMMMM approach

# from a trial with 1000 patients per arm
# Pattern for LoE at trial level
# Pattern for AE in experimental arm
# Pattern for AE in control arm
# Pattern for completers at trial level


# 50% and then use scaling factors to get to 5% the entire range of IE percentages



# Use simulated trial with highest percentage of IE

#View(SimTrial_sm_2000_1_5)
describe(SimTrial_sm_2000_1_5)

# prepare the dataset for reparameterisation to not have any treatment coefficient at baseline

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




head(SimTrial_sm_2000_1_5)

class(SimTrial_sm_2000_1_5$Treat)


fit<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
           Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
         data= SimTrial_sm_2000_1_5,
         correlation = corSymm(form=~1 | id),
         #weights = varIdent(form = ~ 1 | Visit),
         method="REML")

#View(SimTrial_sm_2000_1_5)

summary(fit)
getVarCov(fit, individual = 1)

describe(SimTrial_sm_2000_1_5)


# Pattern for LoE at trial level
fit_LoE_trial<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                     Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
          data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$LoE_YES==1,] ,
          correlation = corSymm(form=~1 | id),
          #weights = varIdent(form = ~ 1 | Visit),
          method="REML")


#View(SimTrial_sm_2000_1_5)


summary(fit_LoE_trial)
getVarCov(fit_LoE_trial, individual = 8)


# Pattern for AE in experimental arm
fit_AE_exp<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42, 
                    data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Exp_Yes==1,] ,
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | Visit),
                    method="REML")



summary(fit_AE_exp)
getVarCov(fit_AE_exp, individual = '2')

#View(SimTrial_1_5[SimTrial_1_5$AE_Exp_Yes==1,])

#corMatrix(fit_AE_exp$modelStruct$corStruct)[[1]]  



# Pattern for AE in control arm
fit_AE_control<-gls(MADRS10 ~ MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                      Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
                    data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Control_Yes==1,] ,
                correlation = corSymm(form=~1 | id),
                weights = varIdent(form = ~ 1 | Visit),
                method="REML")

summary(fit_AE_control)
getVarCov(fit_AE_control, individual = '28')

describe(SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Control_Yes==1,])

#View(SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$AE_Control_Yes==1,])


# Pattern for no intercurrent events at trial level
fit_no_IE<-gls(MADRS10 ~ V7 + V14 + V21 + V28 + V35 + V42 +
                 Treat:V7 + Treat:V14 + Treat:V21 + Treat:V28 + Treat:V35 + Treat:V42, 
               data = SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Behavior== "No IE", ] ,
                    correlation = corSymm(form=~1 | id),
                    weights = varIdent(form = ~ 1 | Visit),
                    method="REML")


summary(fit_no_IE)
getVarCov(fit_no_IE, individual = 3)

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


# code to extract the logit models for LoE at trial level, AE in experimental arm and AE in control arm
### obtain the logistic regression models for LoE at trial level, AE in experimental arm and AE in control arm
class(SimTrial_sm_2000_1_5$AE_Control_Yes)
SimTrial_sm_2000_1_5$AE_Control_Yes <- factor(SimTrial_sm_2000_1_5$AE_Control_Yes)

class(SimTrial_sm_2000_1_5$AE_Exp_Yes)
SimTrial_sm_2000_1_5$AE_Exp_Yes <- factor(SimTrial_sm_2000_1_5$AE_Exp_Yes)


trial_AE_exp <- SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Treat==1,]
trial_AE_control <- SimTrial_sm_2000_1_5[SimTrial_sm_2000_1_5$Treat==0,]



# LoE at trial level logit(Pr(LoE))= fi_LoE * (Y_i0 - Y_i6)
logit_LoE <- glm(LoE_Yes ~ -1 + CfB,
                 data = SimTrial_sm_2000_1_5,
                 family = "binomial")

summary(logit_LoE)

predicted_LoE <- predict(logit_LoE, type="response", newdata=SimTrial_sm_2000_1_5)

SimTrial_sm_2000_1_5$predicted_LoE <- predicted_LoE

View(SimTrial_sm_2000_1_5[,c(5, 20)])

probz_predicted_LoE_Yes <- SimTrial_sm_2000_1_5$predicted_LoE[SimTrial_sm_2000_1_5$LoE_Yes==1]
up_boundary_prob_LoE <- max(probz_predicted_LoE_Yes)  ; up_boundary_prob_LoE
low_boundary_prob_LoE <- min(probz_predicted_LoE_Yes) ; low_boundary_prob_LoE




# AE experimental arm logit(Pr(AE_exp))= fi_AE_exp * (Y_i0 - Y_i2)
logit_AE_exp <- glm(AE_Exp_Yes ~ -1 + CfW2,
                 data = trial_AE_exp,
                 family = "binomial")



summary(logit_AE_exp)

predicted_AE_exp <- predict(logit_AE_exp, type="response", newdata = trial_AE_exp)

trial_AE_exp$predicted_AE_exp <- predicted_AE_exp

#View(trial_AE_exp[,c(6, 20)])
#View(trial_AE_exp)

probz_predicted_AE_exp_Yes <- trial_AE_exp$predicted_AE_exp[trial_AE_exp$AE_Exp_Yes==1]
up_boundary_prob_AE_exp <- max(probz_predicted_AE_exp_Yes)  ; up_boundary_prob_AE_exp
low_boundary_prob_AE_exp <- min(probz_predicted_AE_exp_Yes)  ;  low_boundary_prob_AE_exp





# AE logit(Pr(AE_control))= fi_AE_control * (Y_i0 - Y_i2)


logit_AE_control <- glm(AE_Control_Yes ~ -1 + CfW2,
                    data = trial_AE_control,
                    family = "binomial")


summary(logit_AE_control)

predicted_AE_control <- predict(logit_AE_control, type="response", newdata = trial_AE_control)

trial_AE_control$predicted_AE_control <- predicted_AE_control

#View(trial_AE_control[,c(7, 20)])

probz_predicted_AE_control_Yes <- trial_AE_control$predicted_AE_control[trial_AE_control$AE_Control_Yes==1]
up_boundary_prob_AE_control <- max(probz_predicted_AE_control_Yes)  ; up_boundary_prob_AE_control
low_boundary_prob_AE_control <- min(probz_predicted_AE_control_Yes)  ; low_boundary_prob_AE_control



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


# code for SPM 
# fit SPM model on the simulated trial
# fit the glm to get the intercurrent events models with random effects in the linear predictor


#View(SimTrial_sm_190_1_5)

class(SimTrial_sm_2000_1_5$Visit)
SimTrial_sm_2000_1_5$Visit <- as.numeric(SimTrial_sm_2000_1_5$Visit)-1


fit_lmer <- lmer(MADRS10 ~ Visit + Visit : Treat + (1 + Visit|id), data = SimTrial_sm_2000_1_5, REML = T)
summary(fit_lmer)


fit_lme <- lme(fixed = MADRS10 ~ Visit + Visit : Treat, 
              random = ~ 1 + Visit| id,
              method ="REML", 
               data = SimTrial_sm_2000_1_5)

summary(fit_lme)

vcov(fit_lme)
VarCorr(fit_lme)
random.effects(fit_lme)[,1]
var(random.effects(fit_lme))

# logit(Pr(IE_ij)) = (beta_0 + b0_ij) [aka baseline MADRS] + beta_1 * Treat
# baseline MADRS10 is made of the general intercept + each individual random intercept



d_re <- d_mis_w

View(d_re)

fit_LoE_spm <- glm(LoE_Yes ~ -1 + Baseline + Treat,
                   data = d_re, 
                   family = "binomial")

summary(fit_LoE_spm)


fit_AE_exp_spm <- glm(AE_Exp_Yes ~ -1 + Baseline,
                   data = d_re[d_re$Treat==1,], 
                   family = "binomial")

summary(fit_AE_exp_spm)



fit_AE_control_spm <- glm(AE_Control_Yes ~ -1 + Baseline,
                      data = d_re[d_re$Treat==0,], 
                      family = "binomial")

summary(fit_AE_control_spm)


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~













# Selection model via marginal model for outcomes-generating model and FULLY STOCHASTIC models for generation of intercurrent events


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
m.iterations <- 1# number of generated datasets # number of trials per scaling factor
scaling_factor <-  c(0.5, 1.0, 1.5, 2.0, 2.5)
# total number of simulated trials = m.iterations * length(scaling_factor)
# try with c(0.4, 1.1, 1.8, 2.3, 3)
#c(0.25, 0.5, 1, 2, 2.5) steps for scaling factor
# to multiply for LoE, AE_control and AE_exp ### value 0.5 should yield around 13.5% IEs, 1 ~ %, 2 ~ %, 2.5 ~



n <- 190# number of patients

CFE <- matrix(ncol=4,nrow=length(scaling_factor)*m.iterations)
colnames(CFE) <-c("N ceiled_floored", "% ceiled_floored", "scaling factor", "simulated trial n")


p_LoE_sample <-c(0.31, 0.33, 0.34, 0.35, 0.37); mean(p_LoE_sample)
p_AE_Exp_sample <- c(0.055, 0.065, 0.075, 0.085, 0.095); mean(p_AE_Exp_sample)
p_AE_Control_sample <- c(0.05, 0.10, 0.15, 0.20, 0.25); mean(p_AE_Control_sample)

visits <- as.numeric(c(0, 7, 14, 21, 28, 35, 42)) # number of measurements, baseline + follow-up measurements
delta <- matrix(ncol=1,nrow=m.iterations) # treatment effect estimate at 6 weeks based on MMRM models fitted on each generated dataset
colnames(delta) <-c("TreatmentEffect")
betas <- matrix(ncol=2,nrow=m.iterations)
colnames(betas) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3)

confint_fit <- matrix(ncol=2,nrow=m.iterations) 
colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
delta_errorz <- matrix(ncol=1,nrow=m.iterations)
colnames(delta_errorz) <- c("SE")

bias_f <- matrix(ncol=1,nrow=m.iterations)
colnames(bias_f) <- c("bias_f")

## randomisation objects, allocation
N_Exp  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Exp) <-c("N randomised Exp")

N_Control  <- matrix(ncol=1,nrow=m.iterations)
colnames(N_Control) <-c("N randomised Control")


# IE objects
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

for (s in 1:length(scaling_factor)) {
  for(m in 1:m.iterations) {
    
    
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
      Treat = rep(cccccc, each = length(visits)),
      MADRS10 = rep(NA, n)); d # mean(Treat)
    
    
    d <- d[order(d$visit, d$id),]; #d
    #re
    
    j<-c(re[, 1], re[, 2], re[, 3], re[, 4], re[, 5], re[, 6], re[, 7])
    d$re <-j; #d
    
    #class(d)
    #View(d)
    
    #head(re)
    
    d <- d[order(d$id, d$visit),]; #d
    
    re
    
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
    
    beta_v6_treatment <- -3.5 # this is roughly 1/2 of variance. If 1:1 then 16 to 16, and if 0.5 to 1, then 64 to 64.
    
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
    
    
    ####################################################################################################  
    # MMRM on full outcome data
    ############################
    
    
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
    # Generate Intercurrent events
    
    # use propensity models derived from original trial or from the dataset above 
    # one model for LoE for all patients
    # one model for AE in experimental arm
    # one model for AE in control arm
    
    
    # see paper and use that wording
    
    # if Difference < 5, then missing outcome from visit 3 onwards for LoE
    
    ## take over the original/raw dataset to use for IEGM/MDGM
    d_mis <-d_orig
    
    d_mis_w <- d_mis %>% spread(visit, MADRS10) # reshape to wide in order to create the CfB variable
    #View(d_mis_w)
    
    colnames(d_mis_w)[3:9] <- c("Baseline", "Week1", "Week2", "Week3", "Week4","Week5" ,"Week6"); head(d_mis_w)
    
    d_mis_w$CfB <- d_mis_w[,3] - d_mis_w[,9]; d_mis_w # create the CfB variable
    
    d_mis_w$CfW2 <- d_mis_w[,3] - d_mis_w[,5]; d_mis_w # create the CfW2 variable
    
    
    
    
    
    
    
    
    
    d_mis_L <- d_mis_w %>% gather(Visit, MADRS10, Baseline:Week6) # reshape to long format
    
    d_mis_L <- d_mis_L[order(d_mis_L$id, d_mis_L$Visit),]; #d # order by subject id and Visit
    
    
    
    
    d_mis_L$predicted_LoE <- predict(logit_LoE, type="response", newdata=d_mis_L)
    
    #describe(predict(logit_LoE, type="response", newdata=d_mis_L))
    
    d_mis_L$LoE_yes <- ifelse(d_mis_L$predicted_LoE < up_boundary_prob_LoE & d_mis_L$predicted_LoE > low_boundary_prob_LoE, 1, 0)
    d_mis_L$CfB
    
    
    
    trial_AE_X <- d_mis_L[d_mis_L$Treat==1,]
    

  
    trial_AE_X$predicted_AE <- predict(logit_AE_exp, type="response", newdata = trial_AE_X)
    
    #describe(predict(logit_LoE, type="response", newdata=d_mis_L))
    
    trial_AE_X$AE_yes <- ifelse(trial_AE_X$predicted_AE  < 0.28 & trial_AE_X$predicted_AE  > low_boundary_prob_AE_exp, 1, 0)
    
    trial_AE_X$CfW2
    
      
    
    
    
    trial_AE_C <- d_mis_L[d_mis_L$Treat==0,]
    
    
    
    trial_AE_C$predicted_AE <- predict(logit_AE_control, type="response", newdata = trial_AE_C)
    
    
    trial_AE_C$AE_yes <- ifelse(trial_AE_C$predicted_AE  < up_boundary_prob_AE_control & trial_AE_C$predicted_AE  > low_boundary_prob_AE_control, 1, 0)
    
    trial_AE_C$CfW2
    
    


  
    # AE and LoE in experimental arm
    trial_AE_X$AE_YES <- ifelse(trial_AE_X[,10]==1, 1, 0)
    trial_AE_X$LoE_YES <- ifelse(trial_AE_X[,10]==0 & trial_AE_X[,8]==1 , 1, 0)
    
    
    # AE and LoE in control arm
    trial_AE_C$LoE_YES <- ifelse(trial_AE_C[,8]==1, 1, 0)
    trial_AE_C$AE_YES <- ifelse(trial_AE_C[,8]==0 & trial_AE_C[,10]==1, 1, 0)

    
    
    
    
    
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
    range(d_mis_LL$MADRS10[d_mis_LL$Treat==1], na.rm = T)
    range(d_mis_LL$MADRS10[d_mis_LL$Treat==0], na.rm = T)
    
    
    
    describe(d_mis_LL$AE_YES[d_mis_LL$Treat==0])
    describe(d_mis_LL$AE_YES[d_mis_LL$Treat==1])
    describe(d_mis_LL$AE_YES)
    
    describe(d_mis_LL$Behavior[d_mis_LL$Treat==0])
    describe(d_mis_LL$Behavior[d_mis_LL$Treat==1])
    describe(d_mis_LL$Behavior)
    
    
    describe(d_mis_LL$LoE_YES[d_mis_LL$Treat==0])
    describe(d_mis_LL$LoE_YES[d_mis_LL$Treat==1])
    describe(d_mis_LL$LoE_YES)
    
    
    View(d_mis_LL)
    
    
    
    

    
    #LoE_Y <- d_mis_L[,c(2, 11)]
    #LoE_Y$LoE_YES <- as.numeric(LoE_Y$LoE_YES)-1
    
    
    #LoE_Y_total <- sum(LoE_Y$LoE_YES)/length(visits)
    #LoE_Y_Exp <- sum(LoE_Y$LoE_YES[LoE_Y$Treat==1])/length(visits)
    LoE_Y_Control <- sum(LoE_Y$LoE_YES[LoE_Y$Treat==0])/length(visits)
    
    
    ## LoE
    #tb_LoE_total <- LoE_Y_total
    #tb_LoE_Exp <- LoE_Y_Exp
    #tb_LoE_Control <- LoE_Y_Control
    
    
    #n_LoE_total[m, ] <-  tb_LoE_total
    
    #n_LoE_Exp[m, ] <- tb_LoE_Exp
    
    #n_LoE_Control[m, ] <- tb_LoE_Control
    
    ###
    
    
    #do with colsums instead of table to avoid Error subscript out of bounds when there is no AE generated
    
    
    #AE_Y <- d_mis_L[,c(2, 10)]
    #AE_Y$AE_Yes <- as.numeric(AE_Y$AE_Yes)-1
    
    
    #AE_Y_total <- sum(AE_Y$AE_Yes)/length(visits)
    #AE_Y_Exp <- sum(AE_Y$AE_Yes[AE_Y$Treat==1])/length(visits)
    #AE_Y_Control <- sum(AE_Y$AE_Yes[AE_Y$Treat==0])/length(visits)
    
    
    
    ## AE
    #tb_AE_total <- AE_Y_total
    #tb_AE_Exp <- AE_Y_Exp
    #tb_AE_Control <- AE_Y_Control
    
    
    #n_AE_total[m, ] <-  tb_AE_total
    
    
    #n_AE_Exp[m, ] <- tb_AE_Exp
    
    
    #n_AE_Control[m, ] <- tb_AE_Control
    
    
    
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
    
    
    
    
    
    
    
    d_mis_L_LoE <- d_mis_L[d_mis_L$LoE_YES==1,] # subset only patients that experienced LoE
    d_mis_L_AE <- d_mis_L[d_mis_L$AE_Yes==1,] # subset only patients that experienced AE
    d_mis_L_NoIE <- d_mis_L[d_mis_L$LoE_YES==0 & d_mis_L$AE_Yes==0,] # subset only patients that did not experience any IE
    
  #  View(d_mis_L)
    
    
    # All patients with true trajectory with different colours by LoE (Y/N)
    p<- ggplot(data = d_mis_LL, aes(x = Visit, y = MADRS10, group = id, color=LoE_YES))
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


# resume here
# find an efficient way to store each generated dataset, outcomes and intercurrent events
# find an efficient way to store all descriptive statistics, betas, deltas, percentages for each intercurrent event


