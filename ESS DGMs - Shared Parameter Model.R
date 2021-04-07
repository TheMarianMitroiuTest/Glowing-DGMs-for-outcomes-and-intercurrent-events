## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"

#

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

    set.seed(2147483629)
    b0 <- 29.5
    b1 <- -0.55
    b2 <- -0.583
    bi_means <- c(0, 0)
    bi_covm <- matrix(c(24.611, 0.5869809, 0.5869809, 1.157), nrow = 2)
    bi <- mvrnorm(n, bi_means, bi_covm)	
    eps.sd <-3.247  
    
    c1 <- 0.1
    
    visits <- as.numeric(c(0, 1, 2, 3, 4, 5, 6))	
    
    d <- data.frame(
      id = rep(1:n, each = length(visits)),
      visit = visits,
      Treat = rep(rbinom(n, 1, 0.5), each = length(visits)),
      bi_0 = rep(bi[,1], each = length(visits)),
      bi_1 = rep(bi[,2], each = length(visits))
    )
    
    
    d$MADRS10.true <- (b0 + d$bi_0) + (b1 + d$bi_1) * d$visit +  b2 * d$visit * d$Treat
    
    d$MADRS10_collected <- d$MADRS10.true + rnorm(nrow(d), 0, eps.sd)
    
    #####################################################
    # Logit model and Probabilities for LoE at trial level
    logit_Pr_LoE <- (b0 + d$bi_1)/1000 +  c1 * d$Treat

    Pr_LoE <- 1/(1+exp(-logit_Pr_LoE))
    
    describe(Pr_LoE[d$Treat==1])
    describe(Pr_LoE[d$Treat==0])
    
    d$Pr_LoE <- Pr_LoE
    
    d$LoE_yes[d$Treat==1] <- ifelse(Pr_LoE[d$Treat==1]>0.5331, 1, 0) * rep(rbinom(length(unique(d$id[d$Treat==1])), 1, 0.7), each = length(visits))
    d$LoE_yes[d$Treat==0] <- ifelse(Pr_LoE[d$Treat==0]>0.5081, 1, 0) * rep(rbinom(length(unique(d$id[d$Treat==0])), 1, 0.7), each = length(visits))
  
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
    
    
    d$AE_yes[d$Treat==1] <- ifelse(Pr_AE_exp<0.2648, 1, 0)
    
    
    
    
    
    #####################################################
    # Logit model and Probabilities for AE in control arm
    logit_Pr_AE_control <- (b1 + d$bi_1[d$Treat==0])*2 # for the subset of control arm patients
    
    Pr_AE_control <- 1/(1+exp(-logit_Pr_AE_control))
    
    describe(Pr_AE_control)
    hist(Pr_AE_control)
    d$Pr_AE_control[d$Treat==0] <- Pr_AE_control
    
    
    d$AE_yes[d$Treat==0] <- ifelse(Pr_AE_control>0.38, 1, 0) # add rbinom probabilities
    
    

    #plot(1/(1+exp(seq(-5, 5, 1))))
    #plot(seq(-5, 5, 1))

    
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
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    
    #just AE
    # AE
    p<- ggplot(data = d[d$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    
    
    
    
    # All patients with TRUE trajectory 
    # LoE
    p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    
    
    # AE
    p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    
    
    # All behaviors
    p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
    p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
      scale_y_continuous(limits = c(-10, 60))
    
    
describe(d)
    
    
    
    
    fit_lmer <- lmer(MADRS10_collected ~ visit + visit:Treat + (1 + visit|id), data = d, REML = T)
    
    summary(fit_lmer)
    
    
    
    
    
    fit_lme <- lme(fixed=MADRS10_collected ~ visit + visit:Treat, 
                   random=~1 + visit| id,
                   method="REML", 
                   data=d)
    
    summary(fit_lme)
    
    getVarCov(fit_lme, type= "marginal")
    
    
    

    
    
    
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
    
    
    # resume here
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
