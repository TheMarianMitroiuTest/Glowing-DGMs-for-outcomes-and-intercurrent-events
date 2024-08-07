## License CC-BY-4.0 
## Creative Commons Attribution 4.0 International
## Mitroiu M, Teerenstra S, Oude Rengerink K, Pétavy F, Roes K. Possible data-generating models of longitudinal continuous outcomes and intercurrent events to investigate estimands. Statistics in Biopharmaceutical Research. 
## https://doi.org/10.1080/19466315.2024.2369266


## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"



# add the ceiling to make sure the time to event are at the next visit.
# change 0 to Inf for those that do not experience the intercurrent event.
### DGM
#rm(list=ls())
library(gmailr)
library(MASS)
library(tidyverse)
library(nlme)
library(lme4)
#library(Hmisc) #install.packages("xfun", type="binary")
library(janitor)
library(JM)
library(survival)
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

set.seed(2147483629) # set seed for reproducibility
n <- 190# number of patients to be simulated (sample size)
# this is based on a t-test to ensure  90% power at alpha level=0.025 one-sided 
m.iterations <- 500 # 416 as per the LMM also used in the SPM DGM to verify the longitudinal outcomes

# proportions of Intercurrent events. These correspond to the qt values for the standardisation of time to IE
# for 0.35 LoE
p_LoE_Control_sample <- 0.21
p_LoE_Exp_sample <- 0.14

# for 0.15 AE
#for 0.1 AE exp arm
p_AE_Exp_sample <- 0.10 
#for 0.05 AE control arm
p_AE_Control_sample <- 0.05

treatmenteffect <- -3.5

adjustment.fctr <- 1.15

visits <- as.numeric(c(0, 1, 2, 3, 4, 5, 6))	# scheduled visits

delta_jm <- matrix(ncol=1,nrow=m.iterations) # object to store treatment effect estimates at 6 weeks based on MMRM model fitted on each generated dataset
colnames(delta_jm) <-c("TreatmentEffect")
betas_jm <- matrix(ncol=2,nrow=m.iterations) # object to store parameters for the treatment effect at week 6 based on the MMRM model fitted on each generated dataset
colnames(betas_jm) <-c("Treat", "visit42:Treat")

pb1 <- txtProgressBar(min = 0,  max=m.iterations, style=3) # progress bar in percentages relative to the total number of m.iterations

#confint_fit <- matrix(ncol=2,nrow=m.iterations) # object to store the 95% confidence interval bounds for the estimated treatment effect
#colnames(confint_fit) <-c("Lower boundary 95% CI", "Upper boundary 95% CI")
#delta_errorz <- matrix(ncol=1,nrow=m.iterations) # standard error of the estimated treatment effect
#colnames(delta_errorz) <- c("SE")

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

start_time <- Sys.time() # timestamp for the start time of the nested for loop below.
# it was used to have an estimate of time needed for different larger number of trials to be simulated upon scaling up the simulation parameters (e.g., m.iterations)



## Begin for loop----

for (m in 1:m.iterations) {
  
  
  # The specification of the model: random intercept and random slope + some random noise
  # d$MADRS10.true <- (b0 + d$bi_0) + (b1 + d$bi_1) * d$visit +  b2 * d$visit * d$Treat
  # d$MADRS10_collected <- d$MADRS10.true + rnorm(nrow(d), 0, eps.sd)
  
  
  ### Generate longitudinal outcomes----
  #### Generate correlated (random) intercepts and slopes for each individual patient----
  
  # We used as inspiration the trial 003-002 (ref data analysis paper). Here we used a simplified scenario with linear group trajectories.
  
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

fit_lme <- lme(fixed=MADRS10_collected ~ visit + visit:Treat, 
               random=~1 + visit| id,
               method="REML", 
               data=d)

#std.e <- 0.1749478*length(visits)

summary(fit_lme)
#m<-1
betas_jm[m, ] <- fixef(fit_lme)[2:3]
delta_jm[m, ] <- fixef(fit_lme)[3] * max(visits)


#m<-1
# store the number of patients in the objects defined a priori
Randomised_Exp <- sum(as.numeric(d[,3]))/7 #number of patients in the experimental arm # make sure it is divided by the number of visits (6/7), without or with the baseline included as response
Randomised_Control <- n-sum(as.numeric(d[,3]))/7  #number of patients in the control arm

N_Exp[m,] <- Randomised_Exp
N_Control[m,] <- Randomised_Control




#################################################################
# Weibull model and survival for LoE at trial level and within arms #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


intercept_LoE <- 2.8789#3.3639
c_LoE <- 0.3111#0.1616 #-0.5# coefficient for the Treatment in the linear predictor that contains also the random effects, to be used in the generation of the time to intercurrent event data
Alpha_LoE <- -0.0279#-0.0447
# parameters of the Weibull distributions
# To fit with the assumptions, a trial-and-error/finetuning process can be employed, to find/generate the distribution that fits with the targeted survival data to be generated.

lambda_LoE 	<- 1.7618#1.8329 #3.5			# scale parameter
nu_LoE 		<- 	1.7618#1.8329 #1.4		# shape parameter
# other distributions could be used (e.g., exponential) to describe the time to intercurrent event distribution.

# linear predictor
LP_LoE <- (intercept_LoE + d$bi_0 * Alpha_LoE) +  c_LoE * d$Treat # this can be used to generate time to intercurrent events with differing durations up to intercurrent event between arms.
# e.g., in the control arm the LoE will appear (slightly) earlier than in the experimental/treatment arm 
    #describe(LP1)

    #hist(d$bi_1)

t.event_LoE <- (-log(rep(runif(length(unique(d$id))), each=length(visits)))/(lambda_LoE * exp(LP_LoE))) ^ (1/nu_LoE) ; t.event_LoE
    #hist(t.event_LoE)
    #describe(t.event_LoE)



#t.event_LoE <- ceiling(t.event_LoE*(5/quantile(t.event_LoE, probs = p_LoE_sample*adjustment.fctr)[[1]]) + 1) # standardize the time to event (to fit with the trial duration)
# the standardisation could be to fit tte to the entire trial duration, or to fit it to be at specific visits. We standardise it to fit the entire trial duration (6 weeks) and to have most of the intercurrent events up to and including week 4
# Other Weibull distributions (parameters) can be used to better describe the wanted time to events
# Depending on the desired percentages of intercurrent events and how they can be achieved. The standardisation can be used such that only a certain percentage of patients experience
# the intercurrent event during the trial (e.g., standardise it to fit to longer than end of trial (week 6), or it can be standardised to week 6 and then use a Binomial distribution to
# achieve a certain percentage

#describe(t.event_LoE)
d$time.LoE <- t.event_LoE
#hist(d$t.LoE)



d$t.LoE[d$Treat==0] <- ceiling(d$time.LoE[d$Treat==0]*(5/quantile(d$time.LoE[d$Treat==0], probs = p_LoE_Control_sample*2*adjustment.fctr)[[1]]) + 1) # standardize the time to event (to fit with the trial duration)
d$t.LoE[d$Treat==1] <- ceiling(d$time.LoE[d$Treat==1]*(5/quantile(d$time.LoE[d$Treat==1], probs = p_LoE_Exp_sample*2*adjustment.fctr)[[1]]) + 1) # standardize the time to event (to fit with the trial duration)
# the probabilities are multiplied by 2 to adjust for the trial size and trial level statistics, the trial is double in size than each arm
# the probabilities are multiplied by the adjustment factor to account for patients that experience also AE and AE have priority vs LoE.

#describe(d$t.LoE)
hist(d$t.LoE)
hist(d$t.LoE[d$Treat==0])
hist(d$t.LoE[d$Treat==1])

#d$t.AE <- d$t.AE * rep(rbinom(n, 1, 0.5), each = length(visits))
# this particular step could be also used piecewise with different Binomial probabilities for time intervals, while keeping the intercurrent event percentage similar at  trial level.


#d$t.LoE <- d$t.LoE * rep(rbinom(n, 1, 0.5), each = length(visits))
# this particular step could be also used piecewise with different Binomial probabilities for time intervals, while keeping the intercurrent event percentage similar at  trial level.

    #cbind(d$t.LoE, rep(rbinom(n, 1, 0.5), each = length(visits)))

d$LoE_yes <- ifelse(d$t.LoE <=6 & d$t.LoE!=0, 1, 0)

#d$LoE_yes <-factor(ifelse(d$t.LoE !=0, 1, 0))
#describe(d$LoE_yes)

#d$t.LoE[d$t.LoE==0] <- c("No LoE")

#describe(d$t.LoE)

#cbind(d$Treat,d$bi_0, d$bi_1, t.event_LoE, d$t.LoE)

    #View(d[,c(1, 3, 4, 5, 8, 9)])
    #View(d)

#describe(t.event_LoE[d$Treat==1])
#describe(t.event_LoE[d$Treat==0])
#View(d)
head(d)


#d$t.AE <- d$AE_yes <- 0

# just LoE
p<- ggplot(data = d[d$LoE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=factor(LoE_yes))) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM-LoE pattern")


# LoE
p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=factor(LoE_yes))) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM-LoE all")

#Weibull model for Adverse events in both arms. Standardisation after generating the survival times.

intercept_AE <- 3.7550#4.0580
c_AE <- -0.5394#-0.7469 #-0.5# coefficient for the Treatment in the linear predictor that contains also the random effects, to be used in the generation of the time to intercurrent event data
Alpha_AE <- 0.0076#-0.0014
# parameters of the Weibull distributions
# To fit with the assumptions, a trial-and-error/finetuning process can be employed, to find/generate the distribution that fits with the targeted survival data to be generated.

lambda_AE 	<- 0.9719#0.9785 #3.5			# scale parameter
nu_AE 		<- 	0.9719#0.9785 #1.4		# shape parameter
# other distributions could be used (e.g., exponential) to describe the time to intercurrent event distribution.

# linear predictor
LP_AE <- (intercept_AE + d$bi_0 * Alpha_AE) +  c_AE * d$Treat # this can be used to generate time to intercurrent events with differing durations up to intercurrent event between arms.
# e.g., in the control arm the LoE will appear (slightly) earlier than in the experimental/treatment arm 
#describe(LP1)

#hist(d$bi_1)

t.event_AE <- (-log(rep(runif(length(unique(d$id))), each=length(visits)))/(lambda_LoE * exp(LP_AE))) ^ (1/nu_AE) ; t.event_AE
#hist(t.event_AE)
#describe(t.event_AE)

#describe(t.event_AE)
d$time.AE <- t.event_AE

d$t.AE[d$Treat==0] <- ceiling(d$time.AE[d$Treat==0]*(5/quantile(d$time.AE[d$Treat==0], probs = p_AE_Control_sample*2)[[1]]) + 1) # standardize the time to event (to fit with the trial duration)
d$t.AE[d$Treat==1] <- ceiling(d$time.AE[d$Treat==1]*(5/quantile(d$time.AE[d$Treat==1], probs = p_AE_Exp_sample*2)[[1]]) + 1) # standardize the time to event (to fit with the trial duration)
# the probabilities are multiplied by 2 to adjust for the trial size and trial level statistics, the trial is double in size than each arm

#describe(d$t.AE)
hist(d$t.AE)
hist(d$t.AE[d$Treat==0])
hist(d$t.AE[d$Treat==1])

#d$t.AE <- d$t.AE * rep(rbinom(n, 1, 0.5), each = length(visits))
# this particular step could be also used piecewise with different Binomial probabilities for time intervals, while keeping the intercurrent event percentage similar at  trial level.

#cbind(d$t.AE, rep(rbinom(n, 1, 0.5), each = length(visits)))

d$AE_yes <- ifelse(d$t.AE <=6 & d$t.AE!=0, 1, 0)

#d$AE_yes <-factor(ifelse(d$t.AE !=0, 1, 0))
#describe(d$AE_yes)

#d$t.AE[d$t.AE==0] <- c("No AE")

#describe(d$t.AE)

#describe(t.event_AE[d$Treat==1])
#describe(t.event_AE[d$Treat==0])
#View(d)
head(d)


# just AE all
p<- ggplot(data = d[d$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=factor(AE_yes))) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM-AE pattern")


# AE all
p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=factor(AE_yes))) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM-AE all")


d_united <- d


class(d_united$visit)
class(d_united$Treat)
d_united$Treat <- factor(d_united$Treat)
d_united$LoE_yes <- factor(d_united$LoE_yes)
d_united$AE_yes <- factor(d_united$AE_yes)

#describe(d_united$AE_yes)
d_united$AE_YES <- ifelse(d_united$AE_yes==1, 1, 0)
d_united$LoE_YES <- ifelse(d_united$AE_yes==0 & d_united$LoE_yes==1, 1, 0)


d_united$Behavior <- ifelse(d_united$AE_YES==1, "AE",
                            ifelse(d_united$AE_YES==0 & d_united$LoE_YES==1, "LoE", "No IE"))

#View(d_united)

#describe(d_united$Behavior)


#### assign and save the generated datasets----
# naming sequence is "SimTrial"_"Method"_"trial sample size"_"iteration number"

assign(paste0("SimTrial_jm", "_", n), d)
#View(SimTrial_sm_1_5)
dataset_name.Rdata <- paste0("SimTrial_jm", "_", n, ".Rdata")
dataset_name <- paste0("SimTrial_jm", "_", n)
save(dataset_name, file = dataset_name.Rdata)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Intercurrent events descriptives needed for the verification step ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LoE_Y <- d_united[ ,c("Treat", "LoE_YES")]
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

AE_Y <- d_united[ ,c("Treat", "AE_YES")]
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




p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
  scale_y_continuous(limits = c(-10, 60))

d_united$visit <- as.factor(d_united$visit)

# Plot trajectories

# just LoE
# LoE
p<- ggplot(data = d_united[d_united$LoE_YES==1,], aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_LoE_JM <- p + geom_line(size=0.5, color='#00BA38') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="dark green") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60)) + ggtitle("JM - LoE pattern"); plot_LoE_JM

#describe(d_united[d_united$LoE_YES==1,])

#just AE
# AE
p<- ggplot(data = d_united[d_united$AE_YES==1,], aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_AE_JM <- p + geom_line(size=0.5, color='#F8766D') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM - AE pattern"); plot_AE_JM

#just No IE
# No IE
p<- ggplot(data = d_united[d_united$Behavior=="No IE",], aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_NoIE_JM <- p + geom_line(size=0.5, color='#619CFF') + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="blue") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM - No IE pattern"); plot_NoIE_JM

#View(d_united[d_united$LoE_YES==1,])

#describe(d_united$Behavior)

# All behaviors
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_all_JM <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="black") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60)) + ggtitle("JM - All patterns"); plot_all_JM

the_plot_JM <- (plot_all_JM / plot_LoE_JM) | (plot_AE_JM / plot_NoIE_JM); the_plot_JM


setTxtProgressBar(pb1, m)
}


end_time <- Sys.time()

end_time-start_time



betas_jm; 
colMeans(delta_jm); treatmenteffect

tolerance_margin <- 0.1 
difference_Verification_jm <- abs(treatmenteffect - colMeans(delta_jm))

# check if the result satisfies the inequality
ifelse(isTRUE(paste(difference_Verification_jm) < tolerance_margin), "Verification JM *SUCCESSFUL*", "Verification JM NOT successful :(") 





# Table for the paper ----

table_AE_JM <- data.frame(
  # descriptives AE  
  n_AE_Control,
  n_AE_Exp); table_AE_JM

mean(n_AE_Control)
mean(n_AE_Exp)


# descriptives LoE  
table_LoE_JM <-data.frame(
  n_LoE_Control,
  n_LoE_Exp); table_LoE_JM

mean(n_LoE_Control)
mean(n_LoE_Exp)


#describe(table_IE_JM)


table_AE_JM |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "AE") |> 
  rename(N_C_arm=N.AE.Control) |> 
  rename(N_E_arm=N.AE.Exp)

table_LoE_JM |> 
  as.data.frame() |> 
  mutate("Intercurrent event" = "LoE") |> 
  rename(N_C_arm=N.LoE.Control) |> 
  rename(N_E_arm=N.LoE.Exp)


tab_JM <- tibble(bind_rows(table_AE_JM |> 
                             as.data.frame() |> 
                             mutate("Intercurrent event" = "AE") |> 
                             rename(N_C_arm=N.AE.Control) |> 
                             rename(N_E_arm=N.AE.Exp), 
                           table_LoE_JM |> 
                             as.data.frame() |> 
                             mutate("Intercurrent event" = "LoE") |> 
                             rename(N_C_arm=N.LoE.Control) |> 
                             rename(N_E_arm=N.LoE.Exp))); tab_JM



tab2_JM <- tab_JM |> group_by(`Intercurrent event`) |>
  summarise("N" = round(mean(N_C_arm), digits=1), 
            "%" = round(mean(N_C_arm/n*100), digits=1),
            "N " = round(mean(N_E_arm), digits=1), 
            "% " = round(mean(N_E_arm/n*100), digits=1),
            " N " = round(mean(N_C_arm + N_E_arm), digits=1),
            " % " = round(mean(N_C_arm + N_E_arm)/n*100, digits = 1)) |> 
  adorn_totals("row"); tab2_JM




gt(tab2_JM) |> 
  tab_header(title = md("Table 8d. Descriptive statistics intercurrent events"), subtitle = md("Joint Model DGM")) |>
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
# 6 February 2022
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




# PLOT EVERYTHING----
## ONLY AFTER ALL four DGMs have been ran
## Plots for the paper after running all DGMs
## LoE patterns side by side by each DGM
(plot_LoE_SMd + plot_LoE_PMMM) / (plot_LoE_SPM + plot_LoE_JM)

### AE patterns side by side by each DGM
(plot_AE_SMd + plot_AE_PMMM) / (plot_AE_SPM + plot_AE_JM)


### No IE patterns side by side by each DGM
(plot_NoIE_SMd + plot_NoIE_PMMM) / (plot_NoIE_SPM + plot_NoIE_JM)


### All patterns in a trial side by side by each DGM
(plot_all_SMd + plot_all_PMMM) / (plot_all_SPM + plot_all_JM)



# determine the number of trials needed to simulate for the verification of the longitudinal outcomes
# the underlying linear mixed effects model used in the JM is the exact same linear mixed effects model from the SPM DGM
# please see details there.
# the verification was successful, obviously, as in SPM.

