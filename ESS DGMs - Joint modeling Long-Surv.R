## License CC-BY-4.0 
## Creative Commons Attribution 4.0 International
## Mitroiu M. et al SMMR


## outcomes generates as responses, and analysed with baseline as fixed-effect (as covariate), not as response
# function for DGM, IEGM and analyses methods
#####


#### Simulation study 
#### Scenario A "early separation and treatment effect maintained"



# add the ceiling to make sure the time to event are at the next visit.
# change 0 to Inf for those that do not experience the intercurrent event.
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
#remotes::install_github("r-lib/styler")
#library(styler)




sessionInfo()
installed.packages()


# Selection model via marginal model for outcomes-generating model and deterministic rules for generation of intercurrent events
# Setup to receive e-mails with results of simulations, very useful when running multiple simulations in parallel.

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


n <- 190# number of patients to be simulated (sample size)
# this is based on a t-test to ensure  90% power at alpha level=0.025 one-sided 



set.seed(2147483629) # set seed for reproducibility

# linear mixed effects model parameters to generate the longitudinal outcomes
b0 <- 29.79
b1 <- -0.55
b2 <- -0.583
bi_means <- c(0, 0)
bi_covm <- matrix(c(24.611, 0.5869809, 0.5869809, 1.157), nrow = 2)
bi <- mvrnorm(n, bi_means, bi_covm)	
eps.sd <-3.247 




visits <- as.numeric(c(0, 1, 2, 3, 4, 5, 6))	# protocolled visits

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





#####################################################
# Cox model and survival for LoE at trial level #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #describe(d$bi_0)
    #describe(d$bi_1)

    #hist(d$bi_0)
    #hist(d$bi_1)

    #par(mfrow=c(1,1))

    #hist(d$bi_0[d$Treat==0])
    #hist(d$bi_0[d$Treat==1])

    #hist(d$bi_1[d$Treat==0])
    #hist(d$bi_1[d$Treat==1])
    

    # This chunk of code is only needed if a specific conditioning is made for the intercurrent events on specific values of the intercepts and/or slopes to define the pattern.  
    # e.g., if for instance, the LoE will be experienced by patients with positive slopes, or the AE will be experienced by patients with slopes smaller than a certain (negative) value to reflect the patterns.
    # this would be an implementation with more assumptions, not necessarily bad/wrong, but just with more assumptions about the propensity of a certain IE based on random effects and certain values.
    #subset on positive slopes
    #hist(d$bi_1[d$Treat==1 & d$bi_1>0])
    #d_LoE <- d[d$bi_1>0,]
    #d_no_LoE <- d[d$bi_1<=0,]
    #d_no_LoE$t.LoE <- d_no_LoE$LoE_yes <- 0

c1 <- -0.5 # coefficient for the Treatment in the linear predictor that contains also the random effects, to be used in the generation of the time to intercurrent event data

# parameters of the Weibull distributions
# To fit with the assumptions, a trial-and-error/finetuning process can be employed, to find/generate the distribution that fits with the targeted survival data to be generated.

lambda_LoE 	<- 3.5			# scale parameter
nu_LoE 		<- 	1.4		# shape parameter
# other distributions could be used (e.g., exponential) to describe the time to intercurrent event distribution.

# linear predictor
LP1 <- (d$bi_0 + d$bi_1)/100 +  c1 * d$Treat # this can be used to generate time to intercurrent events with differing durations up to intercurrent event between arms.
# e.g., in the control arm the LoE will appear (slightly) earlier than in the experimental/treatment arm 
    #describe(LP1)

    #hist(d$bi_1)

t.event_LoE <- (-log(rep(runif(length(unique(d$id))), each=length(visits)))/(lambda_LoE * exp(LP1))) ^ (1/nu_LoE) ; t.event_LoE
    #hist(t.event_LoE)
    #describe(t.event_LoE)

#d$t.event_LoE <- t.event_LoE
#hist(t.event_LoE)

# Cumulative hazard function 
#H_0 <- lambda_LoE * t.event_LoE^nu_LoE
#hist(1-H_0)
#hist(-log(H_0))


# hazard function
#h_0 <- lambda_LoE*nu_LoE * t.event_LoE^(nu_LoE-1)
#hist(h_0)


describe(t.event_LoE)




min(t.event_LoE)
median(t.event_LoE)
max(t.event_LoE)

t.event_LoE <- round(t.event_LoE*(5/quantile(t.event_LoE, probs = c(0.50))[[1]]) + 1 , digits=0)
describe(t.event_LoE*(5/quantile(t.event_LoE, probs = c(0.50))[[1]]) + 1)
hist(t.event_LoE)
describe(t.event_LoE)

median(t.event_LoE)
quantile(t.event_LoE, probs = c(0.50))[[1]]

quantile(t.event_LoE, probs = c(seq(0,1, 0.1)))

table(t.event_LoE)

t.event_LoE <- round(t.event_LoE*(5/max(t.event_LoE)) + 1 , digits=0) # standardize the time to event (to fit with the trial duration)
# the standardisation could be to fit tte to the entire trial duration, or to fit it to be at specific visits. We standardise it to fit the entire trial duration (6 weeks) and to have most of the intercurrent events up to and including week 4
# Other Weibull distributions (parameters) can be used to better describe the wanted time to events
# Depending on the desired percentages of intercurrent events and how they can be achieved. The standardisation can be used such that only a certain percentage of patients experience
# the intercurrent event during the trial (e.g., standardise it to fit to longer than end of trial (week 6), or it can be standardised to week 6 and then use a Binomial distribution to
# achieve a certain percentage, as we illustrate and use below).

describe(t.event_LoE)
d$t.LoE <- t.event_LoE

d$t.LoE <- d$t.LoE * rep(rbinom(n, 1, 0.5), each = length(visits))
# this particular step could be also used piecewise with different Binomial probabilities for time intervals, while keeping the intercurrent event percentage similar at  trial level.

    #cbind(d$t.LoE, rep(rbinom(n, 1, 0.5), each = length(visits)))

d$LoE_yes <-factor(ifelse(d$t.LoE !=0, 1, 0))
describe(d$LoE_yes)

d$t.LoE[d$t.LoE==0] <- c("No LoE")

describe(d$t.LoE)

cbind(d$Treat,d$bi_0, d$bi_1, t.event_LoE, d$t.LoE)

    #View(d[,c(1, 3, 4, 5, 8, 9)])
    #View(d)

describe(t.event_LoE[d$Treat==1])
describe(t.event_LoE[d$Treat==0])

head(d)

d$t.AE <- d$AE_yes <- 0

# just LoE
p<- ggplot(data = d[d$LoE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
plot_LoE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM-LoE pattern"); plot_LoE



# LoE
p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM-LoE all")




#####################################################
#####################################################
# Cox model and survival for AE in experimental arm #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# subset arms
d_c_exp <- d[d$Treat==1,]
d_c_control <- d[d$Treat==0,]

c2 <- -1 # parameter in the linear predictor

# Weibull distribution used to generate survival data - > time to adverse events in the experimental arm
# As mentioned above for the distribution used to generate time to lack of efficacy, different distributions can be used and finding the right shape needs finetuning.
lambda_AE_exp 	<- 3			# scale parameter
nu_AE_exp 		<- 	1.2		# shape parameter

LP2 <- (d_c_exp$bi_0 + d_c_exp$bi_1)/100 +  c2
describe(LP2)

#hist(d$bi_1)

t.event_AE_exp <- (-log(rep(runif(length(unique(d_c_exp$id))), each=length(visits)))/(lambda_AE_exp * exp(LP2))) ^ (1/nu_AE_exp) ; t.event_AE_exp
# how do I steer the rule here? I need to subset somehow on the positive slopes to "indicate" LoE, And the same for AE in experimental arm and in control arm.
describe(t.event_AE_exp)

hist(t.event_AE_exp)
describe(d_c_exp)


max(t.event_AE_exp)

d_c_exp$t.event_AE_exp <- t.event_AE_exp

t.event_AE_exp <- round(t.event_AE_exp*(5/max(t.event_AE_exp)) + 1 , digits=0) # at this step you can steer more or less extra the timings of AE or any other intercurrent event
# the standardisation can be used to fit the intercurrent events at a specific visit or for a specific visits interval

    #describe(t.event_AE_exp)
d_c_exp$t.AE_exp <- t.event_AE_exp

d_c_exp$t.AE_exp <- d_c_exp$t.AE_exp * rep(rbinom(length(unique(d_c_exp$id)), 1, 0.5), each = length(visits))

    #cbind(d$t.LoE, rep(rbinom(n, 1, 0.5), each = length(visits)))

d_c_exp$AE_yes <-factor(ifelse(d_c_exp$t.AE_exp !=0, 1, 0))
describe(d_c_exp$AE_yes)

d_c_exp$t.AE_exp[d_c_exp$t.AE_exp==0] <- c("No AE")

#describe(d_c_exp$t.AE_exp)
  
    #d_c_exp_AE <- d_c_exp[d_c_exp$bi_1<-1,]
    #d_c_exp_no_AE <- d_c_exp[d_c_exp$bi_1>=-1,]
    #d_c_exp_no_AE$AE_yes <- d_c_exp_no_AE$t.AE <- 0

#t.event_AE <- (-log(rep(runif(length(unique(d_exp$id))), each=length(visits)))/(lambda * exp((d_exp$bi_0 + d_exp$bi_1)/1000))) ^ (1/nu) ; t.event_AE

#t.event_AE <- round(t.event_AE*42/7, digits=0); t.event_AE
#d_exp$t.AE <- t.event_AE


#describe(d_c_exp_AE$bi_1)

#View(d_c_exp)

#head(d_c_exp_no_AE)

#d_c_cAEexp <- rbind(d_c_exp_AE, d_c_exp_no_AE)
#head(d_c_exp_AE)
#d_c_cAEexp


#View(d_c_cAEexp)


# just AE in experimental arm
p<- ggplot(data = d_c_exp[d_c_exp$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60)) + ggtitle("JM - AE exp pattern")



# AE in experimental arm and the others
p<- ggplot(data = d_c_exp, aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60)) + ggtitle("JM - AE all")



#####################################################
#####################################################
# Cox model and survival for AE in control arm


# subset on slightly-positive slopes in control arm
#resume here
#describe(d_c_control$bi_1)

#d_c_control_AE <- d_c_control[d_c_control$bi_1>0.15,]
#d_c_control_no_AE <- d_c_control[d_c_control$bi_1<=0.15,]
#d_c_control_no_AE$AE_yes <- d_c_control_no_AE$t.AE <- 0


#t.event_AE <- (-log(rep(runif(length(unique(d_c_control_AE$id))), each=length(visits)))/(lambda * exp((d_c_control_AE$bi_0 + d_c_control_AE$bi_1)/1000))) ^ (1/nu) ; t.event_AE

#t.event_AE <- round(t.event_AE*42/7, digits=0); t.event_AE
#d_c_control_AE$t.AE <- t.event_AE

#d_c_control_AE$AE_yes <- factor(ifelse(d_c_control_AE$t.AE == 2, 1,0))



#head(d_c_control_AE)
#head(d_c_control_no_AE)

#d_c_cAEcontrol <- rbind(d_c_control_AE, d_c_control_no_AE)

#d_c_cAEcontrol

#View(d_c_cAEcontrol)


# just AE in control arm
#p<- ggplot(data = d_c_cAEcontrol[d_c_cAEcontrol$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
#p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
#  scale_y_continuous(limits = c(0, 50))



# AE in control arm and the others
#p<- ggplot(data = d_c_cAEcontrol, aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
#p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
 # scale_y_continuous(limits = c(0, 50))


#head(d_c)
#head(d_c_cAEexp)
#head(d_c_cAEcontrol)


  
  
    
  c3 <- 1
  
  lambda_AE_control 	<- 2			# scale parameter
  nu_AE_control 		<- 	3.5	# shape parameter
  
  
  LP3 <- (d_c_control$bi_0 + d_c_control$bi_1)/100 +  c3
  describe(LP3)
  

  t.event_AE_control <- (-log(rep(runif(length(unique(d_c_control$id))), each=length(visits)))/(lambda_AE_control * exp(LP3))) ^ (1/nu_AE_control) ; t.event_AE_control
  # how do I steer the rule here? I need to subset somehow on the positive slopes to "indicate" LoE, And the same for AE in experimental arm and in control arm.
  describe(t.event_AE_control)
  
  hist(t.event_AE_control)
  describe(d_c_control)
  
  
  
  d_c_control$t.event_AE_control <- t.event_AE_control
  
  t.event_AE_control <- round(t.event_AE_control*(1) + 2 , digits=0) 
  
  describe(t.event_AE_control)
  d_c_control$t.AE_control <- t.event_AE_control
  
  d_c_control$t.AE_control <- d_c_control$t.AE_control * rep(rbinom(length(unique(d_c_control$id)), 1, 0.5), each = length(visits))
  
  #cbind(d$t.LoE, rep(rbinom(n, 1, 0.5), each = length(visits)))
  
  d_c_control$AE_yes <-factor(ifelse(d_c_control$t.AE_control !=0, 1, 0))
  describe(d_c_control$AE_yes)
  
  d_c_control$t.AE_control[d_c_control$t.AE_control==0] <- c("No AE")
  
  describe(d_c_control$t.AE_control)
  
  
  
  #d_c_exp_AE <- d_c_exp[d_c_exp$bi_1<-1,]
  #d_c_exp_no_AE <- d_c_exp[d_c_exp$bi_1>=-1,]
  #d_c_exp_no_AE$AE_yes <- d_c_exp_no_AE$t.AE <- 0
  
  
  
  #t.event_AE <- (-log(rep(runif(length(unique(d_exp$id))), each=length(visits)))/(lambda * exp((d_exp$bi_0 + d_exp$bi_1)/1000))) ^ (1/nu) ; t.event_AE
  
  #t.event_AE <- round(t.event_AE*42/7, digits=0); t.event_AE
  #d_exp$t.AE <- t.event_AE
  
  
  #describe(d_c_exp_AE$bi_1)
  
  
  
  #head(d_c_exp_no_AE)
  
  #d_c_cAEexp <- rbind(d_c_exp_AE, d_c_exp_no_AE)
  #head(d_c_exp_AE)
  #d_c_cAEexp
  
  
  #View(d_c_cAEexp)
  
  
  # just AE in control arm
  p<- ggplot(data = d_c_control[d_c_control$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
  p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
    scale_y_continuous(limits = c(-10, 60))
  
  
  
  # AE in control arm and the others
  p<- ggplot(data = d_c_control, aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
  p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
    scale_y_continuous(limits = c(-10, 60))
  
  
  # change parameters to get more at week 2
  
  d_c_exp$t.event_AE_control <-0
  d_c_exp$t.AE_control <-0
  
  
  d_c_control$t.event_AE_exp <-0
  d_c_control$t.AE_exp <-0
  
  head(d_c_exp)
  head(d_c_control)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


d_c_c_all <- rbind(d_c_exp, d_c_control)

describe(d_c_c_all)

#View(d_c_c_all)

d_united <- d_c_c_all



class(d_united$visit)
class(d_united$Treat)
d_united$Treat <- factor(d_united$Treat)
d_united$LoE_yes <- factor(d_united$LoE_yes)
d_united$AE_yes <- factor(d_united$AE_yes)


d_united$AE_YES <- ifelse(d_united$AE_yes==1, 1, 0)
d_united$LoE_YES <- ifelse(d_united$AE_yes==0 & d_united$LoE_yes==1, 1, 0)

View(d_united)

d_united$Behavior <- ifelse(d_united$AE_YES==1, "AE",
                            ifelse(d_united$AE_YES==0 & d_united$LoE_YES==1, "LoE", "No IE"))

#View(d_united)

#describe(d_united$Behavior)

p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
  scale_y_continuous(limits = c(-10, 60))



# Plot trajectories

# just LoE
p<- ggplot(data = d_united[d_united$LoE_YES==1,], aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_LoE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60)) + ggtitle("JM - LoE pattern"); plot_LoE


#just AE
# AE
p<- ggplot(data = d_united[d_united$AE_YES==1,], aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_AE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM - AE pattern"); plot_AE

#just No IE
# AE
p<- ggplot(data = d_united[d_united$Behavior=="No IE",], aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_NoIE <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))+ ggtitle("JM - No IE pattern"); plot_NoIE


# All patients with their trajectory 
# LoE
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))


# AE
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))

# All behaviors
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
plot_all <- p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60)) + ggtitle("JM - All patterns"); plot_all

(plot_all / plot_LoE) | (plot_AE / plot_NoIE)



describe(d_united)
View(d_united)
