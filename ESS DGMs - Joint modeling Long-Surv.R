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


n <- 190 # number of patients





set.seed(2147483629)
b0 <- 29.5
b1 <- -0.55
b2 <- -0.583
bi_means <- c(0, 0)
bi_covm <- matrix(c(24.611, 0.5869809, 0.5869809, 1.157), nrow = 2)
bi <- mvrnorm(n, bi_means, bi_covm)	
eps.sd <-3.247 




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




# if I make the Weibull distribution


#
#####################################################
# Cox model and survival for LoE at trial level

describe(d$bi_0)
describe(d$bi_1)

hist(d$bi_0)
hist(d$bi_1)


#par(mfrow=c(1,1))


hist(d$bi_0[d$Treat==0])
hist(d$bi_0[d$Treat==1])

hist(d$bi_1[d$Treat==0])
hist(d$bi_1[d$Treat==1])

#subset on positive slopes
#hist(d$bi_1[d$Treat==1 & d$bi_1>0])

#d_LoE <- d[d$bi_1>0,]
#d_no_LoE <- d[d$bi_1<=0,]


#d_no_LoE$t.LoE <- d_no_LoE$LoE_yes <- 0

c1 <- -0.5

lambda_LoE 	<- 3.5			# scale parameter
nu_LoE 		<- 	0.9		# shape parameter


LP1 <- (d$bi_0 + d$bi_1)/100 +  c1 * (as.numeric(d$Treat)-1)
describe(LP1)

#hist(d$bi_1)

t.event_LoE <- (-log(rep(runif(length(unique(d$id))), each=length(visits)))/(lambda_LoE * exp(LP1))) ^ (1/nu_LoE) ; t.event_LoE
# how do I steer the rule here? I need to subset somehow on the positive slopes to "indicate" LoE, And the same for AE in experimental arm and in control arm.

#1. add here the probability to rbinom(independent) 
#2. both the occurrence and the timing of IE

hist(t.event_LoE)



d$t.event_LoE <- t.event_LoE

t.event_LoE <- round(t.event_LoE*(4.16) + 1 , digits=0) 

describe(t.event_LoE)
d$t.LoE <- t.event_LoE

d$t.LoE <- d$t.LoE * rep(rbinom(n, 1, 0.5), each = length(visits))

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
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))



# LoE
p<- ggplot(data = d, aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))








#####################################################
#####################################################
# Cox model and survival for AE in experimental arm



d_c_exp <- d[d$Treat==1,]
d_c_control <- d[d$Treat==0,]

c2 <- -1


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




d_c_exp$t.event_AE_exp <- t.event_AE_exp

t.event_AE_exp <- round(t.event_AE_exp*(0.83) + 2 , digits=0) # at this step you can steer more or less extra the timings of AE or any other intercurrent event
# stricter like the above or more loosely distributed with high fidelity to the Weibull distribution


describe(t.event_AE_exp)
d_c_exp$t.AE_exp <- t.event_AE_exp

d_c_exp$t.AE_exp <- d_c_exp$t.AE_exp * rep(rbinom(length(unique(d_c_exp$id)), 1, 0.5), each = length(visits))

#cbind(d$t.LoE, rep(rbinom(n, 1, 0.5), each = length(visits)))

d_c_exp$AE_yes <-factor(ifelse(d_c_exp$t.AE_exp !=0, 1, 0))
describe(d_c_exp$AE_yes)

d_c_exp$t.AE_exp[d_c_exp$t.AE_exp==0] <- c("No LoE")

describe(d_c_exp$t.AE_exp)



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


# just AE in experimental arm
p<- ggplot(data = d_c_exp[d_c_exp$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))



# AE in experimental arm and the others
p<- ggplot(data = d_c_exp, aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))


# change parameters to get more at week 2



# tweak the distributions for each IE, how to make sure that people with the profiles I think they should have an IE. The conclusion could be that this is a way to get IEs
# but not easy to get the right profiles there.





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
  
  d_c_control$t.AE_control[d_c_control$t.AE_control==0] <- c("No LoE")
  
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


p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat) +
  scale_y_continuous(limits = c(-10, 60))



# Plot trajectories

# just LoE
p<- ggplot(data = d_united[d_united$LoE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))


#just AE
# AE
p<- ggplot(data = d_united[d_united$AE_yes==1,], aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))





# All patients with TRUE trajectory 
# LoE
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=LoE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))



# AE
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=AE_yes)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))



# All behaviors
p<- ggplot(data = d_united, aes(x = visit, y = MADRS10_collected, group = id, color=Behavior)) 
p + geom_line() + stat_summary(aes(group = 1), geom = "point", fun = mean, shape = 18, size = 3, col="red") + facet_wrap(~ Treat)+
  scale_y_continuous(limits = c(-10, 60))



describe(d_united)






lp <- seq(-2, 2, 0.1)
plot(lp)
lambda <-10
nu <- 0.5
plot(exp(lp))


time_ie<- - ( log(runif(length(lp)))/(lambda * exp(lp)) )^1/nu
plot(time_ie)



