library(mgcv)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gratia)


## number of observations -- currently using 5
## use measurement error of sd as 5 (rough estimate)

#Data are simulated from Fig. 3 (C) of Vishwanath et. al:
#Vishwanath, Karthik, et al. "Using optical spectroscopy to
#longitudinally monitor physiological changes within solid tumors."
#Neoplasia 11.9 (2009): 889-900.
#https://www.sciencedirect.com/science/article/pii/S1476558609800420

dat <- read_csv(here::here("data", "Vishwanath StO2.csv"))
dat$Group<-as.factor(dat$Group)
## plot the mean response
ggplot(dat, aes(x = Day, y = StO2, color = Group)) +
    geom_line()+
    geom_point()+
    scale_x_continuous(breaks=c(0,2,5,7,10))


#This function simulates data for the tumor data using a dafault number of
#10 observations per time point,
#and a standard deviation of 5% StO2
#Because physiologically StO2 cannot go below 0%, data is generated with a cutoff value of
#0.0001 (the "observation_tobit")

simulate_data <- function(dat, n = 10, sd = 5) {
    dat_sim <- dat %>%
        slice(rep(1:n(), each = n)) %>%
        group_by(Group, Day) %>%
        mutate(observation = rnorm(n, StO2, sd),
               observation_tobit = pmax(rnorm(n, StO2, sd), 0.0001)) %>%
        ungroup()

    return(dat_sim)
}



n <- 10
sd <- 5
set.seed(11) #setting seed
df <- 6
dat_sim <- simulate_data(dat, n, sd)

#plotting simulated data
ggplot(dat_sim, aes(x = Day, y = observation_tobit, color = Group)) +
    geom_point()+
    scale_x_continuous(breaks=c(0,2,5,7,10))


n_time <- length(unique(dat_sim$Day))
n_treatment <- length(unique(dat_sim$Group))



mod1 <- gam(observation_tobit ~ s(Day, by = Group, k = 4), data  = dat_sim)
appraise(mod1) #ALSO, THE CHECK ON THIS MODEL LOOKS WEIRD! THE RESIDUALS ARE SHOWING A TREND

#mod2 <- gam(observation ~ s(Day, by = Group, k = 5), data  = dat_missing)
#if you have n timepoints can you do n knots or should you do n-1 knots? (I read somewhere that
#the n-1 option was preferred for computational issues)

summary(mod1)
#summary(mod2)

par(ask=F)
layout(matrix(1:4, 2, 2))
#plot(mod1, main = "Missing")
plot(mod1, main = "Full")

par(mfrow=c(1,1))

gam.check(mod1)
#gam.check(mod2)

#plotting smooths and simulated data
#full data
#I am doing something wrong here but I can't figure it out!!!!
f_predict<-with(dat_sim,expand.grid(StO2=seq(min(observation),
                                                       max(observation),length=100),
                                                    Group=Group,Day=Day))


f_predict<-cbind(f_predict,
                         predict(mod1,f_predict,
                                 se.fit = TRUE,type='response'))

ggplot(data=dat_sim, aes(x=Day, y=observation, group=Group)) +
    facet_wrap(~Group)+
    geom_point(colour='black',size=1,alpha=0.5)+
    geom_ribbon(aes(ymin=(fit - 2*se.fit), ymax=(fit + 2*se.fit), x=Day),
                data=f_predict, alpha=0.3,
                inherit.aes=FALSE) +
    geom_line(aes(y=(fit),color=factor(Group)), size=1,data=f_predict,show.legend = FALSE)


#This plot looks weird, why is it showing a non-smooth line as the fit? Is it because of the number
#of knots?  I have been able to do this in the past but for some reason this time it doesn't work
##################################################################
############################################################
#missing data
missing_predict<-with(dat_missing,expand.grid(observation=seq(min(observation),
                                                       max(observation),length=400),
                                       Group=levels(Group),Day=Day))

missing_predict<-cbind(missing_predict,
                    predict(mod1,missing_predict,
                            se.fit = TRUE,type='response'))


ggplot(data=dat_missing, aes(x=Day, y=observation, group=Group)) +
    facet_wrap(~Group) +
    geom_point(colour='black',size=1,alpha=0.5)+
    geom_ribbon(aes(ymin=(fit - 2*se.fit), ymax=(fit + 2*se.fit), x=Day),
                data=missing_predict, alpha=0.3,
                inherit.aes=FALSE) +
    geom_line(aes(y=(fit),color=factor(Group)), size=1,data=missing_predict,show.legend = FALSE)
