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
#and a standard deviation of 10% StO2
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
sd <- 10
set.seed(11) #setting seed
df <- 6
dat_sim <- simulate_data(dat, n, sd)

#plotting simulated data
ggplot(dat_sim, aes(x = Day, y = observation_tobit, color = Group)) +
    geom_point()+
    scale_x_continuous(breaks=c(0,2,5,7,10))


n_time <- length(unique(dat_sim$Day))
n_treatment <- length(unique(dat_sim$Group))


#if you have n timepoints can you do n knots or should you do n-1 knots? (I read somewhere that
#the n-1 option was preferred for computational issues)?


mod1 <- gam(observation_tobit ~ Group+s(Day, by = Group, k = 5), data  = dat_sim)
appraise(mod1) #SIMPLEST MODEL: A SMOOTH FOR TIME AS THE COVARIATE, WHICH CAN VARY BY GROUP BUT... APPRAISING THE MODEL SHOWS A WEIRD TREND
#ON THE RESIDUALS

mod2 <- gam(observation_tobit ~ Group+s(Day,by=Group,k=5),data  = dat_sim)
appraise(mod2)
#THIS MODEL is accounting for a linear effect for each group, and a smooth for time by each group
#Why is it fitting better than the "simple" model? What does this "linear" Group term mean?
#The idea for this came from
#https://stats.stackexchange.com/questions/486118/model-building-in-generalized-additive-mixed-models-gamms

#mod2 <- gam(observation ~ s(Day, by = Group, k = 5), data  = dat_missing)


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
#f_predict<-with(dat_sim,expand.grid(StO2=seq(min(observation_tobit),
#                                                        max(observation_tobit),
#                                                         length=10),
#                                                     Group=Group,Day=Day))

 #f_predict<-cbind(f_predict,
   #                       predict(mod2,f_predict,
  #                                se.fit = TRUE,type='response'))
#

#creates a dataframe wusing the lenght of the covariates
f_predict <- expand_grid(Group = factor(c("Control", "Treatment")),
                         Day = seq(0, 10, by = 0.1))

#adds the predictions to the grid and creates a confidence interval
f_predict<-f_predict%>%
    mutate(fit = predict(mod2,f_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(mod2, f_predict,se.fit = TRUE,type='response')$se.fit)


ggplot(data=dat_sim, aes(x=Day, y=observation_tobit, group=Group)) +
    facet_wrap(~Group)+
    geom_point(colour='black',size=1,alpha=0.5)+
    geom_ribbon(aes(ymin=(fit - 2*se.fit), ymax=(fit + 2*se.fit), x=Day),
                data=f_predict, alpha=0.5,
                inherit.aes=FALSE) +
    geom_line(aes(y=(fit),color=factor(Group)),
              size=1,data=f_predict,show.legend = FALSE)



#it works now, but why do the smooths still look so rough? (the connections at the knots
#don't look "smooth")

##trying the pairwise comparisons

pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Treatment', 'Control'))
smooth_diff <- function(model, newdata, f1, f2, alpha = 0.05,
                        unconditional = FALSE) {
    xp <- predict(model, newdata = newdata, type = 'lpmatrix')
    c1 <- grepl(f1, colnames(xp))
    c2 <- grepl(f2, colnames(xp))
    #r1 <- newdata[[var]] == f1
    #r2 <- newdata[[var]] == f2
    r1 <- with(newdata, Group == f1)
    r2 <- with(newdata, Group == f2)
    ## difference rows of xp for data from comparison
    X <- xp[r1, ] - xp[r2, ]
    ## zero out cols of X related to splines for other lochs
    X[, ! (c1 | c2)] <- 0
    ## zero out the parametric cols
    X[, !grepl('^s\\(', colnames(xp))] <- 0
    dif <- X %*% coef(model)
    se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
    crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    data.frame(pair = paste(f1, f2, sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)
}
comp1<-smooth_diff(mod2,pdat,'Control','Treatment')
comp_StO2 <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp1))

c1<-ggplot(comp_StO2, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(color='black',size=1) +
    facet_wrap(~ pair) +
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))


c1
    ##################################################################
############################################################
#missing data
#create a sequence of 40 random numbers between 1 and 100, these numbers will
#correspond to the row numbers to be randomly erased from the original dataset
missing <- sample(1:100, 40)

#create a new dataframe from the simulated data with 40 rows randomly removed
dat_missing <- dat_sim[-missing, ]

#Count the number of remaining observations per day (original dataset had 10 per group per day)
dat_missing %>%
    group_by(Day, Group) %>%
    summarize(count = n())

#plot the simulated data
ggplot(dat_missing, aes(x = Day, y = observation_tobit, color = Group)) +
    geom_point() +
    ggtitle("Simulated time-varying responses") +
    scale_color_viridis_d(end = 0.75) +
    geom_line(aes(y = StO2), alpha = 0.75) +
    theme_bw()



#the same model used for the full dataset
mod_m1 <- gam(observation_tobit ~ Group+s(Day,by=Group,k=5), data  = dat_missing)
#appraise the model
appraise(mod_m1)


m_predict <- expand_grid(Group = factor(c("Control", "Treatment")),
                         Day = seq(0, 10, by = 0.1))

#adds the predictions to the grid and creates a confidence interval
m_predict<-m_predict%>%
    mutate(fit = predict(mod_m1,m_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(mod_m1, m_predict,se.fit = TRUE,type='response')$se.fit)


ggplot(data=dat_missing, aes(x=Day, y=observation_tobit, group=Group)) +
    facet_wrap(~Group)+
    geom_point(colour='black',size=1,alpha=0.5)+
    geom_ribbon(aes(ymin=(fit - 2*se.fit), ymax=(fit + 2*se.fit), x=Day),
                data=m_predict, alpha=0.5,
                inherit.aes=FALSE) +
    geom_line(aes(y=(fit),color=factor(Group)),
              size=1,data=m_predict,show.legend = FALSE)


comp2<-smooth_diff(mod_m1,pdat,'Control','Treatment')
comp_StO2_m <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp2))

c2<-ggplot(comp_StO2_m, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(color='black',size=1) +
    facet_wrap(~ pair) +
    labs(x = 'Days', y = expression(paste('Difference in StO'[2],'(missing data)' )))


c2
