library(mgcv)
library(tidyverse)
library(splines)
library(rstan)
library(rstanarm)
library(patchwork)


## number of observations -- currently using 3
## use measurement error of sd as 5 (rough estimate)

dat <- read_csv(here::here("data", "Vishwanath StO2.csv"))

## plot the data
ggplot(dat, aes(x = Day, y = StO2, color = Group)) +
    geom_line()

simulate_data <- function(dat, n = 10, sd = 5) {
    dat_sim <- dat %>%
        slice(rep(1:n(), each = n)) %>%
        group_by(Group, Day) %>%
        mutate(observation = rnorm(n, StO2, sd),
               observation_log = rlnorm(n, log(StO2), log(sd)),
               observation_tobit = pmax(rnorm(n, StO2, sd), 0.0001)) %>%
        ungroup()

    return(dat_sim)
}



n <- 10
sd <- 5
set.seed(11)
dat_sim <- simulate_data(dat, n, sd)
df <- 6

n_time <- length(unique(dat_sim$Day))
n_treatment <- length(unique(dat_sim$Group))



missing <- sample(1:100, 40)
dat_missing <- dat_sim[-missing, ]


ggplot(dat_missing, aes(x = Day, y = observation, color = Group)) +
    geom_point() +
    ggtitle("Simulated time-varying responses") +
    scale_color_viridis_d(end = 0.75) +
    geom_line(aes(y = StO2), alpha = 0.75) +
    theme_bw()

dat_missing %>%
    group_by(Day, Group) %>%
    summarize(count = n())


dat_missing$Group <- as.factor(dat_missing$Group)
dat_sim$Group <- as.factor(dat_sim$Group)
mod1 <- gam(observation ~ s(Day, by = Group, k = 5), data  = dat_missing)
mod2 <- gam(observation ~ s(Day, by = Group, k = 5), data  = dat_sim)

summary(mod1)
summary(mod2)

par(ask=F)
layout(matrix(1:4, 2, 2))
plot(mod1, main = "Missing")
plot(mod2, main = "Full")

par(mfrow=c(1,1))

gam.check(mod1)
gam.check(mod2)

#plotting smooths and simulated data
#full data
full_predict<-with(dat_sim,expand.grid(observation=seq(min(observation),
                                                       max(observation),length=400),
                                                    Group=Group,Day=Day))

full_predict<-cbind(full_predict,
                         predict(mod2,full_predict,
                                 se.fit = TRUE,type='response'))


ggplot(data=full_predict,aes(x=Day,y=fit))+geom_point()+facet_wrap(~Group)

ggplot(data=dat_sim, aes(x=Day, y=observation, group=Group)) +
    facet_wrap(~Group) +
    geom_point(colour='black',size=1,alpha=0.5)+
    geom_ribbon(aes(ymin=(fit - 2*se.fit), ymax=(fit + 2*se.fit), x=Day),
                data=full_predict, alpha=0.3,
                inherit.aes=FALSE) +
    geom_line(aes(x=Day,y=(fit),color=factor(Group)), size=1,data=full_predict,show.legend = FALSE)

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
