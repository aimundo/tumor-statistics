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


gam.check(mod1)
gam.check(mod2)
