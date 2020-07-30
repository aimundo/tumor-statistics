library(tidyverse)
library(splines)

set.seed(11)
n_time      <- 10
n_treatment <- 3
n_reps      <- 5
## spline degrees of freedom
## controls "wigglieness" of the response
df          <- 4

time      <- 1:n_time
treatment <- 1:n_treatment

## basis expansion over time that allows for smooth functional response
X_bs <- bs(time, df = df)
## beta parameters
## each treatment gets a different set of parameters
beta <- matrix(rnorm(df * n_treatment), df, n_treatment)

## mean time-varying response within each group
mu <- X_bs %*% beta

dat_mu <- data.frame(
    mu        = c(mu),
    time      = rep(1:n_time, times = n_treatment),
    treatment = factor(rep(1:n_treatment, each = n_time))
)
dat_mu %>%
    ggplot(aes(x = time, y = mu, group = treatment, color = treatment)) +
    geom_line() +
    scale_color_viridis_d(end = 0.75) +
    ggtitle("Simulated time-varying mean responses")

## measurement/individual response processes
## can make this more specific to the model as needed
sigma <- 0.25

##
## generate the data
##

y <- array(NA, dim = c(n_time, n_treatment, n_reps))
dimnames(y) <- list(
    time       = 1:n_time,
    treatment  = 1:n_treatment,
    individual = 1:n_reps
)

for (i in 1:n_reps) {
    y[, , i] <- mu + matrix(rnorm(n_time * n_treatment, 0, sigma), n_time, n_treatment)
}

## convert to a data.frame for plotting
dat <- as.data.frame.table(y, responseName = "y")
dat$time <- as.numeric(dat$time)

## merge in the mean response mu
dat <- merge(dat, dat_mu)

dat %>%
    ggplot(aes(x = time, y = mu, group = treatment, color = treatment)) +
    geom_line(alpha = 0.35) +
    scale_color_viridis_d(end = 0.75) +
    ggtitle("Simulated time-varying responses") +
    geom_point(aes(x = time, y = y, color = treatment)) +
    theme_bw()


fit <- stan_gamm4(y ~ s(time, treatment), data = dat)

plot_nonlinear(fit, smooths = "s(treatment)")
