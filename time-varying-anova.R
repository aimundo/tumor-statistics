library(tidyverse)
library(splines)
library(rstan)
library(rstanarm)

set.seed(11)
n_time      <- 10 #repeated observations
n_treatment <- 3 #treatment groups
n_reps      <- 5 #subjects per treatment
## spline degrees of freedom
## controls "wigglieness" of the response
df          <- 4

time      <- 1:n_time
treatment <- 1:n_treatment

## basis expansion over time that allows for smooth functional response
X_bs <- bs(time, df = df) #generates matrix for representing the family of piecewise polynomials with specificed knots
#and df
## beta parameters
## each treatment gets a different set of parameters
beta <- matrix(rnorm(df * n_treatment), df, n_treatment) #generates samples from the normal distribution and
#places them in a matrix with size df and n_treatment

## mean time-varying response within each group
mu <- X_bs %*% beta #matrix multiplication to generate the mean time-varying response

dat_mu <- data.frame(
    mu        = c(mu),
    time      = rep(1:n_time, times = n_treatment),
    treatment = factor(rep(1:n_treatment, each = n_time))
)

#this data frame has mean response for 10 observations for 3 treatment groups
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
#generates empty arrays. n_time rows, n_treatment columns, n_reps (number of subjects per treatment group) number of arrays and then puts names on the rows
#and columns of each array
y <- array(NA, dim = c(n_time, n_treatment, n_reps))
dimnames(y) <- list(
    time       = 1:n_time,
    treatment  = 1:n_treatment,
    individual = 1:n_reps
)
#this cycle places mean response per treatment group from mu, adds a newly generated parameter generated from the
#normal distribution and pastes them in the previously made arrays
for (i in 1:n_reps) {
    y[, , i] <- mu + matrix(rnorm(n_time * n_treatment, 0, sigma), n_time, n_treatment)
}

## convert to a data.frame for plotting
dat <- as.data.frame.table(y, responseName = "y")
dat$time <- as.numeric(dat$time)

## merge in the mean response mu with the generated data per treatment per timepoint per subject
dat <- merge(dat, dat_mu)

dat %>%
    ggplot(aes(x = time, y = mu, group = treatment, color = treatment)) +
    geom_line(alpha = 0.35) +
    scale_color_viridis_d(end = 0.75) +
    ggtitle("Simulated time-varying responses") +
    geom_point(aes(x = time, y = y, color = treatment)) +
    theme_bw()


library(rstan)
dat_fit <- list(
    N = length(y), #length (total number of observations) is No. of subjects*timepoints*treatments
    n_time = n_time,
    n_treatment = n_treatment,
    q = ncol(X_bs), #wiggliness of the response
    time_idx = dat$time, #all the timepoints per group per treatment
    treatment_idx = as.numeric(dat$treatment),
    y = dat$y, #all the generated responses
    X = X_bs #the polynomial spline matrix
)

fit <- stan(file = here::here("time-varying.stan"), data = dat_fit)
e <- rstan::extract(fit)

mu_post <- sapply(1:nrow(e$sigma), function(i) X_bs %*% e$beta[i, , ],
                  simplify = "array")
dimnames(mu_post) <- list(
    time       = 1:n_time,
    treatment  = 1:n_treatment,
    iteration = 1:nrow(e$sigma)
)

mu_post_mean <- apply(mu_post, c(1, 2), mean)
mu_post_lower50 <- apply(mu_post, c(1, 2), quantile, prob = 0.25)
mu_post_upper50 <- apply(mu_post, c(1, 2), quantile, prob = 0.75)
mu_post_lower95 <- apply(mu_post, c(1, 2), quantile, prob = 0.025)
mu_post_upper95 <- apply(mu_post, c(1, 2), quantile, prob = 0.975)

dat_plot <- data.frame(
    mean      = c(mu_post_mean),
    lower_50  = c(mu_post_lower50),
    upper_50  = c(mu_post_upper50),
    lower_95  = c(mu_post_lower95),
    upper_95  = c(mu_post_upper95),
    truth     = c(mu),
    time      = rep(1:n_time, times = n_treatment),
    treatment = factor(rep(1:n_treatment, each = n_time))
)

dat_plot %>%
    ggplot(aes(x = time, y = truth, color = treatment, group = treatment)) +
    geom_line(lwd = 1.5) +
    # geom_line(aes(x = time, y = mean, color = treatment, group = treatment)) +
    scale_color_viridis_d(end = 0.75) +
    scale_fill_viridis_d(end = 0.75) +
    ggtitle("Fitted time-varying responses") +
    theme_bw() +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, fill = treatment), alpha = 0.5, color = NA) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = treatment), alpha = 0.25, color = NA) +
    geom_point(data = dat, aes(x = time, y = y))









## Using pre-built software
## need to figure out the plotting better but it should work too.


# library(rstanarm)
# library(tidybayes)
# fit <- stan_gamm4(y ~ s(time, by = treatment), data = dat)
# get_variables(fit)
#
# plot_nonlinear(fit)
#
#
# XZ <- fit$x
#
# B <- as.matrix(fit)[, colnames(XZ), drop = FALSE]
# labels <- sapply(fit$jam$smooth, "[[", "label")
# xnames <- sapply(fit$jam$smooth, "[[", "vn")
# names(fit$jam$smooth) <- labels
# names(xnames) <- labels
# fs <- sapply(fit$jam$smooth, FUN = "inherits", what = "fs.interaction")
# smooths <- 1:length(labels)
# original <- fit$jam$model
#
#
#
#
# prob = 0.9
# scheme <- bayesplot::color_scheme_get()
# size = 0.75
# alpha = 1
# facet_args <- list()
#
# df_list <- lapply(fit$jam$smooth[smooths], FUN = function(s) {
#     incl <- s$first.para:s$last.para
#     b <- B[, incl, drop = FALSE]
#     if (inherits(s, "fs.interaction")) {
#         xx <- original[, s$base$term]
#         fac <- original[, s$fterm]
#         out <- by(data.frame(fac, xx), list(fac), FUN = function(df) {
#             df <- df[order(df[, 2]), ]
#             names(df) <- c(s$fterm, s$base$term)
#             xz <- mgcv::PredictMat(s, df)
#             f <- linear_predictor.matrix(b, xz)
#             data.frame(predictor = df[, 2], lower = apply(f,
#                                                           2, quantile, probs = (1 - prob)/2), upper = apply(f,
#                                                                                                             2, quantile, probs = prob + (1 - prob)/2),
#                        middle = apply(f, 2, median), term = paste(s$label,
#                                                                   df[, 1], sep = "."))
#         })
#         do.call(rbind, args = out)
#     }
#     else {
#         xz <- XZ[, incl, drop = FALSE]
#         x <- original[, s$term]
#         ord <- order(x)
#         x <- x[ord]
#         xz <- xz[ord, , drop = FALSE]
#         if (!is.null(s$by.level)) {
#             fac <- original[, s$by][ord]
#             mark <- fac == s$by.level
#             x <- x[mark]
#             xz <- xz[mark, , drop = FALSE]
#         }
#         f <- rstanarm:::linear_predictor.matrix(b, xz)
#         data.frame(predictor = x, lower = apply(f, 2, quantile,
#                                                 probs = (1 - prob)/2), upper = apply(f, 2, quantile,
#                                                                                      probs = prob + (1 - prob)/2), middle = apply(f,
#                                                                                                                                   2, median), term = s$label)
#     }
# })
#
#
#
#
# plot_data <- do.call(rbind, df_list)
# str(plot_data)
# facet_args[["facets"]] <- ~term
#
# ## recode the factors
# levels(plot_data$term) <- c("1", "2", "3")
#
# dat_full <- plot_data %>%
#     rename(time = predictor,
#            mean = middle,
#            treatment = term) %>%
#     merge(
#         data.frame(
#             time      = rep(1:n_time, times = n_treatment),
#             treatment = factor(rep(1:n_treatment, each = n_time)),
#             truth     = c(mu)
#         )
#     )
#
# fit %>%
#     spread_draws(Intercept)
#
# ## looks like there needs to be the intercept term
#
# dat_full %>%
#     ggplot(aes(x = time, y = truth, color = treatment, group = treatment)) +
#     geom_line(lwd = 1.5) +
#     # geom_line(aes(x = time, y = mean, color = treatment, group = treatment)) +
#     scale_color_viridis_d(end = 0.75) +
#     scale_fill_viridis_d(end = 0.75) +
#     ggtitle("Fitted time-varying responses") +
#     theme_bw() +
#     geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.5, color = NA) +
#     # geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = treatment), alpha = 0.25, color = NA) +
#     geom_point(data = dat, aes(x = time, y = y))
#
#
#
#
#
# ggplot
# # %>%
#     merge(dat_full, dat_truth)
#
# ggplot(plot_data, aes_(x = ~predictor)) +
#     geom_ribbon(aes_(ymin = ~lower,
#                      ymax = ~upper), fill = scheme[[1]], color = scheme[[2]],
#                 alpha = alpha, size = size) +
#     geom_line(aes_(y = ~middle),
#               color = scheme[[5]], size = 0.75 * size, lineend = "round") +
#     labs(y = NULL) + do.call(facet_wrap, facet_args) + bayesplot::theme_default()
#
#
#
#
#
#
#
