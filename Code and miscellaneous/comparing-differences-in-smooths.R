set.seed(2021)
library(readr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(mgcv)
library(patchwork)


## Example with linear response
example <- function(n_time = 6, 
                    fun_type = "linear", 
                    error_type = "correlated") {
  
  if (!(fun_type %in% c("linear", "quadratic")))
    stop('fun_type must be either "linear", or "quadratic"')
  if (!(error_type %in% c("correlated", "independent")))
    stop('fun_type must be either "correlated", or "independent"')
  
  library(tidyverse)
  library(mvnfast)
  library(nlme)
  
  x <- seq(0, 1, length.out = n_time)
  mu <- matrix(0, length(x), 2)
  # linear response
  if (fun_type == "linear") {
    mu[, 1] <- - (x - .5) + 0.25
    mu[, 2] <- (x - .3)
  } else {
    # nonlinear response
    mu[, 1] <- - 4 * (x - .5)^2 + 0.25
    mu[, 2] <- 4 * (x - .5)^2
  }
  # matplot(mu, type = 'l')
  
  y <- array(0, dim = c(length(x), 2, 10))
  errors <- array(0, dim = c(length(x), 2, 10))
  
  if (error_type == "independent") {
    ## independent errors
    for (i in 1:2) {
      for (j in 1:10) {
        errors[, i, j] <- rnorm(6, 0, 0.25)
        y[, i, j] <- mu[, i] + errors[, i, j]
      }
    }
  } else {
    for (i in 1:2) {     # number of treatments
      for (j in 1:10) {  # number of subjects
        # compound symmetry errors
        errors[, i, j] <- rmvn(1, rep(0, length(x)), 0.1 * diag(6) + 0.25 * matrix(1, 6, 6))
        y[, i, j] <- mu[, i] + errors[, i, j]
      }
    }
  }    
  
  
  ## subject random effects
  
  ## visualizing the difference between independent errors and compound symmetry
  ## why do we need to account for this -- overly confident inference
  
  
  dimnames(y) <- list(time = x, treatment = 1:2, subject = 1:10)
  dimnames(errors) <- list(time = x, treatment = 1:2, subject = 1:10)
  dimnames(mu) <- list(time = x, treatment = 1:2)
  dat <- as.data.frame.table(y, responseName = "y")
  dat_errors <- as.data.frame.table(errors, responseName = "errors")
  dat_mu <- as.data.frame.table(mu, responseName = "mu")
  dat <- left_join(dat, dat_errors, by = c("time", "treatment", "subject"))
  dat <- left_join(dat, dat_mu, by = c("time", "treatment"))
  dat$time <- as.numeric(as.character(dat$time))
  dat <- dat %>%
    mutate(subject = factor(paste(subject, treatment, sep = "-")))
  
  
  ## repeated measures ANOVA in R
  fit_lm <- lm(y ~ time + treatment + time * treatment, data = dat)
  dat$preds_lm <- predict(fit_lm)
  
  fit_lme <- lme(y ~ treatment + time + treatment:time,
                 data = dat,
                 random = ~ 1 | subject,
                 correlation = corCompSymm(form = ~ 1 | subject)
  )
  
  
  pred_dat <- expand.grid(
    treatment = factor(1:2), 
    time = unique(dat$time)
  )
  
  dat$y_pred <- predict(fit_lme)
  
  
  return(list(
    dat = dat,
    pred_dat = pred_dat,
    fit_lm = fit_lm,
    fit_lme = fit_lme
    
  ))
}

dat <- example(fun_type = "quadratic", error_type = "independent")$dat

# make the mean shift relatively large

dat <- dat %>%
  mutate(adjust = case_when(
    treatment == "1" ~ 0,
    treatment == "2" ~ 0.25,
  )) %>%
  mutate(y = y + adjust)



ggplot(dat, aes(x = time, y = y, color = treatment)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top')

m <- gam(y ~ treatment + s(time, by = treatment, k = 6), data = dat)
summary(m)

plot(m, shade = TRUE, pages = 1, scale = 0)


pdat <- expand.grid(time = seq(0, 1, length = 100),
                    treatment = c('1', '2'))
xp <- predict(m, newdata = pdat, type = 'lpmatrix')

## which cols of xp relate to splines of interest?
c1 <- grepl('1', colnames(xp))
c2 <- grepl('2', colnames(xp))
## which rows of xp relate to sites of interest?
r1 <- with(pdat, treatment == '1')
r2 <- with(pdat, treatment == '2')

## difference rows of xp for data from comparison
X <- xp[r1, ] - xp[r2, ]
## zero out cols of X related to splines for other lochs
X[, ! (c1 | c2)] <- 0

#
# This modification keeps in the intercept adjustment
#

## zero out the parametric cols
# X[, !grepl('^s\\(', colnames(xp))] <- 0

dif <- X %*% coef(m)

se <- sqrt(rowSums((X %*% vcov(m)) * X))

crit <- qt(.975, df.residual(m))
upr <- dif + (crit * se)
lwr <- dif - (crit * se)



smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        include_mean = FALSE,
                        freq = FALSE) { # freq = FALSE gives Bayesian posterior distribution
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  if (!include_mean) {
    X[, !grepl('^s\\(', colnames(xp))] <- 0  
  }
  
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, freq = freq)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}




#
# Add in option to use frequentist vs. empirical Bayesian intervals
#
# frequentist intervals without mean shift
comp1 <- cbind(time = seq(0, 1, length = 100), 
              smooth_diff(m, pdat, '1', '2', 'treatment', include_mean = FALSE, freq = TRUE))

# frequentist intervals with mean shift
comp2 <- cbind(time = seq(0, 1, length = 100), 
              smooth_diff(m, pdat, '1', '2', 'treatment', include_mean = TRUE, freq = TRUE))

# Bayesian intervals without mean shift
comp3 <- cbind(time = seq(0, 1, length = 100), 
              smooth_diff(m, pdat, '1', '2', 'treatment', include_mean = FALSE, freq = FALSE))

# Bayesian intervals with mean shift
comp4 <- cbind(time = seq(0, 1, length = 100), 
              smooth_diff(m, pdat, '1', '2', 'treatment', include_mean = TRUE, freq = FALSE))

y_upper <- max(comp1$upper, comp2$upper, comp3$upper, comp4$upper)
y_lower = min(comp1$lower, comp2$lower, comp3$lower, comp4$lower)


#
# Bayesian intervals are slightly larger
#

p1 <- ggplot(comp1, aes(x = time, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  facet_wrap(~ pair, ncol = 2) +
  ylim(c(y_lower, y_upper)) +
  labs(x = NULL, y = 'Difference in trend') +
  ggtitle("Frequentist intervals without mean shift")


p2 <- ggplot(comp2, aes(x = time, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  facet_wrap(~ pair, ncol = 2) +
  ylim(c(y_lower, y_upper)) +
  labs(x = NULL, y = 'Difference in trend') +
  ggtitle("Frequentist intervals with mean shift")


p3 <- ggplot(comp3, aes(x = time, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  facet_wrap(~ pair, ncol = 2) +
  ylim(c(y_lower, y_upper)) +
  labs(x = NULL, y = 'Difference in trend') +
  ggtitle("Bayesian intervals without mean shift")


p4 <- ggplot(comp4, aes(x = time, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  facet_wrap(~ pair, ncol = 2) +
  ylim(c(y_lower, y_upper)) +
  labs(x = NULL, y = 'Difference in trend') +
  ggtitle("Bayesian intervals with mean shift")

(p1 + p2) / (p3 + p4)
