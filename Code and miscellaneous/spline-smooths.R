library(gratia)
library(mvnfast)
library(tidyverse)
library(patchwork)
load_mgcv()
set.seed(2021)
##
## re-write two gratia functions to provide pointwise and simultaneous smooths ----
##

calc_difference_R <- function(f1, f2, smooth, by_var, smooth_var, data, Xp, V, coefs, nrep = 1000) {
  ## make sure f1 and f2 are characters
  f1 <-  as.character(f1)
  f2 <-  as.character(f2)
  cnames <- colnames(Xp)
  ## columns of Xp associated with pair of smooths
  c1 <- grepl(gratia:::mgcv_by_smooth_labels(smooth, by_var, f1), cnames, fixed = TRUE)
  c2 <- grepl(gratia:::mgcv_by_smooth_labels(smooth, by_var, f2), cnames, fixed = TRUE)
  ## rows of Xp associated with pair of smooths
  r1 <- data[[by_var]] == f1
  r2 <- data[[by_var]] == f2

  ## difference rows of Xp for pair of smooths
  X <- Xp[r1, ] - Xp[r2, ]

  ## zero the cols related to other splines
  X[, ! (c1 | c2)] <- 0

  ## zero out the parametric cols
  X[, !grepl('^s\\(', cnames)] <- 0

  ## compute difference
  sm_diff <- drop(X %*% coefs)
  se <- sqrt(rowSums((X %*% V) * X))
  nr <- NROW(X)

  ## Calculate posterior simulation for smooths
  coefs_sim <- t(rmvn(nrep, rep(0, nrow(V)), V))
  rownames(coefs_sim) <- rownames(V)
  simDev <- X %*% coefs_sim
  absDev <- abs(sweep(simDev, 1, se, FUN = "/"))
  masd <- apply(absDev, 2, max)
  crit_s <- quantile(masd, prob = 0.95, type = 8)


  out <- list(smooth = rep(smooth, nr), by = rep(by_var, nr),
              level_1 = rep(f1, nr),
              level_2 = rep(f2, nr),
              diff = sm_diff, se = se,
              lower_s = sm_diff - crit_s * se,
              upper_s = sm_diff + crit_s*se)

  out <- new_tibble(out, nrow = NROW(X), class = "difference_smooth")
  ## Only need rows associated with one of the levels
  out <- bind_cols(out, data[r1, smooth_var])

  out
}

#does both ci and si
difference_smooths_R <- function(model,
                                 smooth,
                                 n = 100,
                                 ci_level = 0.95,
                                 newdata = NULL,
                                 partial_match = TRUE,
                                 unconditional = FALSE,
                                 frequentist = FALSE,
                                 nrep = 10000,
                                 include_means = TRUE,
                                 ...) {
  if (missing(smooth)) {
    stop("Must specify a smooth to difference via 'smooth'.")
  }

  # smooths in model
  S <- gratia::smooths(model) # vector of smooth labels - "s(x)"
  # select smooths
  select <-
    gratia:::check_user_select_smooths(smooths = S, select = smooth,
                                       partial_match = partial_match)#,
  # model_name = expr_label(substitute(object)))
  sm_ids <- which(select)
  smooths <- gratia::get_smooths_by_id(model, sm_ids)
  sm_data <- map(sm_ids, gratia:::smooth_data,
                 model = model, n = n, include_all = TRUE)
  sm_data <- bind_rows(sm_data)
  by_var <- by_variable(smooths[[1L]])
  smooth_var <- gratia:::smooth_variable(smooths[[1L]])
  pairs <- as_tibble(as.data.frame(t(combn(levels(sm_data[[by_var]]), 2)),
                                   stringsAsFactor = FALSE))
  names(pairs) <- paste0("f", 1:2)

  Xp <- predict(model, newdata = sm_data, type = "lpmatrix")
  V <- gratia:::get_vcov(model, unconditional = unconditional,
                         frequentist = frequentist)
  coefs <- coef(model)

  out <- pmap(pairs, calc_difference_R, smooth = smooth, by_var = by_var,
              smooth_var = smooth_var, data = sm_data, Xp = Xp, V = V,
              coefs = coefs, nrep = nrep)
  out <- bind_rows(out)
  crit <- qnorm((1 - ci_level) / 2, lower.tail = FALSE)

  out <- add_column(out,
                    lower = out$diff - (crit * out$se),
                    upper = out$diff + (crit * out$se),
                    .after = 6L)
  out
}






##
## Simulate some data like we have in our example ----
##

#  set.seed(2)
# n <- 50
# dat <- data.frame(x=seq(1:n), group = factor(rep(c(1, 2), each = n)))
# # create the grouped smooths
# X <- splines::bs(seq(1:n))
# beta <- matrix(rnorm(2 * ncol(X)), ncol(X), 2)
#
# dat$mu <- c(X %*% beta)
# dat$y <- dat$mu + rnorm(nrow(dat), 0, 0.1)
#
# ## plot the simulated data
# ggplot(dat, aes(x, y, color = group)) +
#     geom_point() +
#     geom_line(aes(x, mu, group = group))


### *******use the data from the paper*****************

dat<-tibble(StO2=c(4,27,3,2,0.5,7,4,50,45,56),
            x=rep(c(0,2,5,7,10),times=2),
            group=as.factor(rep(c("1","2"),each=5))
)



simulate_data <- function(dat, n = 10, sd = 5) {
    dat_sim <- dat %>%
        slice(rep(1:n(), each = n)) %>%
        group_by(group, x) %>%
        mutate(
            y = pmax(rnorm(n, StO2, sd), 0.0001),
            subject=rep(1:10),
            subject=factor(paste(subject, group))
        ) %>%
        ungroup()

    return(dat_sim)
}


#subject = factor(paste(subject, treatment, sep = "-")))

n <- 10 #number of observations
sd <- 10 #approximate sd from paper

dat_sim <- simulate_data(dat, n, sd)

##
## fit the model to the data ----
##

#mod <- gam(y ~ group + s(x, by = group), data = dat, method = "REML")

mod <- gam(y ~ group + s(x, by = group, k=5), data = dat_sim, method = "REML")

##
## difference smooths with pointwise CI ----
##

diff_df <- gratia::difference_smooths(mod, smooth = "s(x)", unconditional = TRUE, frequentist = FALSE)
ggplot(diff_df, aes(x, diff)) +
    geom_line() +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper),
                alpha = 0.5, color = "grey60")



##
## Plot each interval with pointwise or simultaneous CI ----
##


## point-wise interval
ci <- confint(mod, parm = c("s(x):group1", "s(x):group2"), type = "confidence")
## simultaneous interval
si <- confint(mod, parm = "s(x)", type = "simultaneous", partial_match = TRUE)




p1 <- ggplot(ci, aes(x = x, y = est, group = smooth)) +
    geom_line(lwd = 1) +
    geom_ribbon(data = ci, mapping = aes(ymin = lower, ymax = upper, x = x, group = smooth),
                fill = "grey60", inherit.aes = FALSE, alpha = 0.5) +
    geom_ribbon(data = si,
                mapping = aes(ymin = lower, ymax = upper, x = x, group = smooth),
                fill = "grey80", inherit.aes = FALSE, alpha = 0.5)

p1 #this plot here produces confidence intervals that appear to be shifted, why?




##
## Generate Simultaneous CI for difference in smooths using custom function
##


diff_df <- difference_smooths_R(mod, smooth = "s(x)", newdata = newdat,
                     unconditional = TRUE, frequentist = FALSE,
                     n=20, partial_match = TRUE, nrep=10000)

p2 <- ggplot() +
  geom_line(data = diff_df, aes(x = x, y = diff)) +
  geom_ribbon(data = diff_df, aes(x = x, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey40", inherit.aes = FALSE) +
  geom_ribbon(data = diff_df, aes(x = x, ymin = lower_s, ymax = upper_s),
              alpha = 0.5, fill = "grey60", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, lty = 2, color = "red")


p1 + p2
























##
## Older code used for exploration and to figure out the simultaneous CI (DON'T RUN) ----
##

## Trying simultaneous interval "by hand"
Vb <- vcov(mod)
n_pred <- 2000
newdat <- data.frame(x = seq(1, n, length.out = n_pred), group = factor(rep(c(1, 2), each = n_pred)))
pred <- predict(mod, newdat, se.fit = TRUE)
se.fit <- pred$se.fit
diffs <- pred$fit[1:n_pred] - pred$fit[(n_pred+1):(2*n_pred)]

set.seed(42)
N <- 10000


## pointwise CI using simulation
BUdiff <- t(rmvn(N, rep(0, nrow(Vb)), Vb))
rownames(BUdiff) <- rownames(Vb)
Cg <- predict(mod, newdat, type = "lpmatrix")
simDev <- Cg %*% BUdiff
idx_intercept <- grepl("Intercept", rownames(Vb))
idx_sx1 <- grepl("s\\(x\\):group1", rownames(Vb))
idx_sx2 <- grepl("s\\(x\\):group2", rownames(Vb))
Cg1 <- Cg[1:n_pred, idx_sx1]
Cg2 <- Cg[(n_pred+1):nrow(Cg), idx_sx2]
simDev <- BUdiff[idx_intercept, drop = FALSE] + # includes the mean
  Cg1 %*% BUdiff[idx_sx1, ] -
  Cg2 %*% BUdiff[idx_sx2, ]
str(simDev)

diffs <- Cg1 %*% coef(mod)[idx_sx1] - Cg2 %*% coef(mod)[idx_sx2]
bounds <- apply(simDev, 1, function(x) quantile(x, prob = c(0.025, 0.975)))

df_bounds <- data.frame(x = filter(newdat, group == 1)$x,
                        lower = diffs + bounds[1, ],
                        upper = diffs + bounds[2, ])

n_samples <- 50
posterior <- data.frame(x = filter(newdat, group == 1)$x,
                        curve = rep(diffs, n_samples) + c(simDev[, sample(ncol(simDev), n_samples)]),
                        iteration = rep(1:n_samples, each = n_pred))

##
## Plot of simultaneous intervals (missing mean shift)
##

ggplot() +
  geom_line(data = diff_df, aes(x = x, y = diff)) +
  geom_ribbon(data = diff_df, aes(x = x, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey60") +
  geom_ribbon(data = df_bounds, aes(x = x, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey80", inherit.aes = FALSE) +
  geom_line(data = posterior, aes(x = x, y = curve, group = iteration),
            color = "red", alpha = 0.15, inherit.aes = FALSE)




## Simultaneous interval???
diff_df <- gratia::difference_smooths(mod, smooth = "s(x)", newdata = newdat,
                                      unconditional = TRUE, frequentist = FALSE,
                                      n=20)
absDev <- abs(sweep(simDev, 1, diff_df$se, FUN = "/"))
masd <- apply(absDev, 2, max)

crit <- quantile(masd, prob = 0.95, type = 8)


diff_df <- diff_df %>%
  mutate(crit = crit) %>%
  mutate(lower_s = diff - crit * se,
         upper_s = diff + crit * se)
ggplot() +
  geom_line(data = diff_df, aes(x = x, y = diff)) +
  geom_ribbon(data = diff_df, aes(x = x, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey20") +
  geom_ribbon(data = df_bounds, aes(x = x, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey40", inherit.aes = FALSE) +
  geom_ribbon(data = diff_df, aes(x = x, ymin = lower_s, ymax = upper_s),
              alpha = 0.5, fill = "grey60", inherit.aes = FALSE) +
  geom_line(data = posterior, aes(x = x, y = curve, group = iteration),
            color = "red", alpha = 0.15, inherit.aes = FALSE)







