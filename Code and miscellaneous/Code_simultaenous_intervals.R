set.seed(2021)

thm1<-scale_fill_scico_d(palette="tokyo",begin=0.3, end=0.8, direction = -1, aesthetics = c("colour","fill"))

#data
dat<-tibble(StO2=c(4,27,3,2,0.5,7,4,50,45,56),
            Day=rep(c(0,2,5,7,10),times=2),
            Group=as.factor(rep(c("Control","Treatment"),each=5))
)
#alpha for ribbon
al<-0.8

simulate_data <- function(dat, n = 10, sd = 5) {
    dat_sim <- dat %>%
        slice(rep(1:n(), each = n)) %>%
        group_by(Group, Day) %>%
        mutate(
            StO2_sim = pmax(rnorm(n, StO2, sd), 0.0001),
            subject=rep(1:10),
            subject=factor(paste(subject, Group, sep = "-"))
        ) %>%
        ungroup()

    return(dat_sim)
}


#subject = factor(paste(subject, treatment, sep = "-")))

n <- 10 #number of observations
sd <- 10 #approximate sd from paper
df <- 6
dat_sim <- simulate_data(dat, n, sd)

#gam

gam_02 <- gam(StO2_sim ~ Group+s(Day, by = Group, k = 5),
              method='REML',
              data  = dat_sim)


#creates a dataframe using the length of the covariates for the GAM
gam_predict <- expand_grid(Group = factor(c("Control", "Treatment")),
                           Day = seq(0, 10, by = 0.1),
                           subject=factor(rep(1:10)))

#adds the predictions to the grid and creates a confidence interval for GAM
gam_predict<-gam_predict%>%
    mutate(fit = predict(gam_02,gam_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(gam_02, gam_predict,se.fit = TRUE,type='response')$se.fit)




rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}



Vb <- vcov(gam_02)
newd <- expand_grid(Group = factor(c("Control", "Treatment")),
                                   Day = seq(0, 10, by = 0.1),
                                   subject=factor(rep(1:10)))
pred <- predict(gam_02, newd, se.fit = TRUE)
se.fit <- pred$se.fit


N <- 10000


BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

Cg <- predict(gam_02, newd, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))

masd <- apply(absDev, 2L, max)



crit <- quantile(masd, prob = 0.95, type = 8)




pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))

ggplot(pred, aes(x =Day,group=Group,color=Group,fill=Group)) +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
    geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
    labs(y = "Tumor data",
         x = "Days")+
    thm1





sims <- rmvn(N, mu = coef(gam_02), sig = Vb)
fits <- Cg %*% t(sims)



nrnd <- 50
rnd <- sample(N, nrnd)
stackFits <- stack(as.data.frame(fits[, rnd]))
stackFits <- transform(stackFits, Day = rep(newd$Day, length(rnd)))

stackFits$Group<-rep(newd$Group,length(rnd))



ggplot(pred, aes(x = Day, y = fit,group=Group,color=Group,fill=Group)) +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
    geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
    geom_path(lwd = 2) +
    geom_path(data = stackFits, aes(y = values, x = Day),
              alpha = 0.4, colour = "grey20") +
    labs(y = "StO2",
         x = "Day",
         title = "Point-wise & Simultaneous 95% confidence intervals for fitted GAM",
         subtitle = sprintf("Each line is one of %i draws \nfrom the Bayesian posterior distribution of the model", nrnd))+
    thm1


