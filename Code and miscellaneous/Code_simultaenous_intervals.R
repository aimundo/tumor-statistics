set.seed(2021)

thm1<-scale_fill_scico_d(palette="tokyo",begin=0.3, end=0.8, direction = -1, aesthetics = c("colour","fill"))

#data, "the truth"
dat<-tibble(StO2=c(4,27,3,2,0.5,7,4,50,45,56),
            Day=rep(c(0,2,5,7,10),times=2),
            Group=as.factor(rep(c("Control","Treatment"),each=5))
)

#simulate data
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

#GAM, interaction of time and treatment

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


##Here I follow the exact code from G. Simpson's post to generate the simutaneous interval:
#https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

#Generate random values from multivariate normal
rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}


#Bayesian covariance matrix of the model coefficients
Vb <- vcov(gam_02)

#grid of values for confidence band
newd <- expand_grid(Group = factor(c("Control", "Treatment")),
                                   Day = seq(0, 10, by = 0.1),
                                   subject=factor(rep(1:10)))

#Generate predictions and standard errors
pred <- predict(gam_02, newd, se.fit = TRUE)
se.fit <- pred$se.fit

#number of draws
N <- 10000

#bias of the estimated model coefficients
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)

#evaluate the basis functions at the locations of x
Cg <- predict(gam_02, newd, type = "lpmatrix")

#Computes deviations between the fitted and true parameters
simDev <- Cg %*% t(BUdiff)

#compute absolute values of the standardized deviations from the true model
absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))

#find maximum of the absolute standardized deviations at the grid of values of x for each simulation
masd <- apply(absDev, 2L, max)

#Find critical value used to scale the standard errors to yield the simultaneous interval

crit <- quantile(masd, prob = 0.95, type = 8)




pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))

#plot the simultaneous and pointwise intervals
ggplot(pred, aes(x =Day,group=Group,color=Group,fill=Group)) +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
    geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2) +
    labs(y = "Tumor data",
         x = "Days")+
    thm1


#Below is the modified version of G. Simpson's code for the pairwise comparisons (pointwise)
#https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/

##calculate simultaneous diff

    #xp <- predict(model, newdata = newdata, type = 'lpmatrix')
    #Find columns in xp where the name contains "Control" and "Treatment"
    col1 <- grepl("Control", colnames(Cg))
    col2 <- grepl("Treatment", colnames(Cg))
    #Find rows in xp that correspond to each treatment
    row1 <- with(newd, Group == "Control")
    row2 <- with(newd, Group == "Treatment")
    ## difference rows of xp for data from comparison
    X <- Cg[row1, ] - Cg[row2, ]
    ## zero out cols of X not involved in the comparison
    X[, ! (col1 | col2)] <- 0
    ## zero out the parametric cols
    #X[, !grepl('^s\\(', colnames(Cg))] <- 0
    dif <- X %*% coef(gam_02)
    se <- sqrt(rowSums((X %*% vcov(gam_02, unconditional = FALSE)) * X))
    crit <- qt(0.05/2, df.residual(gam_02), lower.tail = FALSE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
   comp1<- data.frame(pair = paste("Control", "Treatment", sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)




comp_StO2_full <- cbind(Day = seq(0, 10, length = 1010),
                        rbind(comp1)) %>%
    mutate(interval=case_when(
        upper>0 & lower<0~"no-diff",
        upper<0~"less",
        lower>0~"greater"
    ))




c1<-ggplot(comp_StO2_full, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.5,
                fill="gray") +
    geom_line(data=comp_StO2_full,aes(y=0),size=0.5)+
    geom_line(color='black',size=1) +

    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )
c1
