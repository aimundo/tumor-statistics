#This post contains how to perform calculations in a gam model, when testing for
#significance
#https://stats.stackexchange.com/questions/376237/correcting-for-multiple-pairwise-comparisons-with-gam-objects-mgcv-in-r


dat<-tibble(StO2=c(4,27,3,2,0.5,7,4,50,45,56),
            Day=rep(c(0,2,5,7,10),times=2),
            Group=as.factor(rep(c("Control","Treatment"),each=5))
)



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
set.seed(1) #set seed for reproducibility
df <- 6
dat_sim <- simulate_data(dat, n, sd)


gam1 <- gam(StO2_sim ~ Group+s(Day, by = Group, k = 5,bs="gp"),
            method='REML',
            data  = dat_sim)

#creates fine grid at the values of Day
pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Control', 'Treatment'))

##matrix that contains the basis functions evaluated at the points in pdat
    xp <- predict(gam1, newdata = pdat, type = 'lpmatrix')

    #remove intercept column and group identifying column
    xp<-xp[,-c(1:2)]

    control<-xp[1:400,]
    treatment<-xp[401:800,]


#Find columns in xp where the name contains "Control"
    c1 <- grepl('Control', colnames(xp))

#Find columns in xp where the name contains 'Treatment'
    c2 <- grepl('Treatment', colnames(xp))

#Find rows in pdat that correspond to either 'Control' or 'Treatment'
    r1 <- with(pdat, Group == 'Control')
    r2 <- with(pdat, Group == 'Treatment')

# In xp: find the rows that correspond to Control or Treatment, those that do not match will be
    #set to zero. Then, substract the values from the rows corresponding to 'Control' from those that correspond
    #to 'Treatment'
    X <- xp[r1, ] - xp[r2, ]

    test<-control-treatment


    ## remove columns that do not contain name 'Control' or 'Treatment'
    X[, ! (c1 | c2)] <- 0
    ## zero out the parametric cols, those that do not contain in the name 's()'
    X[, !grepl('^s\\(', colnames(xp))] <- 0

    #Multiply matrix by model coefficients. X has (p,n) (rows, columns) and the coefficient matrix has
    #dimensions (n,1). The resulting matrix has dimensions (p,1)
    dif <- X %*% coef(gam1)

    comp<-test %*% coef(gam1)[3:10]


    se <- sqrt(rowSums((X %*% vcov(gam1, unconditional = FALSE)) * X))
    crit <- qt(0.05/2, df.residual(gam1), lower.tail = FALSE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    comp1<-data.frame(pair = paste('Control', 'Treatment', sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)




comp_StO2 <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp1))


c1<-ggplot(comp_StO2, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.5,
                fill='black') +
    geom_line(color='black',size=1) +
    geom_line(data=comp_StO2,aes(y=0),size=0.5)+
    geom_ribbon(data=comp_StO2%>%
                    filter(lower>0),
                aes(ymin =0, ymax =lower),
                alpha = 0.5,
                fill='blue') +
    geom_ribbon(data=comp_StO2 %>%
                    filter(upper<0),
                    aes(ymin =0, ymax =upper),
                alpha = 0.5,
                fill='red') +
    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )


c1



#######################

##original code for the comparisons, from Gavin Simpson's blog.
pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Treatment', 'Control'))

smooth_diff <- function(model, newdata, g1, g2, alpha = 0.05,
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

comp1<-smooth_diff(gam1,pdat,'Control','Treatment')

comp_StO2 <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp1)) %>%
    mutate(interval=case_when(
        upper>0 & lower<0~"no-diff",
        upper<0~"less",
        lower>0~"greater"
    ))

c1<-ggplot(comp_StO2, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5,fill='#993336') +
    geom_line(color='black',size=1) +
    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )
