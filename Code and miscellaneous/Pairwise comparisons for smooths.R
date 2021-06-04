#This post contains how to perform calculations in a gam model, when testing for
#significance
#https://stats.stackexchange.com/questions/376237/correcting-for-multiple-pairwise-comparisons-with-gam-objects-mgcv-in-r

set.seed(123)
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


m1 <- gam(StO2_sim ~ Group+s(Day, by = Group, k = 5,bs="gp"),
            method='REML',
            data  = dat_sim)

######missing data###########
#############################
missing <- sample(1:100, 40)

#create a new dataframe from the simulated data with 40 rows randomly removed, keep the missing values as NA

ind <- which(dat_sim$StO2_sim %in% sample(dat_sim$StO2_sim, 40))

#create a new dataframe, remove the StO2 column
dat_missing <- dat_sim[,-1]

#add NAs at the ind positions
dat_missing$StO2_sim[ind]<-NA

#Count the number of remaining observations per day (original dataset had 10 per group per day)
dat_missing %>%
    group_by(Day,Group) %>%
    filter(!is.na(StO2_sim))%>%
    count(Day)


#the same model used for the full dataset
mod_m1 <- gam(StO2_sim ~ Group+s(Day,by=Group,k=5), data  = dat_missing)
#appraise the model
appraise(mod_m1)


m_predict <- expand_grid(Group = factor(c("Control", "Treatment")),
                         Day = seq(0, 10, by = 0.1))

#adds the predictions to the grid and creates a confidence interval
m_predict<-m_predict%>%
    mutate(fit = predict(mod_m1,m_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(mod_m1, m_predict,se.fit = TRUE,type='response')$se.fit)

##pairwise comparisons###
smooth_diff <- function(model, newdata, g1, g2, alpha = 0.05,
                        unconditional = FALSE) {
    xp <- predict(model, newdata = newdata, type = 'lpmatrix')
    #Find columns in xp where the name contains "Control" and "Treatment"
    col1 <- grepl(g1, colnames(xp))
    col2 <- grepl(g2, colnames(xp))
    #Find rows in xp that correspond to each treatment
    row1 <- with(newdata, Group == g1)
    row2 <- with(newdata, Group == g2)
    ## difference rows of xp for data from comparison
    X <- xp[row1, ] - xp[row2, ]
    ## zero out cols of X related to splines for other lochs
    X[, ! (col1 | col2)] <- 0
    ## zero out the parametric cols
    #X[, !grepl('^s\\(', colnames(xp))] <- 0
    dif <- X %*% coef(model)
    se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
    crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    data.frame(pair = paste(g1, g2, sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)
}

pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Control', 'Treatment'))


comp1<-smooth_diff(m1,pdat,'Control','Treatment')

comp_StO2_full <- cbind(Day = seq(0, 10, length = 400),
                        rbind(comp1)) %>%
    mutate(interval=case_when(
        upper>0 & lower<0~"no-diff",
        upper<0~"less",
        lower>0~"greater"
    ))

###plot for full data####
#function to extract limits for the shades
pairwise_limits<-function(dataframe){
    #extract values where the lower limit of the ribbon is greater than zero
    #this is the region where the control group effect is greater
    v1<-dataframe%>%
        filter(lower>0)%>%
        select(Day)
    #get day  initial value
    init1=v1$Day[[1]]
    #get day final value
    final1=v1$Day[[nrow(v1)]]

    #extract values where the value of the upper limit of the ribbon is lower than zero
    #this corresponds to the region where the treatment group effect is greater
    v2<-comp_StO2_full%>%
        filter(upper<0)%>%
        select(Day)

    init2=v2$Day[[1]]
    final2=v2$Day[[nrow(v2)]]
    #store values
    my_list<-list(init1=init1,
                  final1=final1,
                  init2=init2,
                  final2=final2)
return(my_list)
}


my_list<-pairwise_limits(comp_StO2_full)

c2<-ggplot(comp_StO2_full, aes(x = Day, y = diff, group = pair)) +
    annotate("rect",
                xmin =my_list$init1, xmax =my_list$final1,ymin=-Inf,ymax=Inf,
                fill='blue',
                alpha = 0.1) +
    annotate("text",
             x=1.5,
             y=-10,
             label="Control",size=10)+
    annotate("rect",
             xmin =my_list$init2, xmax =my_list$final2,ymin=-Inf,ymax=Inf,
             fill='red',
             alpha = 0.1,
    ) +
    annotate("text",
             x=6,
             y=-10,
             label="Treatment",
             size=10)+
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.9,
                fill='black') +
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

c2


####plot for missing data####
comp2<-smooth_diff(mod_m1,pdat,'Control','Treatment')

comp_StO2_missing <- cbind(Day = seq(0, 10, length = 400),
                        rbind(comp2)) %>%
    mutate(interval=case_when(
        upper>0 & lower<0~"no-diff",
        upper<0~"less",
        lower>0~"greater"
    ))

my_list<-pairwise_limits(comp_StO2_missing)

c3<-ggplot(comp_StO2_missing, aes(x = Day, y = diff, group = pair)) +
    annotate("rect",
             xmin =my_list$init1, xmax =my_list$final1,ymin=-Inf,ymax=Inf,
             fill='blue',
             alpha = 0.1,
    ) +
    annotate("rect",
             xmin =my_list$init2, xmax =my_list$final2,ymin=-Inf,ymax=Inf,
             fill='red',
             alpha = 0.1,
    ) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.9,
                fill='black') +
    geom_line(data=comp_StO2_missing,aes(y=0),size=0.5)+
    geom_line(color='black',size=1) +
    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )

c3
#######################

##original code for the comparisons, from Gavin Simpson's blog.
pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Treatment', 'Control'))

smooth_diff <- function(model, newdata, g1, g2, alpha = 0.05,
                        unconditional = FALSE) {
    xp <- predict(model, newdata = newdata, type = 'lpmatrix')
    c1 <- grepl(g1, colnames(xp))
    c2 <- grepl(g2, colnames(xp))
    #r1 <- newdata[[var]] == f1
    #r2 <- newdata[[var]] == f2
    r1 <- with(newdata, Group == g1)
    r2 <- with(newdata, Group == g2)
    ## difference rows of xp for data from comparison
    X <- xp[r1, ] - xp[r2, ]
    ## zero out cols of X related to splines for other lochs
    X[, ! (c1 | c2)] <- 0
    ## zero out the parametric cols
    #X[, !grepl('^s\\(', colnames(xp))] <- 0
    dif <- X %*% coef(model)
    se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
    crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    data.frame(pair = paste(g1, g2, sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)
}

comp1<-smooth_diff(m1,pdat,'Control','Treatment')

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

