
## This function plots the rm-ANOVA and LMEM for the data simulated in example.R
plot_example <- function(sim_dat) {
    # Plot the simulated data (scatterplot)
    p1 <- sim_dat$dat %>%
        ggplot(aes(x = time,
                   y = y,
                   group = treatment,
                   color = treatment)
        ) +
        geom_point(alpha=0.5,
                   show.legend=FALSE) +
        labs(y='response')+
        geom_line(aes(x = time,
                      y = mu,
                      color = treatment),
                  size=3,
                  show.legend=FALSE) +
        theme_classic() +
        theme(plot.title = element_text(size = 20,
                                        face = "bold"),
              text=element_text(size=20))+
        thm1


    #plot the model predictions for rm-ANOVA
    p2 <- ggplot(sim_dat$dat,
                 aes(x = time,
                     y = y,
                     color = treatment)) +
        geom_point(alpha=0.5,
                   show.legend=FALSE)+
        labs(y='response')+
        geom_line(aes(y = predict(sim_dat$fit_anova),
                      group = subject, size = "Subjects"),
                  show.legend = FALSE) +
        geom_line(data = sim_dat$pred_dat,
                  aes(y = predict(sim_dat$fit_anova,
                                  level = 0,
                                  newdata = sim_dat$pred_dat),
                      size = "Population"),
                  show.legend=FALSE) +
        guides(color = guide_legend(override.aes = list(size = 2)))+
        scale_size_manual(name = "Predictions",
                          values=c("Subjects" = 0.5, "Population" = 3)) +
        theme_classic() +
        theme(plot.title = element_text(size = 20,
                                        face = "bold"),
              text=element_text(size=20))+
        thm1


    #plot the model predictions for LMEM
    p4 <- ggplot(sim_dat$dat,
                 aes(x = time,
                     y = y,
                     color = treatment)) +
        geom_point(alpha=0.5)+
        labs(y='response')+
        geom_line(aes(y = predict(sim_dat$fit_lme),
                      group = subject, size = "Subjects")) +
        geom_line(data = sim_dat$pred_dat,
                  aes(y = predict(sim_dat$fit_lme,
                                  level = 0,
                                  newdata = sim_dat$pred_dat),
                      size = "Population")) +
        guides(color = guide_legend(override.aes = list(size = 2)))+
        scale_size_manual(name = "Predictions",
                          values=c("Subjects" = 0.5, "Population" = 3)) +
        theme_classic() +
        theme(plot.title = element_text(size = 20,
                                        face = "bold"),
              text=element_text(size=20))+
        thm1

    return((p1+p2+p4)+plot_layout(nrow=1)+plot_annotation(tag_levels = 'A'))

}
#A<-plot_example(example(fun_type = "linear", error_type = "correlated"))

#C<-plot_example(example(fun_type = "quadratic", error_type = "correlated"))
