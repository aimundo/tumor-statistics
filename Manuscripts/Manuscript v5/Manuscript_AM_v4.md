---
title: '**The statistical analysis of non-linear longitudinal data in biomedical research using generalized additive models**'
subtitle: _Beyond repeated measures ANOVA and Linear Mixed Models_
author:
- Ariel Mundo^[Department of Biomedical Engineering, University of Arkansas, Fayetteville]
- Timothy J. Muldoon^[Department of Biomedical Engineering, University of Arkansas,
  Fayetteville]
- John R. Tipton^[Department of Mathematical Sciences, University of Arkansas, Fayetteville]
output:
  bookdown::word_document2:
    fig_caption: yes
    keep_md: yes
  bookdown::pdf_document2:
    keep_tex: yes
    fig_caption: yes
    extra_dependencies:
      subfig: null
      breqn: null
      caption: ["font={small}"]
      float: null
  bookdown::html_document2: default
  link-citations: yes
  css: style.css
  '': default
csl: elsevier-with-titles.csl
bibliography: refs.bib
---



# Background

Longitudinal studies are designed to repeatedly measure a variable of interest in a group (or groups) of subjects, with the intention of observing the evolution of effect across time rather than analyzing a single time point (e.g., a cross-sectional study). Biomedical research frequently uses longitudinal studies to analyze the evolution of a "treatment" effect across multiple time points; and in such studies the subjects of analysis range from animals (mice, rats, rabbits), to human patients, cells, or blood samples, among many others. Tumor response [@roblyer2011;@tank2020;@pavlov2018;@demidov2018], antibody expression [@ritter2001;@roth2017], and cell metabolism [@jones2018;@skala2010] are examples of the different situations where researchers have used longitudinal designs to study some  physiological response. Because the frequency of the measurements in a longitudinal study is dependent on the biological phenomena of interest and the experimental design of the study, the frequency of such measurements can range from minute intervals to study a short-term response such as anesthesia effects in animals[@greening2018], to weekly measurements to analyze a mid-term response like the evolution of dermatitis symptoms in breast cancer patients [@sio2016], to monthly measurements to study a long-term response such as mouth opening following RT in neck cancer patients [@kamstra2015]. 


Traditionally, a “frequentist” or "classical" statistical paradigm is used in biomedical research to derive inferences from a longitudinal study. The frequentist paradigm regards probability as a limiting frequency [@wagenmakers2008] by assuming a null hypothesis under a statistical model that is often an _analysis of variance over repeated measures_ (repeated measures ANOVA or rm-ANOVA). The rm-ANOVA model makes three key assumptions regarding longitudinal data: 1) linearity of the response across time, 2) constant correlation across same-subject measurements, and 3) observations from each subject are obtained at all time points through the study (a condition also known as _complete observations_) [@gueorguieva2004;@schober2018]. 

The expected linear behavior of the response through time is a key requisite in rm-ANOVA [@pinheiro2006]. This "linearity assumption" in rm-ANOVA implies that the model is misspecified when the data does not follow a linear trend, which results in unreliable inference. In biomedical research, non-linear trends are the norm rather than the exception in longitudinal studies. A particular example of this non-linear behavior in longitudinal data arises in measurements of tumor response in preclinical and clinical settings [@roblyer2011;@skala2010;@vishwanath2009]. These studies have shown that the collected signal does not follow a linear trend over time, and presents extreme variability at different time points, making the fit of rm-ANOVA model inconsistent with the observed variation. Therefore, when rm-ANOVA is used to draw inference of such highly-variable data the estimates are inevitably biased, because the model is only able to accommodate linear trends that are far from adequately representing the biological phenomenon of interest.


A _post hoc_ analysis is the statistical test used in conjunction with rm-ANOVA to perform repeated comparisons to estimate a _p-value_, which in turn is used as a measure of significance.
Although it is possible that a _post hoc_ analysis of rm-ANOVA is able to find “significant” _p-values_( _p_<0.05) from non-linear data, the validity of such metric is dependent on how adequate the model fits the data. In other words, _p-values_ are valid only if the model and the data have good agreement; if that is not the case, a "Type III" error (known as "model misspecification") occurs[@dennis2019]. For example, model misspecification will occur when a model that is only able to explain linear responses (such as rm-ANOVA) is fitted to data that follows a quadratic trend, thereby causing the resulting _p-values_ and parameter estimates to be invalid [@wang2019].

Additionally, the _p-value_ itself is highly variable, and multiple comparisons can inflate the false positivity rate (Type I error or $\alpha$) [@liu2010;@halsey2015], consequently biasing the conclusions of the study. Corrections exist to address the Type I error issue of multiple comparisons (such as Bonferroni [@abdi2010]), but they in turn reduce statistical power (1-$\beta$)[@nakagawa2004], and lead to increased Type II error (failing to reject the null hypothesis when the null hypothesis is false) [@gelman2012;@albers2019]. Therefore, the tradeoff of _post hoc_ comparisons in rm-ANOVA between Type I, II and III errors might be difficult to balance in a biomedical longitudinal study where a delicate balance exists between statistical power and sample size.


On the other hand, the assumption of constant correlation in rm-ANOVA (often known as the _compound symmetry assumption_) is typically unreasonable because correlation between the measured responses often diminishes as the time interval between the observation increases [@ugrinowitsch2004].  Corrections can be made in rm-ANOVA in the absence of compound symmetry [@huynh1976;@greenhouse1959], but the effectiveness of the correction is limited by the size of the sample, the number of measurements[@haverkamp2017], and group sizes [@keselman2001]. In the case of biomedical research, where living subjects are frequently used, sample sizes are often not "large" due to ethical and budgetary reasons [@charan2013] which might cause the corrections for lack of compound symmetry to be ineffective.

Due to a variety of causes, the number of observations during a study can vary between all subjects. For example, in a clinical trial patients may voluntarily withdraw, whereas attrition  due to injury or weight loss in preclinical animal studies is possible. It is even plausible that unexpected complications with equipment or supplies arise that prevent the researcher from collecting measurements at certain time points. In each of these missing data scenarios,the _complete observations_ assumption of  classical rm-ANOVA is violated. When incomplete observations occur, a rm-ANOVA model is fit by excluding all subjects with missing observations from the analysis [@gueorguieva2004]. This elimination of partially missing data from the analysis can result in increased costs if the desired statistical power is not met with the remaining observations, because it would be necessary to enroll more subjects. At the same time, if the excluded observations contain insightful information that is not used, their elimination from the analysis may limit the demonstration of significant differences between groups. 

During the last decade, the biomedical community has started to recognize the limitations of rm-ANOVA in the analysis of longitudinal information. The recognition on the shortcomings of rm-ANOVA is exemplified by the  use of  linear mixed effects models (LMEMs) by certain groups to analyze longitudinal tumor response data [@skala2010;@vishwanath2009]. Briefly, LMEMs incorporate _fixed effects_, which correspond to the levels of experimental factors in the study (e.g., the different drug regimens in a clinical trial), and _random effects_, which account for random variation within the population (e.g., the individual-level differences not due to treatment such as weight or age). When compared to the traditional rm-ANOVA, LMEMs are more flexible as they can accommodate missing observations for multiple subjects and allow different modeling strategies for the variability within each measure in every subject [@pinheiro2006]. However, LMEMs impose restrictions in the distribution of the errors  of the random effects, which need to be normally distributed and independent [@gueorguieva2004;@barr2013]. And even more importantly, LMEMs also expect a linear relationship between the response and time [@pinheiro2006], making them unsuitable to analyze non-linear data.

As the rm-ANOVA and the more flexible LMEM approaches make overly restrictive assumptions regarding the linearity of the response, there is a need for biomedical researchers to explore the use of additional statistical tools that allow the data (and not an assumption in trend) to determine the trend of the fitted model, to enable appropriate inference.

In this regard, generalized additive models (GAMs) present an alternative approach to analyze longitudinal data. Although not frequently used by the biomedical community, these semi-parametric models are customarily used in other fields to analyze longitudinal data.   Examples of the use of GAMs include the analysis of temporal variations in geochemical and palaeoecological data  [@rose2012;@pedersen2019;@simpson2018], health-environment interactions [@yang2012] and the dynamics of government in political science [@beck1998] . There are several advantages of GAMs over LMEMs and rm-ANOVA models: 1) GAMs can fit a more flexible class of smooth responses that enable the data to dictate the trend in the fit of the model, 2) they can model non-constant correlation between repeated measurements [@wood2017] and 3) can easily accommodate missing observations. Therefore, GAMs can provide a more flexible statistical approach to analyze non-linear biomedical longitudinal data than LMEMs and rm-ANOVA.

The current advances in programming languages designed for statistical analysis (specifically $\textsf{R}$), have eased the computational implementation of traditional models such as rm-ANOVA and more complex approaches such as LMEMs and GAMs. In particular, $\textsf{R}$[@r] has an extensive collection of documentation and functions to fit GAMs in the package _mgcv_ [@wood2016;@wood2017] that not only speed up the initial stages of the analysis but also enable the use of advanced modeling structures (e.g. hierarchical models, confidence interval comparisons) without requiring advanced programming skills from the user. At the same time, $\textsf{R}$ has many tools that simplify data simulation, an emerging strategy used to test statistical models [@haverkamp2017]. Data simulation methods allow the researcher to create and  explore different alternatives for analysis without collecting information in the field, reducing the time window between experiment design and its implementation, and simulation can be also used for power calculations and study design questions.   

This work provides biomedical researchers with a clear understanding of the theory and the practice of using GAMs to analyze longitudinal data using by focusing on four areas. First, the limitations of LMEMs and rm-ANOVA regarding linearity of response, constant correlation structures and missing observations is explained in detail. Second, the key theoretical elements of GAMs are presented using clear and simple mathematical notation while explaining the context and interpretation of the equations. Third, using simulated data that reproduces patterns in previously reported studies [@vishwanath2009] we illustrate the type of non-linear longitudinal data that often occurs in biomedical research. The simulated data experiments highlight the differences in inference between rm-ANOVA, LMEMs and GAMs on data similar to what is commonly observed in biomedical studies. Finally, reproducibility is emphasized by providing the code to generate the simulated data and the implementation of different models in $\textsf{R}$, in conjunction with a step-by-step guide demonstrating how to fit models of increasing complexity.  

In summary, this work will allow biomedical researchers to identify when the use of GAMs instead of rm-ANOVA or LMEMs is appropriate to analyze longitudinal data, and provide guidance on the implementation of these models by improving the standards for reproducibility  in biomedical research.


# Challenges presented by longitudinal studies

## The repeated measures ANOVA

The _repeated measures analysis of variance_ (rm-ANOVA) is the standard statistical analysis for longitudinal data in biomedical research. This statistical methodology requires certain assumptions for the model to be valid. From a practical view, the assumptions can be divided in three areas: 1) linear relationship between covariates and response, 2) a constant correlation between measurements, and, 3) complete observations for all subjects. Each one of these assumptions is discussed below.
  
## Linear relationship 

### The repeated measures ANOVA case

In a biomedical longitudinal study, two or more groups of subjects (e.g., patients, mice, samples) are subject to different treatments (e.g., a "treatment" group receives a novel drug vs. a "control" group that receives a placebo), and measurements from each subject within each group are collected at specific time points. The collected response is modeled with _fixed_ components. The _fixed_ component can be understood as a constant value in the response which the researcher is interested in measuring, i.e., the average effect of the novel drug in the "treatment" group. 


Mathematically speaking, a rm-ANOVA model with an interaction can be written as:

\begin{equation}
y_{ijt} = \beta_0+\beta_1 \times time_{t} +\beta_2 \times treatment_{j} +\beta_3 \times time_{t}\times treatment_{j}+\varepsilon_{ijt}\\ 
(\#eq:linear-model)
\end{equation}

In this model $y_{ijt}$ is the response for subject $i$, in treatment group $j$ at time $t$, which can be decomposed in a mean value $\beta_0$, _fixed effects_ of time ($time_t$), treatment ($treatment_j$) and their interaction $time_t*treatment_j$ which have linear slopes given by $\beta_1, \beta_2$ and $\beta_3$, respectively. Independent errors $\varepsilon_{tij}$ represent random variation not explained by the _fixed_ effects, and are assumed to be $\sim N(0,\sigma^2)$.
In a  biomedical research context, suppose two treatments groups are used in a study (e.g., "placebo" vs. "novel drug" or "saline" vs. "chemotherapy"). Then, the group terms in Equation \@ref(eq:linear-model) can be written as below with $treatment_j=0$ representing the first treatment group (Group A) and $treatment_j=1$ representing the second treatment group (Group B). The linear models then can be expressed as


\begin{equation}
y_{ijt} = \begin{cases}
\beta_0 + \beta_1\times time_{t}+\varepsilon_{ijt}   & \mbox{if Group A}\\
\beta_0 + \beta_1 \times time_{t} + \beta_2+\beta_3 \times time_{t}+\varepsilon_{ijt}  & \mbox{if Group B}\\
\end{cases}
(\#eq:ANOVA-by-group)
\end{equation}

To further simplify the expression, substitute $\widetilde{\beta_{0}}=\beta_0+\beta_{2}$ and $\widetilde{\beta_{1}}=\beta_{1}+\beta_{3}$ in the equation for Group B. This substitution allows for a different intercept and slope for Groups A and B. The model is then written as

\begin{equation}
y_{ijt} = \begin{cases}
\beta_0 + \beta_1\times time_{t}+\varepsilon_{ijt}   & \mbox{if Group A}\\
\widetilde{\beta_{0}} + \widetilde{\beta_1} \times time_{t}+\varepsilon_{ijt}  & \mbox{if Group B}\\
\end{cases}
(\#eq:ANOVA-lines)
\end{equation}

Presenting the model in this manner makes clear that when treating different groups, an rm-ANOVA model is able to accommodate non-parallel lines in each case (different intercepts and slopes per group). In other words, the rm-ANOVA model "expects" a linear relationship between the covariates and the response, this means that either presented as Equation \@ref(eq:linear-model), Equation \@ref(eq:ANOVA-by-group) or Equation \@ref(eq:ANOVA-lines), an rm-ANOVA model is only able to accommodate linear patterns in the data. If the data show non-linear behavior, the rm-ANOVA model will approximate this behavior with non-parallel lines. 

### The Linear Mixed Model Case

A linear mixed model (LMEM) is a class of statistical model that incorporates _fixed effects_ to model the relationship between the covariates and the response, and _random effects_ to model subject variability that is not the primary focus of the study but that might be important to distinguish [@pinheiro2006;@west2014]. A LMEM with interaction between time and treatment for a longitudinal study  can be written as:


\begin{equation}
y_{ijt} = \beta_0+ \beta_1 \times time_{t} + \beta_2 \times treatment_{j} + \beta_3 \times time_{t}\times treatment_{j}+\mu_{ij} +\varepsilon_{ijt}\\ 
(\#eq:LMEM)
\end{equation}

When Equation \@ref(eq:linear-model) and Equation \@ref(eq:LMEM) are compared, it is easily noticeable that LMEM and rm-ANOVA have the same construction regarding the _fixed effects_ of time and treatment, but that the LMEM incorporates an additional source of variation (the term $\mu_{ij}$). This term $\mu_{ij}$ is the one that corresponds to the  _random effect_, accounting for variability in each subject within each group. The _random_ component can also be understood as used to model some "noise" in the response, but that is intended to be analyzed and disentangled from the "global noise" term $\varepsilon_{ijt}$ from Equation \@ref(eq:linear-model). 

For example, if the blood concentration of the drug is measured in certain subjects in the early hours of the morning while other subjects are measured in the afternoon, it is possible that the difference in the collection time introduces some "noise" in the data. As the name suggests, this "random" variability needs to be modeled as a variable rather than as a constant value.  The _random effect_ $\mu_{ij}$ in Equation \@ref(eq:LMEM) is assumed to be independently normally distributed with mean zero and variance $\sigma^2_\mu$, which can be expressed as $\mu_{ij} \sim N(0,\sigma^2_\mu)$. In essence,the _random effect_ in a LMEM enables to fit models with different slopes at the subject-level[@pinheiro2006]. However,the expected linear relationship of the covariates and the response in Equation \@ref(eq:linear-model) and in Equation \@ref(eq:LMEM) is essentially the same, representing a major limitation of LMEMs to fit a non-linear response.

## Covariance in rm-ANOVA and LMEMs

In a longitudinal study there is an expected _covariance_ between repeated measurements on the same subject, and because repeated measures occur in the subjects within each group, there is a _covariance_ between  measurements at each time point within each group. The _covariance matrix_ (also known as the variance-covariance matrix) is a matrix that captures the variation between and within subjects in a longitudinal study[@wolfinger1996] (For an in-depth analysis of the covariance matrix see [@west2014;@weiss2005]). 

In the case of an rm-ANOVA analysis, it is typically assumed that the covariance matrix has a specific construction known as _compound symmetry_ (also known as "sphericity" or "circularity"). Under this assumption, the between-subject variance and within-subject correlation  are constant across time [@weiss2005;@geisser1958;@huynh1976]. However, it has been shown that this condition is frequently not justified because the correlation between measurements tends to change over time [@maxwell2017]; and it is higher between consecutive measurements [@gueorguieva2004;@ugrinowitsch2004]. Although corrections can be made (such as Huyhn-Feldt or Greenhouse-Geisser)[@huynh1976;@greenhouse1959] the effectiveness of each correction is limited because it depends on the size of the sample,the number of repeated measurements[@haverkamp2017], and they are not robust if the group sizes are unbalanced [@keselman2001]. Because biomedical longitudinal studies are often limited in sample size and can have an imbalanced design, the corrections required to use an rm-ANOVA model may not be able to provide a reasonable adjustment that makes the model valid.


In the case of LMEMs, one key advantage over rm-ANOVA is that they allow different structures for the variance-covariance matrix including exponential, autoregressive of order 1, rational quadratic and others [@pinheiro2006]. Nevertheless, the analysis required to determine an appropriate variance-covariance structure for the data can be a long process by itself. Overall, the spherical assumption for rm-ANOVA may not capture the natural variations of the correlation in the data, and can bias the inferences from the analysis. 


## Missing observations

Missing observations are an issue that arises frequently in longitudinal studies. In biomedical research, this situation can be caused by reasons beyond the control of the investigator [@molenberghs2004]. Dropout from patients and attrition or injury in animals are among the reasons for missing observations. Statistically, missing information can be classified as _missing at random_ (MAR), _missing completely at random_ (MCAR), and _missing not at random_ (MNAR) [@weiss2005].  In a MAR scenario, the pattern of the missing information is related to some variable in the data, but it is not related to the variable of interest [@scheffer2002]. If the data are MCAR, this means that the missingness is completely unrelated to the collected information [@potthoff2006], and in the case of MNAR the missing values are dependent on their value. An rm-ANOVA model assumes complete observations for all subjects, and therefore subjects with one or more missing observations are excluded from the analysis. This is inconvenient because the remaining subjects might not accurately represent the population, and statistical power is affected by this reduction in sample size [@ma2012].

In the case of LMEMs, inferences from the model are valid when missing observations in the data exist that are MAR or MCAR [@west2014]. For example, if attrition occurs in all mice that had lower weights at the beginning of a chemotherapy response study, the missing data can be considered MAR because  the missigness is unrelated to other variables of interest.

This section has presented the assumptions for analyzing longitudinal data using rm-ANOVA and LMEMs and compared their differences regarding linearity, the covariance matrix and  missing data. In particular, LMEMs offer a more robust and flexible approach than rm-ANOVA and if the data follows a linear trend, they provide an excellent choice to derive inferences from a repeated measures study. However, when the data presents high a non-linear behavior, LMEMs and rm-ANOVA fail to capture the trend of the data. To better convey the issues of linearity and correlation  in linear models fitted to non-linear data, simulation is used in the next section.

## How does an rm-ANOVA fit looks like? A visual representation using simulation {#simulation}

To demonstrate the limitations of rm-ANOVA an LMEMs for non-linear longitudinal data, this section presents a simulation experiment of a normally distributed response of two groups of 10 subjects each. An rm-ANOVA model (Equation \@ref(eq:linear-model)) is fitted to each group, using $\textsf{R}$[@r] and the package _nlme_[@nlme]. 
Briefly, two cases for the mean responses for each group are considered: in the first case, the mean response in each group is a linear function with different intercepts and slopes; a negative slope is used for Group 1 and a positive slope is used for Group 2 (Figure \@ref(fig:l-q-response), A). In the second case, a second-degree polynomial (quadratic) function is used for the mean response per group: the quadratic function is concave down for Group 1 and it is concave up for Group 2 (Figure \@ref(fig:l-q-response), C). In both the linear and quadratic simulated data, the groups start with the same mean value at the first time point. This is intentional in order to simulate the expected temporal evolution of some physiological quantity.

Specifically, the rationale for the chosen linear and quadratic functions is the likelihood that a measured response in two treatment groups is similar in the initial phase of the study, but as treatment progresses a divergence in the trend of the response indicates a difference in the effect of each treatment. In other words, Group 1 can be thought as a "Control" group and Group 2 as a "Treatment" group. From the mean response per group (linear or quadratic), the variability or "error" of individual responses  within each group is simulated using a covariance matrix with compound symmetry (constant variance across time). Thus, the response per subject in both the linear and quadratic simulation corresponds to the mean response per group plus the error (Figure \@ref(fig:l-q-response) B,D). A more comprehensive exploration  of the fit of rm-ANOVA for linear and non-linear longitudinal data is in Figure \@ref(fig:linear-cases) and Figure \@ref(fig:quadratic-cases) in the Appendix, where simulation with compound symmetry and independent errors (errors generated from a normal distribution that are not constant over time) and the plot of simulated errors, and fitted parameters in presented.







![(\#fig:l-q-response)Simulated linear responses from two groups with correlated (top row) or independent (bottom row) errors using a rm-ANOVA model. A, C:Simulated data with known mean response (linear or quadratic, thin lines) and individual responses (points) showing the dispersion of the data. B,D: Estimates from the rm-ANOVA model for the mean group response (linear of quadratic). Thick lines are the predicted mean response per group, thin lines are the random effects for each subject and points represent the original raw data. The rm-ANOVA model not only fails to pick the trend of the quadratic data but it also incorrectly estimates the initial conditions.](Manuscript_AM_v4_files/figure-docx/l-q-response-1.png){width=75% }

The simulation shows that the fit produced by the rm-ANOVA model is good for linear data, as the predictions for the mean response are reasonably close to the "truth" of the simulated data  (Figure \@ref(fig:l-q-response),B). When the linearity and compound symmetry assumptions are met, the model approximates well the individual trends and the mean trends by group.

However, consider the case when the data follows a non-linear trend, such as the simulated data in Figure \@ref(fig:l-q-response), C. Here, the mean response per group was simulated using a quadratic function but errors, individual responses and the rm-ANOVA model were produced in the same manner as in (Figure \@ref(fig:l-q-response) A and B) . The mean response in the simulated data with quadratic behavior is changing in each group through the timeline, and the mean value is the same as the initial value by the fifth time point for each group. Fitting an rm-ANOVA model \@ref(eq:linear-model) to this data  produces the fit that appears in panel D in Figure \@ref(fig:l-q-response).

A comparison of the fitted mean response of the rm-ANOVA model to the the simulated data in Figure ((\@ref(fig:l-q-response), D) indicates that the model is not capturing the changes within each group in a good way.  Specifically, note that the fitted mean response of the rm-ANOVA model (panel D) shows that the change (increase for Treatment 1 or decrease for Treatment 2) in the response through time points 2 and 4 is not being captured by the model. Moreover, the rm-ANOVA model is not being able to capture the fact that the initial values are the same in each group, and instead fits non-parallel lines that have initial values that are markedly different from the "true" initial values in each case (compare panels C and D). If such change has important physiological implications,the rm-ANOVA model omits it from the fitted mean response. Thus, even though the model correctly detects a difference in treatment groups, the exact nature of this difference is not correctly identified, limiting valuable inferences from the data.  

This section has used simulation to better convey the limitations of linearity and correlation in the response in non-linear data. Although the model fitted to the simulated data was an rm-ANOVA model, the main issue of an expected linear trend in the response is the same in the case of a LMEM. In the following section, we present generalized additive models (GAMs) as a data-driven  alternative method to analyze longitudinal non-linear data.

# GAMs as a special case of Generalized Linear Models

## GAMs and Basis Functions

Generalized linear models (GLMs) are a family of models that fit a linear response function to data that do not have normally distributed errors[@nelder1972]. In contrast, GAMs are a family of regression-based methods for estimating smoothly varying trends and are a broader class of models that contain the GLM family as a special case[@simpson2018;@wood2017;@hastie1987]. A GAM model can be written as:


\begin{equation}
  y_{ijt}=\beta_0+f(x_t\mid \beta_j)+\varepsilon_{ijt}
  (\#eq:GAM)
\end{equation}

Where $y_{ijt}$ is the response at time $t$ of subject  $i$ in group $j$, $\beta_0$ is the expected value at time 0, the change of $y_{ijt}$ over time is represented by the function $f(x_t\mid \beta_j)$ with inputs as the covariates $x_t$ and parameters $\beta_j$, and $\varepsilon_{ijt}$ represents the residual error.

In contrast to the linear functions used to model the relationship between the covariates and the response in rm-ANOVA or LMEM, GAMs use more flexible _smooth functions_. This approach is advantageous as it does not restrict the model to a linear relationship, although a GAM will estimate a linear relationship if the data is consistent with a linear response. One possible set of functions for $f(x_t\mid \beta_j)$ that allow for non-linear responses are polynomials, but a major limitation is that polynomials create a "global" fit as they assume that the same relationship exists everywhere, which can cause problems with the fit [@beck1998]. In particular, polynomial fits are known to show boundary effects because as $t$ goes to $\pm \infty$, $f(x_t \mid \beta_j)$ goes to $\pm \infty$ which is almost always unrealistic, and causes bias at the endpoints of the time period.


The smooth functional relationship between the covariates and the response in GAMs is specified   using a semi-parametric relationship that can be fit within the GLM framework, by using _basis functions_ expansions of the covariates and by estimating random coefficients for these basis functions. A _basis_ is a set of functions that spans the space where the smooths that approximate $f(x_t \mid \beta_j)$ exist [@simpson2018]. For the linear model in Equation \@ref(eq:linear-model), the basis coefficients are $\beta_1$, $\beta_2$ and $\beta_3$ and the basis vectors are $time_t$, $treatment_j$ and $time_t \times treatment_j$. The basis function then, is the combination of basis coefficients and basis vectors that map the possible relationship between the covariates and the response [@hefley2017], which in the case of Equation \@ref(eq:linear-model) is restricted to a linear family of functions.  In the case of Equation \@ref(eq:GAM), the basis function is $f(x_t\mid \beta_j)$, which means that the model allows for non-linear relationships among the covariates.

A commonly used _basis function_ is the cubic spline, which is a smooth curve constructed from cubic polynomials joined together in a manner that enforces smoothness [@wood2017;@simpson2018]. Cubic splines have a long history in solving semi-parametric statistical problems and are often a default choice to fit GAMs as they are a simple, flexible and powerful option to obtain smoothness [@wegman1983]. Therefore, this data-driven flexibility in GAMs overcomes the limitation that occurs in LMEMs and rm-ANOVA when the data is non linear.

To further clarify the concept of basis functions and smooth functions, consider the simulated response for Group 1 in Figure (\@ref(fig:l-q-response), C). The simplest GAM model that can be used to estimate such response is that of a single smooth term for the time effect; i.e., a model that fits a smooth to the trend of the group through time. The timeline can be divided in equally spaced _knots_, each knot being a region where a different basis function will be used. Because there are six timepoints for this group, five knots can be used. The model with five knots to construct the smooth term means that it will have four basis functions (plus one that corresponds to the intercept). The choice of basis functions is already optimized in the package _mgcv_ depending on the number of knots. In Panel A of Figure \@ref(fig:basis-plot), the four basis functions (and the intercept) are shown. Each of the basis functions is composed of six different points (because there are six points on the timeline). To control the "wigliness" of the fit, each of the basis functions of Panel A is penalized by multiplying it by a coefficient according to the penalty matrix of Panel B. The penalty reduces the "wigliness" of the smooth fit to prevent overfitting: A weak penalty estimate will result in wiggly functions whereas a strong penalty estimate provides evidence that a linear response is appropriate.




In other words, the six points of each basis are multiplied by the corresponding coefficient in panel B, thereby increasing or decreasing the original basis functions of Panel A. In Figure \@ref(fig:basis-plot), Panel C shows the resulting penalized basis functions. Note that the penalization for basis 1 has resulted in a decrease of its overall value (because the coefficient for that basis function is negative and less than 1); on the other hand, basis 3 has roughly doubled its value. Finally, the penalized basis functions are added at each timepoint to produce the smooth term. The resulting smooth term for the effect of _time_ is shown in Panel D (orange line) along the simulated values per group, which appear as points.


![(\#fig:basis-plot)Basis functions for a single smoother for time with five knots. A: Basis functions for a single smoother for time for the simulated data of Group 1 from Figure 2, the intercept basis is not shown. B: Penalty matrix for the basis functions. Each basis function is penalized by a coefficient which can be positive or negative. The coefficient determines the overall effect of each basis in the final smoother. C: Penalized basis functions. Each of the four basis functions of panel A has been penalized by the corresponding coefficient shown in Panel B, note the corresponding increase (or decrease) of each basis. D. Smoother for time and original datapoints. The smoother (line) is the result of the sum of each penalized basis function at each time point, simulated values for the group appear as points.](Manuscript_AM_v4_files/figure-docx/basis-plot-1.png){width=75% }

 \newpage
 
## GAMs and covariance

Although the specific methods of how GAMs model correlation structures is a topic beyond the scope of this work, it suffices to say that GAMs can handle correlation structures beyond compound symmetry. A detailed description on basis functions and correlations can be found in [@hefley2017]. In a practical sense, when a GAM is implemented for longitudinal data, a spline can be added to the model for the _time_ effect to account for the repeated measures over time.  An example where this is covered is in the Appendix.



## Determination of significance in GAMs

At the core of a biomedical longitudinal study lies the question of a significant difference between the effect of two or more treatments in different groups. Whereas in rm-ANOVA a _post-hoc_ analysis is required to answer such question by calculating some _p-values_ after multiple comparisons, GAMs use a different approach to estimate significance. In essence, the idea behind the estimation of significance in GAMs across different treatment groups is that if the _difference_ between the confidence intervals of the fitted smooths for such groups is  non-zero, then a significant difference exists at that timepoint (or timepoints). The absence of a _p-value_ in this case might seem odd, but when a confidence interval comparison is put into context the validity of its rationale becomes apparent. The major advantage of the significance estimation in GAMs  is that they allow to determine significance at a specific timepoint, rather than providing a single _p-value_ for the overall treatment effect.  The estimation of siginificant differences in GAMs is covered in detail in the Appendix using simulated longitudinal data that follows previously reported trends of tumor response between multiple treatment groups [@vishwanath2009].



***

# References

<div id="refs"></div>


\newpage

# (APPENDIX) Appendix {-} 

# Simulation

## Compound symmetry and independent errors in linear and quadratic responses

This section simulated linear and quadratic data in the same manner as in Section \@ref(simulation). The linear simulations using Figure \@ref(fig:linear-cases) show in  panels A and D the simulated mean responses and individual datapoints.  Panels C and G show a visual interpretation of "correlation" in the responses: In panel C, subjects that have a value of the random error $\varepsilon$ either above or below the mean group response are more likely to have other observations that follow the same trajectory, thereby demonstrating correlation in the response. In panel G,because the errors are independent, there is no expectation that responses are likely to follow a similar pattern.  Panels D and H show the predictions from the rm-ANOVA model.


```r
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
  library(ggsci)
  
  x <- seq(1,6, length.out = n_time)
  mu <- matrix(0, length(x), 2)
  # linear response
  if (fun_type == "linear") {
    mu[, 1] <- - (0.25*x)+2  
    mu[, 2] <- 0.25*x+2
  } else {
    # nonlinear response
    
    mu[, 1] <-  -(0.25 * x^2) +1.5*x-1.25
    mu[, 2] <- (0.25 * x^2) -1.5*x+1.25
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

plot_example <- function(sim_dat) {
  library(patchwork)
  ## Plot the simulated data
  p1 <- sim_dat$dat %>%
    ggplot(aes(x = time, y = y, group = treatment, color = treatment)) +
    geom_point(show.legend=FALSE) +labs(y='response')+
    geom_line(aes(x = time, y = mu, color = treatment),show.legend=FALSE) +
    theme_classic() +
    #ggtitle("Simulated data with true response function")+
    theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    scale_color_aaas()
  
  p2 <- sim_dat$dat %>%
    ggplot(aes(x = time, y = y, group = subject, color = treatment)) +
    geom_line(aes(size = "Subjects"),show.legend = FALSE) +
    # facet_wrap(~ treatment) +
    geom_line(aes(x = time, y = mu, color = treatment, size = "Simulated Truth"), lty = 1,show.legend = FALSE) +labs(y='response')+
    scale_size_manual(name = "Type", values=c("Subjects" = 0.5, "Simulated Truth" = 3)) +
    #ggtitle("Simulated data\nIndividual responses with population mean") +
    theme_classic()+
     theme(plot.title = element_text(size = 30, 
                                face = "bold"),
     text=element_text(size=30))+
    scale_color_aaas()
  
   p3 <- sim_dat$dat %>%
    ggplot(aes(x = time, y = errors, group = subject, color = treatment)) +
    geom_line(show.legend=FALSE) +labs(y='errors')+
     theme_classic()+
    # facet_wrap(~ treatment) +
    #ggtitle("Simulated errors") +
     theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    scale_color_aaas()
  
  p4 <- ggplot(sim_dat$dat, aes(x = time, y = y, color = treatment)) +
    geom_point()+labs(y='response')+
    geom_line(aes(y = predict(sim_dat$fit_lme), group = subject, size = "Subjects")) +
    geom_line(data = sim_dat$pred_dat, aes(y = predict(sim_dat$fit_lme, level = 0, newdata = sim_dat$pred_dat), size = "Population")) +
    scale_size_manual(name = "Predictions", values=c("Subjects" = 0.5, "Population" = 3)) +
    theme_classic() +
    #ggtitle("Fitted Model")+
    theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    scale_color_aaas()
  
  return((p1+p3+p2+p4)+plot_layout(nrow=1)+plot_annotation(tag_levels = 'A')) 
  
    
}

txt<-18
A1<-plot_example(example(fun_type = "linear", error_type = "correlated")) 

B1<-plot_example(example(fun_type = "linear", error_type = "independent")) 
  
C1<-plot_example(example(fun_type = "quadratic", error_type = "correlated")) 
  
D1<-plot_example(example(fun_type = "quadratic", error_type = "independent")) 
```



![(\#fig:linear-cases)**Simulated linear responses from two groups with correlated (top row) or independent (bottom row) errors using a rm-ANOVA model. A, C:Simulated data with known mean response (linear or quadratic, thin lines) and individual responses (points) showing the dispersion of the data. B,D: Estimations from the rm-ANOVA model for the mean group response (linear of quadratic). Thick lines are the predicted mean response per group, thin lines are the random effects for each subject and points represent the original raw data. The rm-ANOVA model does not pick the trend of the quadratic data**](Manuscript_AM_v4_files/figure-docx/linear-cases-1.png)

For the quadratic response case, Figure \@ref(fig:quadratic-cases) shows the simulated responses using compound symmetry and independent errors. 

![(\#fig:quadratic-cases)**Simulated quadratic responses from two groups with a rm-ANOVA model fitted. A,E:Simulated data with known mean response (lines) and individual responses (points) showing the dispersion of the data. B,F: Generated errors showing the difference in the behavior of correlated and independent errors.  C,G: Simulated known response per group (thick lines) with individual trajectories (thin lines), note that subjects with observations in the area above the mean response tend to stay in that region through the timeline. D,H: Estimations from the rm-ANOVA model for the mean group response. Thick lines are the predicted mean response per group, thin lines are the random effects for each subject and points represent the original raw data.**](Manuscript_AM_v4_files/figure-docx/quadratic-cases-1.png)

