---
title: '**The statistical analysis of non-linear longitudinal data in biomedical research using generalized additive models**'
subtitle: _Beyond repeated measures ANOVA and Linear Mixed Models_
header-includes:
    \usepackage{placeins}
output:  
  bookdown::word_document2:
    fig_caption: yes #figure caption
    keep_md: yes
  bookdown::html_document2:
    css: "style.css" #style for the HTML document
  bookdown::pdf_document2:
    #template: my-template.tex #if a custom template that removes the additional "and" in the author information is desired
    pandoc_args: --listings #calls the listings package to fit code within the page margins
    keep_tex: yes #keep LaTeX file for submission
    fig_caption: yes #allows captions in figures
    extra_dependencies:
      subfig: null #allows for subfigures
      breqn: null #line breaks for long equations
      caption: ["font={small}"] #size of the figure captions
      float: null #allows for control of placement of figures
    includes:
      in_header: preamble.sty #additional LaTeX formatting
csl: elsevier-with-titles.csl #style for references
bibliography: refs.bib #references
link-citations: yes #adds links to the citations
'': default
---


# Abstract

In biomedical research, the outcome of longitudinal studies has been traditionally analyzed using the _repeated measures analysis of variance_ (rm-ANOVA) or more recently, _linear mixed models_ (LMEMs). Although LMEMs are less restrictive than rm-ANOVA in terms of correlation and missing observations, both methodologies share an assumption of linearity in the measured response, which results in biased estimates and unreliable inference when they are used to analyze data where the trends are non-linear. In contrast, generalized additive models (GAMs) relax the linearity assumption, and allow the data to determine the fit of the model while permitting missing observations and different correlation structures. Therefore, GAMs present an excellent choice to analyze non-linear longitudinal data in the context of biomedical research. This paper summarizes the limitations of rm-ANOVA and LMEMs, presents the basic theory of GAMs, and uses simulated data that follows trends reported in the biomedical literature to demonstrate their implementation in $\textsf{R}$ via the package _mgcv_. To make this work reproducible, the code and data used in this paper are available at:________.






# Background

Longitudinal studies are designed to repeatedly measure a variable of interest in a group (or groups) of subjects, with the intention of observing the evolution of effect across time rather than analyzing a single time point (e.g., a cross-sectional study). Biomedical research frequently uses longitudinal studies to analyze the evolution of a "treatment" effect across multiple time points; and in such studies the subjects of analysis range from animals (mice, rats, rabbits), to human patients, cells, or blood samples, among many others. Tumor response [@roblyer2011;@tank2020;@pavlov2018;@demidov2018], antibody expression [@ritter2001;@roth2017], and cell metabolism [@jones2018;@skala2010] are examples of the different situations where researchers have used longitudinal designs to study some  physiological response. Because the frequency of the measurements in a longitudinal study is dependent on the biological phenomena of interest and the experimental design of the study, the frequency of such measurements can range from minute intervals to study a short-term response such as anesthesia effects in animals[@greening2018], to weekly measurements to analyze a mid-term response like the evolution of dermatitis symptoms in breast cancer patients [@sio2016], to monthly measurements to study a long-term response such as mouth opening following radiotherapy (RT) in neck cancer patients [@kamstra2015]. 


Traditionally, a “frequentist” or "classical" statistical paradigm is used in biomedical research to derive inferences from a longitudinal study. The frequentist paradigm regards probability as the limit of the expected outcome when an experiment is repeated a large number of times [@wagenmakers2008], and such view is applied to the analysis of longitudinal data by assuming a null hypothesis under a statistical model that is often an _analysis of variance over repeated measures_ (repeated measures ANOVA or rm-ANOVA). The rm-ANOVA model makes three key assumptions regarding longitudinal data: 1) linearity of the response across time, 2) constant correlation across same-subject measurements, and 3) observations from each subject are obtained at all time points through the study (a condition also known as _complete observations_) [@gueorguieva2004;@schober2018]. 

The expected linear behavior of the response through time is a key requisite in rm-ANOVA [@pinheiro2006]. This "linearity assumption" in rm-ANOVA implies that the model is misspecified when the data does not follow a linear trend, which results in unreliable inference. In biomedical research, non-linear trends are the norm rather than the exception in longitudinal studies. A particular example of this non-linear behavior in longitudinal data arises in measurements of tumor response to chemo and/or radiotherapy in preclinical and clinical settings [@roblyer2011;@skala2010;@vishwanath2009]. These studies have shown that the collected signal does not follow a linear trend over time, and presents extreme variability at different time points, making the fit of rm-ANOVA model inconsistent with the observed variation. Therefore, when rm-ANOVA is used to draw inference of such data the estimates are inevitably biased, because the model is only able to accommodate linear trends that fail to adequately represent the biological phenomenon of interest.


A _post hoc_ analysis is often used in conjunction with rm-ANOVA to perform repeated comparisons to estimate a _p-value_, which in turn is used as a measure of significance.
Although it is possible that a _post hoc_ analysis of rm-ANOVA is able to find “significant” _p-values_( _p_<0.05) from non-linear data, the validity of such metric is dependent on how adequate the model fits the data. In other words, _p-values_ are valid only if the model and the data have good agreement; if that is not the case, a "Type III" error (known as "model misspecification") occurs[@dennis2019]. For example, model misspecification will occur when a model that is only able to explain linear responses (such as rm-ANOVA) is fitted to data that follows a quadratic trend, thereby causing the resulting _p-values_ and parameter estimates to be invalid [@wang2019].

Additionally, the _p-value_ itself is highly variable, and multiple comparisons can inflate the false positivity rate (Type I error or $\alpha$) [@liu2010;@halsey2015], consequently biasing the conclusions of the study. Corrections exist to address the Type I error issue of multiple comparisons (such as Bonferroni [@abdi2010]), but they in turn reduce statistical power (1-$\beta$)[@nakagawa2004], and lead to increased Type II error (failing to reject the null hypothesis when the null hypothesis is false) [@gelman2012;@albers2019]. Therefore, the tradeoff of _post hoc_ comparisons in rm-ANOVA between Type I, II and III errors might be difficult to resolve in a biomedical longitudinal study where a delicate balance exists between statistical power and sample size.


On the other hand, the assumption of constant correlation in rm-ANOVA (often known as the _compound symmetry assumption_) is typically unreasonable because correlation between the measured responses often diminishes as the time interval between the observation increases [@ugrinowitsch2004].  Corrections can be made in rm-ANOVA in the absence of compound symmetry [@huynh1976;@greenhouse1959], but the effectiveness of the correction is limited by the size of the sample, the number of measurements[@haverkamp2017], and group sizes [@keselman2001]. In the case of biomedical research, where living subjects are frequently used, sample sizes are often not "large" due to ethical and budgetary reasons [@charan2013] which might cause the corrections for lack of compound symmetry to be ineffective.

Due to a variety of causes, the number of observations during a study can vary between all subjects. For example, in a clinical trial patients may voluntarily withdraw, whereas attrition  due to injury or weight loss in preclinical animal studies is possible. It is even plausible that unexpected complications with equipment or supplies arise that prevent the researcher from collecting measurements at certain time points. In each of these missing data scenarios, the _complete observations_ assumption of  classical rm-ANOVA is violated. When incomplete observations occur, a rm-ANOVA model is fit by excluding all subjects with missing observations from the analysis [@gueorguieva2004]. This elimination of partially missing data from the analysis can result in increased costs if the desired statistical power is not met with the remaining observations, because it would be necessary to enroll more subjects. At the same time, if the excluded observations contain insightful information that is not used, their elimination from the analysis may limit the demonstration of significant differences between groups. 

During the last decade, the biomedical community has started to recognize the limitations of rm-ANOVA in the analysis of longitudinal data. The recognition on the shortcomings of rm-ANOVA is exemplified by the  use of  linear mixed effects models (LMEMs) by certain groups to analyze longitudinal tumor response data [@skala2010;@vishwanath2009]. Briefly, LMEMs incorporate _fixed effects_, which correspond to the levels of experimental factors in the study (e.g., the different drug regimens in a clinical trial), and _random effects_, which account for random variation within the population (e.g., the individual-level differences not due to treatment such as weight or age). When compared to the traditional rm-ANOVA, LMEMs are more flexible as they can accommodate missing observations for multiple subjects and allow different modeling strategies for the variability within each measure in every subject [@pinheiro2006]. However, LMEMs impose restrictions in the distribution of the errors  of the random effects, which need to be normally distributed and independent [@gueorguieva2004;@barr2013]. And even more importantly, LMEMs also assume a linear relationship between the response and time [@pinheiro2006], making them unsuitable to analyze non-linear data.

As the rm-ANOVA and the more flexible LMEM approaches make overly restrictive assumptions regarding the linearity of the response, there is a need for biomedical researchers to explore the use of additional statistical tools that allow the data (and not an assumption in trend) to determine the trend of the fitted model, to enable appropriate inference.

In this regard, generalized additive models (GAMs) present an alternative approach to analyze longitudinal data. Although not frequently used by the biomedical community, these semi-parametric models are customarily used in other fields to analyze longitudinal data.   Examples of the use of GAMs include the analysis of temporal variations in geochemical and palaeoecological data  [@rose2012;@pedersen2019;@simpson2018], health-environment interactions [@yang2012] and the dynamics of government in political science [@beck1998] . There are several advantages of GAMs over LMEMs and rm-ANOVA models: 1) GAMs can fit a more flexible class of smooth responses that enable the data to dictate the trend in the fit of the model, 2) they can model non-constant correlation between repeated measurements [@wood2017] and 3) can easily accommodate missing observations. Therefore, GAMs can provide a more flexible statistical approach to analyze non-linear biomedical longitudinal data than LMEMs and rm-ANOVA.

The current advances in programming languages designed for statistical analysis (specifically $\textsf{R}$), have eased the computational implementation of traditional models such as rm-ANOVA and more complex approaches such as LMEMs and GAMs. In particular, $\textsf{R}$[@r] has an extensive collection of documentation and functions to fit GAMs in the package _mgcv_ [@wood2016;@wood2017] that not only speed up the initial stages of the analysis but also enable the use of advanced modeling structures (e.g. hierarchical models, confidence interval comparisons) without requiring advanced programming skills from the user. At the same time, $\textsf{R}$ has many tools that simplify data simulation, an emerging strategy used to test statistical models [@haverkamp2017]. Data simulation methods allow the researcher to create and  explore different alternatives for analysis without collecting information in the field, reducing the time window between experiment design and its implementation, and simulation can be also used for power calculations and study design questions.   

This work provides biomedical researchers with a clear understanding of the theory and the practice of using GAMs to analyze longitudinal data using by focusing on four areas. First, the limitations of LMEMs and rm-ANOVA regarding linearity of response, constant correlation structures and missing observations is explained in detail. Second, the key theoretical elements of GAMs are presented using clear and simple mathematical notation while explaining the context and interpretation of the equations. Third, using simulated data that reproduces patterns in previously reported studies [@vishwanath2009] we illustrate the type of non-linear longitudinal data that often occurs in biomedical research. The simulated data experiments highlight the differences in inference between rm-ANOVA, LMEMs and GAMs on data similar to what is commonly observed in biomedical studies. Finally, reproducibility is emphasized by providing the code to generate the simulated data and the implementation of different models in $\textsf{R}$, in conjunction with a step-by-step guide demonstrating how to fit models of increasing complexity.  

In summary, this work will allow biomedical researchers to identify when the use of GAMs instead of rm-ANOVA or LMEMs is appropriate to analyze longitudinal data, and provide guidance on the implementation of these models by improving the standards for reproducibility  in biomedical research.

\FloatBarrier



# Challenges presented by longitudinal studies

## The repeated measures ANOVA

The _repeated measures analysis of variance_ (rm-ANOVA) is the standard statistical analysis for longitudinal data in biomedical research. This statistical methodology requires certain assumptions for the model to be valid. From a practical view, the assumptions can be divided in three areas: 1) linear relationship between covariates and response, 2) a constant correlation between measurements, and, 3) complete observations for all subjects. Each one of these assumptions is discussed below.
  
## Linear relationship 

### The repeated measures ANOVA case

In a longitudinal biomedical study, two or more groups of subjects (e.g., human subject, mice, samples) are subject to different treatments (e.g., a "treatment" group receives a novel drug or intervention vs. a "control" group that receives a placebo), and measurements from each subject within each group are collected at specific time points. The collected response is modeled with _fixed_ components. The _fixed_ component can be understood as a constant value in the response which the researcher is interested in measuring, i.e., the average effect of the novel drug/intervention in the "treatment" group. 


Mathematically speaking, a rm-ANOVA model with an interaction can be written as:

$$\begin{equation}
y_{ijt} = \beta_0+\beta_1 \times time_{t} +\beta_2 \times treatment_{j} +\beta_3 \times time_{t}\times treatment_{j}+\varepsilon_{ijt}\\ 
\qquad(3.1)
\end{equation}$$

In this model $y_{ijt}$ is the response for subject $i$, in treatment group $j$ at time $t$, which can be decomposed in a mean value $\beta_0$, _fixed effects_ of time ($time_t$), treatment ($treatment_j$) and their interaction $time_t*treatment_j$ which have linear slopes given by $\beta_1, \beta_2$ and $\beta_3$, respectively. Independent errors $\varepsilon_{tij}$ represent random variation not explained by the _fixed_ effects, and are assumed to be $\sim N(0,\sigma^2)$ (independently and identically normally distributed with mean zero and variance $\sigma^2$).
In a  biomedical research context, suppose two treatments groups are used in a study (e.g., "placebo" vs. "novel drug" or "saline" vs. "chemotherapy"). Then, the group terms in Equation (3.1) can be written as below with $treatment_j=0$ representing the first treatment group (Group A) and $treatment_j=1$ representing the second treatment group (Group B). The linear models then can be expressed as


$$\begin{equation}
y_{ijt} = \begin{cases}
\beta_0 + \beta_1\times time_{t}+\varepsilon_{ijt}   & \mbox{if Group A}\\
\beta_0 + \beta_2+\beta_1 \times time_{t} +\beta_3 \times time_{t}+\varepsilon_{ijt}  & \mbox{if Group B}\\
\end{cases}
\qquad(3.2)
\end{equation}$$

To further simplify the expression, substitute $\widetilde{\beta_{0}}=\beta_0+\beta_{2}$ and $\widetilde{\beta_{1}}=\beta_{1}+\beta_{3}$ in the equation for Group B. This substitution allows for a different intercept and slope for Groups A and B. The model is then written as

$$\begin{equation}
y_{ijt} = \begin{cases}
\beta_0 + \beta_1\times time_{t}+\varepsilon_{ijt}   & \mbox{if Group A}\\
\widetilde{\beta_{0}} + \widetilde{\beta_1} \times time_{t}+\varepsilon_{ijt}  & \mbox{if Group B}\\
\end{cases}
\qquad(3.3)
\end{equation}$$

Presenting the model in this manner makes clear that when treating different groups, an rm-ANOVA model is able to accommodate non-parallel lines in each case (different intercepts and slopes per group). In other words, the rm-ANOVA model "expects" a linear relationship between the covariates and the response, this means that either presented as Equation (3.1), Equation (3.2) or Equation (3.3), an rm-ANOVA model is only able to accommodate linear patterns in the data. If the data show non-linear behavior, the rm-ANOVA model will approximate this behavior with non-parallel lines. 

### The Linear Mixed Model Case

A linear mixed model (LMEM) is a class of statistical model that incorporates _fixed effects_ to model the relationship between the covariates and the response, and _random effects_ to model subject variability that is not the primary focus of the study but that might be important to distinguish [@pinheiro2006;@west2014]. A LMEM with interaction between time and treatment for a longitudinal study  can be written as:


$$\begin{equation}
y_{ijt} = \beta_0+ \beta_1 \times time_{t} + \beta_2 \times treatment_{j} + \beta_3 \times time_{t}\times treatment_{j}+\mu_{ij} +\varepsilon_{ijt}\\ 
\qquad(3.4)
\end{equation}$$

When Equation (3.1) and Equation (3.4) are compared, it is easily noticeable that LMEM and rm-ANOVA have the same construction regarding the _fixed effects_ of time and treatment, but that the LMEM incorporates an additional source of variation (the term $\mu_{ij}$). This term $\mu_{ij}$ is the one that corresponds to the  _random effect_, accounting for variability in each subject within each group. The _random_ component can also be understood as used to model some "noise" in the response, but that is intended to be analyzed and disentangled from the "global noise" term $\varepsilon_{ijt}$ from Equation (3.1). 

For example, if the blood concentration of the drug is measured in certain subjects in the early hours of the morning while other subjects are measured in the afternoon, it is possible that the difference in the collection time introduces some "noise" in the data. As the name suggests, this "random" variability needs to be modeled as a variable rather than as a constant value.  The _random effect_ $\mu_{ij}$ in Equation (3.4) is assumed to be $\mu_{ij} \sim N(0,\sigma^2_\mu)$. In essence,the _random effect_ in a LMEM enables to fit models with different slopes at the subject-level[@pinheiro2006]. However,the expected linear relationship of the covariates and the response in Equation (3.1) and in Equation (3.4) is essentially the same, representing a major limitation of LMEMs to fit a non-linear response.

## Covariance in rm-ANOVA and LMEMs

In a longitudinal study there is an expected _covariance_ between repeated measurements on the same subject, and because repeated measures occur in the subjects within each group, there is a _covariance_ between  measurements at each time point within each group. The _covariance matrix_ (also known as the variance-covariance matrix) is a matrix that captures the variation between and within subjects in a longitudinal study[@wolfinger1996] (For an in-depth analysis of the covariance matrix see [@west2014;@weiss2005]). 

In the case of an rm-ANOVA analysis, it is typically assumed that the covariance matrix has a specific construction known as _compound symmetry_ (also known as "sphericity" or "circularity"). Under this assumption, the between-subject variance and within-subject correlation  are constant across time [@weiss2005;@geisser1958;@huynh1976]. However, it has been shown that this condition is frequently not justified because the correlation between measurements tends to change over time [@maxwell2017]; and it is higher between consecutive measurements [@gueorguieva2004;@ugrinowitsch2004]. Although corrections can be made (such as Huyhn-Feldt or Greenhouse-Geisser)[@huynh1976;@greenhouse1959] the effectiveness of each correction is limited because it depends on the size of the sample,the number of repeated measurements[@haverkamp2017], and they are not robust if the group sizes are unbalanced [@keselman2001]. Because biomedical longitudinal studies are often limited in sample size and can have an imbalanced design, the corrections required to use an rm-ANOVA model may not be able to provide a reasonable adjustment that makes the model valid.


In the case of LMEMs, one key advantage over rm-ANOVA is that they allow different structures for the variance-covariance matrix including exponential, autoregressive of order 1, rational quadratic and others [@pinheiro2006]. Nevertheless, the analysis required to determine an appropriate variance-covariance structure for the data can be a challenging process by itself. Overall, the spherical assumption for rm-ANOVA may not capture the natural variations of the correlation in the data, and can bias the inferences from the analysis. 


## Missing observations

Missing observations are an issue that arises frequently in longitudinal studies. In biomedical research, this situation can be caused by reasons beyond the control of the investigator [@molenberghs2004]. Dropout from patients and attrition or injury in animals are among the reasons for missing observations. Statistically, missing information can be classified as _missing at random_ (MAR), _missing completely at random_ (MCAR), and _missing not at random_ (MNAR) [@weiss2005].  In a MAR scenario, the pattern of the missing information is related to some variable in the data, but it is not related to the variable of interest [@scheffer2002]. If the data are MCAR, this means that the missingness is completely unrelated to the collected information [@potthoff2006], and in the case of MNAR the missing values are dependent on their value. 

An rm-ANOVA model assumes complete observations for all subjects, and therefore subjects with one or more missing observations are excluded from the analysis. This is inconvenient because the remaining subjects might not accurately represent the population, and statistical power is affected by this reduction in sample size [@ma2012]. In the case of LMEMs, inferences from the model are valid when missing observations in the data exist that are MAR or MCAR [@west2014]. For example, if attrition occurs in all mice that had lower weights at the beginning of a chemotherapy response study, the missing data can be considered MAR because  the missigness is unrelated to other variables of interest.

## What does an rm-ANOVA fit looks like? A visual representation using simulated data{#simulation}

To visually demonstrate the limitations of rm-ANOVA an LMEMs for non-linear longitudinal data, this section presents a simulation experiment of a normally distributed response of two groups of 10 subjects each. An rm-ANOVA model (Equation (3.1)), and a LMEM  (Equation (3.4)) are fitted to each group, using $\textsf{R}$[@r] and the package _nlme_[@nlme]. 

Briefly, two cases for the mean responses for each group are considered: in the first case, the mean response in each group is a linear function with different intercepts and slopes; a negative slope is used for Group 1 and a positive slope is used for Group 2 (Figure <a href="#fig:l-q-response">3.1</a>, A). In the second case, a second-degree polynomial (quadratic) function is used for the mean response per group: the quadratic function is concave down for Group 1 and it is concave up for Group 2 (Figure <a href="#fig:l-q-response">3.1</a>, C). In both the linear and quadratic simulated data, the groups start with the same mean value at the first time point. This is intentional in order to simulate the expected temporal evolution of some physiological quantity, which is typical in biomedical experiments.

Specifically, the rationale for the chosen linear and quadratic functions is the likelihood that a measured response in two treatment groups is similar in the initial phase of the study, but as therapy progresses a divergence in the trend of the response indicates a treatment effect. In other words, Group 1 can be thought as a "Control" group and Group 2 as a "Treatment" group. From the mean response per group (linear or quadratic), the variability or "error" of individual responses  within each group is simulated using a covariance matrix with compound symmetry (constant variance across time). Thus, the response per subject in both the linear and quadratic simulation corresponds to the mean response per group plus the error (Figure <a href="#fig:l-q-response">3.1</a> B,D). 

A more comprehensive exploration  of the fit of rm-ANOVA and LMEMs for linear and non-linear longitudinal appears in Figure <a href="#fig:linear-cases-Appendix">10.1</a> and Figure <a href="#fig:quadratic-cases-Appendix">10.2</a> in the Appendix, where simulation with compound symmetry and independent errors (errors generated from a normal distribution that are not constant over time) and the plot of simulated errors, and fitted parameters in presented. We are aware that the simulated data used in this section present an extreme case that might not occur frequently in biomedical research, but they are used as a representation of the consequences of modeling non-linear data with a linear model such as rm-ANOVA or LMEMs. Of notice, Section <a href="#longitudinal-GAMs">5</a> uses simulated data that does follow reported trends in the biomedical literature.







![Figure 3.1: Simulated linear responses from two groups with correlated (top row) or independent (bottom row) errors using a LMEM and a rm-ANOVA model. A, D: Simulated data with known mean response (linear or quadratic, thin lines) and individual responses (points) showing the dispersion of the data. B,E: Estimates from the rm-ANOVA model for the mean group response (linear of quadratic). Points represent the original raw data. The rm-ANOVA model not only fails to pick the trend of the quadratic data but also assigns a global estimate that does not take between-subject variation. C, F: Estimates from the LMEM model in the linear and quadratic case.The LMEM incorporates a random effect for each subject, but this model and the rm-ANOVA model are unable to follow the trend of the data in each group and grossly bias the initial estimates for each group.](Full_document_files/figure-docx/l-q-response-1.png){width=75% }

The simulation shows that the fit produced by the LMEM and the rm-ANOVA model is good for linear data, as the predictions for the mean response are reasonably close to the "truth" of the simulated data  (Figure <a href="#fig:l-q-response">3.1</a>,B, E). When the linearity and compound symmetry assumptions are met, the model approximates well the individual trends and the mean trends by group.

However, consider the case when the data follows a non-linear trend, such as the simulated data in Figure <a href="#fig:l-q-response">3.1</a>, C. Here, the mean response per group was simulated using a quadratic function but errors, individual responses and the rm-ANOVA model were produced in the same manner as in (Figure <a href="#fig:l-q-response">3.1</a> A and B) . The mean response in the simulated data with quadratic behavior is changing in each group through the timeline, and the mean value is the same as the initial value by the fifth time point for each group. Fitting an rm-ANOVA model (3.1) or a LMEM (3.4) to this data  produces the fit that appears in panels E and F in Figure <a href="#fig:l-q-response">3.1</a>.

A comparison of the fitted mean response of the LMEM and the rm-ANOVA model to the simulated data in Figure ((<a href="#fig:l-q-response">3.1</a>, E, F) indicates that the models are not capturing the changes within each group.  Specifically, note that the fitted mean response of both models (panel E, F) show that the change (increase for Treatment 1 or decrease for Treatment 2) in the response through time points 2 and 4 is not being captured. The LMEM is only able to account for between-subject variation by providing different intercepts to each subject, but both models are not able to capture the fact that the initial values are the same in each group, and instead fit non-parallel lines that have initial values that are markedly different from the "true" initial values in each case (compare panel D with panels E and F). If such a change has important physiological implications,both rm-ANOVA and LMEMs omit it from the fitted mean response. Thus, even though the model correctly detects a divergence between treatment groups, the exact nature of this difference is not correctly identified, limiting valuable inferences from the data.  

This section has used simulation to better convey the limitations of linearity and correlation in the response in non-linear data. Although the model fitted to the simulated data was an rm-ANOVA model, the main issue of an expected linear trend in the response is the same in the case of a LMEM. In the following section, we present generalized additive models (GAMs) as a data-driven  alternative method to analyze longitudinal non-linear data.

\FloatBarrier


# GAMs as a special case of Generalized Linear Models{#GAM-theory}

## GAMs and Basis Functions

Generalized linear models (GLMs) are a family of models that fit a linear response function to data that do not have normally distributed errors[@nelder1972]. In contrast, GAMs are a family of regression-based methods for estimating smoothly varying trends and are a broader class of models that contain the GLM family as a special case[@simpson2018;@wood2017;@hastie1987]. A GAM model can be written as:


$$\begin{equation}
  y_{ijt}=\beta_0+f(x_t\mid \beta_j)+\varepsilon_{ijt}
  \qquad(4.1)
\end{equation}$$

Where $y_{ijt}$ is the response at time $t$ of subject  $i$ in group $j$, $\beta_0$ is the expected value at time 0, the change of $y_{ijt}$ over time is represented by the _smooth function_ $f(x_t\mid \beta_j)$ with inputs as the covariates $x_t$ and parameters $\beta_j$, and $\varepsilon_{ijt}$ represents the residual error.

In contrast to the linear functions used to model the relationship between the covariates and the response in rm-ANOVA or LMEM, GAMs use more flexible _smooth functions_. This approach is advantageous as it does not restrict the model to a linear relationship, although a GAM will estimate a linear relationship if the data is consistent with a linear response. One possible set of functions for $f(x_t\mid \beta_j)$ that allow for non-linear responses are polynomials, but a major limitation is that polynomials create a "global" fit as they assume that the same relationship exists everywhere, which can cause problems with inference [@beck1998]. In particular, polynomial fits are known to show boundary effects because as $t$ goes to $\pm \infty$, $f(x_t \mid \beta_j)$ goes to $\pm \infty$ which is almost always unrealistic, and causes bias at the endpoints of the time period.


The smooth functional relationship between the covariates and the response in GAMs is specified   using a semi-parametric relationship that can be fit within the GLM framework, by using _basis function_ expansions of the covariates and by estimating random coefficients associated with these basis functions. A _basis_ is a set of functions that spans the mathematical space where the smooths that approximate $f(x_t \mid \beta_j)$ exist [@simpson2018]. For the linear model in Equation (3.1), the basis coefficients are $\beta_1$, $\beta_2$ and $\beta_3$ and the basis vectors are $time_t$, $treatment_j$ and $time_t \times treatment_j$. The basis function then, is the combination of basis coefficients and basis vectors that map the possible relationship between the covariates and the response [@hefley2017], which in the case of Equation (3.1) is restricted to a linear family of functions.  In the case of Equation (4.1), the basis function is $f(x_t\mid \beta_j)$, which means that the model allows for non-linear relationships among the covariates.

Commonly used _basis functions_ are splines (cubic, thin plate regression among others). A cubic spline is a smooth curve constructed from cubic polynomials joined together in a manner that enforces smoothness, and thin plate regression splines are an optimized version that work well with noisy data [@wood2017;@simpson2018]. Splines have a long history in solving semi-parametric statistical problems and are often a default choice to fit GAMs as they are a simple, flexible and powerful option to obtain smoothness [@wegman1983]. Therefore, this data-driven flexibility in GAMs overcomes the limitation that occurs in LMEMs and rm-ANOVA when the data is non linear.

To further clarify the concept of basis functions and smooth functions, consider the simulated response for Group 1 in Figure (<a href="#fig:l-q-response">3.1</a>, C). The simplest GAM model that can be used to estimate such response is that of a single smooth term for the time effect; i.e., a model that fits a smooth to the trend of the group through time. The timeline can be divided in equally spaced _knots_, each knot being a region where a different basis function will be used. Because there are six timepoints for this group, five knots can be used. The model with five knots to construct the smooth term means that it will have four basis functions (plus one that corresponds to the intercept). The choice of basis functions is already optimized in the package _mgcv_ depending on the number of knots. In Panel A of Figure <a href="#fig:basis-plot">4.1</a>, the four basis functions (and the intercept) are shown. Each of the basis functions is composed of six different points (because there are six points on the timeline). To control the "wigliness" of the fit, each of the basis functions of Panel A is weighted by multiplying it by a coefficient according to the matrix of Panel B. The parameter estimates are penalized where the penalty reduces the "wigliness" of the smooth fit to prevent overfitting: A weak penalty estimate will result in wiggly functions whereas a strong penalty estimate provides evidence that a linear response is appropriate.




In other words, the six points of each basis are multiplied by the corresponding coefficient in panel B, thereby increasing or decreasing the original basis functions of Panel A. In Figure <a href="#fig:basis-plot">4.1</a>, Panel C shows the resulting penalized basis functions. Note that the penalization for basis 1 has resulted in a decrease of its overall value (because the coefficient for that basis function is negative and less than 1); on the other hand, basis 3 has roughly doubled its value. Finally, the penalized basis functions are added at each timepoint to produce the smooth term. The resulting smooth term for the effect of _time_ is shown in Panel D (orange line) along the simulated values per group, which appear as points.



![Figure 4.1: Basis functions for a single smoother for time with five knots. A: Basis functions for a single smoother for time for the simulated data of Group 1 from Figure 2, the intercept basis is not shown. B: Matrix for basis function penalization. Each basis function is multiplied by a coefficient which can be positive or negative. The coefficient determines the overall effect of each basis in the final smoother. C: Weighted basis functions. Each of the four basis functions of panel A has been weighted by the corresponding coefficient shown in Panel B, note the corresponding increase (or decrease) of each basis. D: Smoother for time and original data points. The smoother (line) is the result of the sum of each penalized basis function at each time point, with simulated values for the group shown as points.](Full_document_files/figure-docx/basis-plot-1.png){width=75% }


\FloatBarrier


\newpage
 
# The analyisis of longitudinal biomedical data using GAMs{#longitudinal-GAMs}

The previous sections provided the basic framework to understand the GAM framework and how these models are more advantageous to analyze non-linear longitudinal data when compared to rm-ANOVA or LMEMs. This section will use simulation to present the practical implementation of GAMs for longitudinal biomedical data using $\textsf{R}$ and the package `mgcv`. The code for the simulated data and figures, and a brief guide for model selection and diagnostics appear in the Appendix.

## Simulated data

The simulated data is based on the reported longitudinal changes in oxygen saturation ($\mbox{StO}_2$) in subcutaneous tumors that appear in Figure 3, C in [@vishwanath2009]. In the paper, diffuse reflectance spectroscopy was used to quantify $\mbox{StO}_2$ changes in both groups at the same time points (days 0, 2, 5, 7 and 10). In the "Treatment" group (chemotherapy) an increase in $\mbox{StO}_2$ is observed through time, while a decrease is seen in the "Control" (saline) group. Following the reported trend, we simulated 10 normally distributed observations at each time point with a standard deviation (SD) of 10% (matching the SD in the original paper). 
The simulated and real data appear in Figure <a href="#fig:sim-smooth-plot">5.1</a>, A and the inlet, respectively.

## An interaction GAM for longitudinal data

An interaction effect is typically the main interest in longitudinal biomedical data, as it takes into account treatment, time, and their combination. In a practical sense, when a GAM is implemented for longitudinal data, a smooth can be added to the model for the _time_ effect to account for the repeated measures over time.  Although specific methods of how GAMs model correlation structures is a topic beyond the scope of this paper, it suffices to say that GAMs are flexible and can handle correlation structures beyond compound symmetry. A detailed description on basis functions and correlations can be found in [@hefley2017]. 

For the data in Figure <a href="#fig:sim-smooth-plot">5.1</a>, A the main effect of interest is how $\mbox{StO}_2$ changes over time for each treatment. To estimate this, the model needs to incorporate independent smooths for _Group_ and _Day_, respectively. The main thing to consider is that model syntax accounts for the fact that one of the variables is numeric ( _Day_ ) and the other is a factor ( _Group_ ). Because the smooths are centered at 0, the factor variable needs to be specified as a parametric term in order to identify any differences between the groups. Using $\textsf{R}$ and the package `mgcv` the model syntax is:

`m1<-gam(StO2_sim~Group+s(Day,by=Group,k=5), method='REML',data=dat_sim)`


This syntax specifies that `m1` will store the model, and that the change in the simulated oxygen saturation (`StO2_sim`) is modeled using independent smooths for  _Group_ and _Day_ (the parenthesis preceded by _s_) using 5 knots. The smooth is constructed by default using thin plate regression splines. Other splines can be used if desired, including gaussian process smooths [@simpson2018]. The parametric term _Group_ is added to quantify differences in the effect of treatment between groups, and the `method` chosen to select the smoothing parameters is the restricted maximum likelihood (REML) [@wood2017]. When the smooths are plotted over the raw data, it is clear that the model has been able to capture the trend of the change of $\mbox{StO}_2$ for each group across time (Figure <a href="#fig:sim-smooth-plot">5.1</a>,B). Model diagnostics can be obtained using the `gam.check` function, and the function `appraise` from the package _gratia_ [@gratia]. A guide for model selection and diagnostics is in the Appendix, and an in-depth analysis can be found in [@wood2017] and [@harezlak2018]. 

One question that might arise at this point is "what is the fit that an rm-ANOVA model produces for the simulated data?". The rm-ANOVA model, which corresponds to Equation (3.1) is presented in Figure <a href="#fig:sim-smooth-plot">5.1</a>,C. This is a typical case of model misspecification: The slopes of each group are different, which would lead to a _p-value_ indicating significance for the treatment and time effects, but the model is not capturing the changes that occur at days 2 and between days 5 and 7, whereas the GAM model is able to do so (Figure <a href="#fig:sim-smooth-plot">5.1</a>,B) .








Because GAMs do not require equally-spaced or complete observations for all subjects, they are advantageous to analyze longitudinal data where missingness exists. The rationale behind this is that GAMs are able to pick the trend in the data even when some observations are missing. However, this usually causes the resulting smooths to have wider confidence intervals and less ability to pick certain trends. Consider the  simulated $\mbox{StO}_2$ values from Figure (<a href="#fig:sim-smooth-plot">5.1</a>, B). If 40% of the total observations are randomly deleted and the same interaction GAM fitted for the complete dataset is used, the resulting smooths are still able to show a different trend for each group, but it can be seen that the smooths overlap during the first 3 days because with less data points, the trend is less pronounced than in the full dataset (<a href="#fig:sim-smooth-plot">5.1</a>, D). Although the confidence intervals have increased for both smooths, the model still shows different trends with as little as 4 observations per group at certain time points.





![Figure 5.1: Simulated data and smooths for oxygen saturation in tumors. A: Simulated data that follows previously reported trends (inset) in tumors under chemotherapy (Treatment) or saline (Control) treatment. Simulated data is from  a normal distribution with standard deviation of 10% with 10 observations per time point. Lines indicate mean oxygen saturation B: Smooths from the GAM model for the full simulated data with interaction of Group and Treatment. Lines represent trends for each group, shaded regions are 95% confidence intervals. C: rm-ANOVA model for the simulated data, the model does not capture the changes in each group over time. D: Smooths for the GAM model for the simulated data with 40% of its observations missing. Lines represent trends for each group, shaded regions are 95% confidence intervals.](Full_document_files/figure-docx/sim-smooth-plot-1.png){width=75% }







![Figure 5.2: Pairwise comparisons for smooth terms. A: Pairwise comparisons for the full dataset.  B: Pairwise comparisons for the dataset with missing observations. Significant differences exist where the interval does not cover 0. In both cases the effect of treatment is significant after day 5.](Full_document_files/figure-docx/plot-pairwise-comp-1.png){width=75% }


## Determination of significance in GAMs for longitudinal data{#GAM-significance}

At the core of a biomedical longitudinal study lies the question of a significant difference between the effect of two or more treatments in different groups. Whereas in rm-ANOVA a _post-hoc_ analysis is required to answer such question by calculating some _p-values_ after multiple comparisons, GAMs can use a different approach to estimate significance. In essence, the idea behind the estimation of significance in GAMs across different treatment groups is that if the _difference_ between the confidence intervals of the fitted smooths for such groups is  non-zero, then a significant difference exists at that time point(s). The absence of a _p-value_ in this case might seem odd, but the confidence interval comparison can be conceptualized in the following manner: Different trends in each group are an indication of an effect by the treatment. This is what happens for the simulated data in Figure <a href="#fig:sim-smooth-plot">5.1</a>, A where the chemotherapy causes $\mbox{StO}_2$ to increase over time. 

With this expectation of different trends in each group, computing the difference between the trends will identify if the observed change is significant. The difference between groups with similar trends is likely to yield zero, which would indicate that the treatment is not causing a change in the response in one of the groups (assuming the other group is a Control or Reference group).

Consider the calculation of pairwise differences for the smooths in Figure <a href="#fig:sim-smooth-plot">5.1</a>, B and D. Figure <a href="#fig:plot-pairwise-comp">5.2</a>, shows the comparison between each treatment group for the full and missing datasets. Here, the  "Control" group is used as the reference to which "Treatment" group is being compared. Of notice, the pairwise comparison has been set on the response scale (see Appendix for code details), because otherwise the comparison appears shifted and is not intuitively easy to relate to the original data.  

With this correction in mind, the shaded region under or above the confidence interval (that does not cover 0) indicates the interval where each group has a higher effect than the other. Notice that the shaded region between days 0 and 3 for the full dataset indicates that through that time, the "Control" group has higher $\mbox{StO}_2$, but as therapy progresses the effect is reversed and by day 4 it is the "Treatment" group the one that has greater $\mbox{StO}_2$. This would suggest that the effect of chemotherapy in the "Treatment" group becomes significant after day 4 for the model used. Moreover, notice that although there is no actual measurement at day 4, the model is capable of providing an estimate of when the shift in $\mbox{StO}_2$ occurs. 

On the data with missing observations (Figure <a href="#fig:sim-smooth-plot">5.1</a>, D), the confidence intervals of the smooths overlap between days 0 and 3. Consequently, the smooth pairwise comparison (Figure <a href="#fig:plot-pairwise-comp">5.2</a>, B) shows that there is not a significant difference between the groups during that period, but is still able to pick the change on day 4 as the full dataset smooth pairwise comparison. 

In a sense, the pairwise smooth comparison is more informative than a _post-hoc_ _p-value_. For biomedical studies, it is able to provide an estimate of _when_ a biological process becomes significant. This is advantageous because it can help researchers gain insight on metabolic changes and other biological processes that can be worth examining, and can help refine the experimental design of future studies in order to obtain measurements at time points where a significant change is expected.


\FloatBarrier


# Discussion

Biomedical longitudinal non-linear data is particularly challenging to analyze due to the likelihood of missing observations and different correlation structures in the data, which limit the use of rm-ANOVA. Although LMEMs have started to replace rm-ANOVA as the choice to analyze biomedical data, both methods yield biased estimates when they are used to fit non-linear data as we have visually demonstrated in Section <a href="#simulation">3.5</a>. This "model misspecification" error, also is known as"" "Type III" error [@dennis2019] is particularly important because although the _p-value_ is the common measure of statistical significance, the validity of its interpretation is determined by the agreement of the data and the model. Guidelines for statistical reporting in biomedical journals exist (the SAMPL guidelines) [@lang2015] but they have not been widely adopted and in the case of longitudinal data, we consider that researchers would benefit from reporting a visual assessment of the correspondence between the model fit and the data, instead of merely relying on a $R^2$ value.

In this paper we have presented GAMs as a suitable method to analyze non-linear longitudinal data. It is interesting to note that although GAMs are a well established method to analyze temporal data in different fields (among which are palaeoecology, geochemistry, and ecology) [@hefley2017; @pedersen2019] they are not routinely used in biomedical research despite an early publication from Hastie and Tibshirani that demonstrated their use in medical research [@hastie1995]. This is possibly due to the fact that the theory behind GAMs can seem very different from that of rm-ANOVA and LMEMs, but the purpose of Section <a href="#GAM-theory">4</a> is to demonstrate that at its core the theory quite simple: Instead of using a linear relationship to model the response (as rm-ANOVA and LMEMs do), GAMs use basis functions to build smooths that are capable of following non-linear trends in the data.

However, from a practical standpoint is equally important to demonstrate how GAMs are computationally implemented. We have provided an example on how GAMs can be fitted using simulated data that follows trends reported in biomedical literature [@vishwanath2009] using $\textsf{R}$ and the package _mgcv_[@wood2017] in Section <a href="#longitudinal-GAMs">5</a>, while a basic workflow for model selection is in the Appendix. One of the features of GAMs is that they go beyond a mere _p-value_ to indicate differences between groups, and in turn provide a time-based estimate of shifts in the response that can be directly tied to biological values as the pairwise smooth comparisons in Figure  <a href="#fig:plot-pairwise-comp">5.2</a> indicate. The model is therefore able to provide an estimate of significant change between the groups at time points were data was not directly measured even with missing data exists ( $\approx$ day 4 in Figure <a href="#fig:plot-pairwise-comp">5.2</a> A, B ), which can be used by researchers as feedback on experiment design and to further evaluate important biological changes in future studies. 

We have used $\textsf{R}$ as the software of choice for this paper because not only provides a fully developed environment to fit GAMs, but also eases simulation (which is becoming increasingly used for exploratory statistical analysis and power calculations) and provides powerful and convenient methods of visualization, which are key aspects that biomedical researchers might need to consider to make their work reproducible. In this regard, reproducibility is still an issue in biomedical research [@begley2015;@weissgerber2018], but it is becoming apparent that what other disciplines have experienced in this aspect is likely to impact soon rather than later this field. Researchers need to plan on how they will make their data, code, and any other materials open and accessible as more journals and funding agencies recognize the importance and benefits of open science in biomedical research. We have made all the data and code used in this paper accessible, and we hope that this will encourage other researchers to do the same with future projects.

\FloatBarrier



# Conclusion

We have presented GAMs as a method to analyze longitudinal biomedical data. Future directions of this work will include simulation-based estimations of statistical power using GAMs, as well as demonstrating the prediction capabilities of these models using large datasets. 
By making the data and code used in this paper accessible, we hope to address the need of creating and sharing reproducible work in biomedical research.

# Acknowledgements

This work was supported by the National Science Foundation Career Award (CBET 1751554) and the Arkansas Biosciences Institute.


\FloatBarrier

***


\newpage
# References

<div id="refs"></div>








# (APPENDIX) Appendix {-} 

# Code for Manuscript data

This section presents the code used to generate figures, models and simulated data from Sections 3 and 4 from the main manuscript.

## Compound symmetry and independent errors in linear and quadratic responses

This section simulated linear and quadratic data in the same manner as in Section <a href="#simulation">3.5</a>. The linear simulations using Figure <a href="#fig:linear-cases-Appendix">10.1</a> show in  panels A and D the simulated mean responses and individual data points.  Panels C and G show a visual interpretation of "correlation" in the responses: In panel C, subjects that have a value of the random error $\varepsilon$ either above or below the mean group response are more likely to have other observations that follow the same trajectory, thereby demonstrating correlation in the response. In panel G,because the errors are independent, there is no expectation that responses are likely to follow a similar pattern.  Panels D and H show the predictions from the rm-ANOVA model.

The following code produces a more comprehensive exploration of Figure <a href="#fig:l-q-response">3.1</a> in the main manuscript.


```r
set.seed(1)
##########Section for calculations###########


## Example with linear response

#This function simulates data using a linear or quadratic mean response and each with correlated
#or uncorrelated errors. Each group has a different slope/concavity.
example <- function(n_time = 6, #number of time points
                    fun_type = "linear", #type of response
                    error_type = "correlated") {
  
  if (!(fun_type %in% c("linear", "quadratic")))
    stop('fun_type must be either "linear", or "quadratic"')
  if (!(error_type %in% c("correlated", "independent")))
    stop('fun_type must be either "correlated", or "independent"')
  
  
  x <- seq(1,6, length.out = n_time)
  
  #Create mean response matrix: linear or quadratic
  mu <- matrix(0, length(x), 2)
  # linear response
  if (fun_type == "linear") {
    mu[, 1] <- - (0.25*x)+2  
    mu[, 2] <- 0.25*x+2
  } else {
    # quadratic response (non-linear)
    
    mu[, 1] <-  -(0.25 * x^2) +1.5*x-1.25
    mu[, 2] <- (0.25 * x^2) -1.5*x+1.25
  }
  
  #create an array where individual observations per each time point for each group are to be stored. Currently using 10 observations per timepoint
  y <- array(0, dim = c(length(x), 2, 10))
  
  #Create array to store the "errors" for each group at each timepoint. The "errors" are the 
  #between-group variability in the response.
  errors <- array(0, dim = c(length(x), 2, 10))
  #create an array where 10 observations per each time point for each group are to be stored
  
  #The following cycles create independent or correlated responses. To each value of mu (mean response per group) a randomly generated error (correlated or uncorrelated) is added and thus the individual response is created.
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
        # compound symmetry errors: variance covariance matrix
        errors[, i, j] <- rmvn(1, rep(0, length(x)), 0.1 * diag(6) + 0.25 * matrix(1, 6, 6))
        y[, i, j] <- mu[, i] + errors[, i, j]
      }
    }
  }    
  
  
  ## subject random effects
  
  ## visualizing the difference between independent errors and compound symmetry
  ## why do we need to account for this -- overly confident inference
  
#labeling y and errors  
  dimnames(y) <- list(time = x, 
                      treatment = 1:2, 
                      subject = 1:10)

  dimnames(errors) <- list(time = x, 
                           treatment = 1:2, 
                           subject = 1:10)
  
  #labeling the mean response
  dimnames(mu) <- list(time = x, 
                       treatment = 1:2)
  
  #convert y, mu and errors to  dataframes with time, treatment and subject columns
  dat <- as.data.frame.table(y, 
                             responseName = "y")
  dat_errors <- as.data.frame.table(errors, 
                                    responseName = "errors")
  dat_mu <- as.data.frame.table(mu, 
                                responseName = "mu")
  
  #join the dataframes to show mean response and errors per subject
  dat <- left_join(dat, dat_errors, 
                   by = c("time", "treatment", "subject"))
  dat <- left_join(dat, dat_mu, 
                   by = c("time", "treatment"))
  #add time
  dat$time <- as.numeric(as.character(dat$time))
  #label subjects per group
  dat <- dat %>%
    mutate(subject = factor(paste(subject, 
                                  treatment, 
                                  sep = "-")))
  
  
  ## repeated measures ANOVA 
  
  fit_anova <- lm(y ~ time + treatment + time * treatment, data = dat)
  
#LMEM: time and treatment interaction model, compound symmetry 
  fit_lme <- lme(y ~ treatment + time + treatment:time,
                 data = dat,
                 random = ~ 1 | subject,
                 correlation = corCompSymm(form = ~ 1 | subject)
  )
  
  #create a prediction frame where the model can be used for plotting purposes
  pred_dat <- expand.grid(
    treatment = factor(1:2), 
    time = unique(dat$time)
  )
  
  #add model predictions to the dataframe that has the simulated data
  dat$pred_anova <- predict(fit_anova)
  dat$pred_lmem <- predict(fit_lme)

  #return everything in a list
  return(list(
    dat = dat,
    pred_dat = pred_dat,
    fit_anova=fit_anova,
    fit_lme = fit_lme
  ))
}
##################Section for plotting#################################
#######################################################################
#This function will create the plots for either a "linear" or "quadratic" response

plot_example <- function(sim_dat) {
  ## Plot the simulated data (scatterplot)
  
  p1 <- sim_dat$dat %>%
    ggplot(aes(x = time, 
               y = y, 
               group = treatment, 
               color = treatment)
           ) +
    geom_point(show.legend=FALSE) +
    labs(y='response')+
    geom_line(aes(x = time, 
                  y = mu, 
                  color = treatment),
              show.legend=FALSE) +
    theme_classic() +
    theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    thm
  
  #plot the simulated data with trajectories per each subject
  p2 <- sim_dat$dat %>%
    ggplot(aes(x = time, 
               y = y, 
               group = subject, 
               color = treatment)
           ) +
    geom_line(aes(size = "Subjects"),
              show.legend = FALSE) +
    # facet_wrap(~ treatment) +
    geom_line(aes(x = time, 
                  y = mu, 
                  color = treatment,
                  size = "Simulated Truth"), 
              lty = 1,show.legend = FALSE) +
    labs(y='response')+
    scale_size_manual(name = "Type", values=c("Subjects" = 0.5, "Simulated Truth" = 3)) +
    theme_classic()+
     theme(plot.title = element_text(size = 30, 
                                face = "bold"),
     text=element_text(size=30))+
    thm
  
  #plot the errors
   p3 <- sim_dat$dat %>%
    ggplot(aes(x = time, 
               y = errors,
               group = subject, 
               color = treatment)) +
    geom_line(show.legend=FALSE) +
     labs(y='errors')+
     theme_classic()+
     theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    thm
  
   #plot the model predictions for rm-ANOVA
  p4 <- ggplot(sim_dat$dat, 
               aes(x = time, 
                   y = y, 
                   color = treatment)) +
    geom_point(show.legend=FALSE)+
    labs(y='response')+
    geom_line(aes(y = predict(sim_dat$fit_anova), 
                  group = subject, size = "Subjects"),show.legend = FALSE) +
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
    theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    thm
   
   
   
   #plot the LMEM predictions
  p5 <- ggplot(sim_dat$dat, 
               aes(x = time, 
                   y = y, 
                   color = treatment)) +
    geom_point()+
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
    theme(plot.title = element_text(size = 30, 
                                  face = "bold"),
        text=element_text(size=30))+
    thm
  
  return((p1+p3+p2+p4+p5)+plot_layout(nrow=1)+plot_annotation(tag_levels = 'A')) 
  
    
}

txt<-18

#Store each plot in a separate object
A1<-plot_example(example(fun_type = "linear", error_type = "correlated")) 

B1<-plot_example(example(fun_type = "linear", error_type = "independent")) 
  
C1<-plot_example(example(fun_type = "quadratic", error_type = "correlated")) 
  
D1<-plot_example(example(fun_type = "quadratic", error_type = "independent")) 
```





![Figure 10.1: Simulated linear responses from two groups with correlated (top row) or independent (bottom row) errors using a rm-ANOVA model and a LMEM. A, F:Simulated data with known mean response and individual responses (points) showing the dispersion of the data. B,G: Generated errors showing the difference in the behavior of correlated and independent errors. C,H: Simulated data with thin lines representing individual trajectories. D,I: Estimations from the rm-ANOVA model for the mean group response. E, J: Estimations from the LMEM for the mean group response and individual responses (thin lines). In all panels, thick lines are the predicted mean response per group, thin lines are the random effects for each subject and points represent the original raw data. Both rm-ANOVA and the LMEM are able to capture the trend of the data.](Full_document_files/figure-docx/linear-cases-Appendix-1.png)

For the quadratic response case, Figure <a href="#fig:quadratic-cases-Appendix">10.2</a> shows the simulated responses using compound symmetry and independent errors. 





![Figure 10.2: Simulated quadratic responses from two groups with correlated (top row) or independent (bottom row) errors using a rm-ANOVA model and a LMEM. A, F:Simulated data with known mean response and individual responses (points) showing the dispersion of the data. B,G: Generated errors showing the difference in the behavior of correlated and independent errors. C,H: Simulated data with thin lines representing individual trajectories. D,I: Estimations from the rm-ANOVA model for the mean group response. E, J: Estimations from the LMEM for the mean group response and individual responses (thin lines). In all panels, thick lines are the predicted mean response per group, thin lines are the random effects for each subject and points represent the original raw data. Both rm-ANOVA and the LMEM are not able to capture the changes in each group over time.](Full_document_files/figure-docx/quadratic-cases-Appendix-1.png)


## Basis functions and GAMs

This code produces Figure <a href="#fig:basis-plot">4.1</a> from the main manuscript. Briefly, a non-linear (quadratic) response is simulated a gam model is fitted and the basis are extracted in order to explain how the smooth is constructed. The code for data simulation is used again here for the sake of keeping the same structure, although the data can be simulated in a more simple fashion.


```r
#generate the response: the same initial procedure from the previous section to simulate
#the response
set.seed(1)
n_time = 6
 x <- seq(1,6, length.out = n_time)
 mu <- matrix(0, length(x), 2)
 mu[, 1] <-  -(0.25 * x^2) +1.5*x-1.25 #mean response
 mu[, 2] <- (0.25 * x^2) -1.5*x+1.25 #mean response
 y <- array(0, dim = c(length(x), 2, 10))
 errors <- array(0, dim = c(length(x), 2, 10))
 for (i in 1:2) {     # number of treatments
     for (j in 1:10) {  # number of subjects
         # compound symmetry errors
         errors[, i, j] <- rmvn(1, rep(0, length(x)), 0.1 * diag(6) + 0.25 * matrix(1, 6, 6))
         y[, i, j] <- mu[, i] + errors[, i, j]
     }
 }
 
 #label each table
  dimnames(y) <- list(time = x, treatment = 1:2, subject = 1:10)
 dimnames(errors) <- list(time = x, treatment = 1:2, subject = 1:10)
 dimnames(mu) <- list(time = x, treatment = 1:2)
 
 #Convert to dataframes with subject, time and group columns
 dat <- as.data.frame.table(y, responseName = "y")
 dat_errors <- as.data.frame.table(errors, responseName = "errors")
 dat_mu <- as.data.frame.table(mu, responseName = "mu")
 dat <- left_join(dat, dat_errors, by = c("time", "treatment", "subject"))
 dat <- left_join(dat, dat_mu, by = c("time", "treatment"))
 dat$time <- as.numeric(as.character(dat$time))
 
 #label subject per group
 dat <- dat %>%
     mutate(subject = factor(paste(subject, treatment, sep = "-")))
  
 #extract  "Group 1" to fit the GAM
  dat<-subset(dat,treatment==1)
 #keep just the response and timepoint columns
   dat<-dat[,c('y','time')]

   #GAM model of time, 5 knots
gm<-gam(y~s(time,k=5),data=dat)

#model_matrix (also known as) 'design matrix'
#will contain the smooths used to create  model 'gm'
model_matrix<-as.data.frame(predict(gm,type='lpmatrix'))


time<-c(1:6)

basis<-model_matrix[1:6,] #extracting basis (because the values are repeated after every 6 rows)
#basis<-model_matrix[1:6,-1] #extracting basis
colnames(basis)[colnames(basis)=="(Intercept)"]<-"s(time).0"
basis<-basis %>% #pivoting to long format
  pivot_longer(
    cols=starts_with("s")
  )%>%
  arrange(name) #ordering

#length of dataframe to be created: number of knots by number of timepoints (minus 1 for the intercept that we won't plot)
ln<-6*(length(coef(gm))) 

basis_plot<-data.frame(Basis=integer(ln),
                       value_orig=double(ln),
                       time=integer(ln),
                       cof=double(ln)
)

basis_plot$time<-rep(time) #pasting timepoints
basis_plot$Basis<-factor(rep(c(1:5),each=6)) #pasting basis number values
basis_plot$value_orig<-basis$value #pasting basis values
basis_plot$cof<-rep(coef(gm)[1:5],each=6) #pasting coefficients
basis_plot<-basis_plot%>%
  mutate(mod_val=value_orig*cof) #the create the predicted values the bases need to be 
#multiplied by the coefficients

#creating labeller to change the labels in the basis plots

basis_names<-c(
  `1`="Intercept",
  `2`="1",
  `3`="2",
  `4`="3",
  `5`="4"
)

#calculating the final smooth by aggregating the basis functions

smooth<-basis_plot%>% 
  group_by(time)%>%
  summarize(smooth=sum(mod_val))


#original basis
sz<-1
p11<-ggplot(basis_plot,
            aes(x=time,
                y=value_orig,
                colour=as.factor(Basis)
                )
            )+
  geom_line(size=sz,
            show.legend=FALSE)+
  geom_point(size=sz+1,
             show.legend = FALSE)+
  labs(y='Basis functions')+
  facet_wrap(~Basis,
             labeller = as_labeller(basis_names)
             )+
  theme_classic()+
  thm
  

#penalized basis
p12<-ggplot(basis_plot,
            aes(x=time,
                y=mod_val,
                colour=as.factor(Basis)
                )
            )+
  geom_line(show.legend = FALSE,
            size=sz)+
  geom_point(show.legend = FALSE,
             size=sz+1)+
  labs(y='Penalized \n basis functions')+
  scale_y_continuous(breaks=seq(-1,1,1))+
  facet_wrap(~Basis,
             labeller=as_labeller(basis_names)
             )+
  theme_classic()+
  thm

#heatmap of the coefficients
x_labels<-c("Intercept","1","2","3","4")
p13<-ggplot(basis_plot,
            aes(x=Basis,
                y=Basis))+
  geom_tile(aes(fill = cof), 
            colour = "black") +
    scale_fill_gradient(low = "white",
                        high = "#B50A2AFF")+ #color picked from KikiMedium
  labs(x='Basis',
       y='Basis')+
  scale_x_discrete(labels=x_labels)+
  geom_text(aes(label=round(cof,2)),
            size=7,
            show.legend = FALSE)+
  theme_classic()+
  theme(legend.title = element_blank())
  
#plotting simulated datapoints and smooth term
p14<-ggplot(data=dat,
            aes(x=time,y=y))+
  geom_point(size=sz+1)+
  labs(y='Simulated \n response')+
  geom_line(data=smooth,
            aes(x=time,
                y=smooth),
            color="#6C581DFF",
            size=sz+1)+
  theme_classic()
  

#Combining all
b_plot<-p11+p13+p12+p14+plot_annotation(tag_levels='A')&
  theme(
     text=element_text(size=18)
     )
```



![Figure 10.3: Basis functions for a single smoother for time with five knots. A: Basis functions for a single smoother for time for the simulated data of Group 1 from Figure 2, the intercept basis is not shown. B: Penalty matrix for the basis functions. Each basis function is penalized by a coefficient which can be positive or negative. The coefficient determines the overall effect of each basis in the final smoother. C: Penalized basis functions. Each of the four basis functions of panel A has been penalized by the corresponding coefficient shown in Panel B, note the corresponding increase (or decrease) of each basis. D: Smoother for time and original data points. The smoother (line) is the result of the sum of each penalized basis function at each time point, simulated values for the group appear as points.](Full_document_files/figure-docx/basis-plot-appendix-1.png){width=75% }

# Longitudinal biomedical data simulation and GAMs{#tumor-data-simulation}

This section describes how to fit GAMs to longitudinal data using simulated data. First, data is simulated according to Section <a href="#longitudinal-GAMs">5</a>, where reported data of oxygen saturation ($\mbox{StO}_2$) in tumors under either chemotherapy or saline control is used as a starting point to generate individual responses in each group.


```r
set.seed(1)
#Dataframe that contains the original reported trends
dat<-tibble(StO2=c(4,27,3,2,0.5,7,4,50,45,56),
            Day=rep(c(0,2,5,7,10),times=2),
            Group=as.factor(rep(c("Control","Treatment"),each=5))
            )


## plot the mean response
f1<-ggplot(dat, 
           aes(x = Day, 
               y = StO2, 
               color = Group)) +
    geom_line(size=1,
              show.legend = FALSE)+
    geom_point(show.legend = FALSE,
               size=1.5,
               alpha=0.5)+
  labs(y=expression(paste(StO[2],
                          ' (real)')))+
  theme_classic()+
  thm+
    scale_x_continuous(breaks=c(0,5,10))+
    scale_y_continuous(breaks=c(0,40))+
  plot_layout(tag_level = 'new')+
  theme(
    plot.background = element_rect(fill = "transparent", 
                                   color = NA),
    axis.text=element_text(size=14)
  )


#This function simulates data for the tumor data using default parameters of 10 observations per time point,and Standard deviation (sd) of 5%.
#Because physiologically StO2 cannot go below 0%, data is  generated with a cutoff value of 0.0001 (the "StO2_sim")

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

#plotting simulated data
f2<-ggplot(dat_sim, 
           aes(x = Day, 
               y = StO2_sim, 
               color = Group)) +
    geom_point(show.legend=FALSE,
               size=1.5,
               alpha=0.5)+
    stat_summary(aes(y = StO2_sim,
                     group=Group), 
                 fun=mean, geom="line",
                 size=1,
                 show.legend = FALSE)+
  labs(y=expression(atop(StO[2],
                         '(simulated)')))+
  theme_classic()+
  theme(
    axis.text=element_text(size=22)
  )+
  thm+
    scale_x_continuous(breaks=c(0,2,5,7,10))
```

## A basic Workflow for GAMs

This section shows a basic workflow to fit a series of increasingly complex GAMs to the simulated data from the previous section. Graphical and parameter diagnostics for goodness of fit are discussed, as well as model comparison via AIC (Aikake Information Criterion). 

### First model

The first model fitted to the data is one that only accounts for different smooths by day. The model syntax specifies that `gam_00` is the object that will contain all the model information, and that the model attempts to explain changes in `StO2_sim` (simulated $\mbox{StO}_2$) using a smooth per `Day`. The model will use 5 knots (`k=5`) for the smooth. The smooth is constructed by default using thin plate regression splines. The smoothing parameter estimation method used is the restricted maximum likelihood (`REML`).



```r
gam_00<-gam(StO2_sim ~ s(Day, k = 5),
            method='REML',
            data  = dat_sim)
```

To obtain model diagnostics, two methodologies are to be used: 1) graphical diagnostics, and 2) a model check. In the first case, the functions `appraise` and `draw` from the package _gratia_ can be used to obtain a single output with all the graphical diagnostics. For model check, the functions `gam.check` and `summary` from  _mgcv_ provide detailed information about the model fit and its parameters. 




```r
#see https://patchwork.data-imaginist.com/reference/area.html
layout1 <- c(
  area(1, 1),
  area( 1, 2),
  area(2, 1),
  area(2, 2),
  area(1, 3, 2)
)


layout2 <- c(
  area(1, 1),
  area( 1, 2),
  area(2, 1),
  area(2, 2),
  area(1,3,2,5)
)

#plot(layout2)
```

#### Graphical diagnostics

![Figure 11.1: Graphical diagnostics for the first GAM model. Left: Graphical diagnostics provided by the function `appraise` from the package _gratia_. Right: Fitted smooth for the model, provided by the function `draw`.](Full_document_files/figure-docx/first-GAM-diag-1.png){width=75% }

From the output of the function `appraise` in Figure <a href="#fig:first-GAM-diag">11.1</a>, the major indicators of concern about the model are the QQ plot of residuals and the histogram of residuals. The QQ plot shows that the errors are not reasonably located along the 45$^{\circ}$ line (which indicates normality), as there are multiple points that deviate from the trend, specially in the tails. The histogram also shows that the variation (residuals) is not following the assumption of a normal distribution. 

The `draw` function permits to plot the smooths as `ggplot2` objects, which eases subsequent manipulation, if desired. Because model `gam_00` specifies only one smooth for the time covariate (Day), the plot only contains only one smooth. Note that the smooth shows an almost linear profile.

#### Model check


```r
#need to add figure number and caption
gam.check(gam_00)
```

![](Full_document_files/figure-docx/first-GAM-check-1.png)

```
## 
## Method: REML   Optimizer: outer newton
## full convergence after 5 iterations.
## Gradient range [-6.11034e-07,-1.169842e-07]
## (score 439.1428 & scale 414.047).
## Hessian positive definite, eigenvalue range [0.05006795,49.00066].
## Model rank =  5 / 5 
## 
## Basis dimension (k) checking results. Low p-value (k-index<1) may
## indicate that k is too low, especially if edf is close to k'.
## 
##          k'  edf k-index p-value    
## s(Day) 4.00 1.36    0.26  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(gam_00)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## StO2_sim ~ s(Day, k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   21.929      2.035   10.78   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##          edf Ref.df     F p-value   
## s(Day) 1.359  1.624 6.695 0.00273 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.106   Deviance explained = 11.8%
## -REML = 439.14  Scale est. = 414.05    n = 100
```

Special attention must be paid to the 'k-index' from `gam.check`. This parameter indicates if the basis dimension of the smooth is adequate, i.e., it checks that the basis used to create the smooth are adequate to capture the trends in the data. If the model is not adequately capturing the trens in the data, this is indicated by a low k-index value (<1). From the output, it can be seen that the k-index is effectively <0.3, which indicates that the model is not capturing the variability in the data. The 'edf' (effective degrees of freedom) is an indicator of the complexity of the smooth. Here the complexity of the smooth is comparable to that of a 4th degree polynomial.

From the `summary` function, information about the assumed distribution of the errors (Gaussian in this case) and the link function can be obtained. The link function is 'identity' as the model does not make any transformation on the predictors. The 'significance of smooth terms' _p-value_ indicates if each smooth is adding significance to the model. Here, the _p-value_ is low but we have seen that there are issues with the model from the previous outputs. Finally, the 'deviance explained' indicates how much of the data the model is able to capture, which in this case corresponds to $\approx$ 12%.

### Second model

The major flaw of `gam_00` is that this model is not taking into account the fact that the data is nested in groups. The next iteration is a model where a different smooth of time (Day) is assigned for each group using `by=Group` in the model syntax.



```r
gam_01<-gam(StO2_sim ~ s(Day, by=Group,k = 5),
            method='REML',
            data  = dat_sim)

gam.check(gam_01)
```

![](Full_document_files/figure-docx/second-GAM-1.png)

```
## 
## Method: REML   Optimizer: outer newton
## full convergence after 11 iterations.
## Gradient range [-3.57434e-06,1.383186e-06]
## (score 413.523 & scale 230.4732).
## Hessian positive definite, eigenvalue range [0.1532335,48.55232].
## Model rank =  9 / 9 
## 
## Basis dimension (k) checking results. Low p-value (k-index<1) may
## indicate that k is too low, especially if edf is close to k'.
## 
##                         k'  edf k-index p-value    
## s(Day):GroupControl   4.00 3.49    0.38  <2e-16 ***
## s(Day):GroupTreatment 4.00 2.96    0.38  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(gam_01)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## StO2_sim ~ s(Day, by = Group, k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   21.929      1.518   14.45   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                         edf Ref.df      F p-value    
## s(Day):GroupControl   3.488  3.851  5.244 0.00442 ** 
## s(Day):GroupTreatment 2.962  3.461 24.272 < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.502   Deviance explained = 53.5%
## -REML = 413.52  Scale est. = 230.47    n = 100
```

Diagnostics for this model indicate that the k-index is still below 1 (0.32 from `gam.check`), and that the residuals are still not following a normal distribution (Figure  <a href="#fig:second-GAM-diag">11.2</a>). Moreover, the smooths (plotted via the `draw()` function) appear with a fairly linear profile, which indicates they are still not capturing the trends observed in the data. From `summary()`, the deviance explained by the model is $\approx$ 43%.



![Figure 11.2: Graphical diagnostics for the second GAM model. Left: Graphical diagnostics provided by the function `appraise` from the package _gratia_. Right: Fitted smooth for the model, provided by the function `draw`.](Full_document_files/figure-docx/second-GAM-diag-1.png){width=75% }

### Third model

Model `gam_00` was built for didactic purposes to cover the simplest case, but it does not account for the nesting of the data by Group, which is apparent from the type of smooth fitted, the model diagnostics, and, the low variance explained by the model. On the other hand, `gam_01` takes into account the nesting within each group and provides better variance explanation, but as indicated in Section <a href="#longitudinal-GAMs">5</a>, in order to differentiate between each group a parametric term needs to be added to the model for the interaction of _Day_ and _Group_. 


```r
#GAM for StO2

m1 <- gam(StO2_sim ~ Group+s(Day, by = Group, k = 5),
            method='REML',
            data  = dat_sim)

gam.check(m1)
```

![](Full_document_files/figure-docx/final-model-Appendix-1.png)

```
## 
## Method: REML   Optimizer: outer newton
## full convergence after 9 iterations.
## Gradient range [-2.780424e-08,2.076237e-08]
## (score 354.6068 & scale 63.7304).
## Hessian positive definite, eigenvalue range [1.095531,48.08644].
## Model rank =  10 / 10 
## 
## Basis dimension (k) checking results. Low p-value (k-index<1) may
## indicate that k is too low, especially if edf is close to k'.
## 
##                         k'  edf k-index p-value
## s(Day):GroupControl   4.00 3.87    1.02    0.52
## s(Day):GroupTreatment 4.00 3.83    1.02    0.58
```

```r
summary(m1)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## StO2_sim ~ Group + s(Day, by = Group, k = 5)
## 
## Parametric coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       9.781      1.129   8.664 1.68e-13 ***
## GroupTreatment   24.296      1.597  15.217  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                         edf Ref.df     F p-value    
## s(Day):GroupControl   3.867  3.989 19.38  <2e-16 ***
## s(Day):GroupTreatment 3.826  3.981 80.29  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.862   Deviance explained = 87.4%
## -REML = 354.61  Scale est. = 63.73     n = 100
```

The resulting model is`m1`, which is the model fitted in the main manuscript. By using `appraise()` and `draw` on this model (Figure  <a href="#fig:final-GAM-diag">11.3</a>) we see that the trend on the QQ plot has improved, the histogram of the residuals appears to be reasonably distributed, and the smooths are capturing the trend of the data within each group . From `gam.check`, the k-index is now at an acceptable value ($\approx$ 1.02), and `summary` now indicates that the model is able to capture 87% of the variance data.




![Figure 11.3: Graphical diagnostics for the final GAM model. Left: Graphical diagnostics provided by the function `appraise` from the package _gratia_. Right: Fitted smooths for the model, provided by the function `draw`.](Full_document_files/figure-docx/final-GAM-diag-1.png){width=75% }


### Comparing models via AIC

One final  comparison that can be made for model selection involves the use of the Aikake Information Criterion (AIC). This metric is used to estimate information loss, which we want to minimize with an appropriate model. Therefore, when 2 or more models are compared, the model with lower AIC is preferred. In R, the comparison is done using the `AIC` function.


```r
AIC(gam_00,gam_01,m1)
```

```
##               df      AIC
## gam_00  3.623938 891.2466
## gam_01  9.312053 838.6825
## m1     10.970436 710.9994
```

The output in this case is expected: model `gam1` has a lower AIC (711.46) whereas the initial two models have higher AICs (891 and 850). The AIC should not be considered as the only estimator of model goodness, instead to be used as complimentary information to the graphical diagnostics and model checks described above.

#### Pairwise comparisons of smooth confidence intervals

The estimation of significant differences between each treatment group can be achieved via pairwise comparisons of the smooth confidence intervals as described in section <a href="#GAM-significance">5.3</a>. In this case, the "design matrix" is used to estimate the pairwise comparisons (see main manuscript for details and associated references). Briefly, the "design matrix" (also known as the "Xp matrix") from the selected model (`m1`) is used to calculate a 95% confidence interval of the difference between the smooth terms for each group. This approach allows to estimate the time intervals where a significant difference exists between the groups (confidence interval above or below 0). **All pairwise comparisons in this paper have been centered at the response scale to ease interpretation **.


```r
##Pairwise comparisons
pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Control', 'Treatment'))

##matrix that contains the basis functions evaluated at the points in pdat
    xp <- predict(m1, newdata = pdat, type = 'lpmatrix')

    
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

    ## remove columns that do not contain name 'Control' or 'Treatment'
    X[, ! (c1 | c2)] <- 0
    ## zero out the parametric cols, those that do not contain in the characters 's('
    #X[, !grepl('^s\\(', colnames(xp))] <- 0

    #Multiply matrix by model coefficients. X has (p,n) (rows, columns) and the coefficient matrix has
    #dimensions (n,1). The resulting matrix has dimensions (p,1)
    dif <- X %*% coef(m1)

    #comp<-test %*% coef(gam1)[3:10]

#Calculate standard error for the computed differences using the variance-covariance matrix
    #of the model
    se <- sqrt(rowSums((X %*% vcov(m1, unconditional = FALSE)) * X))
    crit <- qt(0.05/2, df.residual(m1), lower.tail = FALSE)
    #upper  limits
    upr <- dif + (crit * se)
    #lower limits
    lwr <- dif - (crit * se)
    #put all components in a dataframe for plotting
    comp1<-data.frame(pair = paste('Control', 'Treatment', sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)



#add time point sequence
comp_StO2 <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp1))

#plot the difference
c1<-ggplot(comp_StO2, aes(x = Day, y = diff, group = pair)) +
  #ribbon for difference confidence interval  
  geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.5,
                fill='#DB3A07FF') +
    geom_line(color='black',size=1) +
    geom_line(data=comp_StO2,aes(y=0),size=0.5)+
  #highlight area under the curve where "Control" is higher
  geom_ribbon(data=comp_StO2%>%
                    filter(lower>0),
                aes(ymin =0, ymax =lower),
                alpha = 0.5,
                fill='#30123BFF') +
  #highlight area under the curve where "Treatment" is higher
  geom_ribbon(data=comp_StO2 %>%
                    filter(upper<0),
                    aes(ymin =0, ymax =upper),
                alpha = 0.5,
                fill='#7A0403FF') +
    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )
```

(ref:pairwise-comp-caption) Smooth pairwise comparisons for model `m1` using a 95% confidence interval for the difference between smooths. 

![Figure 11.4: (ref:pairwise-comp-caption)](Full_document_files/figure-docx/pairwise-comp-workflow-fig-1.png){width=75% }

Of notice, a convenient wrapper for the function described above exists in the package `gratia`. In this package, `difference_smooths` is a function that makes the comparisons and produces Figure <a href="#fig:pairwise-comp-workflow-fig">11.4</a> when is used on a fitted model. The function syntax and an example can be found at:

https://cran.r-project.org/web/packages/gratia/gratia.pdf

Keep in mind that this function **does not** center the pairwise comparison at the response scale, so it has to be shifted in order to be compared to the raw data.


# GAM and Linear model plots and Missing data

This section covers the code used to generate Figure <a href="#fig:sim-smooth-plot">5.1</a>, where the simulated data, fit of the "final" GAM (`m1`), linear model  and GAM on data with missing observations are presented. Note that panel A in Figure <a href="#fig:sim-smooth-plot">5.1</a> and the inlet are generated in the code chunk where the data is simulated in Section <a href="#tumor-data-simulation">11</a>, and are called later to build the figure.

## GAM and Linear model plots

This code chunk creates panels B and D in Figure <a href="#fig:sim-smooth-plot">5.1</a>. Note that this code uses the final GAM from the previous section (`m1`), so the simulated data and the model should be generated before running this section. 


```r
#linear model
lm1<-lm(StO2_sim ~ Day + Group + Day * Group, data = dat_sim)


#creates a dataframe using the length of the covariates for the GAM
gam_predict <- expand_grid(Group = factor(c("Control", "Treatment")),
                         Day = seq(0, 10, by = 0.1),
                         subject=factor(rep(1:10)))

#creates a dataframe using the length of the covariates for rm-ANOVA
lm_predict<-expand_grid(Group = factor(c("Control", "Treatment")),
                         Day = c(0:10),
                        subject=factor(rep(1:10)),
                          )
lm_predict$subject<-factor(paste(lm_predict$subject, lm_predict$Group, sep = "-"))

#adds the predictions to the grid and creates a confidence interval for GAM
gam_predict<-gam_predict%>%
    mutate(fit = predict(m1,gam_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(m1, gam_predict,se.fit = TRUE,type='response')$se.fit)

#using lm
lm_predict<-lm_predict%>%
    mutate(fit = predict(lm1,lm_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(lm1, lm_predict,se.fit = TRUE,type='response')$se.fit)

#plot smooths and confidence interval for GAM
f3<-ggplot(data=dat_sim, aes(x=Day, y=StO2_sim, group=Group)) +
    geom_point(aes(color=Group),size=1.5,alpha=0.5,show.legend = FALSE)+
  geom_ribbon(aes( x=Day,ymin=(fit - 2*se.fit), 
                   ymax=(fit + 2*se.fit),
                   fill=Group
                   ),
              alpha=0.3,
              data=gam_predict,
            show.legend=FALSE,
                inherit.aes=FALSE) +
  geom_line(aes(y=fit,
                color=Group),
              size=1,data=gam_predict,
              show.legend = FALSE)+
  #facet_wrap(~Group)+
  labs(y=expression(atop(StO[2],'complete')))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
      theme_classic()+
  theme(
    axis.text=element_text(size=22)
  )+
      thm+
  thm1
 
#plot linear fit for rm-ANOVA
f4<-ggplot(data=dat_sim, aes(x=Day, y=StO2_sim, group=Group)) +
    geom_point(aes(color=Group),size=1.5,alpha=0.5,show.legend = FALSE)+
  geom_ribbon(aes( x=Day,ymin=(fit - 2*se.fit), 
                   ymax=(fit + 2*se.fit),fill=Group),
              alpha=0.3,
              data=lm_predict,
              show.legend = FALSE,
                inherit.aes=FALSE) +
  geom_line(aes(y=fit,
                color=Group),
              size=1,data=lm_predict,
              show.legend = FALSE)+
  #facet_wrap(~Group)+
  labs(y=expression(paste('StO'[2],' (simulated)')))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
      theme_classic()+
  theme(
    axis.text=element_text(size=22)
  )+
      thm+
  thm1
 


#posthoc comparisons for the linear model
#library(multcomp)


#summary(glht(lm1, linfct = mcp(Group = 'Tukey')))
#summary(glht(lm1, linfct=mcp(Group="Tukey", interaction_average=TRUE)))
```

## Working with Missing data in GAMs

This code chunk first randomly deletes 40% of the total observations in the original simulated data, and then an interaction GAM is fitted to the remaining data. Model diagnostics are presented, and an object that stores the fitted smooths is saved to be called in the final code chunk to build the figure.


```r
#missing data
#create a sequence of 40 random numbers between 1 and 100, these numbers will
#correspond to the row numbers to be randomly erased from the original dataset

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
```

```
## # A tibble: 10 x 3
## # Groups:   Day, Group [10]
##      Day Group         n
##    <dbl> <fct>     <int>
##  1     0 Control       2
##  2     0 Treatment     4
##  3     2 Control       6
##  4     2 Treatment     5
##  5     5 Control       6
##  6     5 Treatment     4
##  7     7 Control       3
##  8     7 Treatment     5
##  9    10 Control       3
## 10    10 Treatment     8
```

```r
#the same model used for the full dataset
mod_m1 <- gam(StO2_sim ~ Group+s(Day,by=Group,k=5), data  = dat_missing,family=scat)
#appraise the model
appraise(mod_m1)
```

![](Full_document_files/figure-docx/missing-data-Appendix-1.png)

```r
m_predict <- expand_grid(Group = factor(c("Control", "Treatment")),
                         Day = seq(0, 10, by = 0.1))

#adds the predictions to the grid and creates a confidence interval
m_predict<-m_predict%>%
    mutate(fit = predict(mod_m1,m_predict,se.fit = TRUE,type='response')$fit,
           se.fit = predict(mod_m1, m_predict,se.fit = TRUE,type='response')$se.fit)


f6<-ggplot(data=dat_missing, aes(x=Day, y=StO2_sim, group=Group)) +
    geom_point(aes(color=Group),size=1.5,alpha=0.5,show.legend = FALSE)+
  geom_ribbon(aes( x=Day,ymin=(fit - 2*se.fit), 
                   ymax=(fit + 2*se.fit),
                   fill=Group
                   ),
              alpha=0.3,
              data=m_predict,
            show.legend=FALSE,
                inherit.aes=FALSE) +
  geom_line(aes(y=fit,
                color=Group),
              size=1,data=m_predict,
              show.legend = TRUE)+
  #facet_wrap(~Group)+
  labs(y=expression(atop(StO[2],'missing')))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
      theme_classic()+
  theme(
    axis.text=element_text(size=22)
  )+
      thm+
  thm1
```





```r
mult_plot<-f2+inset_element(
  f1, left = 0.01, 
  bottom = 0.5, 
  right = 0.5, 
  top = 1.0)+
  f3+f4+f6+
   plot_annotation(tag_levels='A')&
   ylim(c(-7,75)) &
  theme(
     text=element_text(size=18)
     )&
  thm

mult_plot
```

![Figure 12.1: Simulated data and smooths for oxygen saturation in tumors. A: Simulated data that follows previously reported trends (inset) in tumors under chemotherapy (Treatment) or saline (Control) treatment. Simulated data is from  a normal distribution with standard deviation of 10% with 10 observations per time point. Lines indicate mean oxygen saturation B: Smooths from the GAM model for the full simulated data with interaction of Group and Treatment. Lines represent trends for each group, shaded regions are 95% confidence intervals. C: rm-ANOVA model for the simulated data, the model does not capture the changes in each group over time. D: Smooths for the GAM model for the simulated data with 40% of its observations missing. Lines represent trends for each group, shaded regions are 95% confidence intervals.](Full_document_files/figure-docx/sim-smooth-plot-Appendix-1.png){width=75% }


## Pairwise comparisons in GAMs: full and missing data cases

The next code chunk reproduces Figure <a href="#fig:plot-pairwise-comp">5.2</a>. Here pairwise comparisons are made for the full and missing datasets.



```r
##Pairwise comparisons

pdat <- expand.grid(Day = seq(0, 10, length = 400),
                    Group = c('Control', 'Treatment'))

#this function takes the model, grid and groups to be compared using the lpmatrix

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

comp1<-smooth_diff(m1,pdat,'Control','Treatment')

comp_StO2_full <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp1)) %>%
  mutate(interval=case_when(
    upper>0 & lower<0~"no-diff",
    upper<0~"less",
    lower>0~"greater"
  ))

c1<-ggplot(comp_StO2_full, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.5,
                fill='#DB3A07FF') +
    geom_line(color='#E75B64FF',size=1) +
    geom_line(data=comp_StO2_full,aes(y=0),size=0.5)+
    geom_ribbon(data=comp_StO2_full%>%
                    filter(lower>0),
                aes(ymin =0, ymax =lower),
                alpha = 0.5,
                fill='#30123BFF') +
    geom_ribbon(data=comp_StO2_full %>%
                    filter(upper<0),
                    aes(ymin =0, ymax =upper),
                alpha = 0.5,
                fill='#7A0403FF') +
    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )



###for missing data
comp2<-smooth_diff(mod_m1,pdat,'Control','Treatment')
comp_StO2_missing <- cbind(Day = seq(0, 10, length = 400),
                   rbind(comp2))

missing_plot<-ggplot(comp_StO2_missing, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(color='black',size=1) +
    facet_wrap(~ pair) +
    labs(x = 'Days', 
         y = expression(paste('Difference in StO'[2],'\n (missing data)' 
                              )))+
  scale_x_continuous(breaks=c(0,2,5,7,10))+
  theme_classic()+
  theme(
     text=element_text(size=18),
     legend.title=element_blank()
     )

c2<-ggplot(comp_StO2_missing, aes(x = Day, y = diff, group = pair)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.5,
                fill='#DB3A07FF') +
    geom_line(color='#E75B64FF',size=1) +
    geom_line(data=comp_StO2_missing,aes(y=0),size=0.5)+
    # geom_ribbon(data=comp_StO2_missing%>%
    #                 filter(lower>0),
    #             aes(ymin =0, ymax =lower),
    #             alpha = 0.5,
    #             fill='#30123BFF') +
    geom_ribbon(data=comp_StO2_missing %>%
                    filter(upper<0),
                    aes(ymin =0, ymax =upper),
                alpha = 0.5,
                fill='#7A0403FF') +
    facet_wrap(~ pair) +
    theme_classic()+
    labs(x = 'Days', y = expression(paste('Difference in StO'[2] )))+
    scale_x_continuous(breaks=c(0,2,5,7,10))+
    theme(
        text=element_text(size=18),
        legend.title=element_blank()
    )

pair_comp<-c1+c2
```



![Figure 12.2: Pairwise comparisons for smooth terms. A: Pairwise comparisons for the full dataset.  B: Pairwise comparisons for the dataset with missing observations. Significant differences exist where the interval does not cover 0. In both cases the effect of treatment is significant after day 5.](Full_document_files/figure-docx/plot-pairwise-comp-Appendix-1.png){width=75% }






