# Model-agnostic variable importance for survival outcomes: a simulation study under additive and mixed hazards

We conducted a simulation study to examine the behavior of Wolock’s 
model-agnostic exclusion variable importance measure (VIM) under a range of 
survival data-generating mechanisms. The objectives were to evaluate its 
performance and inferential validity under additive hazards, mixed 
additive–multiplicative hazards, and proportional hazards settings, including 
high-dimensional scenarios with rare binary covariates and heterogeneous effect 
sizes. We also compared alternative strategies for nuisance function estimation 
and contrasted exclusion-based importance with permutation-based importance 
across all scenarios.

The simulation study was designed to evaluate the performance, robustness, and 
inferential validity of Wolock’s model-agnostic exclusion variable importance 
measure (VIM) across several survival data-generating mechanisms.

First, under an **additive hazards framework**, we assess the accuracy and 
inferential properties of the VIM and compare different strategies for 
estimating the conditional survival functions, namely 
***global survival stacking***, ***survival Super Learner***, the 
***Aalen additive hazards model***, and ***random survival forests***.

Second, under **mixed additive and multiplicative hazard structures**, we 
investigate the behavior of the VIM when both types of effects coexist. In this 
setting, conditional survival functions are estimated using 
***global survival stacking***, ***survival Super Learner***, the 
***Cox–Aalen model***, and ***random survival forests***.

Third, we assess the sensitivity of the VIM when covariates correspond to 
**rare alterations** with **heterogeneous effect sizes**, ranging from weak to 
strong, under a **proportional hazards model**.

Finally, we compare **exclusion-based variable importance** with 
**permutation-based variable importance** to examine differences in 
interpretation, stability, and empirical performance across the considered 
scenarios.

