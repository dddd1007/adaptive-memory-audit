**Supplementary Table S3. Recovery under observation-equation misspecification**

*Panel A. Family and parameter recovery*

| True learner | Condition | Replications | Family recovery | Exact parameter |
|:--|--:|--:|--:|--:|
| RL | Correct | 200 | .965 | .915 |
| RL | Intercept -0.25 | 200 | .430 | .055 |
| RL | Intercept +0.25 | 200 | 1.000 | .010 |
| RL | State slope ×0.80 | 200 | 1.000 | .005 |
| RL | Intercept +0.25; slope ×0.80 | 200 | 1.000 | .000 |
| Bayesian | Correct | 200 | .990 | .585 |
| Bayesian | Intercept -0.25 | 200 | .895 | .010 |
| Bayesian | Intercept +0.25 | 200 | .570 | .000 |
| Bayesian | State slope ×0.80 | 200 | .315 | .005 |
| Bayesian | Intercept +0.25; slope ×0.80 | 200 | .000 | .000 |

*Panel B. Latent-signal recovery*

| True learner | Condition | State RMSE | State r | PE r |
|:--|--:|--:|--:|--:|
| RL | Correct | .251 | .976 | .990 |
| RL | Intercept -0.25 | .358 | .963 | .981 |
| RL | Intercept +0.25 | .362 | .978 | .990 |
| RL | State slope ×0.80 | .414 | .975 | .988 |
| RL | Intercept +0.25; slope ×0.80 | .846 | .956 | .974 |
| Bayesian | Correct | .218 | .962 | .985 |
| Bayesian | Intercept -0.25 | .295 | .932 | .969 |
| Bayesian | Intercept +0.25 | .281 | .959 | .982 |
| Bayesian | State slope ×0.80 | .294 | .957 | .980 |
| Bayesian | Intercept +0.25; slope ×0.80 | .391 | .942 | .970 |

*Note.* Rows retain the candidate family matching the true learner and average over the two scheduler families. PE denotes unsigned prediction error.

