# CR2 coordinate audit

The archived Julia observational estimator linearizes the fitted logistic model in Pearson coordinates. With fitted probability p, it uses Z = diag(sqrt(p(1-p))) X and residual (Y-p)/sqrt(p(1-p)), then applies the symmetric inverse square root of I-Z_g(Z'Z)^(-1)Z_g' to each cluster. This is a CR2 construction with identity working covariance in the transformed coordinates.

clubSandwich 0.7.0's default GLM method uses the response-scale derivative design diag(p(1-p)) X, response residual Y-p, and target diag(p(1-p)). Its symmetric square-root adjustment is constructed in those coordinates. Symmetric matrix square roots are not invariant under a general nonorthogonal change of coordinates, so the two valid working-model adjustments need not give identical finite-sample estimates of uncertainty.

`validate_pearson_cr2.R` independently fits the frozen data with R GLM, forms the Pearson working response Z beta + (Y-p)/sqrt(p(1-p)), and passes the transformed least-squares model to clubSandwich CR2/Satterthwaite. All 51 models and 153 coefficient/SE/df comparisons pass the original strict tolerances. Maximum absolute differences are 1.82e-12, 1.29e-12, and 1.44e-11, respectively. The six linear randomized models already agreed through the standard R interface.

This audit resolves the implementation discrepancy without changing the archived Julia estimator or its simulation outputs. The earlier standard GLM comparison remains preserved, including maximum SE discrepancy 1.0447e-4 and df discrepancy 2.9872e-4. The revised manuscript must identify the Pearson linearization and must not claim exact identity with clubSandwich's default GLM route. The numerical gradient tolerance used for the R audit is stated in its source; no historical DGP or production convergence rule was changed.

Sources inspected on 2026-09-05: the installed clubSandwich 0.7.0 functions `model_matrix.glm`, `residuals_CS.glm`, `targetVariance.glm`, `CR2`, and `vcov_CR`, saved in the local logs; [maintainer documentation](https://jepusto.github.io/clubSandwich/reference/vcovCR.glm.html). Numerical evidence: `pearson_cr2_refits.csv` and `pearson_cr2_comparison.csv`.
