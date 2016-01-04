[![CRAN Status Badge](http://www.r-pkg.org/badges/version/phylolm)](http://cran.r-project.org/package=phylolm)
[![Research software impact](http://depsy.org/api/package/cran/phylolm/badge.svg)](http://depsy.org/package/r/phylolm)

## phylolm: R package for Phylogenetic Linear Regression on very large trees

The package provides functions for fitting phylogenetic linear models and phylogenetic generalized linear models. 
Computations use an algorithm that is linear in the number of tips in the tree. 
The package also provides functions for simulating continuous or binary traits along the tree.
When a new version is stable, it is pushed to CRAN, so the version available here is newer than that on CRAN.

- Lam Si Tung Ho and Cécile Ané (2014). 
Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models. 
*Methods in Ecology and Evolution* 5(11):1133-1146. 
[(link)](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12285/abstract)

- Lam Si Tung Ho and Cécile Ané (2014). 
A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. 
*Systematic Biology* 63(3):397-408.
[(link to pdf)](http://sysbio.oxfordjournals.org/cgi/reprint/syu005?ijkey=bIsHxa2dpqXCplc&keytype=ref)

### Main features

- phylogenetic signal
- phylogenetic linear, logistic and Poisson regression
- stepwise model selection (from v2.3)
- OU shift detection
- continuous and discrete trait simulators with covariates
- bootstrap-based confidence intervals for phylogenetic logistic regression (from v2.3)
- goodness-of-fit test of a population tree with the coalescent (from v.2.4)
