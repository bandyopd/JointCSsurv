# JointCSsurv
JointCSsurv (which stands for <ins>**Joint**</ins> model for <ins>**C**</ins>luster <ins>**S**</ins>ize and <ins>**surv**</ins>ival outcome) is a package that performs semiparametric estimation and inference for clustered interval-censored data with informative cluster size using the method proposed by Lee et al. (2022) <DOI: xx/xxxx>.

**JointCSsurv** relies on the R-packages `splines2`, `numDeriv`, `statmod`, `plyr`, which are hosted on CRAN.

# How to import the Functions
> install.packages("devtools")<br />
> library(devtools) <br /> 
> source_url("https://github.com/lcyjames/WeibullCMs/blob/845d70f953b0ce692536bfd928fc0c72165fb4a8/CoreFunctions.R?raw=TRUE")

# Usage #
The package contains 3 functions:
|Functions  | Description|
|------------- | -------------|
JointCSsurvSIM  | Generate a data set according to the simulation under scenario I of the proposed model in Lee et al. (2022) <DOI: [xx/xx](https://doi.org/10.1002/sim.8910)>
Est.ICScure  |  Perform the semiparametric estimation method of Lee et al. (2022)  <DOI: [xx/xx](https://doi.org/10.1002/sim.8910)>

<ins>**JointCSsurvSIM**</ins>

```
JointCSsurvSIM(seed = NA, n, m, beta, alpha, kappa, sigma)
```
This function generates a data set according to the simulation under scenario I of the proposed model in Lee et al. (2022) <DOI: 10.1002/sim.8910>.

Example:
```
data <- JointCSsurvSIM(seed = 1234.2,n = 50,m = 10,beta = 1,alpha = c(1,log(2)), kappa = -0.5, sigma= 1)
head(data)

#   family        Ci delta x1         x2
# 1      1 4.0000000     0  0  0.5042357
# 2      1 2.2748313     1  0  1.8434784
# 3      1 4.0000000     0  1  0.2093028
# 4      1 0.4409429     0  0 -0.1799730
# 5      1 4.0000000     1  0  0.5579359
# 6      1 1.4709212     0  0  0.4700130

#   id cs       Lij      Rij DL DI           X          Z
# 1  1  3 1.5548184 2.194611  0  1  0.08005964 -1.2070657
# 2  1  3        NA 4.000000  0  0 -0.63140930 -1.2070657
# 3  1  3 2.3650486 3.344980  0  1 -1.51328812 -1.2070657
# 4  2  8 0.4422604 1.384112  0  1  1.84246363  0.3592891
# 5  2  8 1.7613718 2.400782  0  1  1.11236284  0.3592891
# 6  2  8 1.7778747 2.428052  0  1  0.03266396  0.3592891


# Functions
> wmcmEM(Yi, cen, X, Z, trace=FALSE, tolerance=10^{-4}) <br />
This is the estimation procedure of the Weibull Mixture Cure Model based on the EM algorithm.

> wnmcmEM(Yi, cen, X, Z, trace=FALSE, tolerance=10^{-4}) <br />
This is the estimation procedure of the Weibull Non-Mixture Cure Model based on the EM algorithm

Both functions take the arguments below:
>- Yi is a vector of the right censoring times, with size n 
>- cen is the corresponding censoring indicator with 1 being cases and 0 being censored, with size n
>- X is a covariate matrix for the incidence component with size n times the number of covariates
>- Z is a covariate matrix for the latency component with size n times the number of covariates; X and Z do not contain a column of 1 (i.e. no intercept is required); X can be completely, partially, or not different from Z.<br />
>- trace=FALSE by default. For tracking the converging path of the parameter estimation, set trace=TRUE 
>- tolerance is the converging criteria typically assigned to be 10^{-4}
