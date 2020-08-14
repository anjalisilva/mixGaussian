# `MPLNClust`

## Description
`mixGaussian` is a simple R package for performing clustering using mixtures of multivariate Gaussian distributions. Main function __*mixGaussianEM*__ permit to carry out model-based clustering. Information criteria (AIC, BIC, AIC3 and ICL) are offered for model selection. Function __*mixGaussianVisualize*__ (under construction) permit to visualize clustering results. Function __*mixGaussianDataGenerator*__ (under construction) is available to generate simlulation data. The shiny implementation of *mixGaussian* is available as __*runMixGaussian*__ (under construction). For more information, see details below.  

## Installation

To install the latest version of the package:

``` r
require("devtools")
install_github("anjalisilva/mixGaussian", build_vignettes = TRUE)
library("mixGaussian")
```


## Overview

`mixGaussian` contains 8 functions. 

For carrying out clustering of data using mixtures of multivariate Gaussian distributions via expectation-maximization (EM): __*mixGaussianEM*__. 

For the purpose of generating simlulation data via mixtures of multivariate Gaussian distributions: __*mixGaussianDataGenerator*__ (under construction). 

For visualizing clustering results: __*mixGaussianVisualize*__ (under construction). 

Information criteria are offered for model selection: __*AICFunction*__, __*BICFunction*__, __*AIC3Function*__, __*ICLFunction*__. 

The shiny implementation of *mixGaussian*: __*runMixGaussian*__ (under construction). To list all functions available in the package: 

``` r
lsf.str("package:mixGaussian")
```

For tutorials and plot interpretation, refer to the vignette:

``` r
browseVignettes("mixGaussian")
```


## Details


Finite mixture models assume that a population consists of a finite collection of subpopulations and that each subpopulation can be modeled using a statistical distribution. A Gaussian mixture model assumes a Gaussian strcture for each population. Given *d*-dimensional data vectors $y_1$,..,$y_n$, the model density of a Gaussian mixture is

\begin{equation*}
\begin{split}
f(\boldsymbol{y} | \boldsymbol{\vartheta}) = \sum_{g=1}^{G} \pi_g \phi(\boldsymbol{y}|\boldsymbol{\mu}_g, \boldsymbol{\Sigma}_g),
\end{split}
\end{equation*}

where $G$ is the total number of clusters, $\phi(.)$ is the probability density function of Gaussian distribution with mean $\boldsymbol{\mu}_g$ and covariance $\boldsymbol{\Sigma}_g$, and $\pi_g$ is the mixing proportion of the gth component such that $0 < \pi_g \leq 1$ and $\sum_{g=1}^G \pi_g = 1$.

In model-based clustering applications, the cluster membership of all observations is (assumed) unknown and an indicator variable $\mathbf{Z}_i$ is used, where $i = 1, \ldots, n$. The $\mathbf{Z}_i$ are assumed to be independent and identically distributed with a multinomial distribution. Here, $\mathbf{z}_i = (z_{i1}, \ldots, z_{iG})$ is a realization of $\mathbf{Z}_i$. The ${z}_{ig}$ take on a value of $1$ if the $i$th observation belongs to component $g$ and $0$ otherwise. The predicted group memberships at the maximum likelihood estimates of the model parameters are given by the maximum \textit{a posteriori} probability, $\text{MAP}(\hat{z}_{ig})$, such that $\text{MAP}(\hat{z}_{ig}) = 1$ if $\text{max}_h(\hat{z}_{ih}), h = 1, \ldots, G$, occurs at component $g$ and $\text{MAP}(\hat{z}_{ig}) = 0$ otherwise. The complete-data likelihood for model-based clustering can be written as
\begin{equation*}
\mathcal{L}(\boldsymbol{\vartheta}) = \prod_{i=1}^{n} \prod_{g=1}^{G} \big[ \pi_g f_g( \mathbf{y} | \boldsymbol{\theta}_g)\big]^{z_{ig}}.
\end{equation*}

Here, $\boldsymbol{\theta}_g = (\boldsymbol{\mu}_g, \boldsymbol{\Sigma}_g)$ and $\boldsymbol{\vartheta}$ denotes all model parameters such that $\boldsymbol{\vartheta} = (\pi_1,\ldots,\pi_G, \boldsymbol{\theta}_1,\ldots, \boldsymbol{\theta}_G)$.
 
### EM Framework for Parameter Estimation 

Maximum likelihood estimation is one of the most widely used methods of estimation and involves selecting the values of $\hat{\boldsymbol{\vartheta}}$ which maximizes $\mathcal{L}(\boldsymbol{\vartheta}| \mathbf{y})$. Finding the maximum of $\mathcal{L}(\boldsymbol{\vartheta}| \mathbf{y})$ is an optimization problem and, often, an iterative procedure is needed to solve this problem. The expectation-maximization (EM) algorithm (Dempster et al., 1977) is an iterative procedure for the computation of maximum likelihood estimates and is commonly used in problems with missing or incomplete data. In the expectation (E-) step, the conditional expectation of complete-data log-likelihood given observed data $\mathbf{y}$ and current model parameters, termed $\mathcal{Q}$, is calculated. During the maximization (M-) step, $\mathcal{Q}$ is maximized with respect to the model parameters. The E- and M-steps are iterated until some convergence criterion is met. Note, EM algorithm relies heavily on the initial values (Biernacki et al., 2003).

### Convergence

There are several criteria for stopping an EM algorithm. A popular criterion is Aitken's acceleration (Aitken, 1926). At iteration *j*, Aitken's acceleration criterion is 

\begin{equation*}
\begin{split}
a^{(j)} & = \frac{l^{(j+1)} - l^{(j)}}{l^{(j)} - l^{(j-1)}}, 
\end{split}
\end{equation*}

where $l$ represents the posterior log-likelihoods. Convergence is achieved when 
$$|l_{\infty}^{(j+1)} - l_{\infty}^{(j)}| < \epsilon.$$ Here, $l_{\infty}^{(j+1)}$ is an asymptotic estimate of the log-likelihood (B¨ohning et al., 1994) given by
$$l_{\infty}^{(j+1)} = l^{(m)} + \frac{l^{(j+1)} - l^{(j-1)}}{1 - a^{(j)}}.$$

### Model Selection 

Model selection criteria penalizes the log-likelihood for model complexity, as the log-likelihood value favours the model with more parameters. Among the different criteria available, the Bayesian information criterion \citep[BIC;][]{schwarz1978} remains the most popular criterion for model-based clustering applications \citep{mcnicholas2016}. The BIC is defined as $$\text{BIC} = -2 \log \mathcal{L} (\hat{\boldsymbol{\vartheta}} |\boldsymbol{y}) + K \log n.$$ The $\mathcal{L} (\hat{\boldsymbol{\vartheta}} |\mathbf{y})$ represents maximized log-likelihood, $\hat{\boldsymbol{\vartheta}}$ is the maximum likelihood estimate of the model parameters $\boldsymbol{\vartheta}$, $n$ is the number of observations, and $K$ represents the number of free parameters in the model. There is support for the use of BIC in mixture model selection \citep{Campbell1997, fraley1998, Dasgupta1998} and also support for the use of BIC in selecting the number of latent factors in a factor analysis model \citep{lopes2004}. Other criteria include the Akaike information criterion \citep[AIC;][]{akaike1973}, $$\text{AIC} = -2 \log \mathcal{L} (\hat{\boldsymbol{\vartheta}} |\boldsymbol{y}) + 2K;$$ a variation on the AIC used by \citet{Bozdogan1994}, $$\text{AIC3} =  -2 \log \mathcal{L} (\hat{\boldsymbol{\vartheta}}|\boldsymbol{y}) + 3K;$$ and the integrated completed likelihood \citep[ICL;][]{biernacki2000}, $$\text{ICL} \approx \text{BIC} - 2 \sum_{i=1}^n \sum_{g=1}^G \text{MAP}(\hat{z}_{ig}) \log \hat{z}_{ig}.$$ Here, $\text{MAP}(\hat{z}_{ig})$ is the maximum $\textit{a posteriori}$ classification given $\hat{z}_{ig}$.

As the number of clusters in the finite mixture model is almost always unknown, the parameter estimation methods are fitted for a range of possible number of components and the optimal number is selected using a model selection criterion. Varying model selection criteria differ in terms of the log-likelihood is penalized. The AIC and AIC3 penalize log-likelihood only for the number of free parameters in the model and are constant with respect to the sample size. When the number of observations is large, the AIC tends to favor more complex models \citep{shibata1976, katz1981} and it overestimates number of components in the model.  BIC penalizes the log-likelihood based on both the number of  free parameters in the model and the number of observations. The BIC selects the number of mixture components needed to provide a good approximation to the density, rather than the number of clusters. As a result, BIC can assign multiple mixture components to a single cluster. The ICL penalizes BIC for the estimated mean entropy based on the spread of the mixture components. The ICL favors well separated clusters compared to BIC \citep{biernacki2000}.

## `mixGaussian` Algorithm

In `mixGaussian` the parameter and group membership estimation is carried out using the EM algorithm, because the complete-data consists of the unobserved group membership labels. To check the convergence of EM algorithm, a modified version of Aitken's acceleration criterion as outlined by B¨ohning et al., 1994 is used with $\epsilon = 0.001$. Model selection is performed using AIC, BIC, AIC3 and ICL.



## Citation for Package
``` r
citation("mixGaussian")
```

## References for Package

* [Aitken, A. C. (1926). A series formula for the roots of algebraic and transcendental equations. *Proceedings of the Royal Society of Edinburgh*](https://www.cambridge.org/core/journals/proceedings-of-the-royal-society-of-edinburgh/article/iiia-series-formula-for-the-roots-of-algebraic-and-transcendental-equations/0CC96A97C8B634E2730F5208E506E6A9)

* [B¨ohning, D., E. Dietz, R. Schaub, P. Schlattmann, and B. Lindsay (1994). The distribution of the likelihood ratio for mixtures of densities from the one-parameter exponential family. *Annals of the Institute of Statistical Mathematics*](https://link.springer.com/article/10.1007/BF01720593)

* [Dempster, A. P., N. M. Laird, and D. B. Rubin (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society: Series B*](https://www.ece.iastate.edu/~namrata/EE527_Spring08/Dempster77.pdf)

* For others, refer to help page of inidividual functions via `?function` or `help(function)`.


## Maintainer

* Anjali Silva (anjali.silva@uhnresearch.ca). 


## Contributions

`mixGaussian` welcomes issues, enhancement requests, and other contributions. To submit an issue, use the [GitHub issues](https://github.com/anjalisilva/mixGaussian/issues).

