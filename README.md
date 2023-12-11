# CMP
 Code to generate posteriors from Conway-Maxwell Poisson (CMP) distributions with non- and weakly-informative priors is in the file nicompr.R in the Code folder. A description of the method is available on arXiv, [Meyer, Graye, and Sellers (2023)](https://arxiv.org/abs/2311.18053). Descriptions for some of the key functions are below. The remaining folders contain code and data for the simulations and illustrations, additional details are below and in the manuscript.

## Functions
The primary function in the file nicompr.R is cmp().

 ### CMP model
 Fit CMP models with non- and weakly-informative priors
```
cmp(X, priors = list(a = 1, b = 1, c = 1),
    type = c('conjugate', 'flat', 'Jeffreys'),
    chains = 4, iter = 2000, ...)
```
Takes a vector $X$ containing count data. Defaults to the best performing weakly-informative prior decribed in [Meyer, Graye, and Sellers (2023)](https://arxiv.org/abs/2311.18053), the case where $a = b = c = 1$.

## Simulations folder
This folder contains files that were used to conduct the simulation study described by [Meyer, Graye, and Sellers (2023)](https://arxiv.org/abs/2311.18053).

## Illustrations folder
This folder contains a script (Illustrations.R) and datasets to analyze several datasets that exhibit both under and over dispersion including the horseshoe crab dataset (provide as .txt), the textile dataset (provide as .txt), and Hungarian and Slovak word datasets (hard-coded in Illustrations.R). See [Meyer, Graye, and Sellers (2023)](https://arxiv.org/abs/2311.18053) for details.
