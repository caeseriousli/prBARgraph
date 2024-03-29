## Package Landing Page

### Introduction

We have developed a new L<sub>0</sub> regularized sparse Poisson graphical model with applications to gene network inference from RNA-seq gene expression count data. Assuming a pairwise Markov property, we propose to fit a separate broken adaptive ridge (BAR)  regularized log-linear Poisson regression on each node to evaluate the conditional, instead of marginal, association between two genes in the presence of all other genes. The resulting sparse gene networks are generally more accurate than those generated by the $L_1$ regularized Poisson graphical model as demonstrated by our empirical studies. A real data illustration is given on a kidney renal clear cell carcinoma miRNA-seq data from TCGA. 

---

### Insatall Necessary Packages

To install this package, please use R developer package `library(devtools)` to source and install this Github hosted package.

---

Install devtools, [Package Link](https://www.r-project.org/nosvn/pandoc/devtools.html) 

```markdown

# Install devtools in R

`install.packages("devtools")`

```

For more details see [r-lib/devtools](https://github.com/r-lib/devtools).

---

Install prBARgraph

```{r}
library(devtools)
install_github("caeseriousli/prBARgraph")
```

--- 

To simulate, you will need to install `XMRF` package ($L_1$ regularized Poisson graphical model implementation).


```{r}
require(devtools)
install_github("cran/XMRF")

## If CRAN has dropped support for the package, try 
install_github("zhandong/XMRF")

library(XMRF)
```

--- 

### Simulate a Poisson dataset

Use XMRF Package to simulate a Poisson network data set. Let's try n=100 and p=30, with a [scale-free topology](https://en.wikipedia.org/wiki/Scale-free_network). Simulation process involves manipulations of permutation and adjacency matrices. Citation,

Allen, G. I., & Liu, Z. (2012, October). A log-linear graphical model for inferring genetic networks from high-throughput sequencing data. In 2012 IEEE [International Conference on Bioinformatics and Biomedicine (pp. 1-6). IEEE](https://ieeexplore.ieee.org/abstract/document/6392619).


```{r}
## Randomly generate a Poisson (potentially an RNA-seq) dataset.
Xsim <- XMRF.Sim(n=100, p=30, model="LPGM", graph.type="scale-free")
```

---

### Fit and plot L<sub>0</sub> regularized models

```{r}
## Generate a lambda (penalization strength) path, from 18 to 1e-4. The larger the lambda, 
## the stronger the penalization, and the sparse network would be.

my_sequence = 2^seq(logb(18, base = 2), logb(1e-4, base = 2), length = 10)
fit = fitModel(t(Xsim$X), results, sequen, beta=c(1,2), regularization = 'l0', lchosen = my_sequence)
```

Since we gave 10 lambdas, the fitModel function returns a list of 10 [adjacency matrix](https://en.wikipedia.org/wiki/Adjacency_matrix). This adjacency matrix specifies the structure of a network. Then we move on to plot the networks. As we can see, as we increase the lambda (closer to the 1e-4 end the lambda path, the network becomes denser).

```{r}
plot(outputNetwork(fit[[1]]))

plot(outputNetwork(fit[[4]]))

plot(outputNetwork(fit[[7]]))

```
<img src="inst/network1.png?raw=true"/>

---

<img src="inst/network4.png?raw=true"/>

---

<img src="inst/network7.png?raw=true"/>


