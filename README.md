# HCfused: multi-view hierarchical clustering in R 

![HCfusedLogo](https://github.com/pievos101/HC-fused/blob/master/HCfused.png)


## Installation

```{r}
install.packages("devtools")
library(devtools)

install_github("pievos101/HC-fused")
library(HCfused)
```

## Basic usage

Loading example data views

```{r}
data(view1)
data(view2)
```

Loading the outcome vector

```{r}
data(target)
```

Multi-view clustering using HCfused.
Let's cluster the views using the ward.D method.

```{r}
k   = length(unique(target))
res = HCmv(list(view1, view2), k=k, method="ward.D")
```

The fusion cluster solution can be obtained from 

```{r}
cl = res$cluster
```

Let's check the performance based on the Adjusted R Index (ARI)

```{r}
require(aricode)
ARI(cl, target)
NMI(cl, target)
```

The fused affinity matrix can be accesed via

```{r}
affinityMatrix = res$P
```

which can be clustered by any clustering algorithm

```{r}
distanceMatrix = 1 - affinityMatrix/max(affinityMatrix)
fused = hclust(as.dist(distanceMatrix), method="average")
```

Let's check the performance based on the Adjusted R Index (ARI)

```{r}
cl = cutree(fused, k=k)
require(aricode)
ARI(cl, target)
NMI(cl, target)
```

## The fusion algorithm HCfuse

You may want to use your own clustering algorithm and just employ the hierarchical fusion algorithm.

For instance, lets assume we have two cluster solutions cl1 and cl2.

```{r}
cl1 = c(1,1,1,2,2,2,3,3,3)
cl2 = c(1,1,2,2,2,3,3,3,3)
```

Now, we need to create the co-association matrices

```{r}
ass1 = association(cl1)
ass2 = association(cl2)
```

These two binary matrices can now be fused

```{r}
res = HCfuse(list(ass1, ass2))
affinityMatrix = res$P
```

## Hierarchical ensemble clustering

... will be available soon ...

## References
Please cie the following work in case you find the package useful.


Pfeifer, Bastian, and Michael G. Schimek. "A hierarchical clustering and data fusion approach for disease subtype discovery." Journal of Biomedical Informatics 113 (2021): 103636.

https://www.sciencedirect.com/science/article/pii/S1532046420302641

Please see also

Pfeifer, Bastian, et al. "Integrative hierarchical ensemble clustering for improved disease subtype discovery." 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM). IEEE, 2021.

https://ieeexplore.ieee.org/abstract/document/9669608

for integrative hierarchical ensemble clustering with HC-fused.
