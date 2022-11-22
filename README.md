# Multi-view hierarchical clustering in R
<p align="center">
<img src="https://github.com/pievos101/HC-fused/blob/master/logo.png" width="500">
</p>

## Installation

```{r}
install.packages("devtools")
library(devtools)

install_github("pievos101/HC-fused")
library(HCfused)
```
## Approach
![HCfusedLogo2](https://github.com/pievos101/HC-fused/blob/master/HCfused.png)


## Basic usage

Loading some data views

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

The fused cluster solution can be obtained from 

```{r}
cl = res$cluster
```
Note, when no k (number of clusters) is set, the optimal number 
of clusters is inferred by the silhouette coefficient.

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
distanceMatrix = 1 - affinityMatrix
fused = hclust(as.dist(distanceMatrix), method="average")
cl = cutree(fused, k=k)
```

Let's check the performance based on the Adjusted R Index (ARI)

```{r}
require(aricode)
ARI(cl, target)
NMI(cl, target)
```

## The fusion algorithm HCfuse

You may want to use your own clustering algorithm and employ laste fusion using the hierarchical fusion algorithm HCfuse.

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

These two binary matrices can be fused.

```{r}
res = HCfuse(list(ass1, ass2))
affinityMatrix = res$NETWORK
```

The resulting affinity matrix can then be clustered by any clustering algorithm.

## Parea: multi-view hierarchical ensemble clustering

![Parea1Logo](https://github.com/pievos101/HC-fused/blob/master/Parea1.png)


```{r}
require(GA)
```

First, we need to define the fitness function for the genetic algorithm.

```{r}
# Fitness function for genetic algorithm
check_ensemble <- function(x, methods=FALSE, omics_in=FALSE, fix.k=NaN){
	ens  <- round(x)
	ens  <- methods[ens]
	res  <- Parea(omics=omics_in, this_method=ens, fix.k=fix.k, type=1)
	return(res$SIL) # silhouette cofficient
} # end of fitness function
```

The following methods are available.

```{r}
# Available hierarchical clustering methods
methods = c("single", "complete", "average", "mcquitty", "ward.D",
"ward.D2", "centroid", "median")
```

Starting the genetic algorithm

```{r}
fix.k = 5

# Perform the genetic algorithm
res <- ga(
	type = "real-valued", 
	fitness = check_ensemble, methods, list(view1, view2), fix.k, lower = c(1,1), upper = c(8,8),  
	elitism = 20, maxiter = 20, popSize = 20, 
	run = 20, parallel=FALSE)
```

The inferred methods are within the 'solution' slot

```{r}
sel       <- methods[round(res@solution)]
```

Now, we can cluster the data with the inferred methods

```{r}
res       <-  Parea(list(view1, view2), 
				fix.k=k,
				this_method=sel,
				HC.iter=30, 
				type=1)

cl_ensemble  <- res$cluster
```

Let's check the performance based on ARI and NMI.

```{r}
require(aricode)
ARI(cl_ensemble, target)
NMI(cl_ensemble, target)
```


# Miscellaneous

Logo made by Adobe Express Logo Maker: <https://www.adobe.com/express/create/logo>


## References
Please cite the following work in case you find the package useful.


Pfeifer, Bastian, and Michael G. Schimek. "A hierarchical clustering and data fusion approach for disease subtype discovery." Journal of Biomedical Informatics 113 (2021): 103636.

https://www.sciencedirect.com/science/article/pii/S1532046420302641

Please see also

Pfeifer, Bastian, et al. "Integrative hierarchical ensemble clustering for improved disease subtype discovery." 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM). IEEE, 2021.

https://ieeexplore.ieee.org/abstract/document/9669608

for integrative hierarchical ensemble clustering with HC-fused.
