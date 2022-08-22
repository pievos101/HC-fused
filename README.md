# HCfused: multi-view hierarchical clustering in R 

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

Multi-view clustering using HCfused

```{r}
res = HCmv(list(view1, view2))
```

The fused matrix can be accesed via

```{r}
res$P
```

which can be clustered again

```{r}
fused = hclust(as.dist(res$P), method="ward.D")
```

Let's check the performance based on the Adjusted R Index (ARI)

```{r}
cl = cutree(fused, k=length(unique(target)))
require(aricode)
ARI(cl, target)
NMI(cl, target)
```


## References

Pfeifer, Bastian, and Michael G. Schimek. "A hierarchical clustering and data fusion approach for disease subtype discovery." Journal of Biomedical Informatics 113 (2021): 103636.

https://www.sciencedirect.com/science/article/pii/S1532046420302641

Please see also

Pfeifer, Bastian, et al. "Integrative hierarchical ensemble clustering for improved disease subtype discovery." 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM). IEEE, 2021.

https://ieeexplore.ieee.org/abstract/document/9669608

for integrative hierarchical ensemble clustering with HC-fused.
