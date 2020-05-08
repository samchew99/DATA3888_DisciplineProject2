# Kidney KNN Explorer

## Motivation 

In a K nearest neighbours model, changing the value of K can affect the model performance. This is because a small value for K will result in noise having a higher influence on the result. On the other hand, a larger K value makes model construction more computationally expensive. Furthermore, the number of predictors included in the model will also affect the model performance. Hence, this app enables users to observe the change in KNN model performance as the value for K and the number of predictors change. The KNN model is constructed using highly differential genes as explanatory variables for graft rejection outcome and evaluated using its F1 score obtained through repeated 5-fold cross validation. 


## Dataset Description

The data set used in this app was obtained from the public database GEO, which is accessible at GSE131179 [1]. The data set was from the "Landscape of innate immune system transcriptome and acute T cell-mediated rejection of human kidney allografts." study done by Mueller et al. The data set consisted of 60466 genes and 34 samples, which were collected via RNA sequencing of kidney allograft biopsy specimens from 34 adult kidney transplant recipients. Among the 34 specimens, 16 were categorized as Banff acute T cell-mediated rejection (TCMR) and 18 were normal.

## References
<ol>
<li> NCBI (26 July 2019) Series GSE131179. URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131179 </li>
</ol>
