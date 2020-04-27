# Article:
A. Khan and P. Maji, "Approximate Graph Laplacians for Multimodal Data Clustering," in *IEEE Transactions on Pattern Analysis and Machine Intelligence*, pp. 1--16, 2019.
doi: 10.1109/TPAMI.2019.2945574

URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8859233&isnumber=4359286


Inorder to execute the R code for the **Lower Grade Glioma (LGG) data set **(a type of brain cancer),  within the R environment execute:
>source("LGGdemo.R")



Joint subspace or Eigenvectors of the Joint Laplacian  is written to file : *LrStar_Eigenvectors.txt*
LrStar_Eigenvectors.txt contains a *(n x r)* matrix.
Here *n* is the number of samples in the data set and *r* is the optimal/required rank of the joint subspace.

The eigenvalues of the Joint Laplacian is written to file : *LrStar_Eigenvalues.txt*

*k*-means clustering can be performed on the rows of LrStar_Eigenvectors matrix to get the clusters in the data set. The cluster assignments are written to the file* LGG-ClusterAssignment.txt* for the LGG data set	.

The file *ConvexLaplacian.R* contains the R implementation of the CoALa algorithm as a function `ConvexLaplacian`. 
Details of the fuctions is as follows:

Function Name: `ConvexLaplacian`

###### #Usage 
`ConvexLaplacian(Data,K,rank=NULL,alpha=NULL,beta=NULL,simFromFile=FALSE,mod=NULL)`


Arguments

*Data*:  A list object containing *M* data matrices representing *M* different omic data types measured in a set of *n* samples. 
For each matrix, the rows represent samples, and the columns represent genomic features.
The matrices in the list can have variable numbe of columns(features), but all must have the same number of *n* rows(samples).

*K*: The number of clusters in the data set.

*rank*: The rank of the individual and joint Laplacian. 
Default value: NULL.
if *rank=NULL*, the algorithm varies the rank between ceiling(K/M) to 50 and selects the optimal rank of the joint subspace.
M is the number of modalities in the data set.

*alpha*: Weight of the convex combination in joint Laplacian.
Default value: NULL.
if *alpha=NULL*, the algorithm computes alpha as proposed in the paper.

*beta*: Value of Damping factor during weight selection.
Default value: NULL.
if *beta=NULL*, value of beta is set to *1.25* as in paper.

*simFromFile*: Boolean value (TRUE/FALSE)
if FALSE, algorithm considers the matrices in the 'Data' list as feature-based representation, i.e., as *(n x d_m)* data matrices,
and computes the graph Laplacian from data matrices.
if TRUE, algorithm considers the matrices in the 'Data' list as graph-based representation, i.e., as *(n x n)* similarity matrices,
and computes the graph Laplacian from similarity matrices.
Default value: FALSE.

*mod*: Array containing names of modalities
Default value: NULL.
if mod=NULL, the algorithm names the modalities as* 1,2,...M*.




# Example call:

```r
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/BRCA/mDNA",n),sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/BRCA/RNA",n),sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/BRCA/miRNA",n),sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/BRCA/RPPA",n),sep=" ",header=TRUE,row.names=1))
K=4

#Log Transform of Sequence based Gene and miRNA modality
LogData=Data
LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)


source("ConvexLaplacian.R")
ConvexLaplacian(Data=LogData,K=K)
```
