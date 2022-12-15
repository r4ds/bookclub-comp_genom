expFile=system.file("extdata",
                    "leukemiaExpressionSubset.rds",
                    package="compGenomRData")
mat=readRDS(expFile)
library(pheatmap)
library(cluster)
library(fastICA)

#1
x <- scale(mat)
y <- scale(log2(mat))

boxplot(mat)
boxplot(x)
boxplot(y)

#2
annotation_col = data.frame(LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)

pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D",
         clustering_distance_cols="euclidean")

### IMO this one looks the bet to me
pheatmap(x,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D",
         clustering_distance_cols="euclidean")

pheatmap(y,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D",
         clustering_distance_cols="euclidean")

pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "column",clustering_method="ward.D",
         clustering_distance_cols="euclidean")

#3 number of clusters
set.seed(101)
pamclu=cluster::pam(t(mat),k=5)
plot(silhouette(pamclu),main=NULL)

##even when mat is changed with scaled data it looks the same
Ks=sapply(2:7,
          function(i)
            summary(silhouette(pam(t(mat),k=i)))$avg.width)
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)

#4
library(cluster)
set.seed(101)
# define the clustering function
pam1 <- function(x,k)
  list(cluster = pam(x,k, cluster.only=TRUE))

# calculate the gap statistic
pam.gap= clusGap(t(mat), FUN = pam1, K.max = 8,B=50)

# plot the gap statistic accross k values
plot(pam.gap, main = "Gap statistic for the 'Leukemia' data")

#5
## data all together doesn't tell you anything useful
plot(mat, pch = 19, col=as.factor(annotation_col$LeukemiaType))
princomp(x)
screeplot(princomp(x))




# create the subset of the data with two genes only. notice that we transpose the matrix so samples are on the columns
par(mfrow=c(1,2))
sub.mat=t(mat[rownames(mat) %in% c("ENSG00000100504","ENSG00000105383"),])

# ploting our genes of interest as scatter plots
plot(scale(mat[rownames(mat)=="ENSG00000100504",]),
     scale(mat[rownames(mat)=="ENSG00000105383",]),
     pch=19,
     ylab="CD33 (ENSG00000105383)",
     xlab="PYGL (ENSG00000100504)",
     col=as.factor(annotation_col$LeukemiaType),
     xlim=c(-2,2),ylim=c(-2,2))

# create the legend for the Leukemia types
legend("bottomright",
       legend=unique(annotation_col$LeukemiaType),
       fill =palette("default"),
       border=NA,box.col=NA)

# calculate the PCA only for our genes and all the samples
pr=princomp(scale(sub.mat))

#screeplot
screeplot(pr)

# plot the direction of eigenvectors
# pr$loadings returned by princomp has the eigenvectors
arrows(x0=0, y0=0, x1 = pr$loadings[1,1],
       y1 = pr$loadings[2,1],col="pink",lwd=3)
arrows(x0=0, y0=0, x1 = pr$loadings[1,2],
       y1 = pr$loadings[2,2],col="gray",lwd=3)


# plot the samples in the new coordinate system
plot(-pr$scores,pch=19,
     col=as.factor(annotation_col$LeukemiaType),
     ylim=c(-2,2),xlim=c(-4,4))

# plot the new coordinate basis vectors
arrows(x0=0, y0=0, x1 =-2,
       y1 = 0,col="pink",lwd=3)
arrows(x0=0, y0=0, x1 = 0,
       y1 = -1,col="gray",lwd=3)

#6
par(mfrow=c(1,2))
d=svd(scale(mat)) # apply SVD
assays=t(d$u) %*% scale(mat) # projection on eigenassays
plot(assays[1,],assays[2,],pch=19,
     col=as.factor(annotation_col$LeukemiaType))
#plot(d$v[,1],d$v[,2],pch=19, col=annotation_col$LeukemiaType)

pr=prcomp(t(mat),center=TRUE,scale=TRUE) # apply PCA on transposed matrix

# plot new coordinates from PCA, projections on eigenvectorssince the matrix is transposed eigenvectors represent
plot(pr$x[,1],pr$x[,2],col=as.factor(annotation_col$LeukemiaType))

#7
ica.res=fastICA(t(mat),n.comp=2) # apply ICA
plot(ica.res$S[,1],ica.res$S[,2],col=as.factor(annotation_col$LeukemiaType)) # plot reduced dimensions

#8
library("Rtsne")
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(t(mat),perplexity = 10) # Run TSNE
                                          # image(t(as.matrix(dist(tsne_out$Y))))
                                          # Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
legend("bottomleft", # create the legend for the Leukemia types
       legend=unique(annotation_col$LeukemiaType),
       fill =palette("default"),
       border=NA,box.col=NA)
