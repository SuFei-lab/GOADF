#'Hierarchical clustering was performed to look for co-expression relationships between miRNAs
#'
#' @param data A dataframe or matrix of expression profile with case-sample and control-sample.The list of gene names is the rowname of the dataframe or matrix.
#' @param h A numeric variable with heights where the tree should be cut.
#' @param k A numeric variable with the number of clusters that data is divided.K can be left out.
#'
#' @return Split results of features and cluster dendrogram.
#' @export
#'
#' @examples P_cluster(data,h=0.4) or cluster(data,h=0.4,k=3)

P_cluster<-function(data,h,k){


      a<-apply(data,c(1,2),as.numeric)
      b<-t(scale(t(a)))

      if(dim(b)[1]>=3)
      {

        hr<-hclust(as.dist(1-abs(cor(t(b)))),method="complete")
        hc<-hclust(as.dist(1-abs(cor(b))),method="complete")
       if(missing(k)){
        result<- cutree(hr, h=h)
        plot(hr,hang=-1)
        return(result)
       }else{
         result<- cutree(hr, k=k,h=h)
         plot(hr,hang=-1)
         rect.hclust(hr,k=k)
         return(result)
     }
  }

}
