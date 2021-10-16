#' Build 100 trees
#'
#' @param data A data frame in which the rows are samples, the columns are features, and the first column is the sample class.
#' @param cp A numeric variable with a complexity parameter. The larger the value of cp, the smaller the size of the split.
#' @param prune A logical variable that is choosen between 'T' and 'F' to determine whether to prune.
#'
#' @return a list of result and a plot of decision tree.
#' @export
#'
#' @examples decision_tree(data,cp=0.01,prune=T)
decision_tree<-function(data,cp,prune){
  library(rpart)
  library(rpart.plot)

  fit1<- rpart(class~., data=data,method="class",control=rpart.control(minbucket=1,cp=cp))
  rpart.plot(fit1,type=1,branch=1,fallen.leaves = T,cex = 0.8,sub="Before pruning")
  if(prune==T){
    fit<-prune(fit1,cp=fit1$cptable[which.min(fit1$cptable[,"rel error"]),"CP"])
    rpart.plot(fit,type=1,branch=1,fallen.leaves = T,cex = 0.8,sub="After pruning")
    return(fit)
  }else{
    return(fit1)
  }

}
