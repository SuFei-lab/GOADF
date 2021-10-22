#' Compute the value of sam between case and control,and find the significant miRNAs in train set
#'
#' @param datacase A dataframe of expression profile which the samples are case-sample.The list of gene names is the first column of the dataframe.
#' @param datacontrol A dataframe of expression profile which the samples are control-sample.The list of gene names is the first column of the dataframe.
#' @param delta A numeric variable that adjusts the false discovery rate(FDR)
#' @param p A numeric variable that controls the threshold value of differentially expression of genes.
#'
#' @return A dataframe of feature selection set which is consisted of genes with p value less than or equal to p and the expression value of them.
#' @export
#'
#' @examples SAM(datacase,datacontrol,p=0.05)
SAM<-function(datacase,datacontrol,p){
  library(siggenes)
     gene_name<-data.frame(datacase[,1])
    colnames(gene_name)[1]<-"mir_name"
    case<-data.frame(datacase[,-1])
    control<-data.frame(datacontrol[,-1])
    Feature_selection_set<-as.matrix(cbind(case,control))
    cl<-c(rep(1,dim(case)[2]),rep(0,dim(control)[2]))
    sam.out<-siggenes::sam(as.matrix(Feature_selection_set),cl,gene.names = gene_name[,1])
    p_value<-data.frame(mirname=gene_name,p.value=sam.out@p.value)
    p_value<-p_value[p_value$p.value<=0.05,]
  
    Feature_selection_set<-cbind.data.frame(gene_name,Feature_selection_set)
    select<-Feature_selection_set[match(p_value$mir_name,Feature_selection_set[,1]),]
    select<-select[complete.cases(select),]

      return(select)


}
