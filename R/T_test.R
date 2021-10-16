#' Compute the value of t-test between case and control,and find the significant miRNAs in train set
#'
#' @param datacase A dataframe of expression profile which the samples are case-sample.The list of gene names is the first column of the dataframe.
#' @param datacontrol A dataframe of expression profile which the samples are control-sample.The list of gene names is the first column of the dataframe.
#' @param dataclass A dataframe with two column.The first column is the list of gene names. The second column is the class of genes after clustering.
#' @param p A numeric variable with the threshold value of t-test
#'
#' @return 15*count t-test results
#' @export
#'
#' @examples T_test(datacase,datacontrol,dataclass,0.05)

T_test<-function(datacase,datacontrol,dataclass,p){
  ttest <- c()
        case<-datacase
        control<-datacontrol
        result<-dataclass

        for(m in unique(result[,2])){
          X_set<-data.frame(result[which(result[,2]== m),1])
          # X_set<- data.frame(X_set[,1])
          colnames(control)[1]<-"X"
          colnames(X_set)[1]<-"X"
          colnames(case)[1]<-"X"
          X_case_input_mir<-data.frame(merge(X_set,case,all = FALSE))
          X_control_input_mir <- data.frame(merge(X_set,control,all= FALSE))

          if(length(rownames(X_case_input_mir))>1){

            for(n in 1:length(rownames(X_case_input_mir))){
              tt <- t.test(as.matrix(X_case_input_mir[n,-1]),as.matrix(X_control_input_mir[n,-1]),paired = F)
              pvalue <- tt[["p.value"]]
              tt1 <- cbind(class=m,mir_name=X_case_input_mir[1,1],pvalue)
              ttest <- rbind(ttest,tt1)
            }
          }

          if(length(rownames(X_case_input_mir))==1){
            tt <- t.test(as.matrix(X_case_input_mir[1,-1]),as.matrix(X_control_input_mir[1,-1]),paired = F)
            pvalue <- tt[["p.value"]]
            tt1 <- cbind(class=m,mir_name=X_case_input_mir[1,1],pvalue)
            ttest <- rbind(ttest,tt1)
          }
        }



  ttest<-data.frame(ttest)
  ttest_p<-ttest[ttest$pvalue<=p,]

  result_ttest<-list(ttest,ttest_p)
  names(result_ttest)<-c("t.test","t.test_p")

  return(result_ttest)
}
