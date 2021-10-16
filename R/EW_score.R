#' Score each decision tree
#'
#' @param data  Data is a dataframe with at least 11 columns including names of miRNAs, k, i, frequency of miRNAs, the value of single ROC, combined ROC, prerat, sensitivity, specificity, precision and F_measure of evaluate set and test set.K and i are the ordinal number of the feature selection set.
#'
#'
#' @return The value of score of every decision tree.
#' @export
#'
#' @examples EW_score(data)
#'
EW_score<-function(data){




 dafenn<-data.frame(data[1:nrow(data),2:ncol(data)])


dafen1<-as.numeric(as.character(data$frequency))
dafen3<-as.numeric(as.character(data$roc_cob_value))
dafen4<-as.numeric(as.character(data$prerat_conb))
dafen5<-as.numeric(as.character(data$senrat_conb))
dafenn_1<-data.frame(dafen1,dafen3,dafen4,dafen5)
sum<-apply(dafenn_1,2,sum)
dafen11<-dafen1/sum[1]
dafen13<-dafen3/sum[2]
dafen14<-dafen4/sum[3]
dafen15<-dafen5/sum[4]
dafen_2<-data.frame(dafen11,dafen13,dafen14,dafen15)

log<-log2(dafen_2)
plogp<-dafen_2*log
plogp[which(plogp[,3]=='NaN'),3]=0
plogp[which(plogp[,4]=='NaN'),4]=0
sum_plogp<-apply(plogp,2,sum)
E<--sum_plogp


sumE<-E[1]+E[2]+E[3]+E[4]
w=c()
for(i in 1:4){
  w[i]<-(1-E[i])/(4-sumE)
}


result1<-0.001*dafen1*w[1]
result3<-dafen3*w[2]
result4<-dafen4*w[3]
result5<-dafen5*w[4]
result<-data.frame(result1,result3,result4,result5)
Z<-apply(result,1,sum)

score<-Z
mir_name<-data$mir_name
dafenn2<-cbind.data.frame(mir_name,dafenn,score)

t<-c()

b<-c()


result_fenzu<-c()
mir_name2<-c()
for(p in 1:20){
  data_1<-dafenn2[dafenn2$ids==p,]
  for(j in 1:100){
    data_2<-data_1[data_1$kds==j,]
    a0<-length(data_2[,1])
    a<-as.numeric(a0)
    b<-apply(data_2[,c(-1)],2,mean)
    t<-rbind(t,b)
    mir_name<-paste(data_2$mir_name)
    mir_name1<-paste(mir_name[1:length(mir_name)],collapse=";")
    mir_name2<-rbind(mir_name2,mir_name1)
  }
}
result_fenzu<-data.frame(cbind(mir_name2,t))
result_fenzu <- result_fenzu[-which(result_fenzu[,1]=='NA'),]
colnames(result_fenzu)[1]<-c("miRNA_Name")
rownames(result_fenzu) <- NULL

return(result_fenzu)


}
