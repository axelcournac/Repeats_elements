# Script to do a ROC: significant colocalisation versus enrichment in TFBS

library(ROCR)

# repeats that colocalize set:
# d=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/repeats_pvalues_hesc.dat");

d=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/repeats_all_pvalues_cover.dat2");
# d=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/repeats_all_scores_cover.dat2");

#d=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/test2_GM12978_cover.dat");

# repeats associated with TF : 
s=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/TF_repeats_Bourque_S4_2013.dat");

#s=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/TF_repeats_smail.dat");
#s=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/results_repeats_smail.dat3");
#s=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/TF_repeats_ctcf.dat");

#s1=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/list_repeats_CTCF_nanog_oct4_bourque2006.dat");

# do not forget to have homogenous repeats names between different files 
s1=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/list_TF_nanog_oct4_ctcf_bourque2010.dat");

#s1=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/list_repeats_wang2014.dat")

# selection of hesc cells type exp (for Bourque exp):
s2 =   subset(s,s$V2=="H1hesc" | s$V2== "H9es" | s$V2=="Rep1H1es" | s$V2== "Rep1H7es");
#s2 =   subset(s,s$V2=="Gm12878")


# order by pvalues :
ds <- d[order(d$V2),];

j=0;
labels=vector();
for(i in ds$V1)  #  from significant repeat to less significant 
{
  #   print(i);
  j=j+1;
  #   l=subset(s,s$V1==i);
  
  l=subset(s1,s1$V1==i);
  l2=subset(s2,s2$V3==i);
  
  #   print(dim(l));
  
  if(dim(l)[1]>0 || dim(l2)[1]>0) {labels[j]=1;}
  else {labels[j]=0;}
}


fr <- data.frame(score = - ds$V2 ,    #  the score is a king of age 
                 label =  labels)

pred <- prediction(fr$score, fr$label);
perf <- performance(pred, "tpr", "fpr");

plot(perf@x.values[[1]]*1024,perf@y.values[[1]]*226,col="black",lwd=4,type='l',
     xlab="False positive rate",ylab="True positive rate (Repeat enriched with TF Binding Site)",
     main="ROC curve for repeats ordered by CS for hESC");
abline(a=0,b=226/1024,lw=4,col="red");

# Calculus of area :
perf <- performance(pred, "auc");
perf@y.values[[1]];

# calculus of pvalue :
# install.packages("verification")
library("verification")

roc.area(fr$label,fr$score);

text(800,50,"AUC = 0.5676162",cex=1,col="red");
text(800,40,"p-value = 0.0007220605",cex=1,col="red");
