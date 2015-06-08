# Script to do a ROC: age of TE versus significant colocalisation 

d=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/order_chrono_360repeats");
s=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/sig_hESC")

j=0;
labels=vector();
for(i in d$V1)  #  from older to younger TE 
{
print(i);  
j=j+1;
l=subset(s,s$V1==i)

print(dim(l));

if(dim(l)[1]>0) {labels[j]=1;}
else {labels[j]=0;}
}  

library(ROCR)
fr <- data.frame(score = (360:1),    #  the score is a king of age 
                 label =  labels )

pred <- prediction(fr$score, fr$label);
perf <- performance(pred, "tpr", "fpr");

plot(perf@x.values[[1]]*284,perf@y.values[[1]]*76,col="black",lwd=4,type='l',
     xlab="False positive rate",ylab="True positive rate (that present a significant Colocalization Score)",
     main="ROC curve for repeats ordered by Chronological order");
abline(a=0,b=76/284,lw=4,col="red");

# calculus of area :
perf <- performance(pred, "auc");
perf@y.values[[1]];

# calculus of pvalue :
# install.packages("verification")
# library("verification")

roc.area(fr$label,fr$score)

text(0.7,0.2,"AUC = 0.67",cex=2);
text(0.7,0.1,"p-value = 2.284908e-06",cex=2);


