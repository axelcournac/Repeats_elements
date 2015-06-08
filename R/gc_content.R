#   R script to calculate GC content of bins

library(BSgenome.Hsapiens.UCSC.hg19)

d<-read.table("/bounds_10_corrected.dat2")

chr23<-Hsapiens$chrX
chr24<-Hsapiens$chrY

gc=c(1:dim(d)[1])

for (i in 1:dim(d)[1]) 
{
chro<-d[i,1];
st  <-d[i,2];
ed  <-d[i,3];

if(chro=="chr23") {chro="chrX";}
if(chro=="chr24") {chro="chrY";} 

if(st<ed)
{
s<-getSeq(Hsapiens,chro,start=st,end=ed,as.character=FALSE);
alp<-alphabetFrequency(s);
gc[i]=  (alp['C']+alp['G'])  /  ( alp['A']+alp['T']+alp['C']+alp['G'] );

print(cat(i,gc[i]));
}

else
{
gc[i]= 0;
}

}

write.table(gc,file="gc_content_bornes_10.dat",row.names = FALSE,col.names  = FALSE)
