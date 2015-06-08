library('fitdistrplus')

Args <- commandArgs();


png(file = paste(Args[4],".png" ) )

d=read.table(Args[4])
r=read.table(Args[5])
#hist(r$V1/r$V2,100,col="grey",xlim=c(min(r$V1/r$V2),max(r$V1/r$V2)*1.5 ),main=Args[4] ,xlab="Colocalisation Score")
#points(d$V1/d$V2,30,t="h",col="red",lw=3)

x=r$V1/r$V2;
fitx <- fitdist(x, "lnorm")
p=plnorm(d$V1/d$V2,fitx[[1]][1],fitx[[1]][2],lower.tail=FALSE,log =TRUE)
p=p/log(10);


# if(d$V1/d$V2 > max(r$V1/r$V2) )  
  cat(as.character(d$V4),d$V5,d$V7,d$V1/d$V2, max(r$V1/r$V2),p, "\n", file="test2.dat",append=T)


# graph : 

#hist(r$V1/r$V2,50,col="grey",xlim=c(min(r$V1/r$V2),max(r$V1/r$V2)*1.5),xlab="Colocalisation Score",density=TRUE)

minVal<- min(r$V1/r$V2)
maxVal<- max(r$V1/r$V2*1.5)
numItv<-  100;
v <- seq(minVal, maxVal, length.out=numItv+1);


h =hist(r$V1/r$V2,50,col="grey",xlim=c(min(r$V1/r$V2),max(r$V1/r$V2)*1.5 ),main=Args[4] ,xlab="Colocalisation Score");

#plot(h$breaks[1:(length(h$breaks)-1)],h$counts, log="xy", type='h', lwd=10, lend=2,
 #    col="grey",xlim=c(min(r$V1/r$V2),max(r$V1/r$V2)*1.5 ),main=Args[4] ,xlab="Colocalisation Score")


v2=dlnorm(v,fitx[[1]][1],fitx[[1]][2])

lines(v, v2/max(v2)*max(h$counts)  ,type="l" ,col="deepskyblue4" , lw=5);
points(d$V1/d$V2,0,pch=19,col="red",lw=6);
