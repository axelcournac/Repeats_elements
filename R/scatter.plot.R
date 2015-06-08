# script to plot scatter plots 

require('calibrate')
require('qvalue')

jpeg(file="REPEATS_COVER_SUM_hESC_1000realizations.jpeg", width = 19, height = 4, units = "cm", quality=100, res = 550 , pointsize = 8)


d=read.table("/home/axel/Bureau/REPEATS_FEV2015/Coverage_1000realizations/output.dat")

# Conversion of pvalues in qvalues: 
#qobj <- qvalue(10**d$V6)
#d$V6 <- log(qobj$qvalues)/log(10)

list<-read.table("list_hum")
# list<-read.table("list_mouse")
#list<-read.table("/media/03b8b079-9d7a-4162-8201-6dd5d9923f62/HiC_other/droso_hou/list_family_dm3.dat")

l1<-subset(list,list$V2=="LTR" | list$V2=="LTR?"  );              d1<-subset(d,d$V1 %in% l1$V1 );
l2<-subset(list,list$V2=="SINE" | list$V2=="SINE?" );             d2<-subset(d,d$V1 %in% l2$V1 );
l3<-subset(list,list$V2=="DNA" | list$V2=="DNA?");                d3<-subset(d,d$V1 %in% l3$V1 );
l4<-subset(list,list$V2=="LINE" |  list$V2=="LINE?");             d4<-subset(d,d$V1 %in% l4$V1 );
l5<-subset(list,list$V2=="Unknown" | list$V2=="Unknown?" );       d5<-subset(d,d$V1 %in% l5$V1 );
l6<-subset(list,list$V2=="Low_complexity"  );                     d6<-subset(d,d$V1 %in% l6$V1 );
l7<-subset(list,list$V2=="RC"  );                                 d7<-subset(d,d$V1 %in% l7$V1 );
l8<-subset(list,list$V2=="Satellite" );                           d8<-subset(d,d$V1 %in% l8$V1 );
l9<-subset(list,list$V2=="RNA" );                                 d9<-subset(d,d$V1 %in% l9$V1 );
l10<-subset(list,list$V2=="Simple_repeat" );                      d10<-subset(d,d$V1 %in% l10$V1 );
l11<-subset(list,list$V2=="scRNA" );                              d11<-subset(d,d$V1 %in% l11$V1 );
l12<-subset(list,list$V2=="snRNA" );                              d12<-subset(d,d$V1 %in% l12$V1 );
l13<-subset(list,list$V2=="srpRNA" );                             d13<-subset(d,d$V1 %in% l13$V1 );
l14<-subset(list,list$V2=="tRNA" );                               d14<-subset(d,d$V1 %in% l14$V1 );
l15<-subset(list,list$V2=="dRNA" );                               d15<-subset(d,d$V1 %in% l15$V1 );
l16<-subset(list,list$V2=="Other" );                              d16<-subset(d,d$V1 %in% l16$V1 );
l17<-subset(list,list$V2=="rRNA" );                               d17<-subset(d,d$V1 %in% l17$V1 );

par(mfrow = c(1, 6) )
par(mar=c(0.1,1.7,3,0.5))
par(oma=c(3,5,0.1,1.0) )

#--------------------------------------------------- 
dd= rbind(d8,d6,d10);

plot(-dd$V6, dd$V4 ,log="xy", pch=20, xlab="-Pvalue" ,ylab="CS",main="LC and Satellite",col="purple" , xaxt='n' )
rect(5,1e-10,10000, 1e-01, border = "red", col = "gray90")
points(-dd$V6, dd$V4 , pch=20 ,col="purple")

points(-d6$V6, d6$V4 ,log="xy", pch=20 ,col="red")
points(-d10$V6, d10$V4 ,log="xy", pch=20, col="purple" )

maxv=floor(max(-dd$V6))
min=min(subset(-dd$V6,-dd$V6>0))
axis(1, at=min, labels=expression(1) )
e = bquote(expression(10^.(-maxv)) )
# axis(1, at=max(-dd$V6), labels= eval(e) )
axis(1, at=5, labels=expression(bold(10^-5)) )
box()



plot(-d2$V6, d2$V4 ,log="xy", pch=20, xlab="-Pvalue" ,ylab="CS" ,
     xaxt='n' , main ="SINE",col="green" )
rect(5,1e-10,20000, 1e-01, border = "red", col = "gray90")
points(-d2$V6, d2$V4 , pch=20, xlab="-Pvalue" ,ylab="CS"   ,xaxt='n' , main ="SINE",col="green")

oo =read.table('/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/order_chrono_360repeats_25MYA')
d2_old=subset(d2,d2$V1 %in% oo$V1)

points(-d2_old$V6, d2_old$V4 , pch=20, xlab="-Pvalue" ,ylab="CS"   ,xaxt='n' ,col="green4")

maxv=floor(max(-d2$V6))
min=min(subset(-d2$V6,-d2$V6>0))
axis(1, at=min, labels=expression(1) )
e = bquote(expression(10^.(-maxv)) )
# axis(1, at=max(-d2$V6), labels= eval(e) )
axis(1, at=5, labels=expression(10^-5) )
box()


plot(-d4$V6, d4$V4 ,log="xy", pch=20, 
     xlab="-Pvalue" ,ylab="CS" ,xaxt='n', main ="LINE",col="lightblue" )
rect(5,1e-10,10000, 1e-01, border = "red", col = "gray90")

points(-d4$V6, d4$V4, pch=20, xlab="-Pvalue" ,ylab="CS" ,
       xaxt='n' , main ="LINE",col="lightblue")
maxv=floor(max(-d4$V6))
min=min(subset(-d4$V6,-d4$V6>0))
axis(1, at=min, labels=expression(1) )
e = bquote(expression(10^.(-maxv)) )
# axis(1, at=max(-d4$V6), labels= eval(e) )
axis(1, at=5, labels=expression(10^-5) )
box()

plot(-d1$V6, d1$V4  ,log="xy", pch=20, xlab="-Pvalue" ,
     ylab="CS" ,xaxt='n' , main ="LTR",col="blue" )
rect(5,1e-10,10000, 1e-01, border = "red", col = "gray90")
points(-d1$V6, d1$V4, pch=20, xlab="-Pvalue" ,ylab="CS" ,xaxt='n' , main ="LTR",col="blue")
maxv=floor(max(-d1$V6))
min=min(subset(-d1$V6,-d1$V6>0))
axis(1, at=min, labels=expression(1) )
e = bquote(expression(10^.(-maxv)) )
# axis(1, at=max(-d1$V6), labels= eval(e) )
axis(1, at=5, labels=expression(10^-5) )
box()

plot(-d3$V6, d3$V4,log="xy", pch=20, xlab="-Pvalue" ,ylab="CS" ,
     xaxt='n' , main ="DNA \n transposons",col="yellow2"   )
rect(5,1e-10,10000, 1e-01, border = "red", col = "gray90")
points(-d3$V6, d3$V4 ,log="xy", pch=20, xlab="-Pvalue" ,ylab="CS",
       xaxt='n', main ="DNA \n transposons",col="yellow2")
maxv=floor(max(-d3$V6))
min=min(subset(-d3$V6,-d3$V6>0))
axis(1, at=min, labels=expression(1) )
e = bquote(expression(10^.(-maxv)) )
# axis(1, at=max(-d3$V6), labels= eval(e) )
axis(1, at=5, labels=expression(10^-5) )
box()


dothers= rbind(d11,d12,d13,d14,d15,d16,d17);     

plot(-dothers$V6, dothers$V4,log="xy", pch=20, xlab="-Pvalue" ,
     ylab="CS" ,xaxt='n' , main ="RNA",col="orange"  )
rect(5,1e-10,10000, 1e-01, border = "red", col = "gray90")
points(-dothers$V6, dothers$V4,col="orange", pch=20)
points(-d9$V6, d9$V4, pch=20 ,col="orange")
maxv=floor(max(-dothers$V6))
min=min(subset(-dothers$V6,-dothers$V6>0))
axis(1, at=min, labels=expression(1) )
e = bquote(expression(10^.(-maxv)) )
# axis(1, at=max(-dothers$V6), labels= eval(e) )
axis(1, at=5, labels=expression(10^-5) )
box()


# mtext("Pvalue", side = 1, line = 2, at=0.00000000005)
# mtext("CS", side = 2, line = 62)


dev.off()
