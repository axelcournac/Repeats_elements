# To make a box plot to compare repeats enriched with CTCF, NANOG, OCT4 in hESC and IMR90 cells
# December 2014:

# Repeats enriched with CTCF, NANOG or OCT4 binding sites, (Bourque Lab Results):
# d=read.table("/home/axel/Bureau/NAR/repeat_roc/list_repeats_CTCF_nanog_oct4_bourque2006.dat");
d=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/repeat_roc/list_TF_nanog_oct4_ctcf_bourque2010.dat");

# CS data (ratio of CS between hESC and IMR90) : 
#cs=read.table("/media/axel/RSG3/data/REPEATS/OUT2/OUTPUT_COVER_SUM_hESC/test3.dat");
#cs=read.table("/run/media/axel/RSG3/BACK_UP/data/REPEATS/OUT2/OUTPUT_COVER_SUM_hESC/test3.dat");
#cs=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/figure_boxplot/heasc_imr90_cover.dat");
cs=read.table("/home/axel/Bureau/REPEATS_FEV2015/NAR/figure_boxplot/hesc_imr90_cover.dat2");

cs =subset(cs,cs$V4>0 & cs$V10 >0)
cs=subset(cs,cs$V2>500 & cs$V8>500)

# which repeats are enriched in TFBS
ctcf = subset(d,d$V2=="CTCF");
nanog = subset(d,d$V2=="NANOG");
oct4 = subset(d,d$V2=="OCT4");

cs_ctcf=subset(cs,cs$V1 %in% ctcf$V1);
cs_nanog=subset(cs,cs$V1 %in% nanog$V1);
cs_oct4=subset(cs,cs$V1 %in%  oct4$V1);
cs_nanog_oct4=subset(cs,cs$V1 %in% oct4$V1 & cs$V1 %in% nanog$V1);

#------------------------------------------------------- 

boxplot(log(cs$V4/cs$V10),log(cs_ctcf$V4/cs_ctcf$V10),log(cs_nanog$V4/cs_nanog$V10),log(cs_oct4$V4/cs_oct4$V10),
        log(cs_nanog_oct4$V4/cs_nanog_oct4$V10),
        names=c("All","CTCF","Nanog","Oct4","Nanog and Oct4") ,
        ylim= c(-0.2,0.2) , 
        ylab ="Ratio of CS between hESC and IMR90", 
        xlab = "Group of repeat containing binding sites for TF",
        col='cyan')

# Probality tests:

tes <- wilcox.test(cs$V4/cs$V10,cs_oct4$V4/cs_oct4$V10, alternative = 'less')
tes

tes <- wilcox.test(cs$V4/cs$V10,cs_nanog_oct4$V4/cs_nanog_oct4$V10, alternative = 'less')
tes


