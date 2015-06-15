#include <stdio.h>
#define MATHLIB_STANDALONE 1
#include <Rmath.h>
#include <stdlib.h>
#include <string.h>

// November 2013 :
// To keep bins with a significant number of occurences according to binomial distribution (given by pbinom of R) for the human genome. 
// The library libRmath must be installed to use R functions in C. 
// To compile : gcc binnage_distrib.c -lRmath  -lm 

int main(int argc, char* argv[])
{
int NTF=837603;  // nbre total de frag 
int NBINS;
FILE* fic1 = NULL;FILE* fic2 = NULL;
FILE* fic3 = NULL;FILE* fic5 = NULL;FILE* fic6 = NULL;FILE* fic7 = NULL;
int chr1,chr2;
int locus1,locus2;
int sens1,sens2;
int mq1,mq2;
int match1,match2;
int indice1,indice2;
int d1,d2;
int p;
int left,right;
int j=0,b=0;
int NPOS_repeat=0;
int NBIN_repeat=0;
int NSEUIL=0;
unsigned int i=0;
int c,c_a=0;
int label;
double pval;
int np=0,sum=0;

int *tab_left;tab_left = calloc (837604 , sizeof(int));
int *tab_right;tab_right= calloc (837604 , sizeof(int));
int *IND_REPEAT;IND_REPEAT = calloc (837604 , sizeof(int));
int *CHR_REPEAT;CHR_REPEAT= calloc (837604 , sizeof(int));
int *POS_REPEAT;POS_REPEAT= calloc (837604 , sizeof(int));
int *IND_RANDOM;IND_RANDOM= calloc (837604 , sizeof(int));
int *chro_ind;chro_ind = calloc (837604 , sizeof(int));

char type_repeat[250];
char type_repeat_avant[250];

strcpy(type_repeat, "INIT");strcpy(type_repeat_avant, "INIT");
unsigned int chr1_repeat,position_repeat,position_repeat1,position_repeat2,position_repeat_avant,indice_repeat;

char cmd[300];
int seuil;

int** tab_pos_indice = malloc(24*sizeof(int*));

tab_pos_indice[0]=malloc(249250621*sizeof(int));
tab_pos_indice[1]=malloc(243199373*sizeof(int));
tab_pos_indice[2]=malloc(198022430*sizeof(int));
tab_pos_indice[3]=malloc(191154276*sizeof(int));
tab_pos_indice[4]=malloc(180915260*sizeof(int));
tab_pos_indice[5]=malloc(171115067*sizeof(int));
tab_pos_indice[6]=malloc(159138663*sizeof(int));
tab_pos_indice[7]=malloc(146364022*sizeof(int));
tab_pos_indice[8]=malloc(141213431*sizeof(int));
tab_pos_indice[9]=malloc(135534747*sizeof(int));
tab_pos_indice[10]=malloc(135006516*sizeof(int));
tab_pos_indice[11]=malloc(133851895*sizeof(int));
tab_pos_indice[12]=malloc(115169878*sizeof(int));
tab_pos_indice[13]=malloc(107349540*sizeof(int));
tab_pos_indice[14]=malloc(102531392*sizeof(int));
tab_pos_indice[15]=malloc(90354753*sizeof(int));
tab_pos_indice[16]=malloc(81195210*sizeof(int));
tab_pos_indice[17]=malloc(78077248*sizeof(int));
tab_pos_indice[18]=malloc(59128983*sizeof(int));
tab_pos_indice[19]=malloc(63025520*sizeof(int));
tab_pos_indice[20]=malloc(48129895*sizeof(int));
tab_pos_indice[21]=malloc(51304566*sizeof(int));
tab_pos_indice[22]=malloc(155270560*sizeof(int));
tab_pos_indice[23]=malloc(59373566*sizeof(int));

char chr[100];
char type[100];
char family[100];
int indice,indice_avant;

fic1 = fopen("bornes_10_corrected.dat","r");
if(!fic1) return 1;
j=0;
	      while(fscanf(fic1,"%d %d %d %d",&c,&left,&right,&indice) != EOF)
	      {       
	      tab_left[j]=left;
	      tab_right[j]=right;
	      
		for(p=left;p<right;p++)
		{
		tab_pos_indice[c-1][p-1] = j;
		}
		if(c!=c_a && c_a!=0){printf("%d  %d\n",c-1,j-1);}
		c_a=c;
	      j++;
	      }
fclose(fic1);
NBINS=j;
printf("ALL chromosomes in stock  avec %d bins.\n",NBINS);


fic3 = fopen("/home/axel/Bureau/REPEATS/repeat_masker_2013.dat2","r");
fic6 = fopen("repeat_masker_2013.dat22","w");

while(fscanf(fic3,"%s %s %d %d",type_repeat,chr,&position_repeat1,&position_repeat2) != EOF)
{
	  if (strcmp(type_repeat,type_repeat_avant) != 0  && strcmp(type_repeat_avant,"INIT") != 0 )
	  {  
				fic5 = fopen("hist.dat","w");
				for(i=0;i<NBINS;i++)	  
				{
				if(IND_REPEAT[i] >0) {NBIN_repeat++;fprintf(fic5,"%d\n",IND_REPEAT[i]);sum=sum+IND_REPEAT[i];}  
			        else	              {fprintf(fic5,"%d\n",0);}  
				}
				fclose(fic5);
				printf("There are %d positions and %d bins for the repeat type %s.\n",NPOS_repeat,NBIN_repeat,type_repeat_avant);

				if(NBIN_repeat>=0)
				{  
	  // 		        Fabrication Histogramme et détermination du seuil par R : 
// 				sprintf(cmd, "R --slave --args 10 < hist.r");
//                             system(cmd);
			       
				np=sum;
				pval=1.0;
				int s=0;
				
				while(pval>0.05) {pval = pbinom(s-1,np,1.0/NBINS,0,0);s=s+1;}
				seuil=s-1;
			       
// 				fic7 = fopen("seuil.txt","r");
// 				while(fscanf(fic7,"%d",&seuil)!= EOF) {printf("Seuil de selection pour les bins avec repeat : %d\n",seuil);}
// 				fclose(fic7);
				
				printf("Threshold of selection for bins with repeat: %d\n",seuil);
				
// 				Ecriture new list of repeat positions 
				for(i=0;i<NBINS;i++)	{if(IND_REPEAT[i] >= seuil) {fprintf(fic6,"%s %d %d %d\n",type_repeat_avant,CHR_REPEAT[i],POS_REPEAT[i],i);NSEUIL++;} } 
				printf("There are %d bins above the threshold (%d) for the repeat %s.\n",NSEUIL,seuil,type_repeat_avant);	
				}
			
	  // 		        REINITIALISATION  :
				for(i=0;i<= NBINS;i++)	{IND_REPEAT[i]=0;CHR_REPEAT[i]=0;POS_REPEAT[i]=0;}
				NPOS_repeat=0;NBIN_repeat=0;NSEUIL=0;sum=0;
	  }

  
     if (strcmp(chr,"chr1") == 0) {chr1=1;}
else if (strcmp(chr,"chr2") == 0) {chr1=2;}
else if (strcmp(chr,"chr3") == 0) {chr1=3;}
else if (strcmp(chr,"chr4") == 0) {chr1=4;}
else if (strcmp(chr,"chr5") == 0) {chr1=5;}
else if (strcmp(chr,"chr6") == 0) {chr1=6;}
else if (strcmp(chr,"chr7") == 0) {chr1=7;}
else if (strcmp(chr,"chr8") == 0) {chr1=8;}
else if (strcmp(chr,"chr9") == 0) {chr1=9;}
else if (strcmp(chr,"chr10") == 0) {chr1=10;}
else if (strcmp(chr,"chr11") == 0) {chr1=11;}
else if (strcmp(chr,"chr12") == 0) {chr1=12;}
else if (strcmp(chr,"chr13") == 0) {chr1=13;}
else if (strcmp(chr,"chr14") == 0) {chr1=14;}
else if (strcmp(chr,"chr15") == 0) {chr1=15;}
else if (strcmp(chr,"chr16") == 0) {chr1=16;}
else if (strcmp(chr,"chr17") == 0) {chr1=17;}
else if (strcmp(chr,"chr18") == 0) {chr1=18;}
else if (strcmp(chr,"chr19") == 0) {chr1=19;}
else if (strcmp(chr,"chr20") == 0) {chr1=20;}
else if (strcmp(chr,"chr21") == 0) {chr1=21;}
else if (strcmp(chr,"chr22") == 0) {chr1=22;}
else if (strcmp(chr,"chrX") == 0) {chr1=23;}
else if (strcmp(chr,"chrY") == 0) {chr1=24;}
else  				  {chr1=0;}
  
if(chr1 >0)
{  

// indice_repeat=tab_pos_indice[chr1-1][position_repeat1-1];    
// NPOS_repeat++;
// CHR_REPEAT[indice_repeat]=chr1;
// POS_REPEAT[indice_repeat]=position_repeat;
// IND_REPEAT[indice_repeat]++;    
//   
  
// To take into account sequences that are on two bins: 
  
for(position_repeat=position_repeat1;position_repeat<=position_repeat2;position_repeat++)
{      
indice_repeat=tab_pos_indice[chr1-1][position_repeat-1];      	 

	if(indice_repeat != indice_avant || position_repeat1 != position_repeat_avant )  
	{
	NPOS_repeat++;
	CHR_REPEAT[indice_repeat]=chr1;
	POS_REPEAT[indice_repeat]=position_repeat;
	IND_REPEAT[indice_repeat]++;  
	}
indice_avant =indice_repeat;
position_repeat_avant=position_repeat1;
}

}

strcpy(type_repeat_avant,type_repeat);
}

fclose(fic3);
fclose(fic6);

return 0;
}
