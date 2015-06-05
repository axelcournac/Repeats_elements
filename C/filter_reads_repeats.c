#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// To filter reads or other from repeats sequences 
// 3 : for data of independent alignments 

int main(int argc, char *argv[])
{
FILE* fic1 = NULL;
FILE* fic2 = NULL;
FILE* fic3 = NULL;

int chr1,chr2;
int locus1,locus2;
int sens1,sens2;
int mq1,mq2;
int match1,match2;
int indice1,indice2;
int d1,d2;
int p;
int left,right;
int j=0,jj=0;
unsigned int i=0;

int tab_left[837604];
int tab_right[837604];

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
char c[100];
char type[100];
char family[100];
int indice;
int somme1,somme2; 
long int n_bases_repeat=0;

fic1 = fopen("/data/axel/Hi_seq_Dixon/repeat_masker_hg19.dat1","r");
if(!fic1) return 1;
j=0;
	      while(fscanf(fic1,"%s %d %d",c,&left,&right) != EOF)
	      {
	     if (strcmp(c,"chr1") == 0) {chr1=1;}
	else if (strcmp(c,"chr2") == 0) {chr1=2;}
	else if (strcmp(c,"chr3") == 0) {chr1=3;}
	else if (strcmp(c,"chr4") == 0) {chr1=4;}
	else if (strcmp(c,"chr5") == 0) {chr1=5;}
	else if (strcmp(c,"chr6") == 0) {chr1=6;}
	else if (strcmp(c,"chr7") == 0) {chr1=7;}
	else if (strcmp(c,"chr8") == 0) {chr1=8;}
	else if (strcmp(c,"chr9") == 0) {chr1=9;}
	else if (strcmp(c,"chr10") == 0) {chr1=10;}
	else if (strcmp(c,"chr11") == 0) {chr1=11;}
	else if (strcmp(c,"chr12") == 0) {chr1=12;}
	else if (strcmp(c,"chr13") == 0) {chr1=13;}
	else if (strcmp(c,"chr14") == 0) {chr1=14;}
	else if (strcmp(c,"chr15") == 0) {chr1=15;}
	else if (strcmp(c,"chr16") == 0) {chr1=16;}
	else if (strcmp(c,"chr17") == 0) {chr1=17;}
	else if (strcmp(c,"chr18") == 0) {chr1=18;}
	else if (strcmp(c,"chr19") == 0) {chr1=19;}
	else if (strcmp(c,"chr20") == 0) {chr1=20;}
	else if (strcmp(c,"chr21") == 0) {chr1=21;}
	else if (strcmp(c,"chr22") == 0) {chr1=22;}
	else if (strcmp(c,"chrX") == 0) {chr1=23;}
	else if (strcmp(c,"chrY") == 0) {chr1=24;}
	else  				{chr1=0;}
       
	if(chr1>0)
	{	
		for(p=left;p<=right;p++)
		{
		if(tab_pos_indice[chr1-1][p-1] == 0 ) {n_bases_repeat++;}   
		tab_pos_indice[chr1-1][p-1] = 1;   // we mark the repeat 
		}
	      }
	}

fclose(fic1);
printf("ALL repeats marked on the human genome\n");
printf("There are %ld bases marked with a repeat\n",n_bases_repeat);

//----------------------------------------------------------------------------------------
fic2 = fopen("/data/axel/Hi_seq_TCC/tcc.ai.d120.pcr", "r");
fic3 = fopen("/data/axel/Hi_seq_TCC/tcc.ai.d120.pcr.fr2", "w");
j=0;
jj=0;
while(fscanf(fic2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d",&chr1,&locus1,&sens1,&mq1,&match1,&indice1,&d1,&chr2,&locus2,&sens2,&mq2,&match2,&indice2,&d2) != EOF)
	{
	j++;
	match1=36;
	match2=36;

	if(j % 1000000 == 0) {fprintf(fic3,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}
        somme1=0;somme2=0;

	if(sens1==0)
	{  
	for(p=locus1;p<=locus1+match1;p++)
		{
		somme1 = somme1 + tab_pos_indice[chr1-1][p-1];	
		}
	}
	
	if(sens2==0)
	{  
	for(p=locus2;p<=locus2+match2;p++)
		{
		somme2 = somme2 + tab_pos_indice[chr2-1][p-1];	
		}
	}
	
	if(sens1==16)
	{  
	for(p=locus1;p>=locus1-match1;p--)
		{
		somme1 = somme1 + tab_pos_indice[chr1-1][p-1];	
		}
	}
	
	if(sens2==16)
	{  
	for(p=locus2;p>=locus2-match2;p--)
		{
		somme2 = somme2 + tab_pos_indice[chr2-1][p-1];	
		}
	}
	
if(somme1 == 0  && somme2 == 0)  {jj++;fprintf(fic3,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}

	}

printf("There %d positions from the file.\n",j);
printf("There are %d positions after the repeat filter.\n",jj);

fclose(fic2);
fclose(fic3);

return 0;
}


