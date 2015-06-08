#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// November 2013 :
// To have the cover of the bins from a HiC experiment for human genome (hg19). 


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
int *COVER;COVER = calloc (837604 , sizeof(int));

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

fic1 = fopen("bounds_corrected.dat","r");
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

fic2 = fopen("/data/axel/Hi_seq_TCC/tcc.ai.d120.pcr.fr2", "r");
fic6 = fopen("covering_tcc_fr.dat","w");

    while(fscanf(fic2,"%d %d %d %d %d %d %d    %d %d %d %d %d %d %d",&chr1,&locus1,&sens1,&mq1,&match1,&indice1,&d1,&chr2,&locus2,&sens2,&mq2,&match2,&indice2,&d2) != EOF)
    {
    j++;
//     if(j % 1000000 == 0){printf("%d %d %d %d %d %d %d %d   %d %d %d %d %d %d %d\n",j,chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}
    
	if(chr1 > 0)
	{  
	indice1 =tab_pos_indice[chr1-1][locus1-1]; 
	indice2 =tab_pos_indice[chr2-1][locus2-1]; 

        if(indice1 ==0 ) {printf("%d %d %d %d %d %d %d %d   %d %d %d %d %d %d %d\n",j,chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}
    
	COVER[indice1]++;
	COVER[indice2]++;
	}
    }
    
printf("There are %d interactions 3D from the HiC file.\n",j);
printf("MATRICE OF RAW INTERACTIONS CONSTRUCTED\n");
fclose(fic2);


for(i=0;i<NBINS;i++)	  
{
fprintf(fic6,"%d %d\n",i,COVER[i]); 
}  
 
fclose(fic6);

return 0;
}
