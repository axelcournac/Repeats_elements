#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  To assign a fragment to reads 
//  and make several filters according to the distance between beginning of the mapping and next restriction site 

int main(int argc, char *argv[])
{
FILE* fic1 = NULL;
FILE* fic2 = NULL;
FILE* fic3 = NULL;

int th1 = 200;  //  minimum size of the DNA segment to be valid
int th2 = 500;  //  maximum size of the DNA segment to be valid
int th3 = 20;   //  minimum distance betweem start of the read and next restriction site to be considered valid

int chr1,chr2;
int locus1,locus2;
int sens1,sens2;
int mq1,mq2;
int match1,match2;
int indice1,indice2;
int d1,d2;
int p;
int left,right;
int j,jj;
unsigned int i=0;
int c,c_a=0;
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


int tab_fin[25];

tab_fin[0]=249250621;
tab_fin[1]=243199373;
tab_fin[2]=198022430;
tab_fin[3]=191154276;
tab_fin[4]=180915260;
tab_fin[5]=171115067;
tab_fin[6]=159138663;
tab_fin[7]=146364022;
tab_fin[8]=141213431;
tab_fin[9]=135534747;
tab_fin[10]=135006516;
tab_fin[11]=133851895;
tab_fin[12]=115169878;
tab_fin[13]=107349540;
tab_fin[14]=102531392;
tab_fin[15]=90354753;
tab_fin[16]=81195210;
tab_fin[17]=78077248;
tab_fin[18]=59128983;
tab_fin[19]=63025520;
tab_fin[20]=48129895;
tab_fin[21]=51304566;
tab_fin[22]=155270560;
tab_fin[23]=59373566;


fic1 = fopen("/data/axel/frag_hindiii_chrALL.dat1","r");
if(!fic1) return 1;
j=0;
	      while(fscanf(fic1,"%d %d %d",&c,&left,&right) != EOF)
	      {       
	      j++;  //indice du frag 
	      tab_left[j]=left;
	      tab_right[j]=right;
	      
		for(p=left;p<right;p++)
		{
		tab_pos_indice[c-1][p-1] = j;
		}
		if(c!=c_a && c_a!=0){printf("%d  %d\n",c-1,j-1);}
		c_a=c;
	      }
fclose(fic1);
printf("ALL chromosomes in memory  avec %d frag.\n",j);

//----------------------------------------------------------------------------------------
fic2 = fopen(argv[1], "r");
char name[100];
strcpy(name,argv[1]);
strcat(name,".d120.inter");
fic3 = fopen(name, "w");

j=0;
jj=0;

 while(fscanf(fic2,"%d %d %d %d %d %d %d %d \n",&chr1,&locus1,&sens1,&chr2,&locus2,&sens2,&indice1,&indice2) != EOF)
	{
		      j++;
		      //  conversions from the output of Mirny Lib:
		      chr1=chr1+1;
		      chr2=chr2+1;
		      locus1=locus1+1;
		      locus2=locus2+1;

		      if(sens1==0){sens1=16;}       if(sens1==1){sens1=0;}
          if(sens2==0){sens2=16;}       if(sens2==1){sens2=0;}
	      
		      mq1=00;match1=00;
		      mq2=00;match2=00;	
		      if(chr1 >0  && chr2 >0 && locus1 < tab_fin[chr1-1] && locus2 < tab_fin[chr2-1])	
		      {	
		      //correction pour les négatifs car le logiciel de mapping donne la borne inférieure en position absolue : 
		      //if(sens1 == 16) {locus1=locus1+ match1;}  //position mapping corrigé pour les - 
		      //if(sens2 == 16) {locus2=locus2+ match2;}  //position mapping corrigé pour les - 

		      if(tab_pos_indice[chr1-1][locus1-1] >0 ) {indice1 = tab_pos_indice[chr1-1][locus1-1];}
		      if(tab_pos_indice[chr2-1][locus2-1] >0 ) {indice2 = tab_pos_indice[chr2-1][locus2-1];}

		      //distance between the beginning of the mapping and the next restriction site  : 
		      if(tab_pos_indice[chr1-1][locus1-1] >0 ) 
		      {
		      if(sens1 == 0) {d1=tab_right[indice1]-locus1;}
		      else  	     {d1=locus1-tab_left[indice1];}
		      }

		      if(tab_pos_indice[chr2-1][locus2-1] >0 ) 
		      {
		      if(sens2 == 0) {d2=tab_right[indice2]-locus2;}
		      else  	     {d2=locus2-tab_left[indice2];}
		      }
		      
if(j % 1000000 == 0)
{printf("%d\n",j);printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}  

if (th1 < d1+d2 && d1+d2 < th2  && d1 > th3 && d2 > th3 &&  chr1 != chr2)
{jj++;fprintf(fic3,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);} 

}
	}

printf("There are %d interactions 3D from the file.\n",j);
printf("There are %d interactions 3D from the file after filters.\n",jj);

fclose(fic2);
fclose(fic3);

return 0;
}
