#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//  THis code makes the matrices from the files of reads, normalizes it, 
//  calculates the Colocalization Scores and generate random sets 
//  with the conservation of the distribution in GC content or coverage. 

int main(int argc, char *argv[])
{
#define foreach(idxtype, idxpvar, col, colsiz ) idxtype* idxpvar; for( idxpvar=col ; idxpvar < (col + (colsiz)) ; idxpvar++)
#define arraylen( ary ) ( sizeof(ary)/sizeof(ary[0]) )

int BIN= atoi(argv[1]);   //   number of fragments that contains the bin   

int NTF=837603;  // nbre total de frags
int NBINS;
FILE* fic1 = NULL;FILE* fic2 = NULL;FILE* fic7 = NULL;FILE* fic55 = NULL;
FILE* fic3 = NULL;FILE* fic5 = NULL;FILE* fic6 = NULL;FILE* fic4 = NULL;
int chr1,chr2;
int locus1,locus2;
int sens1,sens2;
int mq1,mq2;
int match1,match2;
int indice1,indice2;
int d1,d2;
int left,right;
int j=0,b=0,NPOS_repeat=0;
unsigned int i=0;
short int c,c_a=0;
int label;

int *IND_REPEAT;IND_REPEAT = calloc (837604 , sizeof(int));
int *POS_REPEAT;POS_REPEAT = calloc (837604 , sizeof(int));
int *keep;keep             = calloc (837604 , sizeof(int));
int Npos_random;
char type_repeat[100];
char type_repeat_avant[100];
strcpy(type_repeat, "INIT");strcpy(type_repeat_avant, "INIT");
unsigned int chr1_repeat,position_repeat,indice_repeat,indice_repeat_before;
int NBIN_repeat=0, NBIN_repeat2=0,NBIN_DETECT_repeat=0;
char cmd[300],file1[300],file2[300];
int ind,nber_pol,seuil,sig;
int groupe00=0,groupe01=0,groupe11=0;
float groupe00d=0,groupe01d=0,groupe11d=0; 
int n_bin_detect = 0;
const gsl_rng_type * T;
gsl_rng * r;      // create a generator chosen by the environment variable GSL_RNG_TYPE 
gsl_rng_env_setup();
T = gsl_rng_default;
r = gsl_rng_alloc (T);
int rr;

ind=0;
int keep1,keep2;
float norm[83761];
int sum[83761];
int chromo[83761];
int chromo_area[83761];

float component,GC;
int cover;
int polII,ctcf,dnase;
float R_time;

int TAB_REPEAT[10000];
int n,ii;
float moy_GC_repeat=0,moy_GC_random=0;

fic2 = fopen(argv[2], "r");

// fic3 = fopen("/home/romain/Desktop/divers/CELL/hesc_repeat/repeat_masker_hg19_10.dat22","r");
// fic3 = fopen("/home/romain/Desktop/divers/CELL/hesc_repeat/ITS-T2AG3-10-7.bed.dat.correct.ind.binomial.format","r");
fic3 = fopen("hesc_repeat/repeat_masker_2013.dat22","r");
// fic3 = fopen("nads_10.dat22","r");
// fic3 = fopen("TF_and_others_trna_nads_10.dat23","r");
// fic3 = fopen("/home/romain/Desktop/divers/CELL/TF_and_others_trna_nads_10.dat22","r");

// fic3 = fopen("/home/romain/Desktop/divers/CELL/syn_seq_repeats.dat22","r");
// fic3 = fopen("/home/romain/Desktop/divers/CELL/out_s20_masked.dat22","r");

fic7 = fopen(argv[3],"r");

fic4 =fopen("hist_hesc_10.dat","w");

for(i=0;i<=(int) NTF/BIN;i++)	  {IND_REPEAT[i]=0;keep[i]=0;}

// -----------------------------------------MATRICE INT CONSTRUCTION-------------------------------------------------------------------
float** MAT_INT = malloc(83761*sizeof(float*));
for(i=0;i<83761;i++){MAT_INT[i] = malloc(83761*sizeof(float));}

j=0;
    while(fscanf(fic2,"%d %d %d %d %d %d %d    %d %d %d %d %d %d %d",&chr1,&locus1,&sens1,&mq1,&match1,&indice1,&d1,&chr2,&locus2,&sens2,&mq2,&match2,&indice2,&d2) != EOF)
    {
    j++;
    if(j % 1000000 == 0){printf("%d %d %d %d %d %d %d %d   %d %d %d %d %d %d %d\n",j,chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}
 
    if(chr1 > 0 && chr2 > 0)
    {
     //  Binnage : 
    indice1 = (int) indice1/BIN;
    indice2 = (int) indice2/BIN;
	    
    MAT_INT[indice1][indice2]++;
    MAT_INT[indice2][indice1]++;
    
    chromo[indice1]=chr1;
    chromo[indice2]=chr2;
    }
    }
    
printf("There are %d interactions 3D from the HiC file.\n",j);
printf("MATRICE OF RAW INTERACTIONS CONSTRUCTED\n");
fclose(fic2);

// -------------------------------------------------------------------------------------------------------------------------------------------------
//  NORMALISATION MATRICE :  SCN Norme 2 


float seuil_bin = atoi(argv[4]);
float seuil_bin2 = atoi(argv[5]);
int i_bin=0;

if(1==1)
{

for(j=0;j<= (int) NTF/BIN;j++)  // colonne 
{
  norm[j] = 0; 
  for(i=0;i<= (int) NTF/BIN;i++)	  
  { 
  norm[j] = norm[j] + MAT_INT[i][j];  
  sum[j] = sum[j] + MAT_INT[i][j]; 
  }

  fprintf(fic4,"%f\n",norm[j]);
  
  if(norm[j] >=seuil_bin && norm[j] < seuil_bin2 ){keep[j]=1;i_bin++;}
  if(norm[j] <seuil_bin){norm[j]=0;keep[j]=0;} 
}

printf("We keep %d bins with norm > %f .\n",i_bin,seuil_bin);
fclose(fic4);


for(i=0;i<= (int) NTF/BIN;i++)	  // ligne 
{
    for(j=0;j<= (int) NTF/BIN;j++)	  // colonne 
    {
    if( keep[i]==1 &&  keep[j]==1) {MAT_INT[i][j] = MAT_INT[i][j]/norm[j];}
    else 		 {MAT_INT[i][j] = 0;}
    }
}

printf("1\n");
// ----------------------------------------------------
// ----------------------------------------------------
for(i=0;i<= (int) NTF/BIN;i++)  // ligne
{
  norm[i] = 0; 
  for(j=0;j<= (int) NTF/BIN;j++)	  	  
  { 
  norm[i] = norm[i] + MAT_INT[i][j];  // 2
  }

}

for(i=0;i<= (int) NTF/BIN;i++)	  // ligne 
{
for(j=0;j<= (int) NTF/BIN;j++)	  // colonne 
    {  
    if( keep[i]==1 &&  keep[j]==1)  {MAT_INT[i][j] = MAT_INT[i][j]/norm[i];}
    else 		 {MAT_INT[i][j] = 0;}
    }
}
printf("2\n");
// ----------------------------------------------------
// ------------------------------------------------
for(j=0;j<= (int) NTF/BIN;j++)  // colonne 
{
  norm[j] = 0; 
  for(i=0;i<= (int) NTF/BIN;i++)	  
  { 
  norm[j] = norm[j] + MAT_INT[i][j];  // 2
  }

}


for(i=0;i<= (int) NTF/BIN;i++)	  // ligne 
{
    for(j=0;j<= (int) NTF/BIN;j++)	  // colonne 
    {
    if( keep[i]==1 &&  keep[j]==1) {MAT_INT[i][j] = MAT_INT[i][j]/norm[j];}
    else 		 {MAT_INT[i][j] = 0;}
    }
}
printf("3\n");
// ----------------------------------------------------
// ----------------------------------------------------
for(i=0;i<= (int) NTF/BIN;i++)  // ligne
{
  norm[i] = 0; 
  for(j=0;j<= (int) NTF/BIN;j++)	  	  
  { 
  norm[i] = norm[i] + MAT_INT[i][j];  // 2
  }

}


for(i=0;i<= (int) NTF/BIN;i++)	  // ligne 
{
    for(j=0;j<= (int) NTF/BIN;j++)	  // colonne 
    {  
    if( keep[i]==1 &&  keep[j]==1)  {MAT_INT[i][j] = MAT_INT[i][j]/norm[i];}
    else 		 {MAT_INT[i][j] = 0;}
    }
}
printf("4\n");

printf("MATRICE INT NORMALISED WITH SCN NORM SUM\n");
}


// -------------------------------------------------------------------------------------------------------------------------------------------------
int boundary =0;
int boundary2 =0;
int ind_area= 0,ind_area2= 0;
int N_bins_parameter= 200;
int GC_bin;
int *GC_freq_genome;GC_freq_genome = calloc (N_bins_parameter, sizeof(int));

int** TAB_GC; 
TAB_GC= calloc(N_bins_parameter,sizeof(int*));
for(i=0;i<N_bins_parameter;i++) {TAB_GC[i] = (int *) calloc(10000,sizeof(int));}

float bin_GC = 10;
int GC_indice[120000];  // give for one bin indice the GC bin indice

// -------------------------------------------------------------------------------------------------------------------------------------------------
    while(fscanf(fic7,"%d %d",&ind,&cover) != EOF)
    { 
     if(cover >= seuil_bin && cover < seuil_bin2 )   { keep[ind]=1; } 
//      keep[ind]=1;
      
     if(keep[ind]==1) 
     {     
     ind_area++;
     GC_bin = (int) (cover / bin_GC);
     
	if(GC_bin < N_bins_parameter)
	{
	GC_indice[ind]=GC_bin;
	TAB_GC[GC_bin][ GC_freq_genome[GC_bin] ] = ind;
	
	GC_freq_genome[GC_bin]++;
	}
     }
    } 
printf("BINS IN THE AREA DETECTABLE : %d \n",ind_area);
printf("Initial Distribution of cover in memory \n");




// ---------------------------------------------------------------------------------------------------------------------------------------------------   
int *GC_freq_repeat;GC_freq_repeat = calloc (N_bins_parameter, sizeof(int));

indice_repeat_before=-1;
while(fscanf(fic3,"%s %d %d %d",type_repeat,&chr1_repeat,&position_repeat,&indice_repeat) != EOF)
{
	if (strcmp(type_repeat,type_repeat_avant ) != 0  && strcmp(type_repeat_avant,"INIT") != 0  )
	{
	//  Output files : 
	sprintf(file1, "OUTPUT_TEMP2/repeat_%s",type_repeat_avant);
	sprintf(file2, "OUTPUT_TEMP2/random_%s",type_repeat_avant);

	fic5 = fopen(file1,"w");    //    NAME FILE   
	fic55 =fopen(file2,"w");
	
	  int * TAB_REPEAT2;
	  TAB_REPEAT2  = (int*)calloc(NBIN_repeat2, sizeof(int));
	  for(n=0;n<NBIN_repeat2;n++) {TAB_REPEAT2[n] = TAB_REPEAT[n];}   
	  
		foreach (int, p, TAB_REPEAT2, NBIN_repeat2)  
		{
		 i=*p; 
		 if(IND_REPEAT[i] > 0 ) {NBIN_repeat++;} 
		 if(IND_REPEAT[i] > 0 && keep[i] == 1) {NBIN_DETECT_repeat++;}
		 
		    foreach (int, p, TAB_REPEAT2, NBIN_repeat2)   
		    { 
		      j=*p;
		      if( chromo[i] != chromo[j] && keep[i] ==1 && keep[j] ==1)  {groupe11++; if(MAT_INT[i][j] > 0) {groupe11d = groupe11d + MAT_INT[i][j];}  }	
		    }
		}
		
		printf("There are %d bins for the repeat type %s.\n",NBIN_repeat2,type_repeat_avant);
		printf("%f %d  %d %s %d %d  %f\n",groupe11d,groupe11,(int) BIN,type_repeat_avant,NBIN_repeat,NBIN_DETECT_repeat,moy_GC_repeat/NBIN_DETECT_repeat);
		
		fprintf(fic5,"%f %d  %d %s %d %d  %f\n",groupe11d,groupe11,(int) BIN,type_repeat_avant,NBIN_repeat,NBIN_DETECT_repeat,moy_GC_repeat/NBIN_DETECT_repeat);
		
		
// -------------------------------------------------------------------------------------------------------------------------------------------------------------
// RANDOM PART 
//--------------------------------------------------------------------------------------------------------------------------------------------------------------

int * vect_good_bin;
vect_good_bin  = (int*)calloc(ind_area, sizeof(int));
 for(n=0;n<NBIN_repeat2;n++) {TAB_REPEAT2[n] = 0;}  

      for(rr=1;rr<=1000;rr++) 
      {  
	groupe00=0;groupe00d=0;groupe01=0;groupe01d=0;groupe11=0;groupe11d=0;
      
		      // ramdom generation of detectables bins Uniform model 
// 		      int* sample = NULL;
// 		      sample = malloc(Npos_random * sizeof(int));
// 		      gsl_ran_choose (r, sample, Npos_random, vect_good_bin, ind_area, sizeof (int)); 
		      
// 		      Null Model with the same GC distribution : 
		      moy_GC_random =0;
		      ii=0;	
		      
		      for(i=0;i<N_bins_parameter;i++)
		      {  
		      int * sample; sample = calloc(GC_freq_repeat[i], sizeof(int));
		      gsl_ran_choose (r, sample, GC_freq_repeat[i],TAB_GC[i],GC_freq_genome[i], sizeof (int));
		      foreach (int, p, sample,GC_freq_repeat[i] ) {TAB_REPEAT2[ii] = *p;ii++;}
		      free(sample);
		      }	
		      
		      foreach (int, p, TAB_REPEAT2, NBIN_repeat2) 
		      {	
		      i=*p;	
		      moy_GC_random = moy_GC_random + GC_indice[i]*bin_GC;
 
			  foreach (int, p, TAB_REPEAT2, NBIN_repeat2)   
			  {
			  j=*p;  
			  if( chromo[i] != chromo[j] && keep[i] ==1 && keep[j] ==1) {groupe11++;  if(MAT_INT[i][j] > 0) {groupe11d = groupe11d + MAT_INT[i][j];}  }
			  }
		      }
		      fprintf(fic55,"%f %d  %d %d %d %d %f\n",groupe11d,groupe11,(int) BIN,rr,NBIN_repeat,NBIN_DETECT_repeat,moy_GC_random/NBIN_repeat2);
          }		  
	

	  // 	REINITIALISATION  :
		for(i=0;i<= (int) NTF/BIN;i++) {IND_REPEAT[i]=0;}
		NPOS_repeat=0;NBIN_repeat=0;NBIN_repeat2=0;NBIN_DETECT_repeat=0;groupe00=0;groupe00d=0.0;groupe01=0;groupe01d=0.0;groupe11=0;groupe11d=0.0;
		for(i=0;i< N_bins_parameter;i++) {GC_freq_repeat[i]=0;}
		indice_repeat_before=-1;
		moy_GC_random =0;moy_GC_repeat =0;
		
		sprintf(cmd, "Files repeat_%s and random_%s sent to the OUTPUT directory.",type_repeat_avant,type_repeat_avant);
		printf("%s\n",cmd);

		
		fclose(fic5);
		fclose(fic55);
	}	


  if(keep[indice_repeat] ==1)
  {  
  NPOS_repeat++;
  IND_REPEAT[indice_repeat]++;
      if(indice_repeat != indice_repeat_before) 
      {
	TAB_REPEAT[NBIN_repeat2] = indice_repeat; 
	NBIN_repeat2++; moy_GC_repeat =   moy_GC_repeat + GC_indice[indice_repeat]*bin_GC;
	GC_freq_repeat[ GC_indice[indice_repeat] ]++; 
      }
  
  indice_repeat_before = indice_repeat;

  }
 	
strcpy(type_repeat_avant,type_repeat);
}


free(MAT_INT);
fclose(fic3);

return 0;
}
