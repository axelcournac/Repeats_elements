#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash-1.9.6/src/uthash.h"

// 05/15/2013
// Remove PCR duplicates according to their genomic positions 

typedef struct struct_entry {
    char seq[200];                        
    UT_hash_handle hh;         /* makes this structure hashable */
} struct_entry;

int main(int argc, char *argv[])
{
FILE* fic1 = NULL;FILE* fic2 = NULL;FILE* fic3 = NULL;

unsigned int i=0;
char word12[300];
char word21[300];
char line_before[300];
int bool;
int f=-2,ff=0;
fic1 = fopen(argv[1],"r");
if(!fic1) return 1;

struct_entry *group=NULL,*tmp;
struct_entry *entry1,*entry2,*entry3,*entry4;

char chr1[100],chr2[100];
char locus1[100],locus2[100];
int sens1,sens2;
int mq1,mq2;
int match1,match2;
int indice1,indice2;
int d1,d2;

char name[100];
strcpy(name,argv[1]);
strcat(name,".pcr");
fic3 = fopen(name, "w");

	      while(fscanf(fic1,"%s %s %d %d %d %d %d  %s %s %d %d %d %d %d\n",chr1,locus1,&sens1,&mq1,&match1,&indice1,&d1,chr2,locus2,&sens2,&mq2,&match2,&indice2,&d2) != EOF)
	      {
	      f++;
	      strcpy(word12,"");		    strcpy(word21,""); 
	      strcat(word12,chr1);                 strcat(word21,chr2);
	      strcat(word12," ");                  strcat(word21," ");
	      strcat(word12,locus1);               strcat(word21,locus2);
	      strcat(word12," ");                  strcat(word21," ");
	      strcat(word12,chr2);                 strcat(word21,chr1);
	      strcat(word12," ");                  strcat(word21," ");
	      strcat(word12,locus2);               strcat(word21,locus1);
	      strcat(word12," ");                  strcat(word21," "); 

	      if ( (entry1 =  (struct_entry*)malloc(sizeof(struct_entry))) == NULL) exit(-1);
// 	      if ( (entry2 =  (struct_entry*)malloc(sizeof(struct_entry))) == NULL) exit(-1);
	      
	      strcpy(entry1->seq, word12);
// 	      strcpy(entry2->seq, word21);
	      
	      HASH_FIND_STR(group,word12,entry3);
	      HASH_FIND_STR(group,word21,entry4);
	      
	      if(entry3 || entry4) {}  //  if these sequences were already found 
	      else {HASH_ADD_STR(group,seq, entry1);bool=1;ff++;fprintf(fic3,"%s %s %d %d %d %d %d    %s %s %d %d %d %d %d\n",chr1,locus1,sens1,mq1,match1,indice1,d1,chr2,locus2,sens2,mq2,match2,indice2,d2);}
	      
	      }	
	      
fclose(fic1);

printf("There are %d interactions 3D from the file.\n",f);
printf("There are %d interactions 3D from the file after the PCR duplicates filter.\n",ff);

/* free the hash table contents */
HASH_ITER(hh, group, entry1, tmp) 
{
HASH_DEL(group, entry1);
free(entry1);  
}

return 0;
}

