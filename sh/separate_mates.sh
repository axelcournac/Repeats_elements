#!/bin/bash

# To separate the two ends of the reads for bowtie2. 

for i in *.fastq ;
do echo $i;


echo "end1";

awk 'BEGIN {FS = "";c=0;}  {c++;   \
if(c%4==1) {print $0;}   \
if(c%4==2) {for(i=1;i<=NF/2;i++) {printf "%s",$i;} printf "\n";} \
if(c%4==3) {print $0;} \
if(c%4==0) {for(i=1;i<=NF/2;i++) {printf "%s",$i;} printf "\n";} 
}'  $i > $i.end1 ;


echo "end2";


awk 'BEGIN {FS = "";c=0;}  {c++;   \
if(c%4==1) {print $0;}   \
if(c%4==2) {for(i=NF/2+1;i<=NF;i++) {printf "%s",$i;} printf "\n";} \
if(c%4==3) {print $0;} \
if(c%4==0) {for(i=NF/2+1;i<=NF;i++) {printf "%s",$i;} printf "\n";} 
}'  $i > $i.end2 ;



done

