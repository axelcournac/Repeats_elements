#!/bin/bash


# R --slave --args 'repeat_AluJb' 'random_AluJb' < null_model.R


for i in random_*; do j=$(echo $i | sed "/random_/s///" ); echo "$j" ; 
R --slave --args "repeat_$j" "random_$j" < null_model.R ; done 
