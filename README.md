# Repeats_elements
This repository contains scripts to recreate figures in paper The 3D folding of metazoan genomes correlates with the association of similar repetitive elements (2015). For queries or help getting these running, you can contact me on mail or open an issue at the github repository.


# Dependencies

Scripts will run on OS X and other Unix-based systems. External dependencies should be installed somewhere on your `$PATH`.

## R packages

Lots of commonly-installed R packages are also used, including but not limited to: 

### CRAN

* `calibrate`
* `verification`
* `ROCR`

## External programs

* `bedtools`
* `bigWigAverageOverBed`<sup>*</sup>
* `ICE` / [hiclib](http://mirnylab.bitbucket.org/hiclib)<sup>*</sup>
* `GNU Scientific Library` / [GSL](http://www.gnu.org/software/gsl/) 
* `SRA tool` / [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)


# Raw data
We downloaded raw data from Short Read Archive server at the following address :

ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR400/SRR400260/SRR400260.sra

We convert the files in fastq files using sra tool fastq-dump:   

./fastq-dump /data/human/*



## Filtering of the data
## Normalization of the data
We used the normalization procedure called SCN that we implemented in C. 
With dynamic allocation of memory, C language allows us to allocate big matrices (100000 x 100000) in a station with 50G of ram. 

## Calcul of Colocalisation Scores and random generations

Colocalisation Scores were calculated. 
For random generetion and vectors manipulation, we used the generators from the GNU Scientific Library. 

## Plots of the results
To plot scatter plots of the different repeats elements. we use the R script: scatter_plot.R


## Scripts used for Figures 3
We used R scripts:



## sessionInfo()

Below is an output of sessionInfo() for troubleshooting purposes, some loaded packages may not be required and likewise, some required packages may not be loaded. An exception caused by attached packages is likely due to version issues.
```r
R version 3.2.0 (2015-04-16)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

locale:
 [1] LC_CTYPE=fr_FR.utf8       LC_NUMERIC=C             
 [3] LC_TIME=fr_FR.utf8        LC_COLLATE=fr_FR.utf8    
 [5] LC_MONETARY=fr_FR.utf8    LC_MESSAGES=fr_FR.utf8   
 [7] LC_PAPER=fr_FR.utf8       LC_NAME=C                
 [9] LC_ADDRESS=C              LC_TELEPHONE=C           
[11] LC_MEASUREMENT=fr_FR.utf8 LC_IDENTIFICATION=C      

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base  
