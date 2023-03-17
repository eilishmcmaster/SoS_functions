# SoS_functions
Functions used for my SOS analyses 


### Notes

* **Average expected heterozygosity (He)** is calculated by BasicStats by getting the reference allele frequency for each locus in a population, calculating He using HW, and then averaging across all loci
* **Average observed heterozygosity (Ho)** is calculated by finding all heterozygous loci (==1) and then dividing by all loci in a population 

* **FST** only uses common loci with a specified MAF when comparing groups, so it is not ideal for genetic groups that are not similar
