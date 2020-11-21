# variant_counter

Created: Summer 2020 (year of COVID-19)

Prerequisites: python, os, argparse

Description: command line tool for counting variants

Usage: countVars.py [-h] [-input INPUT] [-output OUTPUT] [-option OPTION] [-k]
                    [-outputDir OUTPUTDIR]
                    
 -input user created text file conaining paths to the vcf files that the user would like counted
 -output where the tables containing variants are outputted to
 -option rawCalls,mafCouting,vepCounting which describe different options for counting variants.
 They go as follows: rawCalls simply counts the variants
 mafCounting outputs all the variants with MAF <= 0.01 bcfGnomAD_AF tag and all the variants
 vepCounting outputs all the variants with MAF <=0.01, bcfGnomAD_AF tag, all the variants, the number
 of variants with greater than 20 cadd, number greater than 15 but less than 20 cadd, and the same but
 with 0.01 MAF as well
 -k flag that tells VEP to count MAF using VEP notations (experimental - do not use)
 -outputDir output directory for vepCounting option, outputs all the variants with metrics described in the option option.
 
 countAll.sh allows you to quickly count a large cohort of files. Say you have seperate directories that
 hold vcfs that contain those variants from haplotype caller, those variants that have been annotated 
 with MAF, and those variants with VEP. Then you would make text files containing the paths to those 
 files, and place each of them in their respective directories located within this project (maf, vep, raw) and run ./countAll.sh

pros/cons:
- very good on memory
- Count variants/ filter them on large(ish) scale.
- accurate
- slower than if I had used python standard libraries such as numpy/pandas
