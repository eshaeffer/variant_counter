#! /bin/bash

raws=$(ls raw)
mafs=$(ls maf)
veppy=$(ls vepped)

##Replace with your directory to store variants filtered variants/ make tables
rareDir=../../vcfs/rareAndDeleterious
tables=../../tables

for raw in $raws
do
./countVars.py -input raw/$raw -output $tables/$raw
done

for af in $mafs
do
./countVars.py -input maf/$af -option mafCounting -output $tables/$af
done

for vep in $veppy
do
./countVars.py -input vepped/$vep -option vepCounting -output $tables/$vep -outputDir $rareDir
done
