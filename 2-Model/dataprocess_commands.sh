#!/bin/bash

#cd DATA/

unzip -qq -o matrice.zip

## targeted data
NBCELLS=$(($(($(wc -l matrice/R2/Matrice_0* | tail -n1 | cut -f1 -d't')-$(ls matrice/R2/ | wc -l)))+$(($(wc -l matrice/R3/Matrice_0* | tail -n1 | cut -f1 -d't')-$(ls matrice/R3/ | wc -l)))))
echo "nb cells = "$NBCELLS

echo "sample,STIM,PND" > final_metadata_targeted.csv
head -n1 matrice/R2/Matrice_006F.csv > final_targeted.csv 
touch gene_order.csv

for FILE in $(ls matrice/R2/Matrice_0*)
do
    FNAME=$(echo $FILE | cut -d'/' -f3)
    head -n1 $FILE >> gene_order.csv
    sed '1d' $FILE | sed "s/^/R2-/g" >> final_targeted.csv
    COND=$(grep "R2,$FNAME" metadata_targeted.csv | cut -d"," -f3,4) ## metadata_targeted.ods -> metadata_targeted.csv
    for SAMP in $(cat $FILE | sed "1d" | cut -d',' -f1)
    do
        echo "R2-"$SAMP","$COND >> final_metadata_targeted.csv
    done
done

for FILE in $(ls matrice/R3/Matrice_0*)
do
    FNAME=$(echo $FILE | cut -d'/' -f3)
    head -n1 $FILE >> gene_order.csv
    sed '1d' $FILE | sed "s/^/R3-/g" >> final_targeted.csv
    COND=$(grep "R3,$FNAME" metadata_targeted.csv | cut -d"," -f3,4 | sed 's/pns/pbs/g') ## metadata_targeted.ods -> metadata_targeted.csv 
    for SAMP in $(cat $FILE | sed "1d" | cut -d',' -f1)
    do
        echo "R3-"$SAMP","$COND >> final_metadata_targeted.csv
    done
done

echo $((NBCELLS+1))" == "$(wc -l final_targeted.csv | cut -d' ' -f1)
echo $((NBCELLS+1))" == "$(wc -l final_metadata_targeted.csv | cut -d' ' -f1)
echo "is the same gene order 1 == "$(cat gene_order.csv | uniq -u | wc -l)
rm gene_order.csv metadata_targeted.csv
rm -rf matrice/
