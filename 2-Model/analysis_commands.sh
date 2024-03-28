#!/bin/bash

cp params.py params_save.py

## whole standard
sed -i 's/^pred_data_type.*/pred_data_type = "whole"/' params.py
sed -i 's/^normalisation_type.*/normalisation_type = "standard"/' params.py
echo "---- Prediction on whole data with standard normalization"
python3 -m cluster_alignment ## fit model
python3 -m cluster_curves ## get model performance
python3 -m cluster_prediction # get prediction
python3 -m cluster_barplot ## get barplots

## targeted standard 
echo ""
sed -i 's/^pred_data_type.*/pred_data_type = "targeted"/' params.py
sed -i 's/^normalisation_type.*/normalisation_type = "standard"/' params.py
echo "---- Prediction on targeted data with standard normalization"
python3 -m cluster_alignment ## fit model
python3 -m cluster_curves ## get model performance
python3 -m cluster_prediction # get prediction
python3 -m cluster_barplot ## get barplots

mv params_save.py params.py
