#!/bin/bash

# *** Confirm Directory Structure *************************************************
if [ ! -d "./raw_data" ]; then
  echo "Missing raw_data directory. Please download the zip file per the README"
fi

if [ ! -d "./inputs" ]; then
  mkdir -p inputs;
fi

if [ ! -d "./output" ]; then
  mkdir -p output;
  cd output;
  mkdir -p baseline;
  mkdir -p JJ_eff;
  mkdir -p no_hes;
  mkdir -p rt_13;
  cd ..;
fi

if [ ! -d "./exhibits" ]; then
  mkdir -p exhibits;
  cd exhibits;
  mkdir -p baseline;
  mkdir -p JJ_eff;
  mkdir -p no_hes;
  mkdir -p rt_13;
  mkdir -p 15_mil;
  cd ..;
fi

# *** Run R Scripts to Populate /inputs *******************************************
Rscript -e 'source("./code/R/main.R")';

# *** Run MATLAB Code to Populate /ouputs, /exhibits/ *****************************
cd sir_code
matlab -r "generate_plot_data;exit;"; 
matlab -r "flags.no_hes   = true; generate_plot_data; exit;"; 
matlab -r "flags.JJ_eff   = true; generate_plot_data; exit;"; 
matlab -r "flags.high_mit = true; generate_plot_data; exit;"; 
matlab -r "make_plots;";
matlab -r "compute_paper_numbers;";
 
