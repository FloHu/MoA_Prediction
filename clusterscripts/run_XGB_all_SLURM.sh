#!/bin/bash
#SBATCH -J ALL_XBG_test  
#SBATCH -A typas              
#SBATCH -N 1                        # number of nodes
#SBATCH -n 24                      # number of cores
#SBATCH --mem 20G                # memory pool for all cores
#SBATCH -t 0-12:00                  # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=leonard.dubois@embl.de


module load R/3.4.3-foss-2017b-X11-20171023

Rscript /home/dubois/MoA_Prediction/clusterscripts/cluster_run_testAll_XGB.R