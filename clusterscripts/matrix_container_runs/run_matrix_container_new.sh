#!/bin/bash
#SBATCH -J new_container_run_correctedResamp
#SBATCH -A typas              
#SBATCH -N 1                        # number of nodes
#SBATCH -n 24                      # number of cores
#SBATCH --mem 20G                # memory pool for all cores
#SBATCH -t 7-00:00                  # runtime limit (D-HH:MM:SS)
#SBATCH -o new_container_run_correctedResamp.out         # STDOUT
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=leonard.dubois@embl.de


module load R/3.4.3-foss-2017b-X11-20171023

Rscript /home/dubois/MoA_Prediction/clusterscripts/matrix_container_runs/run_matrix_container_new.R
