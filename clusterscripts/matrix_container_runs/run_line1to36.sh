#!/bin/bash
#SBATCH -J lines1to36_container
#SBATCH -A typas              
#SBATCH -N 1                        # number of nodes
#SBATCH -n 24                      # number of cores
#SBATCH --mem 20G                # memory pool for all cores
#SBATCH -t 3-00:00                  # runtime limit (D-HH:MM:SS)
#SBATCH -o lines1to36_container.out          # STDOUT
#SBATCH -e lines1to36_container.err          # STDERR
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=leonard.dubois@embl.de


module load R/3.4.3-foss-2017b-X11-20171023

Rscript /home/dubois/MoA_Prediction/clusterscripts/matrix_container_runs/run_matrix_lines.R 1 36