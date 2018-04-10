#!/bin/bash
#SBATCH -J MoA_pred_ML_BR_RF      
#SBATCH -A typas              
#SBATCH -N 1                        # number of nodes
#SBATCH -n 24                       # number of cores
#SBATCH --mem 20G                # memory pool for all cores
#SBATCH -t 7-00:00                  # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=leonard.dubois@embl.de


Rscript cluster_run_trees.R