#!/bin/bash

# Script name: rsync_from_server.sh

# Description: This script synchronises files from the directory /home/dubois/MoA_Prediction/run_results
# to the local (Github) MoA_Prediction project into run_results_from_server/

SOURCEFILES="/scratch/typas/cluster_run_results/matrix_container_result/*"
SOURCEHOST="login.cluster.embl.de"
USER=$(whoami)
TARGETDIR="./run_results_from_server/matrix_container_result"

# echo "Sending from ${FROM_HOSTNAME}@${USER} to ${TO_HOSTNAME}@${USER}"

# options:
# -v = increase verbosity
# -z = use compression
# --progress: show progress during transfer

rsync -vuz --progress ${USER}@${SOURCEHOST}:${SOURCEFILES} ${TARGETDIR}
