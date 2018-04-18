#!/bin/bash

# Script name: rsync_from_server.sh

# Description: This script synchronises files from the directory /home/dubois/MoA_Prediction/run_results
# to the local (Github) MoA_Prediction project into run_results_from_server/

SOURCEFILES="/home/dubois/MoA_Prediction/run_results/*.RData"
SOURCEHOST="login.cluster.embl.de"
USER=$(whoami)
TARGETDIR="./run_results_from_server/"

# echo "Sending from ${FROM_HOSTNAME}@${USER} to ${TO_HOSTNAME}@${USER}"

# options:
# -a = archive mode
# -v = increase verbosity
# -z = use compression
# --progress: show progress during transfer

rsync -avz --progress ${USER}@${SOURCEHOST}:${SOURCEFILES} ${TARGETDIR}
