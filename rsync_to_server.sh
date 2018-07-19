#!/bin/bash

# Script name: rsync_to_server.sh

# Description: This script synchronises files from the directories R/ and clusterscripts/
# in the 'MoA_Prediction' project (Github) to the EMBL server (/home/dubois/MoA_Prediction)
# so that jobs can be executed on the cluster.

SOURCEFILES="clusterscripts/*.sh clusterscripts/*.R R/*.R data/*.RData data/*.rds"
TO_HOSTNAME="login.cluster.embl.de"
USER=$(whoami)
TARGETDIR="/home/dubois/MoA_Prediction"

# echo "Sending from ${FROM_HOSTNAME}@${USER} to ${TO_HOSTNAME}@${USER}"

# options:
# -v = increase verbosity
# -z = use compression
# --progress: show progress during transfer
# --relative: to preserve the full path (so that directories are reproduced at the target)

rsync -vz --progress --relative ${SOURCEFILES} ${USER}@${TO_HOSTNAME}:${TARGETDIR}
