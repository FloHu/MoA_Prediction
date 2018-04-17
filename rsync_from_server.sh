#!/bin/bash

# Script name: rsync_to_server.sh

# Description: This script synchronises files from 'MoA_Prediction' to the EMBL
# server so that jobs can be executed on the cluster.

SOURCEFILES="/home/dubois/results/*.txt"
SOURCEHOST="login.cluster.embl.de"
USER=$(whoami)
TARGETDIR="./results_from_server/"

# echo "Sending from ${FROM_HOSTNAME}@${USER} to ${TO_HOSTNAME}@${USER}"

# options:
# -a = archive mode
# -v = increase verbosity
# -z = use compression
# --progress: show progress during transfer
# --relative: to preserve the full path (so that directories are reproduced at the target)

rsync -avz --progress ${USER}@${SOURCEHOST}:${SOURCEFILES} ${TARGETDIR}
