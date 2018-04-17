#!/bin/bash

# Script name: rsync_to_server.sh

# Description: This script synchronises files from 'MoA_Prediction' to the EMBL
# server so that jobs can be executed on the cluster.

SOURCEFILES="yetanothertestfile.txt testfolder/*.txt testfolder2/*.txt"
TO_HOSTNAME="login.cluster.embl.de"
USER=$(whoami)
TARGETDIR="/home/dubois"

# echo "Sending from ${FROM_HOSTNAME}@${USER} to ${TO_HOSTNAME}@${USER}"

# options:
# -a = archive mode
# -v = increase verbosity
# -z = use compression
# --progress: show progress during transfer
# --relative: to preserve the full path (so that directories are reproduced at the target)

rsync -avz --progress --relative ${SOURCEFILES} ${USER}@${TO_HOSTNAME}:${TARGETDIR}
