#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 27, 2024
# Description : Wrapper to directly run the tools in singularity image. 

# exit when any command fails
set -e

#Get present working directory
pwd=$(pwd)

if [ -z "$MARSSIF" ] ; then
    sif="${pwd}/mars_latest.sif"
else
    sif="${MARSSIF%/}/mars_latest.sif"
fi
sif=$(realpath $sif);

if [ -z "$sif" ] ; then
    echo "The sif file ${file} does not exists. Please specify the path by MARSSIF env variable."
    exit 1;
fi

#Check for at least one argument
if [ -z "$1" ] ; then
    echo -e "Usage: mars-run.sh <command> [options]"
    echo -e "Please specify the command !"
    exit 1;
fi

singularity exec -e -B ${pwd} $sif "$@";