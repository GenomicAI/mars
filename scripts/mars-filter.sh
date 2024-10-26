#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 25, 2024
# Description : This command filter a specified chromosome from a fasta file
# Usage : mars-filter.sh -f <fasta,fa file> -n <chromosome name to filter> -r <rename chromosome (optional)>
# ./mars-filter.sh -f GRCh38_full_analysis_set_plus_decoy_hla.fa.gz -n chr20 -r 20
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
# GCF_000001405.40_GRCh38.p14_genomic.fna -> NC_000020.11 Homo sapiens chromosome 20, GRCh38.p14 Primary Assembly
# ./mars-filter.sh -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz -n NC_000020.11 -r 20

SHORT=f:,n:,r:,w:,h
LONG=file:,name:,rename:,write:,help
OPTS=$(getopt -a -n mars-filter.sh --options $SHORT --longoptions $LONG -- "$@")

help_text="Usage: mars-filter.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa sequence file\n"
help_text+="-n | --name STR Chromosome name to be filtered.\n"
help_text+="-r | --rename STR rename the chromosome in output sequence file (optional)\n"
help_text+="-w | --write STR write logs to this file (optional, default 'mars.log')\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
    case "$1" in
        -f | --file )
            file="$2"
            shift 2
        ;;
        -n | --name )
            name="$2"
            shift 2
        ;;
        -r | --rename )
            rename="$2"
            shift 2
        ;;
        -w | --write )
            write="$2"
            shift 2
        ;;
        -h | --help )
            echo "Program : mars-filter.sh"
            echo "Version : 1.0"
            echo -e $help_text
            exit 2
        ;;
        --)
            shift;
            break
        ;;
        *)
            echo "Unexpected option: $1";
            exit 2
        ;;
    esac
done

# exit when any command fails
set -e

#Get present working directory
pwd=$(pwd)
#Log file
if [ -z "$write" ] ; then
    write='mars.log'
    elif ! [[ $write =~ ^[0-9a-zA-Z._-]+$ ]]; then
    echo "Invalid log file name !"
    echo -e $help_text
    exit 1;
fi
log="${pwd}/${write}"

#Function to print and log messages
function mlog(){
    echo $1;
    echo $1 >> $log
}
mlog " "
#Get Date
d=$(date)

mlog ">>> Starting mars-filter workflow on ${d} ..."
mlog ">>> Checking for mars_latest.sif ..."
if [ -z "$MARSSIF" ] ; then
    sif="${pwd}/mars_latest.sif"
else
    sif="${MARSSIF%/}/mars_latest.sif"
fi
sif=$(realpath $sif);

if [ -z "$sif" ] ; then
    mlog "The sif file ${file} does not exists. Please specify the path by MARSSIF env variable."
    exit 1;
fi
mlog ">>> Checking for fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    mlog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

if [ -z $name ]; then
    mlog "Chromosome name must be specified by -n|--name <filename> !"
    echo -e $help_text
    exit 1;
fi

mlog ">>> Checking for index file ..."
if [ ! -f "${file}.fai" ]; then
    mlog ">>> Index file does not exists. Creating it..."
    singularity exec -e -B ${pwd} $sif samtools faidx $file;
fi

mlog ">>> Filtering the Chromosome and creating the file ${name}.fa ..."
if [ -z $rename ]; then
    singularity exec -e -B ${pwd} $sif samtools faidx $file $name | singularity exec -e -B ${pwd} $sif tr [:lower:] [:upper:] > "${pwd}/${name}.fa"
else
    singularity exec -e -B ${pwd} $sif samtools faidx $file $name | singularity exec -e -B ${pwd} $sif tr [:lower:] [:upper:] | singularity exec -e -B ${pwd} $sif sed "s/^>${name}/>${rename}/i" > "${pwd}/${name}.fa"
fi

mlog ">>> Creating the index for ${name}.fa ..."
singularity exec -e -B ${pwd} $sif samtools faidx "${pwd}/${name}.fa"
singularity exec -e -B ${pwd} $sif bwa index "${pwd}/${name}.fa";
d=$(date)
mlog ">>> Done mars-filter on ${d} !";
