#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 20, 2024
# Description : Silmulate reads using wgsim, ngsngs or pbsim
# Usage : mars-reads.sh -f <sample fasta,fa file> -l <read length> -d <depth> -s <simulator>
# ./mars-reads.sh -f HG00096.fa -s n -l 150 -d 60

SHORT=f:,l:,d:,s:,w:,h
LONG=file:,length:,depth:,sim:,write:,help
OPTS=$(getopt -a -n mars-reads.sh --options $SHORT --longoptions $LONG -- "$@")

help_text="Usage: mars-reads.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa sequence file to read from\n"
help_text+="-l | --length INT read length (Default 100)\n"
help_text+="-d | --depth INT read coverage depth (Default 30)\n"
help_text+="-s | --sim STR read simulator. 'n' for NGSNGS, 'w' for wgsim,'p' pbsim and 'b' for badread (Default 'n')\n"
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
        -l | --length )
            length="$2"
            shift 2
        ;;
        -d | --depth )
            depth="$2"
            shift 2
        ;;
        -s | --sim )
            sim="$2"
            shift 2
        ;;
        -w | --write )
            write="$2"
            shift 2
        ;;
        -h | --help )
            echo "Program : mars-reads.sh"
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

mlog ">>> Starting mars-reads workflow on ${d} ..."
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

mlog ">>> Checking for sequence fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    mlog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

mlog ">>> Checking the index file for $(basename ${file}) ..."
if [ ! -f "${file}.fai" ]; then
    mlog ">>> Index file does not exists. Creating it..."
    singularity exec -e -B ${pwd} $sif samtools faidx $file;
fi

mlog ">>> Finding the length of the sample sequence file ..."
sample_len=$(singularity exec -e -B ${pwd} $sif awk 'BEGIN {t=0} {t+=$2} END {print t}' ${file}.fai)
mlog "Length is ${sample_len}"

re='^[0-9]+$'
if [ -z $depth ]; then
    depth=30;
    mlog "No coverage depth is specified. Assigning default 30";
fi
if ! [[ $depth =~ $re ]] ; then
    echo "--depth must be an integer !"; exit 1;
fi

if [ -z $length ]; then
    length=100;
    mlog "No reads length is specified. Assigning default 100";
fi

if ! [[ $length =~ $re ]] ; then
    echo "--length must be an integer !"; exit 1;
fi

if [ -z $sim ]; then
    sim="n";
    mlog "No reads simulator is specified. Assigning default NGSNGS";
fi

prefix=$(basename ${file%.*})
if [[ $sim == 'w' ]] ; then
    mlog ">>> Calculating number of reads required ..."
    reads=$(($sample_len*$depth/(2*$length))) # Divided by 2 is for paired read senario
    mlog "Need ${reads} reads for ${depth} coverage depth with ${length} reads";
    mlog ">>> Simulating reads with wgsim ..."
    singularity exec -e -B ${pwd} $sif wgsim -r 0 -e 0.001 -N $reads -1 $length -2 $length $file "${pwd}/${prefix}_reads_R1.fq" "${pwd}/${prefix}_reads_R2.fq";
    singularity exec -e -B ${pwd} $sif bgzip -f "${pwd}/${prefix}_reads_R1.fq";
    singularity exec -e -B ${pwd} $sif bgzip -f "${pwd}/${prefix}_reads_R2.fq";
elif [[ $sim == 'n' ]] ; then
    mlog ">>> Simulating reads with NGSNGS ..."
    cd $pwd;
    #ngsngs -i "../${sample}" -r $reads -l $length -seq PE -qs 20 -f fq.gz -o "${sample%.*}_reads";
    singularity exec -e -B ${pwd} $sif ngsngs -i $file -c $depth -l $length -seq PE -qs 30 -f fq.gz -o "${prefix}_reads";
elif [[ $sim == 'p' ]] ; then
    mlog ">>> Simulating long reads with pbsim (strategy wgs and method qshmm). Length parameter will be ignored ..."
    cd $pwd
    singularity exec -e -B ${pwd} $sif pbsim --strategy wgs --method qshmm --qshmm /opt/pbsim3/data/QSHMM-RSII.model --depth ${depth} --genome $file
    mv "${pwd}/sd_0001.fastq" "${pwd}/${prefix}_long_reads.fq"
    singularity exec -e -B ${pwd} $sif bgzip -f "${pwd}/${prefix}_long_reads.fq";
elif [[ $sim == 'b' ]] ; then
    mlog ">>> Simulating long reads with badread (Oxford Nanopore R10.4.1). Length parameter will be ignored ..."
    singularity exec -e -B ${pwd} $sif badread simulate --reference $file --quantity "${depth}x" > "${pwd}/${prefix}_long_reads.fq";
    singularity exec -e -B ${pwd} $sif bgzip -f "${pwd}/${prefix}_long_reads.fq";
else
    mlog "Unknown option \"${sim}\" for --sim!"
    echo -e $help_text
    exit 1;
fi

d=$(date)
mlog ">>> Done mars-reads on ${d} !";
