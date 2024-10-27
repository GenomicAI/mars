#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 25, 2024
# Description : Construct a graph using the reference fasta file and givem samples. 
# Usage : mars-graph.sh -f <reference fasta,fa file> -v <VCF file> -s <comma seperated sample names> -r <region>
# ./mars-graph.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00097,HG00099,HG00100,HG00101,HG00102

SHORT=f:,v:,s:,r:,w:,h
LONG=file:,name:,rename:,write:,help
OPTS=$(getopt -a -n mars-graph.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: mars-prepare.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa reference sequence file\n"
help_text+="-v | --vcf STR the main VCF file\n"
help_text+="-s | --sample STR Sample names (at least 4 other than seleted sample for the simulation) to be considered from the VCF file\n"
help_text+="-r | --region STR chromosome name and region in chr:from-to format (Optional)\n"
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
	-v | --vcf )
      vcf="$2"
      shift 2
      ;;
	-s | --sample )
      sample="$2"
      shift 2
      ;;
	-r | --region )
      region="$2"
      shift 2
      ;;
	-w | --write )
      write="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : mars-graph.sh"
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

mlog ">>> Starting mars-graph workflow on ${d} ..."
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
mlog ">>> Checking for reference fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    mlog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

mlog ">>> Checking for reference fasta file ..."

if [ -z "$file" ] || [ ! -f "$file" ]; then
    mlog "The fasta file ${file} does not exists or not specified by -f|--file <filename> !"
    echo -e $help_text
    exit 1;
fi

file=$(realpath $file)

if [ -z "$vcf" ] || [ ! -f "$vcf" ]; then
    mlog "The vcf file ${vcf} does not exists or not specified by -v|--vcf <filename> !"
    echo -e $help_text
    exit 1;
fi

vcf=$(realpath $vcf)

if [ -z $sample ]; then
    mlog "Sample name must be specified by -s|--sample <filename> !"
    echo -e $help_text
    exit 1;
fi

mlog ">>> Checking the index file for $(basename ${file}) ..."
if [ ! -f "${file}.fai" ]; then
	mlog ">>> Index file does not exists. Creating it..."
	samtools faidx $file;
fi

mlog ">>> Checking the index file for $(basename ${vcf}) ..."
if [ ! -f "${vcf}.csi" ]; then
	mlog ">>> Index file does not exists. Creating it..."
	singularity exec -e -B ${pwd} $sif bcftools index $vcf;
fi

mlog ">>> Generating the vcf file for the samples ${sample} ..."
if [ -z $region ]; then
	singularity exec -e -B ${pwd} $sif bcftools view --min-ac=1 -s $sample $vcf | singularity exec -e -B ${pwd} $sif bcftools norm -d all - | singularity exec -e -B ${pwd} $sif bgzip > "${pwd}/vg.vcf.gz"
else
	singularity exec -e -B ${pwd} $sif bcftools view -r $region --min-ac=1 -s $sample $vcf | singularity exec -e -B ${pwd} $sif bcftools norm -d all - | singularity exec -e -B ${pwd} $sif bgzip > "${pwd}/vg.vcf.gz"
fi

mlog ">>> construct the graph and indexes ..."
singularity exec -e -B ${pwd} $sif tabix "${pwd}/vg.vcf.gz"
singularity exec -e -B ${pwd} $sif vg autoindex --workflow giraffe -r $file -v "${pwd}/vg.vcf.gz" -p vgindex

d=$(date)
mlog ">>> Done mars-graph on ${d} !";