#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 25, 2024
# Description : This command will compare 2 vcf files and generate a report
# Usage : mars-compare.sh -g <ground truth vcf file> -v <mars generated vcf file> -f <fa file of the sample>
# ./mars-compare.sh -g HG00096.vcf.gz -v HG00096.mars.d.vcf.gz -f HG00096.fa


SHORT=g:,v:,f:,w:,h
LONG=gtvcf:,vfvcf:,file:,write:,help
OPTS=$(getopt -a -n mars-compare.sh --options $SHORT --longoptions $LONG -- "$@")

 help_text="Usage: mars-compare.sh [options]\n"
help_text+="-g | --gtvcf STR ground truth vcf file file\n"
help_text+="-v | --vfvcf STR marser call generated vcf file.\n"
help_text+="-f | --file STR .fasta or .fa sequence file of the sample\n"
help_text+="-w | --write STR write logs to this file (optional, default 'mars.log')\n"
help_text+="-h | --help Display this help message\n"

eval set -- "$OPTS"
while :
do
  case "$1" in
	-g | --gtvcf )
      gtvcf="$2"
      shift 2
      ;;
	-v | --vfvcf )
      vfvcf="$2"
      shift 2
      ;;
	-f | --file )
      file="$2"
      shift 2
      ;;
	-w | --write )
      write="$2"
      shift 2
      ;;
    -h | --help )
      echo "Program : mars-compare.sh"
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

mlog ">>> Starting mars-compare workflow on ${d} ..."
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
mlog ">>> Checking for ground truth vcf file ..."

if [ -z "$gtvcf" ] || [ ! -f "$gtvcf" ]; then
    mlog "The ground truth vcf file ${gtvcf} does not exists or not specified by -g|--gtvcf <filename> !"
    echo -e $help_text
    exit 1;
fi

gtvcf=$(realpath $gtvcf)

mlog ">>> Checking for marser vcf file ..."

if [ -z "$vfvcf" ] || [ ! -f "$vfvcf" ]; then
    mlog "The marser vcf file ${vfvcf} does not exists or not specified by -v|--vfvcf <filename> !"
    echo -e $help_text
    exit 1;
fi

vfvcf=$(realpath $vfvcf)

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

mlog ">>> Finding ground truth SNPs and INDELs from ${vcf} ...."

snp=$(singularity exec -e -B ${pwd} $sif bcftools stats "${gtvcf}" | grep "number of SNPs:" | cut -f 4)
indel=$(singularity exec -e -B ${pwd} $sif bcftools stats "${gtvcf}" | grep "number of indels:" | cut -f 4)

mlog ">>> Comparing the ground truth vcf and mars vcf using bcftools isec ..."
singularity exec -e -B ${pwd} $sif bcftools isec -c none -p mars_compare "${gtvcf}" "${vfvcf}"

mlog ">>> Generating stats";
v_snp=$(singularity exec -e -B ${pwd} $sif bcftools stats "${vfvcf}" | grep "number of SNPs:" | cut -f 4)
v_indel=$(singularity exec -e -B ${pwd} $sif bcftools stats "${vfvcf}" | grep "number of indels:" | cut -f 4)
p_snp=$(singularity exec -e -B ${pwd} $sif bcftools stats mars_compare/0001.vcf | grep "number of SNPs:" | cut -f 4)
p_indel=$(singularity exec -e -B ${pwd} $sif bcftools stats mars_compare/0001.vcf | grep "number of indels:" | cut -f 4)
m_snp=$(singularity exec -e -B ${pwd} $sif bcftools stats mars_compare/0002.vcf | grep "number of SNPs:" | cut -f 4)
m_indel=$(singularity exec -e -B ${pwd} $sif bcftools stats mars_compare/0002.vcf | grep "number of indels:" | cut -f 4)

#INDEL Stats
indel_tp=$m_indel
indel_fp=$p_indel
indel_tn=$(($sample_len-$indel-$p_indel))
indel_fn=$(($indel-$indel_tp))
indel_sensitivity=$(bc <<< "scale=4; (${indel_tp}*100/(${indel_tp}+${indel_fn}));")
indel_specificity=$(bc <<< "scale=4; (${indel_tn}*100/(${indel_tn}+${indel_fp}));")
indel_f1=$(bc <<< "scale=4; (${indel_tp}*100/(${indel_tp}+(0.5*(${indel_fn}+${indel_fp}))));")

#SNP stats
snp_tp=$m_snp
snp_fp=$p_snp
snp_tn=$(($sample_len-$snp-$p_snp))
snp_fn=$(($snp-$snp_tp))
snp_sensitivity=$(bc <<< "scale=4; (${snp_tp}*100/(${snp_tp}+${snp_fn}));")
snp_specificity=$(bc <<< "scale=4; (${snp_tn}*100/(${snp_tn}+${snp_fp}));")
snp_f1=$(bc <<< "scale=4; (${snp_tp}*100/(${snp_tp}+(0.5*(${snp_fn}+${snp_fp}))));")

tp=$(($m_snp+$m_indel))
fp=$(($p_snp+$p_indel))
tn=$(($sample_len-$snp-$indel-$p_snp-$p_indel))
fn=$(($snp+$indel-$tp))
sensitivity=$(bc <<< "scale=4; (${tp}*100/(${tp}+${fn}));")
specificity=$(bc <<< "scale=4; (${tn}*100/(${tn}+${fp}));")
f1=$(bc <<< "scale=4; (${tp}*100/(${tp}+(0.5*(${fn}+${fp}))));")

mlog  "|  Description                      | Stats      |"
mlog  "|:----------------------------------|-----------:|"
mlog  "$(printf "|  Ground Truth SNPs                | %'10d |\n" ${snp})"                  
mlog  "$(printf "|  Ground Truth INDELs              | %'10d |\n" ${indel})"
mlog  " "                
mlog  "$(printf "|  mars SNPs                     | %'10d |\n" ${v_snp})"              
mlog  "$(printf "|  mars INDELs                   | %'10d |\n" ${v_indel})"
mlog  " "             
mlog  "$(printf "|  SNPs Private to mars vcf      | %'10d |\n" ${p_snp})"            
mlog  "$(printf "|  INDELs Private to mars vcf    | %'10d |\n" ${p_indel})"  
mlog  " "         
mlog  "$(printf "|  Exact Matched SNPs               | %'10d |\n" ${m_snp})"            
mlog  "$(printf "|  Exact Matched INDELs             | %'10d |\n" ${m_indel})"
mlog  " "             
mlog  "$(printf "|  True Positive (TP)               | %'10d |\n" ${tp})"
mlog  "$(printf "|  False Positive (FP)              | %'10d |\n" ${fp})"
mlog  "$(printf "|  True Negative (TN)               | %'10d |\n" ${tn})"
mlog  "$(printf "|  False Negative (FN)              | %'10d |\n" ${fn})"
mlog  " " 
mlog  "$(printf "|  SNP Sensitivity                  | %'9.4f%% |\n" ${snp_sensitivity})"
mlog  "$(printf "|  SNP Specificity                  | %'9.4f%% |\n" ${snp_specificity})"
mlog  "$(printf "|  SNP F1 Score                     | %'9.4f%% |\n" ${snp_f1})"
mlog  " " 
mlog  "$(printf "|  INDEL Sensitivity                | %'9.4f%% |\n" ${indel_sensitivity})"
mlog  "$(printf "|  INDEL Specificity                | %'9.4f%% |\n" ${indel_specificity})"
mlog  "$(printf "|  INDEL F1 Score                   | %'9.4f%% |\n" ${indel_f1})"
mlog  " "
mlog  "$(printf "|  Overall Sensitivity              | %'9.4f%% |\n" ${sensitivity})" 
mlog  "$(printf "|  Overall Specificity              | %'9.4f%% |\n" ${specificity})"
mlog  "$(printf "|  Overall F1 Score                 | %'9.4f%% |\n" ${f1})"

mlog  " "

mlog "Ground Truth SNPs,Ground Truth INDELs,Identified SNPs,Identified INDELs,Private SNPs,Private INDELs,Matched SNPs,Matched INDELs,TP,FP,TN,FN,SNP Sensitivity,SNP Specificity,SNP F1 Score,INDEL Sensitivity,INDEL Specificity,INDEL F1 Score,Overall Sensitivity,Overall Specificity,Overall F1 Score"

mlog "${snp},${indel},${v_snp},${v_indel},${p_snp},${p_indel},${m_snp},${m_indel},${tp},${fp},${tn},${fn},${snp_sensitivity},${snp_specificity},${snp_f1},${indel_sensitivity},${indel_specificity},${indel_f1},${sensitivity},${specificity},${f1}"

mlog  " "

d=$(date)
mlog ">>> Done mars-compare.sh on ${d} !";
