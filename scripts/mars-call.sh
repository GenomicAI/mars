#!/usr/bin/bash
# Author : Shanika Amarasoma, Nuzla Ismail
# Date : October 24, 2024
# Description : This command will Call the variants
# Usage : mars-call.sh -f <reference fasta,fa file> -g <reference graph gbz file> -m <bam or gam file> -t <no of threads> -c <variant caller>
# ./mars-call.sh -f NC_000020.11.fa -m HG00096.bam -c b -t 48
# ./mars-call.sh -g vgindex.giraffe.gbz -m HG00096.gam -c v -t 48

SHORT=f:,g:,c:,y:,m:,t:,w:,h
LONG=file:,gbz:,caller:,delly:,map:,threads:,write:,help
OPTS=$(getopt -a -n mars-call.sh --options $SHORT --longoptions $LONG -- "$@")

help_text="Usage: mars-call.sh [options]\n"
help_text+="-f | --file STR .fasta or .fa reference sequence file to map (If the mapper is not 'vg call')\n"
help_text+="-g | --gbz STR reference graph .gbz file (e.g. generated by mars-graph.sh) if the mapper is 'vg call' \n"
help_text+="-m | --map STR .bam or .gam (if the caller is 'vg call') sequence alignment/map file\n"
help_text+="-c | --caller STR variant caller to use. 'b' for 'bcftools', 'f' for 'freebayes', \n
'g' for 'gatk HaplotypeCaller', 'v' for 'vg call' and 'd' for 'delly' (SV detection only) . (Default 'b')\n"
help_text+="-y | --delly STR read simulation type, 's' for short read, 'p' for 'PacBio' and 'o' for 'ONT' (For delly only, Default 's')\n"
help_text+="-t | --threads INT number of threads to use (Default 'nproc')\n"
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
        -g | --gbz )
            gbz="$2"
            shift 2
        ;;
        -m | --map )
            map="$2"
            shift 2
        ;;
        -c | --caller )
            caller="$2"
            shift 2
        ;;
        -y | --delly )
            delly="$2"
            shift 2
        ;;
        -t | --threads )
            threads="$2"
            shift 2
        ;;
        -w | --write )
            write="$2"
            shift 2
        ;;
        -h | --help )
            echo "Program : mars-call.sh"
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

if [ -z $caller ]; then
    caller="b";
fi

if [ -z $delly ]; then
    delly="s";
fi

mlog ">>> Starting mars-call workflow on ${d} ..."
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

if [[ $caller == 'v' ]] ; then
    mlog "Variant caller is 'vg call'"
    mlog ">>> Checking for gbz file ..."
    
    if [ -z "$gbz" ] || [ ! -f "$gbz" ]; then
        mlog "The .gbz file ${gbz} does not exists or not specified by -g|--gbz <filename> !"
        echo -e $help_text
        exit 1;
    fi
    gbz=$(realpath $gbz)
    mlog ">>> Checking for gam file ..."
    if [ -z "$map" ] || [ ! -f "$map" ] || [ ${map##*.} != "gam" ]; then
        mlog "The .gam file ${map} does not exists, wrong or not specified by -m|--map <filename> !"
        echo -e $help_text
        exit 1;
    fi
    map=$(realpath $map)
else
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
    
    mlog ">>> Checking for bam file ..."
    if [ -z "$map" ] || [ ! -f "$map" ] || [ ${map##*.} != "bam" ]; then
        mlog "The .bam file ${map} does not exists, wrong or not specified by -m|--map <filename> !"
        echo -e $help_text
        exit 1;
    fi
    
    map=$(realpath $map)
    
    mlog ">>> Checking the index file for $(basename ${map}) ..."
    if [ ! -f "${map}.csi" ]; then
        mlog ">>> Index file does not exists. Creating it..."
        singularity exec -e -B ${pwd} $sif tabix -f $map;
    fi
fi

re='^[0-9]+$'
if [ -z $threads ]; then
    threads=$(nproc);
    mlog "No number of threads are specified. Assigning default nproc=${threads}";
fi
if ! [[ $threads =~ $re ]] ; then
    echo "--threads must be an integer !"; exit 1;
fi

prefix=$(basename ${map%.*})
outvcf="${pwd}/${prefix}.mars.${caller}.vcf"
if [[ $caller == 'b' ]] ; then
    mlog ">>> Variant call using bcftools call"
    singularity exec -e -B ${pwd} $sif bcftools mpileup -Ou -f $file $map | singularity exec -e -B ${pwd} $sif bcftools call -vmO z -o "${outvcf}.gz"
elif [[ $caller == 'f' ]] ; then
    mlog ">>> Variant call using freebayes call"
    singularity exec -e -B ${pwd} $sif freebayes -f $file $map | singularity exec -e -B ${pwd} $sif bgzip > "${outvcf}.gz"
elif [[ $caller == 'g' ]] ; then
    echo ">>> Variant call using gatk HaplotypeCaller call"
    if [ ! -f "${file%.*}.dict" ]; then
        singularity exec -e -B ${pwd} $sif gatk CreateSequenceDictionary -R "${file}" -O "${file%.*}.dict"
    fi
    singularity exec -e -B ${pwd} $sif gatk HaplotypeCaller -R $file -I $map -O "${outvcf}"
    singularity exec -e -B ${pwd} $sif bgzip -f "${outvcf}"
elif [[ $caller == 'd' ]] ; then
    mlog ">>> SV detect using delly"
    if [[ $delly == 'p' ]] ; then
        singularity exec -e -B ${pwd} $sif delly lr -y pb -g $file $map | singularity exec -e -B ${pwd} $sif bgzip > "${outvcf}.gz"
    elif [[ $delly == 'o' ]] ; then
        singularity exec -e -B ${pwd} $sif delly lr -y ont -g $file $map | singularity exec -e -B ${pwd} $sif bgzip > "${outvcf}.gz"
    else #Assum others are short reads
        singularity exec -e -B ${pwd} $sif delly call -g $file $map | singularity exec -e -B ${pwd} $sif bgzip > "${outvcf}.gz"
    fi
elif [[ $caller == 'v' ]] ; then
    mlog ">>> Variant call using vg call"
    singularity exec -e -B ${pwd} $sif vg pack -x $gbz -g $map -t $threads -o "${pwd}/${prefix}.pack"
    singularity exec -e -B ${pwd} $sif vg call $gbz -k "${pwd}/${prefix}.pack" -t $threads > "${outvcf}"
    singularity exec -e -B ${pwd} $sif bgzip -f "${outvcf}"
else
    mlog "Unknown option \"${caller}\" for --caller!"
    echo -e $help_text
    exit 1;
fi
singularity exec -e -B ${pwd} $sif bcftools index "${outvcf}.gz"

d=$(date)
mlog ">>> Done mars-call on ${d} !";
