# "mars"
### Simplifying Bioinformatics Workflows through a Containerized Approach to Tool Integration and Management
In the rapidly evolving field of bioinformatics, numerous specialized tools have been developed for essential genomic analysis tasks, such as read simulation, mapping, and variant calling. However, the management of these tools often presents significant challenges due to varied dependencies, execution steps, and output formats, complicating the installation and configuration processes. To address these issues, we introduce "mars," a bioinformatics solution encapsulated within a Singularity container that pre-loads a comprehensive suite of widely used genomic tools. Mars not only simplifies the installation of these tools but also automates critical workflow functions, including sequence sample preparation, read simulation, read mapping, variant calling, and result comparison. By streamlining the execution of these workflows, mars enables users to easily manage input-output formats and compare results across different tools, thereby enhancing reproducibility and efficiency. Furthermore, by providing a cohesive environment that integrates tool management with a flexible workflow interface, mars empowers researchers to focus on their analyses rather than the complexities of tool configuration. This integrated solution facilitates the testing of various combinations of tools and algorithms, enabling users to evaluate performance based on different metrics and identify the optimal tools for their specific genomic analysis needs. Through mars, we aim to enhance the accessibility and usability of bioinformatics tools, ultimately advancing research in genomic analysis.

## Installation
Since "mars" run on a [singularity](https://sylabs.io/singularity/) container, first we need to install the singularity. This [guide](https://docs.sylabs.io/guides/main/user-guide/quick_start.html) will explain how to install singularity.

Then pull the "mars" image.
```
singularity pull library://shanuz/genomicai/mars:latest
ls -ltrh
total 1.8G
-rwxr-xr-x 1 shanika shanika 1.8G Oct 27 20:21 mars_latest.sif
```
Copy all "mars' scripts from this [git](https://github.com/GenomicAI/mars/tree/main/scripts) to your environment. All the "mars" scripts assumes the `mars_latest.sif` is in the current folder or you can specify diffenet folder with env variable `MARSSIF`. Run the following command to see everything working fine. 
```
./mars-run.sh cat /etc/os-release
  PRETTY_NAME="Ubuntu 22.04.5 LTS"
  NAME="Ubuntu"
  VERSION_ID="22.04"
  VERSION="22.04.5 LTS (Jammy Jellyfish)"
  VERSION_CODENAME=jammy
  ID=ubuntu
  ID_LIKE=debian
  HOME_URL="https://www.ubuntu.com/"
  SUPPORT_URL="https://help.ubuntu.com/"
  BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
  PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
  UBUNTU_CODENAME=jammy
```

## Available tools in "mars"
|#|Tool|Version|Description|Execute|
|-:|----|-------|-----------|-------|
|1| [htslib, tabix, bgzip](https://github.com/samtools/htslib) |1.21|C library for high-throughput sequencing data formats|`./mars-run.sh tabix -h`, `./mars-run.sh bgzip -h`|
|2| [SAMtools](https://github.com/samtools/samtools) |1.21|Tools (written in C using htslib) for manipulating next-generation sequencing data|`./mars-run.sh samtools version`|
|3| [BCFtools](https://github.com/samtools/bcftools) |1.21|Utilities for variant calling and manipulating VCFs and BCFs.|`./mars-run.sh bcftools version`|
|4| [wgsim](https://github.com/lh3/wgsim) |1.21|Short read simulator|`./mars-run.sh wgsim -h`|
|5| [ngsngs](https://github.com/RAHenriksen/NGSNGS) |v0.9.2.2|Next Generation Simulator for Next Generation Sequencing Data|`./mars-run.sh ngsngs -v`|
|6| [bwa](https://github.com/lh3/bwa) |0.7.18-r1243-dirty|Burrow-Wheeler Aligner for short-read alignment|`./mars-run.sh bwa`|
|7| [bowtie2](https://github.com/BenLangmead/bowtie2)|2.5.4|A fast and sensitive gapped read aligner|`./mars-run.sh bowtie2 --version`|
|8| [freebayes](https://github.com/freebayes/freebayes) |v1.3.6|Bayesian haplotype-based genetic polymorphism discovery and genotyping.|`./mars-run.sh freebayes -h`|
|9| [vg](https://github.com/vgteam/vg) |v1.60.0 "Annicco"|Tools for working with genome variation graphs|`./mars-run.sh vg help`|
|10| [GATK](https://github.com/broadinstitute/gatk) |4.6.0.0|A wide variety of tools with a primary focus on variant discovery and genotypin|`./mars-run.sh gatk --list`|
|11| [pbsim](https://github.com/yukiteruono/pbsim3)|v3.0.4|A simulator for all types of PacBio and ONT long reads|`./mars-run.sh pbsim`|
|12| [minimap2](https://github.com/lh3/minimap2/)|2.28-r1209|A versatile pairwise aligner for genomic and spliced nucleotide sequences|`./mars-run.sh minimap2 --version`|
|13| [wfmash](https://github.com/waveygang/wfmash/)|v0.21.0|Base-accurate DNA sequence alignments using WFA and mashmap3|`./mars-run.sh wfmash -h`|
|14| [seqwish](https://github.com/ekg/seqwish)|v0.7.11-0-g0eb6468|Alignment to variation graph inducer|`./mars-run.sh seqwish -v`|
|15| [smoothxg](https://github.com/pangenome/smoothxg)|v0.8.0-0-g66b17ae|Linearize and simplify variation graphs using blocked partial order alignment|`./mars-run.sh smoothxg -v`|
|16| [gfaffix](https://github.com/marschall-lab/GFAffix)|0.1.5b|GFAffix identifies walk-preserving shared affixes in variation graphs and collapses them into a non-redundant graph structure|`./mars-run.sh gfaffix -V`|
|17| [odgi](https://github.com/pangenome/odgi)|v0.9.0-3-g237fc1b0|Optimized Dynamic Genome/Graph Implementation: understanding pangenome graphs|`./mars-run.sh odgi`|
|18| [badread](https://github.com/rrwick/Badread)|v0.4.1|A long read simulator that can imitate many types of read problems|`./mars-run.sh badread --version`|
|19| [delly](https://github.com/dellytools/delly)|1.2.9|Structural variant discovery by integrated paired-end and split-read analysis|`./mars-run.sh delly`|
|20| [diamond](https://github.com/bbuchfink/diamond)|v2.1.9.163|A sequence aligner for protein and translated DNA searches|`./mars-run.sh diamond help`|
|21| [nextflow](https://github.com/nextflow-io/nextflow)|24.09.2-edge|A workflow system for creating scalable, portable, and reproducible workflows|`./mars-run.sh nextflow info`|

## "mars" utilities
#### `mars-run.sh`
This works as a wrapper to run the command in the singulairy container. As explained in the above table, all the can executed through this. 

#### `mars-reads.sh`
All the read simulation tools can be executed through this utility in a convienient manner. However this supports limited number of options and you can use adbavced options through `mars-run.sh`.
```
./mars-reads.sh -h
Program : mars-reads.sh
Version : 1.0
Usage: mars-reads.sh [options]
-f | --file STR .fasta or .fa sequence file to read from
-l | --length INT read length (Default 100)
-d | --depth INT read coverage depth (Default 30)
-s | --sim STR read simulator. 'n' for 'NGSNGS', 'w' for 'wgsim','p' 'pbsim' and 'b' for 'badread' (Default 'n')
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

#### `mars-map.sh`
Most of the tools which are used for read mapping to the reference can be executed through this utility. For advanced options check the documentton of required tool and excute it through `mars-run.sh`.
```
./mars-map.sh -h
Program : mars-map.sh
Version : 1.0
Usage: mars-map.sh [options]
-f | --file STR .fasta or .fa reference sequence file to map (If the mapper is not 'vg giraffe')
-g | --gbz STR reference graph .gbz file (e.g. generated by mars-graph.sh) if the mapper is 'vg giraffe' 
-m | --mapper STR mapper/Aligner to use. 'm' for 'bwa mem', '2' for 'minimap2',
 'mp' for 'bwa mem -x pacbio', 'mo' for 'bwa mem -x ont2d' (For single ended reads only),
 's' for 'bwa sampe (Pared ended read), bwa samse (Single ended read)',
 '2p' for 'minimap2 -ax map-pb', '2o' for 'minimap2 -ax map-ont', '2i' for 'minimap2 -ax map-iclr' (Single ended read only),
 'b' for 'bowtie2' and 'g' for 'vg giraffe'. (Default 'm')
-1 | --read1 STR pared read .fastq or .fq file 1 
-2 | --read2 STR pared read .fastq or .fq file 2 
-t | --threads INT number of threads to use (Default 'nproc')
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

#### `mars-call.sh`
This utility will support to execute variant call tools in a convenient manner. However this also supports limited number of options and you can use adbavced options through `mars-run.sh` by reffering the documentations.
```
./mars-call.sh -h
Program : mars-call.sh
Version : 1.0
Usage: mars-call.sh [options]
-f | --file STR .fasta or .fa reference sequence file to map (If the mapper is not 'vg call')
-g | --gbz STR reference graph .gbz file (e.g. generated by mars-graph.sh) if the mapper is 'vg call' 
-m | --map STR .bam or .gam (if the caller is 'vg call') sequence alignment/map file
-c | --caller STR variant caller to use. 'b' for 'bcftools', 'f' for 'freebayes', 
 'g' for 'gatk HaplotypeCaller', 'd' for 'delly' and 'v' for 'vg call'. (Default 'b')
-t | --threads INT number of threads to use (Default 'nproc')
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

#### `mars-compare.sh`
This utility will compare 2 vcf files (ground truth vs mars generated) using `bcftools` and will produced comprehensive report on SNPs and INDELs. 
```
./mars-compare.sh -h
Program : mars-compare.sh
Version : 1.0
Usage: mars-compare.sh [options]
-g | --gvcf STR ground truth vcf file
-v | --mvcf STR mars-call generated vcf file.
-f | --file STR .fasta or .fa sequence file of the sample
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

#### `mars-filter.sh`, `mars-prepare.sh` and `mars-graph.sh`
These utilities are used for preparing input files for testing. A detailed example is decribed in the next section.
```
./mars-filter.sh -h
Program : mars-filter.sh
Version : 1.0
Usage: mars-filter.sh [options]
-f | --file STR .fasta or .fa sequence file
-n | --name STR Chromosome name to be filtered.
-r | --rename STR rename the chromosome in output sequence file (optional)
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```
```
./mars-prepare.sh -h
Program : mars-prepare.sh
Version : 1.0
Usage: mars-prepare.sh [options]
-f | --file STR .fasta or .fa reference sequence file
-v | --vcf STR ground truth VCF file
-s | --sample STR Sample name to be considered from the VCF file
-r | --region STR chromosome name and region in chr:from-to format (Optional)
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```
```
./mars-graph.sh -h
Program : mars-graph.sh
Version : 1.0
Usage: mars-prepare.sh [options]
-f | --file STR .fasta or .fa reference sequence file
-v | --vcf STR the main VCF file
-s | --sample STR Sample names (at least 4 other than seleted sample for the simulation) to be considered from the VCF file
-r | --region STR chromosome name and region in chr:from-to format (Optional)
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

#### `mars-pipe.sh`
This utility will execte `mars-reads.sh`, `mars-map.sh` and `mars-call.sh` in a pipeline. 
```
./mars-pipe.sh -h
Program : mars-pipe.sh
Version : 1.0
Usage: mars-pipe.sh [options]
-r | --ref STR reference .fasta sequence file or reference .gbz graph file
-f | --file STR .fasta file of the sample considered
-v | --vcf STR Ground truth vcf file of the sample considered
-l | --length INT read length (Default 100)
-d | --depth INT read coverage depth (Default 30)
-s | --sim STR read simulator. 'n' for 'NGSNGS', 'w' for 'wgsim','p' 'pbsim' and 'b' for badread (Default 'n')
-m | --mapper STR mapper/Aligner to use. 'm' for 'bwa mem', '2' for 'minimap2',
 'mp' for 'bwa mem -x pacbio', 'mo' for 'bwa mem -x ont2d' (For single ended reads only),
 's' for 'bwa sampe (Pared ended read), bwa samse (Single ended read)',
 '2p' for 'minimap2 -ax map-pb', '2o' for 'minimap2 -ax map-ont', '2i' for 'minimap2 -ax map-iclr' (Single ended read only),
 'b' for 'bowtie2' and 'g' for 'vg giraffe'. (Default 'm')
-c | --caller STR variant caller to use. 'b' for 'bcftools', 'f' for 'freebayes',
 'g' for 'gatk HaplotypeCaller', 'd' for 'delly' and 'v' for 'vg call'. (Default 'b' or 'v' if the --ref is .gbz graph file)
-t | --threads INT number of threads to use (Default 'nproc')
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

------
## Example Using Genome assembly GRCh38.p14 Chromosome 20  
### Preparing Input files

Here, we use reference [Genome assembly GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). Download it and get the file `GCF_000001405.40_GRCh38.p14_genomic.fna`. 
Filtering the chromosome 20 can be done using the following command. 

```bash
./mars-filter.sh -f GCF_000001405.40_GRCh38.p14_genomic.fna -n NC_000020.11 -r 20
```
This will produce the files `NC_000020.11.fa` and `NC_000020.11.fa.fai`
Next, we'll prepare a sequence file and a ground truth VCF file for a random sample (HG00096) of the [GRCh38](http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz) VCF file. Download that file and use the below command. We'll restrict the example for the region 30000000-32000000

```bash
./mars-prepare.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00096
```
This will produce the sequence file `HG00096.fa`, index `HG00096.fa.fai`, and the ground truth VCF file `HG00096.vcf.gz`.
If we need to consider the graph-based method as well, we have to prepare a reference graph using the below command. Here, we use five or more random samples other than the selected sample for the analysis. 

```bash
./mars-graph.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00097,HG00099,HG00100,HG00101,HG00102
```
This will produce `vg.vcf.gz`, its index file `vg.vcf.gz.tbi` and the graph file `vgindex.giraffe.gbz` with its other supporting files `vgindex.dist` and `vgindex.min`. 

### 1. Executing the workflow (wgsim -> bwa mem -> bcftools call)
First, we'll simulate the reading using the command `mars-reads.sh` with a read length of 100 and a coverage of 60.  

```bash
./mars-reads.sh -f HG00096.fa -s w -l 100 -d 60
```
This will produce paired read files named `HG00096_reads_R1.fq.gz` and `HG00096_reads_R2.fq.gz`. 

The next step is to map the reads to the reference sequence using the command `mars-map.sh`

```bash
./mars-map.sh -f NC_000020.11.fa -1 HG00096_reads_R1.fq.gz -2 HG00096_reads_R2.fq.gz -t 48 -m m
```
This will produce a bam file called `HG00096.bam`. 

Now we can do the variant calling using `mars-call.sh`. 

```bash
./mars-call.sh -f NC_000020.11.fa -m HG00096.bam -c b -t 48
```
This step will create a VCF file called `HG00096.mars.b.vcf.gz`. 

The final step is to compare the 2 VCF files `HG00096.mars.b.vcf.gz` and `HG00096.vcf.gz` using `mars-compare.sh`.

```bash
./mars-compare.sh -g HG00096.vcf.gz -v HG00096.mars.b.vcf.gz -f HG00096.fa
```
This command will produce the below report with all the stats. 

| Description | Stats |
|:----------------------------------|-----------:|
| Ground Truth SNPs | 1,886 |
| Ground Truth INDELs | 271 |
| mars SNPs | 2,001 |
| mars INDELs | 249 |
| SNPs Private to mars vcf | 121 |
| INDELs Private to mars vcf | 209 |
| Exact Matched SNPs | 1,880 |
| Exact Matched INDELs | 40 |
| True Positive (TP) | 1,920 |
| False Positive (FP) | 330 |
| True Negative (TN) | 1,997,225 |
| False Negative (FN) | 237 |
| SNP Sensitivity | 99.6818% |
| SNP Specificity | 99.9939% |
| SNP F1 Score | 96.7326% |
| INDEL Sensitivity | 14.7601% |
| INDEL Specificity | 99.9895% |
| INDEL F1 Score | 15.3846% |
| Overall Sensitivity | 89.0125% |
| Overall Specificity | 99.9834% |
| Overall F1 Score | 87.1341% |

### 2. Executing the graph-based workflow (wgsim -> vg giraffe -> vg call)

In this scenario, the mapping step and the calling step will be different below. 

```bash
./mars-map.sh -g vgindex.giraffe.gbz -1 HG00096_reads_R1.fq.gz -2 HG00096_reads_R2.fq.gz -t 48 -m g
```
This will produce a file called `HG00096.fa`, which can be used in the next step to call the variants using `mars-call.sh`.

```bash
./mars-call.sh -g vgindex.giraffe.gbz -m HG00096.gam -c v -t 48
```
```
./mars-compare.sh -g HG00096.vcf.gz -v HG00096.mars.v.vcf.gz -f HG00096.fa
```

| Description | Stats |
|:----------------------------------|-----------:|
| Ground Truth SNPs | 1,886 |
| Ground Truth INDELs | 271 |
| mars SNPs | 1,577 |
| mars INDELs | 267 |
| SNPs Private to mars vcf | 9 |
| INDELs Private to mars vcf | 26 |
| Exact Matched SNPs | 1,568 |
| Exact Matched INDELs | 241 |
| True Positive (TP) | 1,809 |
| False Positive (FP) | 35 |
| True Negative (TN) | 1,997,520 |
| False Negative (FN) | 348 |
| SNP Sensitivity | 83.1389% |
| SNP Specificity | 99.9995% |
| SNP F1 Score | 90.5573% |
| INDEL Sensitivity | 88.9298% |
| INDEL Specificity | 99.9986% |
| INDEL F1 Score | 89.5910% |
| Overall Sensitivity | 83.8664% |
| Overall Specificity | 99.9982% |
| Overall F1 Score | 90.4273% |

### 3. Executing the whole workflow using `mars-pipe.sh`
It is possible to execute the above-discussed workflow using the command `mars-pipe.sh`. e.g:

```
./mars-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c b -t 48
#For the graph approach
./mars-pipe.sh -r vgindex.giraffe.gbz -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m g -c v -t 48
```

### 4. Executing Workflow with NextFlow
Refer to the file [example.nf](https://github.com/GenomicAI/mars/blob/main/nextflow/example.nf).

```
./mars-run.sh nextflow run example.nf

 N E X T F L O W   ~  version 24.09.2-edge

Launching `example.nf` [marvelous_euler] DSL2 - revision: 6ecd4b5d6a

executor >  local (3)
[7e/b63d83] process > SIMULATE_READS [100%] 1 of 1 ✔
[91/98a530] process > ALIGN_READS    [100%] 1 of 1 ✔
[10/b0a259] process > CALL_VARIANTS  [100%] 1 of 1 ✔
Completed at: 30-Oct-2024 18:02:33
Duration    : 2m 41s
CPU hours   : (a few seconds)
Succeeded   : 3
```

## Building the Singulairty Image
Refer to the singularity build definition file [mars.def](https://github.com/GenomicAI/mars/blob/main/singularity/mars.def). 
```
sudo singularity build -F mars.sif mars.def
```

