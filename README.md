# "Mars"
### Simplifying Bioinformatics Workflows through a Containerized Approach to Tool Integration and Management
In the rapidly advancing field of bioinformatics, numerous specialized tools have been developed for critical genomic analysis tasks such as read simulation, mapping, and variant calling. However, managing these tools often presents significant challenges due to varied dependencies, execution steps, and output formats, making installation and configuration complex. To address these issues, we introduce "mars," a bioinformatics solution packaged in a Singularity container that preloads a comprehensive suite of widely used genomic tools.

Mars simplifies tool installation and automates key workflow functions, including sequence sample preparation, read simulation, read mapping, variant calling, and result comparison. By streamlining these workflows, mars enables users to easily manage input-output formats and compare results across different tools, enhancing reproducibility and efficiency. Additionally, by providing an integrated environment that combines tool management with a flexible workflow interface, mars empowers researchers to focus on analysis rather than tool configuration. This cohesive solution allows users to test different combinations of tools and algorithms, evaluate performance based on various metrics, and identify the best tools for their specific genomic analysis needs. Through mars, we aim to make bioinformatics tools more accessible and user-friendly, advancing research in genomic analysis.

## Installation
Since "mars" runs on a [Singularity](https://sylabs.io/singularity/) container, you’ll need to install Singularity first. This guide explains how to install it.

Afterward, pull the "mars" image:
```
singularity pull library://shanuz/genomicai/mars:latest
ls -ltrh
total 1.8G
-rwxr-xr-x 1 shanika shanika 1.8G Oct 27 20:21 mars_latest.sif
```
Copy all "mars" scripts from the [repository]((https://github.com/GenomicAI/mars/tree/main/scripts)) to your environment. All "mars" scripts assume that `mars_latest.sif` is in the current directory, or you can specify a different directory with the environment variable `MARSSIF`. Run the following command to ensure everything is working:

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
|4| [BEDtools](https://github.com/arq5x/bedtools2) |v2.31.1|The swiss army knife for genome arithmetic|`./mars-run.sh bedtools --version`|
|5| [wgsim](https://github.com/lh3/wgsim) |1.21|Short read simulator|`./mars-run.sh wgsim -h`|
|6| [ngsngs](https://github.com/RAHenriksen/NGSNGS) |v0.9.2.2|Next Generation Simulator for Next Generation Sequencing Data|`./mars-run.sh ngsngs -v`|
|7| [bwa](https://github.com/lh3/bwa) |0.7.18-r1243-dirty|Burrow-Wheeler Aligner for short-read alignment|`./mars-run.sh bwa`|
|8| [bowtie2](https://github.com/BenLangmead/bowtie2)|2.5.4|A fast and sensitive gapped read aligner|`./mars-run.sh bowtie2 --version`|
|9| [freebayes](https://github.com/freebayes/freebayes) |v1.3.6|Bayesian haplotype-based genetic polymorphism discovery and genotyping.|`./mars-run.sh freebayes -h`|
|10| [vg](https://github.com/vgteam/vg) |v1.60.0 "Annicco"|Tools for working with genome variation graphs|`./mars-run.sh vg help`|
|11| [GATK](https://github.com/broadinstitute/gatk) |4.6.0.0|A wide variety of tools with a primary focus on variant discovery and genotypin|`./mars-run.sh gatk --list`|
|12| [pbsim](https://github.com/yukiteruono/pbsim3)|v3.0.4|A simulator for all types of PacBio and ONT long reads|`./mars-run.sh pbsim`|
|13| [minimap2](https://github.com/lh3/minimap2/)|2.28-r1209|A versatile pairwise aligner for genomic and spliced nucleotide sequences|`./mars-run.sh minimap2 --version`|
|14| [wfmash](https://github.com/waveygang/wfmash/)|v0.21.0|Base-accurate DNA sequence alignments using WFA and mashmap3|`./mars-run.sh wfmash -h`|
|15| [seqwish](https://github.com/ekg/seqwish)|v0.7.11-0-g0eb6468|Alignment to variation graph inducer|`./mars-run.sh seqwish -v`|
|16| [smoothxg](https://github.com/pangenome/smoothxg)|v0.8.0-0-g66b17ae|Linearize and simplify variation graphs using blocked partial order alignment|`./mars-run.sh smoothxg -v`|
|17| [gfaffix](https://github.com/marschall-lab/GFAffix)|0.1.5b|GFAffix identifies walk-preserving shared affixes in variation graphs and collapses them into a non-redundant graph structure|`./mars-run.sh gfaffix -V`|
|18| [odgi](https://github.com/pangenome/odgi)|v0.9.0-3-g237fc1b0|Optimized Dynamic Genome/Graph Implementation: understanding pangenome graphs|`./mars-run.sh odgi`|
|19| [badread](https://github.com/rrwick/Badread)|v0.4.1|A long read simulator that can imitate many types of read problems|`./mars-run.sh badread --version`|
|20| [delly](https://github.com/dellytools/delly)|1.2.9|Structural variant discovery by integrated paired-end and split-read analysis|`./mars-run.sh delly`|
|21| [diamond](https://github.com/bbuchfink/diamond)|v2.1.9.163|A sequence aligner for protein and translated DNA searches|`./mars-run.sh diamond help`|
|22| [nextflow](https://github.com/nextflow-io/nextflow)|24.09.2-edge|A workflow system for creating scalable, portable, and reproducible workflows|`./mars-run.sh nextflow info`|

## "mars" utilities
#### `mars-run.sh`
This wrapper allows commands to be run within the Singularity container. Refer to the table above for examples of how each tool can be executed.

#### `mars-reads.sh`
This utility allows the execution of read simulation tools in a streamlined way. Note that this script supports a limited set of options; for advanced options, use `mars-run.sh`.
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
Most tools for mapping reads to reference sequences can be executed using this utility. For more options, consult the required tool’s documentation and execute it with `mars-run.sh`.

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
This utility facilitates running variant-calling tools. It also supports limited options; use `mars-run.sh` for advanced options, following the tool’s documentation.

```
 ./mars-call.sh -h
Program : mars-call.sh
Version : 1.0
Usage: mars-call.sh [options]
-f | --file STR .fasta or .fa reference sequence file to map (If the mapper is not 'vg call')
-g | --gbz STR reference graph .gbz file (e.g. generated by mars-graph.sh) if the mapper is 'vg call' 
-m | --map STR .bam or .gam (if the caller is 'vg call') sequence alignment/map file
-c | --caller STR variant caller to use. 'b' for 'bcftools', 'f' for 'freebayes', 
 'g' for 'gatk HaplotypeCaller', 'v' for 'vg call' and 'd' for 'delly' (SV detection only) . (Default 'b')
-y | --delly STR read simulation type, 's' for short read, 'p' for 'PacBio' and 'o' for 'ONT' (For delly only, Default 's')
-t | --threads INT number of threads to use (Default 'nproc')
-w | --write STR write logs to this file (optional, default 'mars.log')
-h | --help Display this help message
```

#### `mars-compare.sh`
This utility compares two VCF files (ground truth vs. mars-generated) using `bcftools` and produces a comprehensive report on SNPs and INDELs.
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
These utilities help prepare input files for testing. Detailed examples are provided in the following section.
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
This utility exectes `mars-reads.sh`, `mars-map.sh` and `mars-call.sh` in a pipeline. 
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

In this example, we use the reference [Genome assembly GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). Download the file `GCF_000001405.40_GRCh38.p14_genomic.fna`. To filter for chromosome 20, use the following command:

```bash
./mars-filter.sh -f GCF_000001405.40_GRCh38.p14_genomic.fna -n NC_000020.11 -r 20
```
This command produces the files `NC_000020.11.fa` and `NC_000020.11.fa.fai`.

Next, prepare a sequence file and a ground truth VCF file for a sample (HG00096) from the [GRCh38](http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz) VCF file. Download the file and use the command below, restricting the example to the region 30000000-32000000:

```bash
./mars-prepare.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00096
```
This command generates the sequence file `HG00096.fa`, its index `HG00096.fa.fai`, and the ground truth VCF file `HG00096.vcf.gz`.
If using a graph-based method, prepare a reference graph with the following command. Here, select five or more random samples other than the sample chosen for analysis:

```bash
./mars-graph.sh -r 20:30000000-32000000 -f NC_000020.11.fa -v ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -s HG00097,HG00099,HG00100,HG00101,HG00102
```
This produces `vg.vcf.gz`, its index file `vg.vcf.gz.tbi`, and the graph file `vgindex.giraffe.gbz` with supporting files `vgindex.dist` and `vgindex.min`.

### 1. Executing the workflow (wgsim -> bwa mem -> bcftools call)
First, simulate the reads using the `mars-reads.sh` command with a read length of 100 and a coverage of 60:

```bash
./mars-reads.sh -f HG00096.fa -s w -l 100 -d 60
```
This command generates paired read files named `HG00096_reads_R1.fq.gz` and `HG00096_reads_R2.fq.gz`.

Next, map the reads to the reference sequence using `mars-map.sh`:

```bash
./mars-map.sh -f NC_000020.11.fa -1 HG00096_reads_R1.fq.gz -2 HG00096_reads_R2.fq.gz -t 48 -m m
```
This produces a BAM file named `HG00096.bam`.

Now, perform variant calling with `mars-call.sh`:

```bash
./mars-call.sh -f NC_000020.11.fa -m HG00096.bam -c b -t 48
```
This step generates a VCF file called `HG00096.mars.b.vcf.gz`.

Finally, compare the two VCF files (`HG00096.mars.b.vcf.gz` and `HG00096.vcf.gz`) using `mars-compare.sh`:

```bash
./mars-compare.sh -g HG00096.vcf.gz -v HG00096.mars.b.vcf.gz -f HG00096.fa
```
This command outputs a report with statistics:

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

In the graph-based workflow, the mapping and calling steps differ. First, map the reads as follows:

```bash
./mars-map.sh -g vgindex.giraffe.gbz -1 HG00096_reads_R1.fq.gz -2 HG00096_reads_R2.fq.gz -t 48 -m g
```
This step produces a file `HG00096.fa`, which can be used in the next step to call variants using `mars-call.sh`:

```bash
./mars-call.sh -g vgindex.giraffe.gbz -m HG00096.gam -c v -t 48
```
Finally, compare the results with:
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

### 3. Executing the Entire Workflow with `mars-pipe.sh`
To run the complete workflow, use `mars-pipe.sh` as shown below:

```
./mars-pipe.sh -r NC_000020.11.fa -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s w -m m -c b -t 48
#For the graph approach
./mars-pipe.sh -r vgindex.giraffe.gbz -f HG00096.fa -v HG00096.vcf.gz -l 100 -d 60 -s n -m g -c v -t 48
```

### 4. Executing Workflow with NextFlow
To run the workflow in NextFlow, refer to [example.nf](https://github.com/GenomicAI/mars/blob/main/nextflow/example.nf).

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

## Structural Variant (SV) Detection Workflow.
For this experiment please use the data available in [svdata](/svdata/) folder. 

Ground truth VCF file with only one SV. 
```
less sample_sv.vcf.gz
```
![image](https://github.com/user-attachments/assets/1c6e77b5-9f88-4c1f-8d89-acfd63191652)

#### Workflow (pbsim -> minimap2 -ax map-pb -> delly lr -y pb) 
```
./mars-pipe.sh -r NC_000020.11.fa -f sample_sv.fa -s p -m 2p -c d
```
Output VCF file with detected SV 
```
less sample_sv.mars.d.vcf.gz
```
![image](https://github.com/user-attachments/assets/d2bc253c-c6fa-4efb-931b-c87855561233)

### Workflow (badread (ont) -> minimap2 -ax map-ont -> delly lr -y ont) 
```
./mars-pipe.sh -r NC_000020.11.fa -f sample_sv.fa -s b -m 2o -c d
```
Output VCF file with detected SV

![image](https://github.com/user-attachments/assets/3d7d4929-e31f-4019-ab0d-9b5c8ea3934d)

## Appendix
### How to Build the Singulairty Image
Refer to the singularity build definition file [mars.def](https://github.com/GenomicAI/mars/blob/main/singularity/mars.def). 
```
sudo singularity build -F mars.sif mars.def
```

