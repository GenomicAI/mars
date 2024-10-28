# mars
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
|21| [nextflow](https://github.com/nextflow-io/nextflow)|24.09.2-edge|A workflow system for creating scalable, portable, and reproducible workflows|./mars-run.sh nextflow info|
