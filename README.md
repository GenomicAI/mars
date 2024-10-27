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
|#|Tool|version|Description|Execute|
|-|----|-------|-----------|-------|
|1| [htslib, tabix, bgzip](https://github.com/samtools/htslib) |(Tested with htslib 1.17/18)|d||x|
|2| [SAMtools](https://github.com/samtools/samtools) |(Tested with Version: 1.17/18 (using htslib 1.17/18))|d||x|
|3| [BCFtools](https://github.com/samtools/bcftools) |(Tested with Version: 1.17/18 (using htslib 1.17/18))|d||x|
|4| [wgsim](https://github.com/lh3/wgsim) |(Tested with Version: 1.17/18)|d||x|
|5| [ngsngs](https://github.com/RAHenriksen/NGSNGS) |(Tested with current master branch 2023-Oct-18)|d||x|
|6| [bwa](https://github.com/lh3/bwa) |(Tested with version 0.7.17-r1198-dirty)|d||x|
|7| [bowtie2](https://github.com/BenLangmead/bowtie2)| (Tested with version 2.5.1)|d||x|
|8| [freebayes](https://github.com/freebayes/freebayes) |(Tested with version 1.3.7, master branch)|d||x|
|9| [vg](https://github.com/vgteam/vg) |(Tested with version v1.51.0 "Quellenhof")|d||x|
|10| [GATK](https://github.com/broadinstitute/gatk) |(Tested version singularity image [gatk_4.4.0.0.sif](https://hub.docker.com/r/broadinstitute/gatk))|d||x|
|11| [pbsim](https://github.com/yukiteruono/pbsim3)|v|d||x|
|12| minimap2|v|d||x|
|13| wfmash|v|d||x|
|14| seqwish|v|d||x|
|15| smoothxg|v|d||x|
|16| gfaffix|v|d||x|
|17| odgi|v|d||x|
|18| badread|v|d||x|
|19| delly|v|d||x|
|20| [diamond](https://github.com/bbuchfink/diamond)|v|d||x|
|21| nextflow||v|d||x|
