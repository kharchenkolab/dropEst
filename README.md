# dropEst - Pipeline
Pipeline for estimating molecular count matrices for droplet-based single-cell RNA-seq measurements. Implements methods, described in [this paper](https://doi.org/10.1101/171496).

## News
### [0.8.0] - 2018-03-20
**Important changes:**
* Algorithm of UMI correction now uses UMI quality. To provide this information you have to extract it on the dropTag 
phase with `-s` option and pass it to dropEst with `-r filename.gz`.
* Now, information about reads is kept in separate file instead (*.reads.gz), which should be passed to dropEst
* New format of *reads_per_umi_per_cell* in cell.counts.rds

See [CHANGELOG.md](CHANGELOG.md) for full list.

## Table of contents
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [dropEst - Pipeline](#dropest-pipeline)
    - [News](#news)
	- [Table of contents](#table-of-contents)
	- [General processing steps](#general-processing-steps)
	- [Setup](#setup)
		- [System requirements](#system-requirements)
		- [Installation](#installation)
		- [Troubleshooting](#troubleshooting)
		- [Dockers](#dockers)
	- [dropTag](#droptag)
		- [Protocols](#protocols)
			- [inDrop v1 & v2](#indrop-v1--v2)
			- [inDrop v3](#indrop-v3)
			- [10x](#10x)
			- [iCLIP](#iCLIP)
		- [Command line arguments for dropTag](#command-line-arguments-for-droptag)
	- [Alignment](#alignment)
	    - [Alignment with TopHat](#alignment-with-tophat)
	    - [Alignment with Kallisto](#alignment-with-kallisto)
	- [dropEst](#dropest)
		- [Usage of tagged bam files (e.g. 10x, Drop-seq) as input](#usage-of-tagged-bam-files-eg-10x-drop-seq-as-input)
		- [Usage of pseudoaligners](#usage-of-pseudoaligners)
		- [Count intronic / exonic reads only](#count-intronic--exonic-reads-only)
		- [Command line arguments for dropEst](#command-line-arguments-for-dropest)
		- [Output](#output)
	- [dropReport](#dropreport)
		- [Troubleshooting](#troubleshooting)
	- [dropEstR package](#dropestr-package)
	- [Additional notes](#additional-notes)

<!-- /TOC -->

## General processing steps
1. **dropTag**: extraction of cell barcodes and UMIs from the library. Result: demultiplexed .fastq.gz files, which should be aligned to the reference.
2. Alignment of the demultiplexed files to reference genome. Result: .bam files with the alignment.
3. **dropEst**: building count matrix and estimation of some statistics, necessary for quality control. Result: .rds file with the count matrix and statistics. *Optionally: count matrix in MatrixMarket format.*
4. **dropReport** - Generating report on library quality.
5. **dropEstR** - UMI count corrections, cell quality classification

## Setup
### System requirements
* Boost >= 1.54
* BamTools library
	* Note that some linux distributions have separate packages for the library and the executable, i.e. on Ubuntu you need `libbamtools-dev`, but not `bamtools`.
	* or you can [build it locally](https://github.com/pezmaster31/bamtools/wiki/Building-and-installing) and then specify the location of the build when running cmake (e.g. `cmake -D BAMTOOLS_ROOT=/home/username/bamtools .`)
* Zlib
* R >= 3.2.2 with packages:
  * Rcpp
  * RcppEigen
  * RInside
  * Matrix
* Compiler with c++11 support *(was tested with gcc >= 4.8.5 and CLang 3.9.1)*

### Installation
Install R packages:

```R
install.packages(c("Rcpp","RcppEigen", "RInside", "Matrix"))
```

Clone this repository:

```bash
git clone https://github.com/hms-dbmi/dropEst.git
```

Build:

```bash
cd dropEst
cmake . && make
```

### Troubleshooting
#### Local libraries installation
If `cmake` can't find one of the libraries, or you want to use some specific versions, which are currently not in the default path, use corresponding cmake variables:
* Boost: BOOST_ROOT.
* BamTools: BAMTOOLS_ROOT.
* R: R_ROOT. Can be found by running the `cat(R.home())` in R.

These variables should be set to the path to the installed library. It can be done either by using command line options: `cmake -D R_ROOT="path_to_r"` or by adding the variable declaration to the beginning of CMakeLists.txt: `set(R_ROOT path_to_r)`.


In case you have some issues with the linker for specific library, please build this library manually with the version of compiler, which you're going to use for dropEst build.

#### Problems with std::__cxx11::string
If you have messages like "*(path to some library)*: undefined reference to *(some name)* for std::__cxx11::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >", it means that you're trying to link a library, built with gcc < 5.0, while dropEst is built with gcc >= 5.0. It's a compiler issue, and you have to guarantee consistency of compiler versions by rebuilding either the library or dropEst. For more details see [question on stackoverflow](https://stackoverflow.com/questions/33394934/converting-std-cxx11string-to-stdstring).

If you have several compilers in your system, please use cmake flags `-DCMAKE_CXX_COMPILER=(c++ compiler)` and `-DCMAKE_C_COMPILER=(c compiler)` to choose a compiler. Here, `(c++ compiler)` and `(c compiler)` denotes path to the prefered compiler version.

#### Boost 1.65
CMake < 3.10 has known issues with boost 1.65. If you have such combination, please try either to upgrade cmake or to downgrade boost.

### Dockers
In case you still can't build the project, dockerfiles for the most popular linux distributions are provided (see `dropEst/dockers/`). 
You can either build and run these dockers or just read dockerfiles for the further instructions on dropEst installation for specific distribution.
Manual boost installation is shown in **CentOS 6** docker.

To install docker on your system see [installation instruction](https://github.com/wsargent/docker-cheat-sheet#installation). After installing the docker, use the following commands to start it:
```bash
cd dropEst/dockers/debian9
docker build -t dropest .
docker run --name dropest -it dropest
```

You can find more info about dockers at [Docker Cheat Sheet](https://github.com/wsargent/docker-cheat-sheet)

## dropTag
    droptag -- generate tagged fastq files for alignment

Example command:
```bash
./droptag -c config.xml [-S] [-s] reads1.fastq reads2.fastq [...]
```

Positional arguments of the dropTag phase contain paths to the read files, obtained with a sequencer. These files can be 
either in *.fastq* of *.fastq.gz* format. The file order depends on the type of used protocol. Below is the instruction 
on the usage for currently supported protocols.

**Note.** Algorithm of UMI correction requires information about base-call quality of UMIs. To save this information, use `-s` option.
In this case, in addition to fastq files with reads, the pipeline saves separate gzipped file with read parameters. These files 
must be passed to dropEst with `-r` option. 

### Protocols
#### inDrop v1 & v2
* File 1: barcode reads. Structure:
  * Cell barcode, part 1
  * Spacer
  * Cell barcode, part 2
  * UMI
* File 2: gene reads

Example config file is located at "*dropEst/configs/indrop_v1_2.xml*".  
Example command:
```bash
./droptag -c [-S] [-s] ./configs/indrop_v1_2.xml barcode_reads.fastq gene_reads.fastq
```

#### inDrop v3
* File 1: cell barcode, part 1 *(default length: 8bp)*
* File 2: cell barcode + UMI, part 1 *(default length: >= 14bp)*
* File 3: gene read
* *File 4 (optional): library tag*

If a file with library tags provided, option "-t" is required. ***This option wasn't tested properly, so it's better to avoid using it.***
<!-- To get data from multiple libraries **TODO: understand what happens with barcode in the case of multiple library tags**. -->

Example config file is located at "*dropEst/configs/indrop_v3.xml*".  
Example command:
```bash
./droptag -c dropEst/configs/indrop_v3.xml [-S] [-s] [-t library_tag] barcode1_reads.fastq barcode2_reads.fastq gene_reads.fastq [library_tags.fastq]
```

#### 10x
* File 1: library tag *(default length: 8bp)*
* File 2: cell barcode + UMI, part 1 *(default length: 16+10=26bp)*
* File 3: gene read

Example config file is located at "*dropEst/configs/10x.xml*".  
Example command:
```bash
./droptag -c dropEst/configs/10x.xml [-S] [-s] lib_tag_reads.fastq barcode_reads.fastq gene_reads.fastq
```

While dropTag provides way to demultiplex 10x data, [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) is still recommended tool for this. [dropEst](##dropEst) phase can be ran on the Cell Ranger demultiplexed .bam file to obtain data in the format, optimized for the subsequent analysis.

**NOTE.** Sometimes 10x CellRanger isn't able to determine gene, from which a read originated. In this cases it fills gene info with a list of possible genes, separated by semicolon (e.g. 'ENSG00000255508;ENSG00000254772'). These genes **must be filtered out** prior to further analysis:
```r
holder <- readRDS('./cell.counts.rds')
cm <- holder$cm[grep("^[^;]+$", rownames(holder$cm)),]
```

#### iCLIP
* File 1: Gene reads with barcodes at the beginning of the sequence

Example config file is located at "*dropEst/configs/iclip.xml*".  
Example command:
```bash
./droptag -c dropEst/configs/iclip.xml [-S] [-s] data.fastq
```
**NOTE.** Implementation of iCLIP wasn't tested properly. Please, be careful using it. Anyone who used it is very welcome to comment it either in Issues or by e-mail.


### Command line arguments for dropTag
*  -c, --config filename: xml file with droptag parameters  
*  -l, --log-prefix prefix: logs prefix
*  -n, --name name: alternative output base name
*  -p, --parallel number: number of threads (usage of more than 6 threads should lead to significant speed up)
*  -S, --save-stats : save stats to rds file. This data is used on the dropReport phase.
*  -t, --lib-tag library tag : (for IndropV3 with library tag only)
*  -q, --quiet : disable logs

Please, use `./droptag -h` for additional help.

## Alignment
dropTag writes the tagged reads into multiple files. All these files must be aligned to reference, and all bam files with the alignments must be provided as input for the dropEst stage. In the paper we used [TopHat2](https://ccb.jhu.edu/software/tophat/tutorial.shtml) aligner, however any RNA-seq aligners (i.e. [Kallisto](https://pachterlab.github.io/kallisto/) or [STAR](https://github.com/alexdobin/STAR)) can be used.

### Alignment with TopHat

1. Install [bowtie index](http://bowtie-bio.sourceforge.net/tutorial.shtml#preb) for the sequenced organism.
2. Download corresponding [gene annotation](http://genome.ucsc.edu/cgi-bin/hgTables) in the gtf format.
3. Run TopHat2 for each file:
```bash
tophat2 -p number_of_threads --no-coverage-search -g 1 -G genes.gtf -o output_dir Bowtie2Index/genome reads.fastq.gz
```
4. The result needed for the count estimation is *./output_dir/accepted_hits.bam*.

### Alignment with Kallisto

**Prior to v0.42.4, Kallisto didn't use BAM flags to mark primary/secondary alignments. Only versions $\geq$ 0.42.4 are supported.**

1. Download transcript sequences in .fasta format (i.e. [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) genes).
2. Build Kallisto index: `kallisto index -i genes.fa.gz`.
3. Run `kallisto quant --pseudobam --single -i genes.index -o out -l mean_length -s std_length reads.fastq.gz`. Here, *mean_length* is mean length of RNA fragment (not read length) and *std_length* is standard deviataion of RNA fragment length. You should specify values, according to the experiment design.

## dropEst

    dropest: estimate molecular counts per cell

This phase requires aligned .bam files as input and uses them to estimate count matrix. These files must contain information about cell barcode and UMI for each read (reads, which don't contain such information are ignored). Three possible ways to encode this information are acceptable:
1. Encode it in read names (as a result of dropTag phase).
2. Use .bam tags (i.e. output of 10x Cell Ranger). Tag names can be configured in [*config.xml* file](##additional-notes).
3. Save this information in a separate file using `-s` option on dropTag phase and pass it to dropEst with `-r` option.

**Note.** To run Bayesian algorithm of UMI error correction with dropestr, you need to pass information about UMI base call quality
to dropEst. You can do it either using .bam tags or `-s` dropTag option (i.e. option 1 isn't available in this case). You still can run
simpler algorithms (i.e. *cluster* or *directional*, see paper) without this information. 

Count matrix estimation also requires information about the source gene for the reads. It can be provided in two ways:
1. Use gene annotation in either *.bad* or *.gtf* format. To provide such file, "*-g*" option should be used.
2. Use .bam tags. Tag name can be configured in [*config.xml* file](##additional-notes).

Another crucial moment in estimation of count matrix is correction of cell barcode errors. Most protocols provide the list of real barcodes, which simplifies the task of correction. If such file is available, path to the file **should be specified in the *config.xml* file** (*Estimation/Merge/barcodes_file*). This can dramatically increase quality of the result. Lists for inDrop protocols can be found at *dropEst/data/barcodes/*. Two algorithms of barcode correction are available:
1. Simple, "*-m*" option. This algorithm is recommended in the case, where barcodes list is supplied.
2. Precise, "*-M*" option. Doesn't requires list of real barcodes to obtain high-quality results, however has a lower performance for large datasets.

Example command:
```bash
dropest [options] [-m] [-r pipeline_res.params.gz] -g ./hg38/genes.gtf -c ./config.xml ./alignment.*/accepted_hits.bam
```

### Usage of tagged bam files (e.g. 10x, Drop-seq) as input
Some protocols provide pipelines, which create .bam files with information about CB, UMI and gene.
To use these files as input, specify "*-f*" option. Example:
```bash
dropest [options] -f [-g ./genes.gtf] -c ./config.xml ./pipeline_res_*.bam
```

If "*-g*" option is provided, genes are parsed from the gtf, and information about genes from the bam file is ignored.

To specify corresponding .bam tag names, use "*Estimation/BamTags*" section in the config (see *configs/config_desc.xml*).  

### Usage of pseudoaligners
Pseudoaligners, such as Kallisto store gene / transcript names in the field of chromosome name. To parse such files, 
use "*-P*" option. Example:
 ```bash
 dropest [options] -P [-r pipeline_res.params.gz] -c ./config.xml ./kallisto_res_*.bam
 ```

### Count intronic / exonic reads only
One feature of the pipeline is the ability to count only UMIs, which reads touch only specific parts of a genome. 
Option *"-L"* allows to specify all acceptable types of regions:
* e: UMIs with exonic reads only
* i: UMIs with intronic reads only
* E: UMIs, which have both exonic and not annotated reads
* I: UMIs, which have both intronic and not annotated reads
* B: UMIs, which have both exonic and intronic reads
* A: UMIs, which have exonic, intronic and not annotated reads

Thus, to count all UMIs with exonic **or** not annotated reads, use *"-L eE"*. Default value: *"-L eEBA"*.

Example commands:
* Intronic reads only:
    ```bash
    dropest [-f] [-g ./genes.gtf] [-r pipeline_res.params.gz] -L i -c ./config.xml ./alignment_*.bam
    ```
* Exonic reads only:
    ```bash
    dropest [-f] [-g ./genes.gtf] [-r pipeline_res.params.gz] -L i -c ./config.xml ./alignment_*.bam
    ```
* Exon/intron spanning reads:
    ```bash
    dropest [-f] [-g ./genes.gtf] [-r pipeline_res.params.gz] -L BA -c ./config.xml ./alignment_*.bam
    ```

The pipeline can determine genome regions either using .gtf annotation file or using .bam tags, i.e. for CellRanger 
output (see *Estimation/BamTags/Type* in *configs/config_desc.xml*). If .gtf file isn't provided and .bam file doesn't containt 
annotation tags, all reads with not empty gene tag are considered as exonic. 

#### Velocyto integration
For some purposes (i.e. [velocyto](http://velocyto.org/)) it can be useful to look separately at the fraction of intronic and exonic UMIs.
Option *"-V"* allows to output three separate count matrices, each of which contains only UMIs of a specific type: 
intronic, exonic or exon/intron spanning. These matrices are stored in the separate file *"cell.counts.matrices.rds"*. 

### Command line arguments for dropEst
*  -b, --bam-output: print tagged bam files  
*  -c, --config filename: xml file with estimation parameters  
*  -C, --cells num: maximal number of output cells  
*  -f, --filled-bam: bam file already contains genes/barcodes tags  
*  -F, --filtered-bam: print tagged bam file after the merge and filtration  
*  -g, --genes filename: file with genes annotations (.bed or .gtf)  
*  -G, --genes-min num: minimal number of genes in output cells  
*  -l, --log-prefix : logs prefix  
*  -m, --merge-barcodes : merge linked cell tags  
*  -M, --merge-barcodes-precise : use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available  
*  -o, --output-file filename : output file name
*  -r, --read-params filenames: file or files with serialized params from tags search step. If there are several files, they should be provided in quotes, separated by space: "file1.reads.gz file2.reads.gz file3.reads.gz"  
*  -R, --reads-output: print count matrix for reads and don't use UMI statistics  
*  -q, --quiet : disable logs  
*  -V, --velocyto : save separate count matrices for exons, introns and exon/intron spanning reads
*  -w, --write-mtx : write out matrix in MatrixMarket format  

### Output
<!-- TODO: add that output has data of two types: all cells and filtered cells -->
Result of this phase is cell.counts.rds file with the next fields:
* **cm** (sparse matrix): count matrix in sparse format
* **reads_per_chr_per_cell** (list of data.frame): number of reads per cell (row) for each chromosome (column):
  * **Exon** (data.frame): exonic reads
  * **Intron** (data.frame): intronic reads
  * **Intergenic** (data.frame): intergenic reads
* **mean_reads_per_umi** (vector): mean number of reads per UMI for each cell
* some additional info <!-- TODO: describe it -->

To additionaly print the file with count matrix in MatrixMarket format use "*-w*" option.

## dropReport
To run the report you have to install:
* [dropestr](#dropestr-package) R package with all dependencies (`dependencies = T`).
* [pandoc](https://pandoc.org/installing.html)
* R packages:
    ```R
    install.packages(c("optparse","rmarkdown"))
    ```

To run the report use:
```bash
./dropReport.Rsc cell.counts.rds
```

To see full list of options use:
```bash
./dropReport.Rsc -h
```

### Troubleshooting
If you get the error *"pandoc version 1.12.3 or higher is required and was not found"*, try to set path in the corresponding environment variable:  `export RSTUDIO_PANDOC=/path/to/pandoc`.

## dropEstR package

This package implements UMI errors corrections and low-quality cells filtration, which are not performed during dropEst phase.

To install the package, use:

```r
devtools::install_github('hms-dbmi/dropEst/dropestr', dependencies = T)
```

Package content:
  * filtration of low-quality cells (see vignette "*low-quality-cells*")
  * correction of UMI errors (see vignette "*umi-correction*")
  * quality control (see *dropReport.Rsc*)  <!-- TODO: create vignette -->

### Troubleshooting
*dropestr* depends on *ks* package, which requires installed X Server. If you don't want to install it, use [this fork](https://github.com/VPetukhov/ks), which has almost the same functionality without X Server requirenment.

## Additional notes
Description of the fields for the config file is provided in *dropEst/configs/config_desc.xml*.
