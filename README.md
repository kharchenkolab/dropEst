# dropEst - Pipeline
Pipeline for estimating molecular count matrices for droplet-based single-cell RNA-seq measurements. Implements methods, described at [this paper](https://doi.org/10.1101/171496).

## Table of contents
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [dropEst - Pipeline](#dropest-pipeline)
	- [Table of contents](#table-of-contents)
	- [General processing steps](#general-processing-steps)
	- [Setup](#setup)
		- [System requirements](#system-requirements)
		- [Installation](#installation)
		- [Troubleshooting](#troubleshooting)
	- [dropTag](#droptag)
		- [Protocols](#protocols)
			- [inDrop v1 & v2](#indrop-v1-v2)
			- [inDrop v3](#indrop-v3)
			- [10x](#10x)
		- [Command line arguments for dropTag](#command-line-arguments-for-droptag)
	- [Alignment](#alignment)
	- [dropEst](#dropest)
		- [Command line arguments for dropEst](#command-line-arguments-for-dropest)
		- [Output](#output)
	- [dropReport](#dropreport)
	- [dropEstR package](#dropestr-package)
	- [Additional notes](#additional-notes)

<!-- /TOC -->


## General processing steps
1. **dropTag**: extraction of cell barcodes and UMIs from the library. Result: demultiplexed .fastq.gz files, which should be aligned to the reference.
2. Alignment of the demultiplexed files to reference genome. Result: .bam files with the alignment.
3. **dropEst**: building count matrix and estimation of some statistics, necessary for quality control. Result: .rds file with the count matrix and statistics. *Optionally: count matrix in MatrixMarket format.*
4. **dropReport** - Generating report on library quality.

## Setup
### System requirements
* Boost >= 1.54
* BamTools
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
If `cmake` can't find one of the libraries, or you want to use some specific versions, which are currently not in the default path, use corresponding cmake variables:
* Boost: BOOST_ROOT.
* BamTools: BAMTOOLS_ROOT.
* R: R_ROOT. Can be found by running the `cat(R.home())` in R.

These variables should be set to the path to the installed library. It can be done either by using command line options: `cmake -D R_ROOT="path_to_r"` or by adding the variable declaration to CMakeLists.txt: `set(R_ROOT path_to_r)`.


In case you have some issues with the linker for specific library, please try to build this library manually with the version of compiler, which you're going to use for dropEst build.

## dropTag
    droptag -- generate tagged fastq files for alignment

Example command:
```bash
./droptag -c config.xml [-S] reads1.fastq reads2.fastq [...]
```

Positional arguments of the dropTag phase contain paths to the read files, obtained with a sequencer. These files can be either in *.fastq* of *.fastq.gz* format. The file order depends on the type of used protocol. Below is the instruction on the usage for currently supported protocols.

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
./droptag -c ./configs/indrop_v1_2.xml barcode_reads.fastq gene_reads.fastq
```

#### inDrop v3
* Reads to provide: (can be fastq or fastq.gz)
  * Read 1: Cell barcode part 1 (length 8)
  * Read 2: Cell barcode part 2 + UMI (length 14)
  * Read 3: Gene read (length 61)
  * Read 4: Library barcode (length 8) *(optional)*
* File 1: cell barcode, part 1 *(default length: 8bp)*
* File 2: cell barcode + UMI, part 1 *(default length: >= 14bp)*
* File 3: gene read
* *File 4 (optional): library tag*

If a file with library tags provided, option "-t" is required.
<!-- To get data from multiple libraries **TODO: understand what happens with barcode in the case of multiple library tags**. -->

Example config file is located at "*dropEst/configs/indrop_v1_2.xml*".  
Example command:
```bash
./droptag -c dropEst/configs/indrop_v3.xml [-S] [-t library_tag] barcode1_reads.fastq barcode2_reads.fastq gene_reads.fastq [library_tags.fastq]
```

#### 10x
[Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) is recommended tool for demultiplexing of 10x data. However, [dropEst](##dropEst) can be ran on the demultiplexed .bam file to obtain data in the format, optimized for the subsequent analysis.

### Command line arguments for dropTag
*  -c, --config filename: xml file with droptag parameters  
*  -l, --log-prefix prefix: logs prefix  
*  -n, --name name: alternative output base name  
*  -S, --save-stats : save stats to rds file. This data is used on the dropReport phase.
*  -t, --lib-tag library tag : (for IndropV3 with library tag only)  
*  -q, --quiet : disable logs  

Please, use `./droptag -h` for additional help.

## Alignment
dropTag writes the tagged reads into multiple files. All these files must be aligned to reference, and all bam files with the alignments must be provided as input for the dropEst stage. In the paper we used [TopHat2](https://ccb.jhu.edu/software/tophat/tutorial.shtml) aligner, however any RNA-seq aligners (i.e. [Kallisto](https://pachterlab.github.io/kallisto/) or [STAR](https://github.com/alexdobin/STAR)) can be used.

Alignment with TopHat2:
1. Install [bowtie index](http://bowtie-bio.sourceforge.net/tutorial.shtml#preb) for the sequenced organism.
2. Download corresponding [gene annotation](http://genome.ucsc.edu/cgi-bin/hgTables) in the gtf format.
3. Run TopHat2 for each file:
```bash
tophat2 -p number_of_threads --no-coverage-search -g 1 -G genes.gtf -o output_dir Bowtie2Index/genome reads.fastq.gz
```
4. The result needed for the count estimation is *./output_dir/accepted_hits.bam*.

## dropEst

    dropest: estimate molecular counts per cell

This phase requires aligned .bam files as input and uses them to estimate count matrix. These files must contain information about cell barcode and UMI for each read (reads, which don't contain such information are ignored). Two possible ways to encode this information are acceptable:
1. Encode it in read names (as a result of dropTag phase).
2. Use .bam tags (i.e. output of 10x Cell Ranger). Tag names can be configured in [*config.xml* file](##additional-notes).

Count matrix estimation also requires information about the source gene for the reads. It can be provided in two ways:
1. Use gene annotation in either *.bad* or *.gtf* format. To provide such file, "*-g*" option should be used.
2. Use .bam tags. Tag name can be configured in [*config.xml* file](##additional-notes).

Another crucial moment in estimation of count matrix is correction of cell barcode errors. Most protocols provide the list of real barcodes, which simplifies the task of correction. If such file is available, path to the file **should be specified in the *config.xml* file** (*Estimation/Merge/barcodes_file*). This can dramatically increase quality of the result. Lists for inDrop protocols can be found at *dropEst/data/barcodes/*. Two algorithms of barcode correction are available:
1. Simple, "*-m*" option. This algorithm is recommended in the case, where barcodes list is supplied.
2. Precise, "*-M*" option. Doesn't requires list of real barcodes to obtain high-quality results, however has significantly lower performance.

Example command:
```bash
dropest [options] [-m] -g ./hg38/genes.gtf -c ./config.xml ./alignment.*/accepted_hits.bam
```

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
*  -R, --reads-output: print count matrix for reads and don't use UMI statistics  
*  -q, --quiet : disable logs  
*  -w, --write-mtx : write out matrix in MatrixMarket format  

### Output
<!-- TODO: add that output has data of two types: all cells and filtered cells -->
Result of this phase is .rds file with the next fields:
* **cm** (sparse matrix): count matrix in sparse format
* **reads_per_chr_per_cell** (list of data.frame): number of reads per cell (row) for each chromosome (column):
  * **Exon** (data.frame): exonic reads
  * **Intron** (data.frame): intronic reads
  * **Intergenic** (data.frame): intergenic reads
* **mean_reads_per_umi** (vector): mean number of reads per UMI for each cell
* some additional info <!-- TODO: describe it -->

To additionaly print the file with count matrix in MatrixMarket format use "*-w*" option.

## dropReport
To run the report you have to install [dropestr](#dropestr-package) R package.

Required R packages for the report script:

```R
install.packages(c("rmarkdown","preseqR"))
```
You need pandoc for the creation of the report.html.

The Report can be called with

```bash
Rscript dropReport.Rsc cell.counts.rds
```

## dropEstR package

This package implements UMI errors corrections and low-quality cells filtration, which are not performed during dropEst phase.

To install the package, use:

```r
devtools::install_github('hms-dbmi/dropEst/dropestr')
```

Package content:
  * filtration of low-quality cells (see vignette "*low-quality-cells*")
  * correction of UMI errors (see vignette "*umi-correction*")
  * quality control (see *dropReport.Rsc*)  <!-- TODO: create vignette -->

## Additional notes
Description of the fields for the config file is provided in *dropEst/configs/config_desc.xml*.
