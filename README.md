# dropEst - Pipeline
Pipeline for estimating molecular count matrices for droplet-based single-cell RNA-seq measurements

## General processing steps
1. **dropTag** - Tagging and merging of the reads
2. Alignment to reference genome
3. **dropEst** - Counting of Reads, Merging of cells by barcode and export to .rds count matrix.
4. **dropReport** - Generating report on library quality.

## Setup & Installation
### System requirements
* Boost >= 1.54
* BamTools
* Zlib
* R >= 3.2.2 with pacakges:
  * Rcpp
  * RcppAramdillo
  * RInside
  * Matrix
  * fitdistrplus *(optional, used only with "-M" option)*
* Compiler with c++11 support *(was tested with gcc >= 4.8.5 and CLang 3.9.1)*

### Installation
Install R packages:

```R
install.packages(c("Rcpp","RcppArmadillo", "RInside", "Matrix", "fitdistrplus"))
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

## dropTag
    droptag -- generate tagged fastq files for alignment

Currently supported Indrop Library versions for Tag stage:  
* Indrop v1
  * Read structure:
    * Biological Read (gene)
    * Metadata Read (barcode)
* Indrop v2
  * Read structure:
    * Biological Read (gene)
    * Metadata Read (barcode)
* Indrop v3
  * Read structure (You need 4 different read-files from your sequencing facility)
    * Biological Read (~ 60 bp)
    * Barcode + UMI (14 bp)
    * Cell Barcode (8 bp)
    * Library Barcode (8 bp)

### inDrop v1 & v2
* Reads to provide: (can be fastq or fastq.gz)
  * Read_1: Index reads (length?)
  * Read_2: Biological read (length longer?)
* Provide the config.xml file *(Explanation of the Paramters in "config-desc.xml"
(write more about important parameters here?)*
* Call `droptag`:
```bash
./droptag -c ./config.xml barcode_reads gene_reads
```

### inDrop v3
* Reads to provide: (can be fastq or fastq.gz)
  * Read 1: Cell barcode part 1 (length 8)
  * Read 2: Cell barcode part 2 + UMI (length 14)
  * Read 3: Gene read (length 61)
  * Read 4: Library barcode (length 8) *(optional)*
* Provide the config.xml file *(Explanation of the Paramters in "config-desc.xml" (write more about important parameters here?))*
* Call `droptag`:
```bash
./droptag -c ./config.xml [-S] [-t library_tag] barcode1_reads barcode2_umi_reads gene_reads [library_tags]
```

#### further options / flags:
*  -c, --config filename: xml file with droptag parameters  
*  -l, --log-prefix prefix: logs prefix  
*  -n, --name name: alternative output base name  
*  -S, --save-stats : save stats to rds file  
*  -t, --lib-tag library tag : (for IndropV3 with library tag only)  
*  -q, --quiet : disable logs  

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

**You need to provide the right barcode file in the config.xml file!** (How to find the right file?)

After running tophat you can run the dropest command on all ./align-output/Sample_1.\*/accepted_hits.bam at once.

You have to provide a genome annotation (`genomes.gtf`) file with the genome information to run dropest.

For indrop-v3 you should use the option -m which fixes barcode errors and improves estimation per cell.

```bash
dropest [options] -g ./hg38/genes.gtf -c ./config.xml ./align-output/Sample_1.\*/accepted_hits.bam
```

### Options for dropest  
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

## dropReport
To run the report you have to install [dropestr](#dropest-r-package) R package.



Required R packages for the report script:


```R
install.packages(c("rmarkdown","preseqR"))
```
You need pandoc for the creation of the report.html.

The Report can be called with

```bash
Rscript dropReport.Rsc cell.counts.rds
````

# dropEst R package
I haven't used this yet! I'd like to try it and help you with documentation.

To install the package, use

```R
devtools::install_github('hms-dbmi/dropEst/dropestr').
```
