# dropEst - Pipeline
Pipeline for initial analysis of droplet-based single-cell RNA-seq data

General outline of the pipeline.

## Prerequisites: 
### Sequencing output of an Indrop Experiment.

Currently supported Indrop Library versions:  
- Indrop v1
    - Read structure: 
        - Biological Read (gene)
        - Metadata Read (barcode)
- Indrop v2
    - Read structure: 
        - Biological Read (gene)
        - Metadata Read (barcode)
- Indrop v3
    - Read structure (You need 4 different read-files from your sequencing facility)
        - Biological Read (~ 60 bp)
        - Barcode + UMI (14 bp)
        - Cell Barcode (8 bp)
        - Library Barcode (8 bp)


# Pipeline Overview
## General processing steps:
1. **dropTag** - Tagging and merging of the reads
2. Alignment to reference genome (**tophat**)
3. **dropEst** - Counting of Reads, Merging of cells by barcode and export to .rds count matrix. 


## Setup & Installation
>> It is quite hard to install
>> I tried it quickly on my macbook (yeah, unrealistic, but still) and I couldn't even finish running cmake, as eg R is quite different if installed "precompiled"

>> Installing it on orchestra is straightforward, but can we assume others have the same setups? 

You need:
- gcc (which versions work? Only 4.x?) 
- boost (probably all versions work?)
- R (all versions?) with Rcpp and RcppArmadillo installed
to install them: 
> R
> install.packages(c("Rcpp","RcppArmadillo"))

Download this repository:
git clone git@github.com/hms-dbmi/dropEst.git

You need to change these variables in the CMakeLists.txt file:
You have to locate your bamtools installation, your R and the R-package-library:

line 22
set(BAMTOOLS_ROOT "`$path to your bamtools` /opt/bamtools-2.2.3/")
line 25
set(R_ROOT "`path to your R: ` /opt/R-3.3.1/lib64/R/")
line 26
set(R_PACKAGES "`path to your packages` ~/R/x86_64-pc-linux-gnu-library/3.3/")


Then run in the directory:
> cmake . && make 


# **dropTag** 

	droptag -- generate tagged fastq files for alignment

Understanding your read files to submit them: Extracting the actual reads and count the occuring barcodes:

zcat `indropread_Readnumber1.fastq.gz` | head -n 20000 | awk ‘(NR+2)%4==0{print}’ | sort -S 4G | uniq -c | sort -k1,1n

Two reads should contain a clear barcode of length 8. One of these contains a smaller number of different barcodes, this is / (might be) the library id. 

## inDrop v3
- Reads to provide: (can be fastq or fastq.gz)
    - Read_1: Cell barcode part 1 (length 8)
    - Read_2: Cell barcode part 2 + UMI (length 14)
    - Read_3: Gene read (length 61)
    - Read_4: Library barcode (length 8)
 (optional)

- Provide the config.xml file  
Explanation of the Paramters in "config-desc.xml"
(write more about important parameters here?)

Call the `droptag` function:

> ./droptag -c ./config.xml Read_1 (Barcode) Read_2 (barcode + UMI) Read_3 (gene read, long) -t Read_4 (library barcode)

## inDrop v1 & v2 
- Reads to provide: (can be fastq or fastq.gz)
    - Read_1: Index reads (length?)
    - Read_2: Biological read (length longer?)

- Provide the config.xml file  
Explanation of the Paramters in "config-desc.xml"
(write more about important parameters here?)

Call the `droptag` function:

> ./droptag -c ./config.xml Read_1 (Barcodes) Read_2 (gene Read)

### further options / flags:

OPTIONS:  
-	-c, --config filename: xml file with droptag parameters  
-	-l, --log-prefix : logs prefix  
-	-n, --name BASE_NAME : specify alternative output base name  
-	-s, --save-reads-names : serialize reads parameters to save names  
-	-S, --save-stats filename : save stats to rds file  
-	-t, --lib-tag library tag : (for IndropV3 with library tag only)  
-	-q, --quiet : disable logs  

  
Running droptag results in


Running droptag:
### Example Reads of a indropV3:

#### Barcode Read (Read_1)
@NS500144:857:HWK3FBGX2:1:11101:23748:1039 2:N:0:1  
`AAAGAGG`  
+  
AAAAAEE/  

#### Barcode Read + UMI (Read_2)
@NS500144:857:HWK3FBGX2:1:11101:23748:1039 4:N:0:1  
GTCATATT`ACGTTC`  
+  
AAAAAEEEEEAEAA  

#### Biological Read (Read_3)
@NS500144:857:HWK3FBGX2:1:11101:23748:1039 1:N:0:1  
*CCTACNACGACAACTAAAATTTCACTNCACATNAAAACATCACTTCGGATTGCAAGCCGCA*  
+  
AA/AA#E<EAEEEEEEEEEEEEEEEE#AEEEE#EEEEAAEEAEE/<EEE///////E////  

#### Optional Library Tag (Read_4)
@NS500144:857:HWK3FBGX2:1:11101:23748:1039 3:N:0:1  
TACTCCTT  
+  
AA//AEEE  

### Example: Resulting fastq file
@1150996099X2!`AAAGAGGT`GTCATATT#`ACGTTC` 
*CCTACNACGACAACTAAAATTTCACTNCACATNAAAACATCACTTCGGATTGCAAGCCGC*  
+  
AA/AA#E<EAEEEEEEEEEEEEEEEE#AEEEE#EEEEAAEEAEE/<EEE///////E///  

# **Alignment**
    Here you are required to have experience with alignment of RNAseq data. We can't provide a full tutorial.
    This tasks takes a while for large sequencing libraries, you will need a high performance cluster


[Link to Tophat `Getting Started`](https://ccb.jhu.edu/software/tophat/tutorial.shtml)

To align the resulting fastq files to a reference genome call tophat on all files. droptag exports the tagged reads into multiple files. 

Run tophat on each resulting fastq file
The files are written by droptag in the scheme: "Samplename".n.fastq.gz

You need a Bowtie2 - indexed genome file. Please refer to the tophat manual for this.  

> tophat -o ./align-output/Sample_1.1 ./hg19/Bowtie2Index/genome ./tag-output/Sample_1.1.fastq.gz

The result needed for the count estimation is the ./align-output/Sample_1.*/accepted_hits.bam


## dropest
    dropest: estimate molecular counts per cell

**You need to provide the right barcode file in the config.xml file!** (How to find the right file?)

After running tophat you can run the dropest command on all ./align-output/Sample_1.*/accepted_hits.bam at once. 

You have to provide a `genomes.gtf` file with the genome information to run dropest.

For indrop-v3 you should use the option -m which fixes barcode errors and improves estimation per cell. 


> dropest [options] -g ./hg19/genes.gtf -c ./config.xml ./align-output/Sample_1.*/accepted_hits.bam 

  
### Options for dropest  
-	-b, --bam-output: print tagged bam files  
-	-c, --config filename: xml file with estimation parameters  
-	-C, --cells num: maximal number of output cells  
-	-f, --filled-bam: bam file already contains genes/barcodes tags  
-	-F, --filtered-bam: print tagged bam file after the merge and filtration  
-	-g, --genes filename: file with genes annotations (.bed or .gtf)  
-	-G, --genes-min num: minimal number of genes in output cells  
-	-l, --log-prefix : logs prefix  
-	-L, --gene-match-level :  
	-	e: count UMIs with exonic reads only;  
	-	i: count UMIs with intronic reads only;  
	-	E: count UMIs, which have both exonic and not annotated reads;  
	-	I: count UMIs, which have both intronic and not annotated reads;  
	-	B: count UMIs, which have both exonic and intronic reads;  
	-	A: count UMIs, which have exonic, intronic and not annotated reads.  
	-	Default: -L eEBA.  
-	-m, --merge-barcodes : merge linked cell tags  
-	-M, --merge-barcodes-precise : use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available  
-	-o, --output-file filename : output file name  
- **(does not work!)**	-p, --parallel number_of_threads : number of threads  
-	-R, --reads-output: print count matrix for reads and don't use UMI statistics  
-	-q, --quiet : disable logs  
-	-w, --write-mtx : write out matrix in MatrixMarket format  
          
# dropEst R package
I haven't used this yet! I'd like to try it and help you with documentation.

To install the package, use 
> devtools::install_github('hms-dbmi/dropEst/dropestr').
