dropEst
-------

::

    dropest: estimate molecular counts per cell

This phase requires aligned .bam files as input and uses them to
estimate count matrix. These files must contain information about cell
barcode and UMI for each read (reads, which don't contain such
information are ignored). Three possible ways to encode this information
are acceptable: 1. Encode it in read names (as a result of dropTag
phase). 2. Use .bam tags (i.e. output of 10x Cell Ranger). Tag names can
be configured in `*config.xml* file <##additional-notes>`__. 3. Save
this information in a separate file using ``-s`` option on dropTag phase
and pass it to dropEst with ``-r`` option.

**Note.** To run Bayesian algorithm of UMI error correction with
dropestr, you need to pass information about UMI base call quality to
dropEst. You can do it either using .bam tags or ``-s`` dropTag option
(i.e. option 1 isn't available in this case). You still can run simpler
algorithms (i.e. *cluster* or *directional*, see paper) without this
information.

Count matrix estimation also requires information about the source gene
for the reads. It can be provided in three ways: 1. Use gene annotation
in either *.bad* or *.gtf* format. To provide such file, "*-g*" option
should be used. 2. Use .bam tags. Tag name can be configured in
`*config.xml* file <##additional-notes>`__. 3. Pseudoaligners encode
gene names as chromosome names

Another crucial moment in estimation of count matrix is correction of
cell barcode errors. Most protocols provide the list of real barcodes,
which simplifies the task of correction. If such file is available, path
to the file **should be specified in the *config.xml* file**
(*Estimation/Merge/barcodes\_file*). This can dramatically increase
quality of the result. Lists for inDrop protocols can be found at
*dropEst/data/barcodes/*. Two algorithms of barcode correction are
available: 1. Simple, "*-m*" option. This algorithm is recommended in
the case, where barcodes list is supplied. 2. Precise, "*-M*" option.
Doesn't requires list of real barcodes to obtain high-quality results,
however has lower performance for large datasets.

Example command:

.. code:: bash

    dropest [options] [-m] [-r pipeline_res.params.gz] -g ./hg38/genes.gtf -c ./config.xml ./alignment.*/accepted_hits.bam

Usage of tagged bam files (e.g. 10x, Drop-seq) as input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some protocols provide pipelines, which create .bam files with
information about CB, UMI and gene. To use these files as input, specify
"*-f*" option. Example:

.. code:: bash

    dropest [options] -f [-g ./genes.gtf] -c ./config.xml ./pipeline_res_*.bam

If "*-g*" option is provided, genes are parsed from the gtf, and
information about genes from the bam file is ignored.

To specify corresponding .bam tag names, use "*Estimation/BamTags*"
section in the config (see *configs/config\_desc.xml*).

Usage of pseudoaligners
~~~~~~~~~~~~~~~~~~~~~~~

Pseudoaligners, such as Kallisto store gene / transcript names in the
field of chromosome name. To parse such files, use "*-P*" option.
Example:
``bash  dropest [options] -P [-r pipeline_res.params.gz] -c ./config.xml ./kallisto_res_*.bam``

Count intronic / exonic reads only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One feature of the pipeline is the ability to count only UMIs, which
reads touch only specific parts of a genome. Option *"-L"* allows to
specify all acceptable types of regions: \* e: UMIs with exonic reads
only \* i: UMIs with intronic reads only \* E: UMIs, which have both
exonic and not annotated reads \* I: UMIs, which have both intronic and
not annotated reads \* B: UMIs, which have both exonic and intronic
reads \* A: UMIs, which have exonic, intronic and not annotated reads

Thus, to count all UMIs with exonic **or** not annotated reads, use *"-L
eE"*. Default value: *"-L eEBA"*, i.e. to count all UMIs, which have at
least one reads, touching exon.

Example commands: \* Intronic reads only:
``bash     dropest [-f] [-g ./genes.gtf] [-r pipeline_res.params.gz] -L i -c ./config.xml ./alignment_*.bam``
\* Exonic reads only:
``bash     dropest [-f] [-g ./genes.gtf] [-r pipeline_res.params.gz] -L e -c ./config.xml ./alignment_*.bam``
\* Exon/intron spanning reads:
``bash     dropest [-f] [-g ./genes.gtf] [-r pipeline_res.params.gz] -L BA -c ./config.xml ./alignment_*.bam``

The pipeline can determine genome regions either using .gtf annotation
file or using .bam tags, i.e. for CellRanger output (see
*Estimation/BamTags/Type* in *configs/config\_desc.xml*). If .gtf file
isn't provided and .bam file doesn't containt annotation tags, all reads
with not empty gene tag are considered as exonic.

**Note.** There is no way to extract information about intronic reads
from pseudoalignment, so you can't use pseudoaligners at this stage.

Velocyto integration
^^^^^^^^^^^^^^^^^^^^

For some purposes (i.e. `velocyto <http://velocyto.org/>`__) it can be
useful to look separately at the fraction of intronic and exonic UMIs.
Option *"-V"* allows to output three separate count matrices, each of
which contains only UMIs of a specific type: intronic, exonic or
exon/intron spanning. These matrices are stored in the separate file
*"cell.counts.matrices.rds"*.

**Note.** Please **ensure that you provided gtf file with genes** with
``-g`` option.

Command line arguments for dropEst
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  -b, --bam-output: print tagged bam files
-  -c, --config filename: xml file with estimation parameters
-  -C, --cells num: maximal number of output cells
-  -f, --filled-bam: bam file already contains genes/barcodes tags
-  -F, --filtered-bam: print tagged bam file after the merge and
   filtration
-  -g, --genes filename: file with genes annotations (.bed or .gtf)
-  -G, --genes-min num: minimal number of genes in output cells
-  -l, --log-prefix : logs prefix
-  -m, --merge-barcodes : merge linked cell tags
-  -M, --merge-barcodes-precise : use precise merge strategy (can be
   slow), recommended to use when the list of real barcodes is not
   available
-  -o, --output-file filename : output file name
-  -P, --pseudoaligner: use chromosome name as a source of gene id
-  -q, --quiet : disable logs
-  -r, --read-params filenames: file or files with serialized params
   from tags search step. If there are several files, they should be
   provided in quotes, separated by space: "file1.params.gz
   file2.params.gz file3.params.gz"
-  -R, --reads-output: print count matrix for reads and don't use UMI
   statistics
-  -u, --merge-umi: apply 'directional' correction of UMI errors. This
   option prevents output of ``reads_per_umi_per_cell``. If you want to
   apply more advanced UMI correction, don’t use ‘-u’, but use follow up
   R analysis.
-  -V, --velocyto : save separate count matrices for exons, introns and
   exon/intron spanning reads
-  -w, --write-mtx : write out matrix in MatrixMarket format

Output
~~~~~~

.. raw:: html

   <!-- TODO: add that output has data of two types: all cells and filtered cells -->

Result of this phase is cell.counts.rds file with the next fields: \*
**cm** (sparse matrix): count matrix in sparse format \*
**reads\_per\_chr\_per\_cell** (list of data.frame): number of reads per
cell (row) for each chromosome (column): \* **Exon** (data.frame):
exonic reads \* **Intron** (data.frame): intronic reads \*
**Intergenic** (data.frame): intergenic reads \*
**mean\_reads\_per\_umi** (vector): mean number of reads per UMI for
each cell \* some additional info

To additionaly print the file with count matrix in MatrixMarket format
use "*-w*" option.
