Alignment
---------

dropTag writes the tagged reads into multiple files. All these files
must be aligned to reference, and all bam files with the alignments must
be provided as input for the dropEst stage. In the paper we used
`TopHat2 <https://ccb.jhu.edu/software/tophat/tutorial.shtml>`__
aligner, however any RNA-seq aligners (i.e.
`Kallisto <https://pachterlab.github.io/kallisto/>`__ or
`STAR <https://github.com/alexdobin/STAR>`__) can be used.

Alignment with TopHat
~~~~~~~~~~~~~~~~~~~~~

1. Install `bowtie
   index <http://bowtie-bio.sourceforge.net/tutorial.shtml#preb>`__ for
   the sequenced organism.
2. Download corresponding `gene
   annotation <http://genome.ucsc.edu/cgi-bin/hgTables>`__ in the gtf
   format.
3. Run TopHat2 for each file:

   .. code:: bash

       tophat2 -p number_of_threads --no-coverage-search -g 1 -G genes.gtf -o output_dir Bowtie2Index/genome reads.fastq.gz

4. The result needed for the count estimation is
   *./output\_dir/accepted\_hits.bam*.

Alignment with Kallisto
~~~~~~~~~~~~~~~~~~~~~~~

**Prior to v0.42.4, Kallisto didn't use BAM flags to mark
primary/secondary alignments. Only versions :math:`\geq` 0.42.4 are
supported.**

1. Download transcript sequences in .fasta format (i.e.
   `Ensembl <https://www.ensembl.org/info/data/ftp/index.html>`__
   genes).
2. Build Kallisto index: ``kallisto index -i genes.fa.gz``.
3. Run

   .. code:: bash

       kallisto quant --pseudobam --single -i genes.index -o out -l mean_length -s std_length reads.fastq.gz

   Here, *mean\_length* is mean length of RNA fragment (not read length)
   and *std\_length* is standard deviataion of RNA fragment length. You
   should specify values, according to the experiment design.
