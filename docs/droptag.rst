dropTag
-------

::

    droptag -- generate tagged fastq files for alignment

Example command:

.. code:: bash

    ./droptag -c config.xml [-S] [-s] reads1.fastq reads2.fastq [...]

Positional arguments of the dropTag phase contain paths to the read
files, obtained with a sequencer. These files can be either in *.fastq*
of *.fastq.gz* format. The file order depends on the type of used
protocol. Below is the instruction on the usage for currently supported
protocols.

**Note.** Algorithm of UMI correction requires information about
base-call quality of UMIs. To save this information, use ``-s`` option.
In this case, in addition to fastq files with reads, the pipeline saves
separate gzipped file with read parameters. These files must be passed
to dropEst with ``-r`` option.

Protocols
~~~~~~~~~

inDrop v1 & v2
^^^^^^^^^^^^^^

-  File 1: barcode reads. Structure:
-  Cell barcode, part 1
-  Spacer
-  Cell barcode, part 2
-  UMI
-  File 2: gene reads

| Example config file is located at
  "*dropEst/configs/indrop\_v1\_2.xml*".
| Example command:

.. code:: bash

    ./droptag -c [-S] [-s] ~/dropEst/configs/indrop_v1_2.xml barcode_reads.fastq gene_reads.fastq

inDrop v3
^^^^^^^^^

-  File 1: cell barcode, part 1 *(default length: 8bp)*
-  File 2: cell barcode + UMI, part 1 *(default length: >= 14bp)*
-  File 3: gene read
-  *File 4 (optional): library tag*

If a file with library tags provided, option "-t" is required. ***This
option wasn't tested properly, so it's better to avoid using it.***

| Example config file is located at "*dropEst/configs/indrop\_v3.xml*".
| Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/indrop_v3.xml [-S] [-s] [-t library_tag] barcode1_reads.fastq barcode2_reads.fastq gene_reads.fastq [library_tags.fastq]

10x
^^^

-  File 1: library tag *(default length: 8bp)*
-  File 2: cell barcode + UMI, part 1 *(default length: 16+10=26bp)*
-  File 3: gene read

| Example config file is located at "*dropEst/configs/10x.xml*".
| Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/10x.xml [-S] [-s] lib_tag_reads.fastq barcode_reads.fastq gene_reads.fastq

While dropTag provides way to demultiplex 10x data, `Cell
Ranger <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger>`__
is still highly recommended tool for this. `dropEst <##dropEst>`__ phase can be
ran on the Cell Ranger demultiplexed .bam file to obtain data in the
format, optimized for the subsequent analysis.

**Note.** Sometimes 10x CellRanger isn't able to determine gene, from
which a read originated. In this cases it fills gene info with a list of
possible genes, separated by semicolon (e.g.
'ENSG00000255508;ENSG00000254772'). These genes **must be filtered out**
prior to further analysis:

.. code:: r

    holder <- readRDS('./cell.counts.rds')
    cm <- holder$cm[grep("^[^;]+$", rownames(holder$cm)),]

iCLIP
^^^^^

-  File 1: Gene reads with barcodes at the beginning of the sequence

| Example config file is located at "*dropEst/configs/iclip.xml*".
| Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/iclip.xml [-S] [-s] data.fastq

**Note.** Implementation of iCLIP wasn't tested properly. Please, be
careful using it. Anyone who used it is very welcome to comment it
either in Issues or by e-mail.

SPLiT-seq
^^^^^^^^^

-  File 1: UMI + cell barcode (3 parts)
-  File 2: Gene reads

Example config file is located at "*dropEst/configs/split\_seq.xml*".
Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/split_seq.xml [-S] [-s] barcode_reads.fastq gene_reads.fastq

Seq-Well
^^^^^^^^

-  File 1: Cell barcode + UMI
-  File 2: Gene reads

Example config file is located at "*dropEst/configs/seq\_well.xml*".
Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/seq_well.xml [-S] [-s] barcode_reads.fastq gene_reads.fastq

Drop-seq
^^^^^^^^

-  File 1: Cell barcode + UMI
-  File 2: Gene reads

Example config file is located at "*dropEst/configs/drop\_seq.xml*".
Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/drop_seq.xml [-S] [-s] barcode_reads.fastq gene_reads.fastq

**Note.** Implementation of Drop-seq wasn't tested properly. Please, be
careful using it. Anyone who used it is very welcome to comment it
either in Issues or by e-mail.

Alternatively, to run processing on the bam files, obtained with
`Drop-seq tools <https://github.com/broadinstitute/Drop-seq>`__,
see `Usage of tagged bam files <dropest.html#usage-of-tagged-bam-files-e-g-10x-drop-seq-as-input>`__.

CEL-Seq2
^^^^^^^^

-  File 1: Cell barcode + UMI
-  File 2: Gene reads

Example config file is located at "*dropEst/configs/cel\_seq2.xml*".
Example command:

.. code:: bash

    ./droptag -c ~/dropEst/configs/cel_seq2.xml [-S] [-s] barcode_reads.fastq gene_reads.fastq

**Note.** Implementation of CEL-Seq2 wasn't tested properly. Please, be
careful using it. Anyone who used it is very welcome to comment it
either in Issues or by e-mail.

Command line arguments for dropTag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  -c, --config filename: xml file with droptag parameters
-  -l, --log-prefix prefix: logs prefix
-  -n, --name name: alternative output base name
-  -p, --parallel number: number of threads (usage of more than 6
   threads should lead to significant speed up)
-  -r, --reads-per-out-file : maximum number of reads per output file;
   (0 - unlimited). Overrides corresponding xml parameter.
-  -S, --save-stats : save stats to rds file. This data is used on the
   dropReport phase.
-  -t, --lib-tag library tag : (for IndropV3 with library tag only)
-  -q, --quiet : disable logs

Please, use ``./droptag -h`` for additional help.
