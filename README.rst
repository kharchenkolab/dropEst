dropEst - Pipeline
==================

Pipeline for estimating molecular count matrices for droplet-based
single-cell RNA-seq measurements. If you used the pipeline in your
research, please `cite <#citation>`__ the corresponding
`paper <https://doi.org/10.1186/s13059-018-1449-6>`__. To reproduce
results from the paper see `this
repository <https://github.com/VPetukhov/dropEstAnalysis>`__.

Documentation
-------------

For detailed explanations, please see the `documentation <https://dropest.readthedocs.io/en/latest/>`__

Particularly:

- `Installation <https://dropest.readthedocs.io/en/latest/setup.html#installation>`__
- `Integration with Velocyto <https://dropest.readthedocs.io/en/latest/dropest.html#velocyto-integration>`__

If you have problems with installation, please look at the `Troubleshooting <https://dropest.readthedocs.io/en/latest/setup.html#troubleshooting>`__ page and open an `issue <https://github.com/hms-dbmi/dropEst/issues>`__ if there is nothing.

News
----

[0.8.6] - 2019-08-01
~~~~~~~~~~~~~~~~~~~~

-  Added support for Drop-seq and CEL-Seq2

See `Changelog <https://github.com/hms-dbmi/dropEst/blob/develop/CHANGELOG.rst>`__ for the full list.

General processing steps
------------------------

1. **dropTag**: extraction of cell barcodes and UMIs from the library.
   Result: demultiplexed .fastq.gz files, which should be aligned to the
   reference.
2. **Alignment** of the demultiplexed files to reference genome. Result:
   .bam files with the alignment.
3. **dropEst**: building count matrix and estimation of some statistics,
   necessary for quality control. Result: .rds file with the count
   matrix and statistics. *Optionally: count matrix in MatrixMarket
   format.*
4. **dropReport** - Generating report on library quality.
5. **dropEstR** - UMI count corrections, cell quality classification

Examples
--------

Complete examples of the pipeline can be found at
`EXAMPLES.md <examples/EXAMPLES.md>`__.

`Here <http://pklab.med.harvard.edu/viktor/dropest_paper/dropest_0.8.5.zip>`__
are results of processing of
`neurons\_900 <https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neurons_900>`__
10x dataset.

Supported protocols
-------------------

- 10x
- CEL-Seq2
- Drop-seq
- iCLIP
- inDrop (v1-3)
- Seq-Well
- SPLiT-seq

Citation
--------

If you find this pipeline useful for your research, please consider citing the paper:

Petukhov, V., Guo, J., Baryawno, N., Severe, N., Scadden, D. T.,
Samsonova, M. G., & Kharchenko, P. V. (2018). dropEst: pipeline for
accurate estimation of molecular counts in droplet-based single-cell
RNA-seq experiments. Genome biology, 19(1), 78.
doi:10.1186/s13059-018-1449-6
