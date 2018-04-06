# Changelog

## [0.8.1] - 2018-04-06
### Changed
* Fixed some bugs in dropestr
### Added
* Added protocol type 10x (which is alias for indrop3) to dropTag

## [0.8.0] - 2018-03-20
### Changed
* Now, information about reads is kept in separate file instead (*.reads.gz), which should be passed to dropEst
* New format of *reads_per_umi_per_cell* in cell.counts.rds

### Added
* UMI correction algorithm now uses UMI quality (only Illumina 1.8 or later formats are supported)

## [0.7.6] - 2018-02-22
### Added
* Version number output in dropEst
### Changed
* Fixed bug with report generation
* Fixed bugs with estimation of low-quality cells

## [0.7.5] - 2018-01-28
### Added
* Dockers for Centos6, Centos7 and Debian9
* iclip protocol support
* 10x barcodes
* Published code for filtration of multialigned reads from bam files with cell mixture.
### Changed
* Optimized precise merge performance
* Fixed bug with fastq split during dropTag
* **New format of barcode files**
* Optimized memory usage in parsing read params from file
* Algorithm of filtration of low-quality cells was significantly improved

## [0.7.1] - 2017-10-18
### Added
* Integration with velocyto

## [0.7.0] - 2017-10-17
### Changed
* Optimized cmake
* **Secondary alignments are filtered now**

### Added
* Output UMIs with only exonic or only intronic reads

## [0.6.8] - 2017-09-22
### Added
* Filtration of reads by barcode quality ("*TagsSearch/Processing/min_barcode_quality*" and
"*Estimation/Other/min_barcode_quality*" fields in the config)
* dropEst is now able to parse read type (e.g. exonic/intronic) from .bam file (see *config_desc.xml*)

## [0.6.7] - 2017-09-13
### Changed
* Fixed bug, which led to erroneous parsing of incorrect read (e.g. reads without spacer for Indrop V1)

### Added
* Parallelized dropTag ("*-p*" option)

## [0.6.5] - 2017-09-07
### Changed
* Optimized memory usage and performance of dropEst
* Sorting for cells selection (by number of genes) is stable now
* Fixed bug with merge_targets in low-quality cells estimation
* Fixed bug with N's in UMIs after the merge

## [0.6.1] - 2017-09-05
### Added
* Support for pseudoaligners .bam format (usage of chromosome name as a source of gene name)
* Changelog

### Changed
* Check R libraries immediately after dropEst start

## [0.6.0] - 2017-09-04
### Added
* Versioning
