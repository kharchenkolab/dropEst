# Changelog

## [Upcoming]
### Changed
* Now, information about reads is kept in separate file instead (*.reads.gz), which should be passed to dropEst.
* New format of *reads_per_umi_per_cell* in cell.counts.rds.

### Added
* Filtration of reads by barcode quality ("*TagsSearch/Processing/min_barcode_quality*" and 
"*Estimation/Other/min_barcode_quality*" fields in the config)
* UMI correction algorithm now uses UMI quality (only Illumina 1.8 or later formats are supported)

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
