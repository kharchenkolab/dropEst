dropEstR package
----------------

This package implements UMI errors corrections and low-quality cells
filtration, which are not performed during dropEst phase.

To install the package, use:

.. code:: r

    devtools::install_github('hms-dbmi/dropEst/dropestr', dependencies = T)

Package content:

- filtration of low-quality cells (see vignette "`low-quality-cells.Rmd <https://github.com/hms-dbmi/dropEst/blob/master/dropestr/vignettes/low-quality-cells.Rmd>`__")
- correction of UMI errors (see vignette "`umi-correction.Rmd <https://github.com/hms-dbmi/dropEst/blob/master/dropestr/vignettes/umi-correction.Rmd>`__"). *This step isn't needed for most of modern protocols.*
- quality control (see "`report.Rmd <https://github.com/hms-dbmi/dropEst/blob/master/scripts/report.Rmd>`__")

