dropEstR package
----------------

This package implements UMI errors corrections and low-quality cells
filtration, which are not performed during dropEst phase.

To install the package, use:

.. code:: r

    devtools::install_github('hms-dbmi/dropEst/dropestr', dependencies = T)

Package content: \* filtration of low-quality cells (see vignette
"*low-quality-cells*") \* correction of UMI errors (see vignette
"*umi-correction*") \* quality control (see *dropReport.Rsc*)

