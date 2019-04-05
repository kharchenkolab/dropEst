dropReport
----------

To run the report you have to install: \*
`dropestr <#dropestr-package>`__ R package with all dependencies
(``dependencies = T``). \*
`pandoc <https://pandoc.org/installing.html>`__ \* R packages:
``R     install.packages(c("optparse","rmarkdown"))``

To run the report use:

.. code:: bash

    ./dropReport.Rsc cell.counts.rds

To see full list of options use:

.. code:: bash

    ./dropReport.Rsc -h

Troubleshooting
~~~~~~~~~~~~~~~

If you get the error *"pandoc version 1.12.3 or higher is required and
was not found"*, try to set path in the corresponding environment
variable: ``export RSTUDIO_PANDOC=/path/to/pandoc``.
