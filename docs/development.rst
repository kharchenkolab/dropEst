For developers
--------------

Adding new protocols to dropTag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each new protocol you need to write new c++ class, which inherits TagsFinderBase
(*“TagsSearch/TagsFinderBase.cpp”*). There are several examples for
existed protocols (see all classes with suffix “TagsFinder”). All you
really need there is to define function ``parse_fastq_records``. After
that you need to add new protocol type to function ``get_tags_finder``
in *“droptag.cpp”*. It takes only several lines, e.g.:

.. code:: cpp

   if (protocol_type == "ddSEQ")
       return ...

Having this your are able to set “TagsSearch/protocol” to “ddSEQ” in the
config.xml and work with ddSEQ as with any other supported protocol in
very efficient and parallel manner.

General workflow of dropTag under the hood
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Entry point is *"droptag.cpp"* file, ``main()`` function. Most of it is just parsing
CLI parameters and logging. So, there are only two important places:
(1) ``get_tags_finder``: factory-like function, which creates TagsFinder for a
specific protocol based on configs, and (2) ``finder->run``: TagsFinder's method,
which actually does the whole work.

Workflow of ``finder->run`` is the following:

#. Read fastq records from the files, which contain gene and barcode reads.
   Reading is performed synchronously over all files, as records in different files correspond to each other.
#. This set of fastq records is parsed to a single gene record with its parameters
   (i.e. cell barcode, UMI and its quality). There are two options on how to store the
   parameters: it can either be done in the read name of the gene record or in a
   separate data structure. In the later case, this information is stored as a gzipped
   tsv file (see ``-s`` option).
#. Parsed record is converted to a string and gzipped.
#. Gzipped info is written to the output file.

Though code is more complicated then that, because it's implemented in a
multi-threaded way, and there it can't be done by simple "parallel map" style.
So, parallelism style looks more like MPI (though of course it's implemented in C++11
threads), and works as the follows (see ``TagsFinderBase::run_thread`` function):

- All data is stored in concurrent queues, which have limited maximal size.
- Each workers independently iterates over all 4 tasks and check if it can do some
  work on it. If task is single-threaded and another worker is already doing it, or if
  corresponding queue is already full, the worker goes to the next task.

Such scheme allows to achieve ~10 times higher performance.