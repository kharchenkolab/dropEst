For developers
-----

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