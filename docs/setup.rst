Setup
-----

System requirements
~~~~~~~~~~~~~~~~~~~

-  Boost >= 1.54
-  BamTools library >= 2.5.0

   -  Note that some linux distributions have separate packages for the
      library and the executable, i.e. on Ubuntu you need
      ``libbamtools-dev``, but not ``bamtools``.
   -  or you can `build it
      locally <https://github.com/pezmaster31/bamtools/wiki/Building-and-installing>`__
      and then specify the location of the build when running cmake
      (e.g. ``cmake -D BAMTOOLS_ROOT=/home/username/bamtools .``)

-  Zlib *(was tested on 1.2.11 version)*
-  Bzip2 *(was tested on 1.0.5 version)*
-  R >= 3.2.2 with packages:
-  Rcpp
-  RcppEigen
-  RInside
-  Matrix
-  Compiler with c++11 support *(was tested with gcc >= 4.8.5 and CLang
   3.9.1)*
-  CMake >= 3.0

Installation
~~~~~~~~~~~~

Install R packages:

.. code:: r

    install.packages(c("Rcpp","RcppEigen", "RInside", "Matrix"))

Clone this repository:

.. code:: bash

    git clone https://github.com/hms-dbmi/dropEst.git

Build (**replace "/installation/path" with path to the folder, where you want to install the pipeline** or drop this option if you want to keep the build without installation):

.. code:: bash

    mkdir dropEst/build
    cd dropEst/build
    cmake -DCMAKE_INSTALL_PREFIX="/installation/path" .. && make

Install:

.. code:: bash

    make install

After this, ``droptag``, ``dropest`` and ``dropReport.Rsc`` binaries
must be available at ``installation/path/bin/`` folder.

Manual installation of dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the instruction on how to install specific versions of libraries
to the local folder (i.e. ``~/local/``). Let's store this directory in
``LOCAL_LIBS`` variable:

.. code:: bash

    export LOCAL_LIBS=$HOME/local/

To create this directory use:

.. code:: bash

    mkdir $LOCAL_LIBS

To add installed libraries to ``PATH`` use:

.. code:: bash

    export PATH=$LOCAL_LIBS/bin:$LOCAL_LIBS/usr/local/bin/:$PATH

CMake
^^^^^

Download version 3.12:

.. code:: bash

    wget https://cmake.org/files/v3.12/cmake-3.12.0-rc1.tar.gz
    tar xvf cmake-3.12.0-rc1.tar.gz
    cd cmake-3.12.0-rc1

Build and install:

.. code:: bash

    ./bootstrap --prefix=$LOCAL_LIBS
    make
    make install

For the detailed instruction see `instruction
page <https://cmake.org/install/>`__.

Zlib
^^^^

Download version 1.2.11:

.. code:: bash

    wget https://zlib.net/zlib-1.2.11.tar.gz
    tar xvf zlib-1.2.11.tar.gz
    cd zlib-1.2.11

Build and install:

.. code:: bash

    ./configure --prefix=$LOCAL_LIBS
    make
    make install

BamTools
^^^^^^^^

Clone repository, version 2.5.0:

.. code:: bash

    git clone https://github.com/pezmaster31/bamtools.git
    cd bamtools
    git reset --hard 94f072

Build and install:

.. code:: bash

    mkdir build && cd build
    cmake ../
    make
    make install DESTDIR=$LOCAL_LIBS

For the detailed instruction see `instruction
page <https://github.com/pezmaster31/bamtools/wiki/Building-and-installing>`__.

Bzip2
^^^^^

Download version 1.0.6:

.. code:: bash

    wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
    tar xvf bzip2-1.0.6.tar.gz
    cd bzip2-1.0.6

Build and install:

.. code:: bash

    make -f Makefile-libbz2_so
    make install PREFIX=$LOCAL_LIBS
    cp -a libbz2.so* $LOCAL_LIBS/lib/
    ln -s $LOCAL_LIBS/lib/libbz2.so.1.0 $LOCAL_LIBS/lib/libbz2.so

For the detailed instruction see `this
page <http://www.linuxfromscratch.org/lfs/view/stable/chapter06/bzip2.html>`__.

Boost
^^^^^

Download version 1.60:

.. code:: bash

    wget http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
    tar xzf boost_1_60_0.tar.gz
    cd boost_1_60_0

Build and install:

.. code:: bash

    ./bootstrap.sh --with-libraries=filesystem,iostreams,log,system,thread,test
    ./b2 cxxflags="-std=c++11" include="$LOCAL_LIBS/include/" search="$LOCAL_LIBS/lib/" link=shared threading=multi install --prefix=$LOCAL_LIBS

For the detailed instruction see `tutorial
page <https://www.boost.org/doc/libs/1_60_0/tools/build/tutorial.html>`__.

Dockers
~~~~~~~

Alternatively, you can use dropEst through Docker. Dockerfiles for the most
popular linux distributions are provided (see ``dropEst/dockers/``). You
can either build and run these dockers or just read dockerfiles for the
further instructions on dropEst installation for specific distribution.

To install docker on your system see `installation
instruction <https://github.com/wsargent/docker-cheat-sheet#installation>`__.

To pull pre-built CentOS-based docker from DockerHub use:

.. code:: bash

    docker pull vpetukhov/dropest:latest

Dockers for older dropEst versions are available on `DockerHub <https://hub.docker.com/r/sgosline/dropest>`__.

Or, you can build docker by hands, using the following commands:

.. code:: bash

    cd dropEst/dockers/centos7
    docker build -t dropest .

Or, for CentOS 7:

.. code:: bash

    cd dropEst/dockers/centos7
    docker build -t dropest .

To run the docker, use

.. code:: bash

    docker run --name dropest -it dropest

You can find more info about dockers at `Docker Cheat
Sheet <https://github.com/wsargent/docker-cheat-sheet>`__

Updating dropEst inside docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please note, that docker container isn't wired to a specific dropEst
version, it just builds the latest commit from the master branch of the
git repo. To update the code inside a compiled container, you need to
log into it, pull the latest version and rebuild the code:

.. code:: bash

    docker exec -it dropest /bin/bash

    cd /home/user/dropEst/build
    rm -rf ./*
    git pull origin master
    cmake .. && make

Troubleshooting
~~~~~~~~~~~~~~~

CMake can't find installed libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``cmake`` can't find one of the libraries, or you want to use some
specific versions, which are currently not in the default path, use
corresponding cmake variables: \* Boost: BOOST\_ROOT. \* BamTools:
BAMTOOLS\_ROOT. \* R: R\_ROOT. Can be found by running the
``cat(R.home())`` in R.

These variables should be set to the path to the installed library. It
can be done either by using command line options:
``cmake -D R_ROOT="path_to_r"`` or by adding the variable declaration to
the beginning of CMakeLists.txt: ``set(R_ROOT path_to_r)``.

In case you have some issues with the linker for specific library,
please build this library manually with the version of compiler, which
you're going to use for dropEst build.

Problems with std::\_\_cxx11::string
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have messages like "*(path to some library)*: undefined reference to
*(some name)* for std::\_\_cxx11::basic\_ostringstream, std::allocator >", it
means that you're trying to link a library, built with gcc < 5.0, while dropEst
is built with gcc >= 5.0. It's a compiler issue, and you have to guarantee
consistency of compiler versions by rebuilding either the library or dropEst.
For more details see `question on stackoverflow <https://stackoverflow.com/questions/33394934/converting-std-cxx11string-to-stdstring>`__.

If you have several compilers in your system, please use cmake flags
``-DCMAKE_CXX_COMPILER=(c++ compiler)`` and
``-DCMAKE_C_COMPILER=(c compiler)`` to choose a compiler. Here,
``(c++ compiler)`` and ``(c compiler)`` denotes path to the prefered
compiler version.

Boost 1.65
^^^^^^^^^^

CMake < 3.10 has known issues with boost 1.65. If you have such
combination, please try either to upgrade cmake or to downgrade boost.
