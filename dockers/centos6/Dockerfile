FROM library/centos:6.8
MAINTAINER Viktor Petukhov "viktor.s.petuhov@ya.ru"

RUN \
  yum -y install epel-release && \
  yum -y install \
    bzip2-devel \
    centos-release-scl \
    cmake \
    cmake3 \
    git \
    libcurl-devel \
    openssl-devel \
    R \
    vim \
    wget

RUN yum -y install devtoolset-6-gcc*

RUN source /opt/rh/devtoolset-6/enable

RUN \
  cd /root && \
  wget https://github.com/jgm/pandoc/releases/download/2.1.3/pandoc-2.1.3-linux.tar.gz && \
  tar xvzf pandoc-2.1.3-linux.tar.gz --strip-components 1 -C /usr/local/ && \
  git clone git://github.com/pezmaster31/bamtools.git && \
  mkdir bamtools/build && \
  cd bamtools/build && \
  cmake3 .. && make && make install

RUN useradd -m user
USER user

WORKDIR /home/user

RUN echo "source /opt/rh/devtoolset-6/enable" > ~/.bashrc

RUN \
  source ~/.bashrc && \
  mkdir ~/local && \
  wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz && \
  tar -xvzf boost_1_55_0.tar.gz && \
  cd boost_1_55_0 && \
  ./bootstrap.sh --with-libraries=filesystem,iostreams,log,system,thread,test && \
  ./b2 cxxflags=-std=c++11 link=shared threading=multi install --prefix=/home/user/local

RUN \
  source ~/.bashrc && \
  cd && \
  git clone https://github.com/VPetukhov/ks.git && \
  git clone https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/3.4 && \
  echo '.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths()))' > .Rprofile && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools"))' && \
  R -e 'devtools::install_local("~/ks/")' && \
  R -e 'devtools::install_local("~/dropEst/dropestr/", dependencies=T)' && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown"))'

RUN \
  source ~/.bashrc && \
  mkdir dropEst/build && \
  cd dropEst/build && \
  cmake -D BOOST_ROOT=~/local/ .. && \
  make

ENV PATH=/home/user/dropEst/build:$PATH
