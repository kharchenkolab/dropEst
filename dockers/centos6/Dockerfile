FROM library/centos:6.8
MAINTAINER Viktor Petukhov "viktor.s.petuhov@ya.ru"

RUN \
  yum -y install epel-release && \
  yum -y install \
    boost-devel \
    centos-release-scl \
    cmake \
    cmake3 \
    git \
    libcurl-devel \
    libpng-devel \
    libX11-devel \
    mesa-libGL-devel \
    mesa-libGLU-devel \
    openssl-devel \
    R-3.4.1 \
    vim \
    wget \
    xorg-x11-server-Xorg \
    xorg-x11-server-Xorg-devel && \
  yum install -y devtoolset-4-gcc* && \
  scl enable devtoolset-4 bash

RUN \
  cd /root && \
  git clone git://github.com/pezmaster31/bamtools.git && \
  mkdir bamtools/build && \
  cd bamtools/build && \
  cmake3 .. && make && make install

RUN useradd -m user
USER user

ENV PATH="/opt/rh/devtoolset-4/root/usr/bin/:$PATH"

ENTRYPOINT ["scl", "enable", "devtoolset-4", "bash"]
WORKDIR /home/user

RUN \
  wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz && \
  tar -xvzf boost_1_55_0.tar.gz && \
  cd boost_1_55_0 && \
  ./bootstrap.sh --with-libraries=filesystem,iostreams,log,system,thread,test && \
  mkdir ~/local

RUN \
  cd boost_1_55_0 && \
  ./b2 cxxflags=-std=c++11 link=shared threading=multi install --prefix=/home/user/local

RUN \
  cd && \
  git clone https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/3.4 && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools"), lib = "~/R/x86_64-redhat-linux-gnu-library/3.4")' && \
  R -e 'devtools::install_local("~/dropEst/dropestr/", dependencies=T)' && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown"))'

RUN \
  mkdir dropEst/build && \
  cd dropEst/build && \
  cmake -D BOOST_ROOT=~/local/ .. && \
  make

ENV PATH=/home/user/dropEst/build:$PATH
