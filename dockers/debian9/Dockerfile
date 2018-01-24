FROM rocker/r-ver:3.4
MAINTAINER Viktor Petukhov "viktor.s.petuhov@ya.ru"

RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  less \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libz-dev \
  nano

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools"))' && \
  R -e 'devtools::install_github("hms-dbmi/dropEst/dropestr", dependencies = T)' && \
  R -e 'install.packages(c("Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown"))'

RUN useradd -m user
USER user
ENTRYPOINT ["/bin/bash"]
WORKDIR "/home/user"

RUN \
  git clone https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p dropEst/build && \
  cd dropEst/build && \
  cmake ../ && \
  make

ENV PATH=/home/user/dropEst/build:$PATH \
  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH
