FROM library/centos:7.4.1708
MAINTAINER Viktor Petukhov "viktor.s.petuhov@ya.ru"

RUN \
  yum install -y wget && \
  wget http://repo.enetres.net/enetres.repo -O /etc/yum.repos.d/enetres.repo && \
  wget ftp://fr2.rpmfind.net/linux/Mandriva/official/2010.0/x86_64/media/main/release/lib64icu42-4.2.1-1mdv2010.0.x86_64.rpm && \
  rpm -ivh lib64icu42-4.2.1-1mdv2010.0.x86_64.rpm && \
  yum -y install epel-release && \
  yum -y install \
    boost-devel \
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
    xorg-x11-server-Xorg \
    xorg-x11-server-Xorg-devel

RUN \ 
    cd /root && \
    git clone git://github.com/pezmaster31/bamtools.git && \
    mkdir -p bamtools/build && \
    cd bamtools/build && \
    cmake3 .. && make && make install

RUN useradd -m user
USER user

ENTRYPOINT ["/bin/bash"]
WORKDIR /home/user

RUN \
  git clone https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/3.4 && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools"), lib = "~/R/x86_64-redhat-linux-gnu-library/3.4")' && \
  R -e 'devtools::install_local("~/dropEst/dropestr/", dependencies=T)' && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown"))'

RUN \
  cd && \
  mkdir -p ~/dropEst/build && \
  cd ~/dropEst/build && \
  cmake3 .. && make

ENV PATH=/home/user/dropEst/build:$PATH
