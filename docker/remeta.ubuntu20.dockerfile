FROM ubuntu:20.04 AS builder

ENV CMAKE_VERSION 3.10
ENV CMAKE_VERSION_PATCH 0
ENV HTSLIB_VERSION 1.20
ENV TMP_DIR /tmp

COPY .git ${TMP_DIR}/remeta/.git
COPY lib ${TMP_DIR}/remeta/lib
COPY src ${TMP_DIR}/remeta/src
COPY CMakeLists.txt ${TMP_DIR}/remeta/CMakeLists.txt
COPY VERSION ${TMP_DIR}/remeta/VERSION

ADD http://cmake.org/files/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}.${CMAKE_VERSION_PATCH}-Linux-x86_64.sh cmake_install.sh
ADD https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB intel_key.PUB

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    g++ \
    make \
    gnupg \
    gpg-agent \
    wget \
    bzip2 \
    apt-transport-https \
    ca-certificates \
    git-all \
    zlib1g-dev \
    libboost-all-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgomp1 \
    libdeflate-dev \
    && sh -c 'cat intel_key.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null' \
    && sh -c 'echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list' \
    && apt-get update \
    && apt-get install -y --no-install-recommends intel-oneapi-mkl-devel \
    && . /opt/intel/oneapi/setvars.sh \
    && echo "MKL_THREADING_LAYER=GNU" >> /etc/environment \
    && sh cmake_install.sh --prefix=/usr/local --skip-license --exclude-subdir \
    && mkdir -p $TMP_DIR \
    && wget -q --no-check-certificate https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2 \
    && tar -xf eigen-3.4.0.tar.bz2 -C /usr/local/lib \
    && rm eigen-3.4.0.tar.bz2 \
    && wget -q --no-check-certificate https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 \
    && tar -xf htslib-$HTSLIB_VERSION.tar.bz2 -C $TMP_DIR \
    && rm htslib-$HTSLIB_VERSION.tar.bz2 \
    && cd $TMP_DIR/htslib-$HTSLIB_VERSION/ \
    && ./configure \
    && make \
    && make install \
    && cp tabix /usr/local/bin/tabix \
    && cd $TMP_DIR/remeta/lib/pgenlib \
    && make clean \
    && cd $TMP_DIR/remeta/lib/faddeeva \
    && make clean \
    && cd $TMP_DIR/remeta/lib/qfc \
    && make clean \
    && cd $TMP_DIR/remeta \
    && cmake -D EIGEN_PATH=/usr/local/lib/eigen-3.4.0 \
             -D CMAKE_CXX_COMPILER=g++ \
             -D MKLROOT=${MKLROOT} \
             . \
    && make remeta \
    && cp remeta /usr/local/bin/remeta \
    && cd \
    && rm -rf $TMP_DIR

FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    g++ \
    make \
    gnupg \
    wget \
    apt-transport-https \
    ca-certificates \
    git-all \
    zlib1g-dev \
    libboost-all-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgomp1 \
    libdeflate-dev

COPY --from=builder /usr/local/bin/remeta /usr/local/bin/
COPY --from=builder /usr/local/bin/tabix /usr/local/bin/