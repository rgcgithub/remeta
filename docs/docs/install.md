# Install

## Precompiled binary
Precompiled binaries are available on the [**remeta** Github](https://github.com/rgcgithub/remeta/releases).

## From source
**remeta** requires compilation using g++>=8.

It uses the following libraries:

* HTSlib: <http://www.htslib.org/download/>
* EIGEN: <https://eigen.tuxfamily.org/>
* Intel MKL (optional) <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>

### Steps to install (Ubuntu 22.04)
#### Install dependencies
```bash
$ apt-get update && apt-get install build-essential gnupg gpg-agent wget bzip2 apt-transport-https  ca-certificates git-all zlib1g-dev libboost-all-dev libz-dev libbz2-dev liblzma-dev  libcurl4-openssl-dev libssl-dev libgomp1 libdeflate-dev
```

#### Install CMake
```bash
$ wget http://cmake.org/files/cmake-v3.10.0-Linux-x86_64.sh
$ ./cmake-v3.10.0-Linux-x86_64.sh
```

#### Install HTSlib
```bash
$ wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
$ tar -xf htslib-1.20.tar.bz2
$ cd htslib-1.20
$ ./configure
$ make
$ make install
```

#### Install IntelMKL
```bash
$ wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
$ cat https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
$ echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
$ apt-get update \
$ apt-get install -y --no-install-recommends intel-oneapi-mkl-devel \
$ . /opt/intel/oneapi/setvars.sh \
$ echo "MKL_THREADING_LAYER=GNU" >> /etc/environment
```

#### Install EIGEN
```bash
$ wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
$ tar -xf eigen-3.4.0.tar.bz2 -C /usr/local/lib
```

#### Build remeta
```bash
$ git clone https://github.com/rgcgithub/remeta.git
$ cd remeta
$ cmake . -D EIGEN_PATH=/usr/local/lib/eigen-3.4.0 \
          -D CMAKE_CXX_COMPILER=g++ \
          -D MKLROOT=${MKLROOT} # set by . /opt/intel/oneapi/setvars.sh above
$ make remeta
```

## With Docker
You can also run **remeta** in docker. A docker image is provided on [Github](https://github.com/rgcgithub/remeta/).