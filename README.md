# SNPLIB


## Introduction
SNPLIB is a MATLAB library for Quantitative Genetics




## Prerequisit

SNPLIB is accelerated by BLAS/LAPACK compatible matrix libraries.  Before build SNPLIB, it is required to install following libraries.

### Intel&reg; Math Kernel Library

>Intel Math Kernel Library (Intel MKL) is a library of optimized math routines for science, engineering, and financial applications. The routines in MKL are hand-optimized specifically for Intel processors. The library supports Intel processors and is available for Windows, Linux and macOS operating systems.

You can use Intel&reg; freely under the [Intel Simplified Software License](https://software.intel.com/en-us/license/intel-simplified-software-license).

1. [Download Intel&reg; MKL](https://software.intel.com/en-us/mkl)
2. Install Intel&reg; MKL following the installation instructions

### OpenBLAS
>OpenBLAS is an optimized BLAS library based on GotoBLAS2 1.13 BSD version.

OpenBLAS is an alternative open-source optimized math library. It is recommended to use OpenBLAS instead of Intel&reg; MKL when you are using non-Intel processors (e.g., AMD Ryzen). 

It is recommend to build OpenBLAS from source following [Installation Guide](https://github.com/xianyi/OpenBLAS/wiki/Installation-Guide). It is also possible to download the [binary packages](https://sourceforge.net/projects/openblas/files/).

## Compiling SNPLIB

SNPLIB is written in C++ and provide the interfaces for Python 3 (not 2) and MATLAB. To compile the core library of SNPLIB, a C/C++ compiler is needed. 

### Python Package

To install Python 3 interfaces, you will need CMake. To install the Python 3 interface, type:

```bash
python setup.py install
```

The script will ask you to install SNPLIB with MKL or OpenBLAS and where they can be found. Here, using Intel MKL for example:

>Use Intel MKL(m) or OpenBLAS(o)?m

>The directory of Intel MKL [C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows]:


### MATLAB Package

To install the SNPLIB in MATLAB&reg;, you need first config the mex:

```matlab
mex -setup c;
mex -setup c++;
```
and then

```matlab
%% using Intel MKL
setup_mkl;
%% using OpenBLAS
setup_openblas;
```

## License

[BSD @ Jiarui Li](LICENSE)

## Contributors

* Jiarui Li (jiarui_li@outlook.com)

Contact me freely if you'd like to contribute and I'll add you to the repository.
