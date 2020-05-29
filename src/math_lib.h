#ifndef _SNPLIB_SRC_MATH_LIB_H_
#define _SNPLIB_SRC_MATH_LIB_H_

#ifdef USE_MKL
#include <mkl.h>
#define set_num_threads mkl_set_num_threads
#endif
#ifdef USE_OPENBLAS
#ifdef _MSC_VER
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif
#include <cblas.h>
#include <lapacke.h>
#define set_num_threads openblas_set_num_threads
#endif

#endif //_SNPLIB_SRC_MATH_LIB_H_
