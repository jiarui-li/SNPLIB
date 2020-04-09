#ifndef _SNPLIB_SRC_MATH_LIB_H_
#define _SNPLIB_SRC_MATH_LIB_H_

#ifdef USE_MKL
#include <mkl.h>
#endif
#ifdef USE_OPENBLAS
#include <cblas.h>
#include <lapacke.h>
#endif

#endif  //_SNPLIB_SRC_MATH_LIB_H_
