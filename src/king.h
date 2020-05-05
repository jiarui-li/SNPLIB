#ifndef _SNPLIB_SRC_KING_H_
#define _SNPLIB_SRC_KING_H_

#include <algorithm>
#include <thread>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include "snp.h"

namespace snplib {
void CalcKINGMatrix(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
}

#endif  //_SNPLIB_SRC_KING_H_
