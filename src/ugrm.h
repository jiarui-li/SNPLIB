#ifndef SNPLIB_SRC_UGRM_H_
#define SNPLIB_SRC_UGRM_H_

#include <algorithm>
#include <thread>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include "transpose_geno.h"

namespace snplib {
void CalcUGRMMatrix(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
}

#endif  // SNPLIB_SRC_UGRM_H_
