#ifndef SNPLIB_SRC_IBS_H_
#define SNPLIB_SRC_IBS_H_

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
void CalcIBSMatrix(const uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *matrix, size_t num_threads);
void CalcIBSConnection(const uint8_t *src_geno, size_t num_src_samples,
                       uint8_t *dest_geno, size_t num_dest_samples,
                       size_t num_snps, double *connection, size_t num_threads);
}  // namespace snplib

#endif  // SNPLIB_SRC_IBS_H_
