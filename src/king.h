#ifndef _SNPLIB_SRC_KING_H_
#define _SNPLIB_SRC_KING_H_

#include <algorithm>
#include <list>
#include <numeric>
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
std::vector<int32_t> FindUnrelatedGroup(const double *matrix,
                                        size_t num_samples, double threshold);
}  // namespace snplib

#endif  //_SNPLIB_SRC_KING_H_
