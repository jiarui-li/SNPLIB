#ifndef SNPLIB_KING_H
#define SNPLIB_KING_H

#include <algorithm>
#include <thread>
#include <vector>

#include "snp.h"
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

extern "C" {
void CalcKINGMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
};

#endif  // SNPLIB_KING_H
