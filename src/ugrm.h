#ifndef SNPLIB_UGRM_H
#define SNPLIB_UGRM_H

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
void CalcUGRMMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
void UnpackUGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *geno_d);
};

#endif  // SNPLIB_UGRM_H
