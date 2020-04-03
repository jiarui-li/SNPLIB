#ifndef SNPLIB_IBS_H
#define SNPLIB_IBS_H

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
void CalcIBSMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *ibs, size_t num_threads);
void CalcIBSConnection(uint8_t *src_geno, size_t num_src_samples,
                       uint8_t *dest_geno, size_t num_dest_samples,
                       size_t num_snps, double *connection, size_t num_threads);
};

#endif  // SNPLIB_IBS_H
